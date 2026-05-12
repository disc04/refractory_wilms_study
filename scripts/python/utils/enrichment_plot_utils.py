"""
enrichment_plot_utils.py — plotting helpers for pathway enrichment analysis.

Functions
---------
plot_ora_dotplot   — multi-panel ORA dot plot (one column per database)
plot_gsea_nes      — GSEA NES horizontal bar chart (one panel per database)
plot_ora_heatmap   — ORA summary heatmap for any database (Hallmarks, KEGG, etc.)
"""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

from utils.configuration_pathways import (PATHWAY_DATABASES, PATHWAY_DB_SHORT,
                                          PATHWAY_LABEL_LONG, PATHWAY_STEP_LABEL,
                                          ORA_PADJ_THRESH, GSEA_NES_THRESH,
                                          GSEA_FDR_THRESH, TOP_N_PLOT)
from utils.data_utils import savefig


# ── Internal helpers ──────────────────────────────────────────────────────────

def _count_from_genes(genes_val) -> int:
    """
    Count overlapping genes from gseapy's 'genes' column.
    gseapy ≥1.1 stores them as a semicolon-separated string
    (e.g. 'GENE1;GENE2;GENE3') rather than a fraction string ('3/200').
    """
    if not genes_val or (isinstance(genes_val, float) and np.isnan(genes_val)):
        return 0
    return len([g for g in str(genes_val).split(';') if g.strip()])


# ── ORA dot plot ──────────────────────────────────────────────────────────────

def plot_ora_dotplot(ora_all: pd.DataFrame,
                     gene_set_labels: list[str],
                     out_path: str,
                     panel_title: str,
                     top_n: int = TOP_N_PLOT) -> None:
    """
    Multi-panel ORA dot plot.

    Layout  : rows = gene set labels, columns = databases.
    x-axis  : Odds Ratio.
    y-axis  : pathway term name (truncated to 55 chars).
    Colour  : −log10(adj. P-value).
    Size    : overlap gene count.

    Parameters
    ----------
    ora_all         : combined ORA results DataFrame (all gene sets, all databases)
    gene_set_labels : ordered list of gene_set_label values to show as rows
    out_path        : save path without extension (PDF + PNG saved by savefig)
    panel_title     : figure suptitle
    top_n           : maximum terms per panel cell (top by significance)
    """
    sub = ora_all[
        ora_all['gene_set_label'].isin(gene_set_labels) &
        (ora_all['adj_pval'] < ORA_PADJ_THRESH)
    ].copy()

    if sub.empty:
        print(f"  [skip] {panel_title}: no significant ORA terms")
        return

    sub['overlap_count']  = sub['genes'].apply(_count_from_genes)
    sub['neg_log10_padj'] = -np.log10(sub['adj_pval'].clip(lower=1e-300))

    databases = [PATHWAY_DB_SHORT[d] for d in PATHWAY_DATABASES]
    n_labels  = len(gene_set_labels)
    n_dbs     = len(databases)
    fig_h     = max(4, n_labels * top_n * 0.18 + 1.5)
    fig_w     = n_dbs * 4.5

    fig, axes = plt.subplots(n_labels, n_dbs,
                             figsize=(fig_w, fig_h),
                             squeeze=False,
                             constrained_layout=True)
    fig.suptitle(panel_title, fontsize=11, fontweight='bold')

    vmax = sub['neg_log10_padj'].quantile(0.95)
    norm = mcolors.Normalize(vmin=0, vmax=max(vmax, 1))
    cmap = plt.cm.RdYlBu_r

    for row_i, label in enumerate(gene_set_labels):
        for col_j, db_short in enumerate(databases):
            ax = axes[row_i][col_j]
            panel = sub[
                (sub['gene_set_label'] == label) &
                (sub['database_short']  == db_short)
            ].nsmallest(top_n, 'adj_pval')

            if row_i == 0:
                ax.set_title(db_short, fontsize=9, fontweight='bold')
            if col_j == 0:
                ax.set_ylabel(PATHWAY_LABEL_LONG.get(label, label),
                              fontsize=7, labelpad=4)

            if panel.empty:
                ax.text(0.5, 0.5, 'no sig. terms',
                        ha='center', va='center',
                        transform=ax.transAxes, fontsize=7, color='grey')
                ax.set_xticks([])
                ax.set_yticks([])
                continue

            panel = panel.sort_values('adj_pval', ascending=False)
            y_pos = range(len(panel))
            ax.scatter(
                panel['odds_ratio'],
                y_pos,
                c=panel['neg_log10_padj'],
                s=panel['overlap_count'] * 4 + 10,
                cmap=cmap,
                norm=norm,
                edgecolors='white',
                linewidths=0.4,
                zorder=3,
            )
            ax.set_yticks(list(y_pos))
            ax.set_yticklabels(
                [t[:55] + '…' if len(t) > 55 else t for t in panel['term']],
                fontsize=6,
            )
            ax.set_xlabel('Odds Ratio', fontsize=7)
            ax.tick_params(axis='x', labelsize=7)
            ax.grid(axis='x', linestyle='--', alpha=0.3, linewidth=0.5)
            ax.set_axisbelow(True)

    # Shared colorbar — flat list so matplotlib carves space from all axes correctly
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(),
                        fraction=0.015, pad=0.02, shrink=0.6)
    cbar.set_label('−log₁₀(adj. P)', fontsize=8)

    savefig(fig, out_path)


# ── GSEA NES bar chart ────────────────────────────────────────────────────────

def plot_gsea_nes(gsea_all: pd.DataFrame,
                  step: str,
                  out_path: str) -> None:
    """
    Horizontal NES bar chart for one GSEA step, one panel per database.

    Only significant terms are shown (|NES| ≥ GSEA_NES_THRESH, FDR < GSEA_FDR_THRESH).
    Bars are coloured blue (negative NES) or red (positive NES).
    FDR q-values are annotated on each bar.

    Parameters
    ----------
    gsea_all : combined GSEA results DataFrame (all steps, all databases)
    step     : which step to plot, e.g. 'step1', 'step2', 'step3'
    out_path : save path without extension
    """
    sub = gsea_all[
        (gsea_all['step'] == step) &
        (gsea_all['NES'].abs() >= GSEA_NES_THRESH) &
        (gsea_all['fdr'] < GSEA_FDR_THRESH)
    ].copy()

    if sub.empty:
        print(f"  [skip] GSEA NES plot {step}: no significant terms")
        return

    databases = [PATHWAY_DB_SHORT[d] for d in PATHWAY_DATABASES]
    fig, axes = plt.subplots(1, len(databases),
                             figsize=(len(databases) * 4.2, 6),
                             squeeze=False)
    fig.suptitle(
        f'GSEA — {PATHWAY_STEP_LABEL[step]}\n'
        f'Significant: |NES| ≥ {GSEA_NES_THRESH}, FDR < {GSEA_FDR_THRESH}',
        fontsize=10, fontweight='bold',
    )

    for col_j, db_short in enumerate(databases):
        ax = axes[0][col_j]

        panel_all = sub[sub['database_short'] == db_short]
        pos_terms = panel_all[panel_all['NES'] > 0].nlargest(TOP_N_PLOT, 'NES')
        neg_terms = panel_all[panel_all['NES'] < 0].nsmallest(TOP_N_PLOT, 'NES')
        panel = pd.concat([neg_terms, pos_terms]).sort_values('NES')

        ax.set_title(db_short, fontsize=9, fontweight='bold')

        if panel.empty:
            ax.text(0.5, 0.5, 'no sig. terms',
                    ha='center', va='center',
                    transform=ax.transAxes, fontsize=8, color='grey')
            continue

        colors = ['#1f77b4' if n < 0 else '#d62728' for n in panel['NES']]
        y_pos  = range(len(panel))
        bars   = ax.barh(list(y_pos), panel['NES'], color=colors,
                         edgecolor='white', linewidth=0.4)
        ax.set_yticks(list(y_pos))
        ax.set_yticklabels(
            [t[:52] + '…' if len(t) > 52 else t for t in panel['term']],
            fontsize=6,
        )
        ax.axvline(0,                color='black', linewidth=0.6)
        ax.axvline( GSEA_NES_THRESH, color='grey',  linewidth=0.5, linestyle='--')
        ax.axvline(-GSEA_NES_THRESH, color='grey',  linewidth=0.5, linestyle='--')
        ax.set_xlabel('NES', fontsize=8)
        ax.tick_params(axis='x', labelsize=7)
        ax.grid(axis='x', linestyle='--', alpha=0.3, linewidth=0.5)
        ax.set_axisbelow(True)

        for bar, (_, row) in zip(bars, panel.iterrows()):
            fdr_str = f"q={row['fdr']:.3f}"
            x_pos   = row['NES'] + (0.05 if row['NES'] > 0 else -0.05)
            ha      = 'left' if row['NES'] > 0 else 'right'
            ax.text(x_pos, bar.get_y() + bar.get_height() / 2,
                    fdr_str, va='center', ha=ha, fontsize=5, color='#444')

    fig.tight_layout()
    savefig(fig, out_path)


# ── ORA summary heatmap ───────────────────────────────────────────────────────

_DB_TITLES = {
    'Hallmarks': 'MSigDB Hallmarks',
    'KEGG':      'KEGG Pathways',
    'GO_BP':     'GO Biological Process',
    'TRRUST':    'TRRUST Transcription Factors',
}


def plot_ora_heatmap(ora_all: pd.DataFrame,
                     out_path: str,
                     database_short: str = 'Hallmarks',
                     max_terms: int = 40) -> None:
    """
    Heatmap of −log10(adj_pval) for a chosen database across all ORA gene sets.

    Rows  : pathway terms significant in ≥1 gene set, sorted by total signal
            descending, capped at max_terms.
    Columns: gene sets in study-logic order (UP then DOWN for each comparison).
    Values : −log10(adj_pval); 0 (white) = not significant in that gene set.

    Parameters
    ----------
    database_short : short label as stored in ora_all['database_short'];
                     e.g. 'Hallmarks', 'KEGG', 'GO_BP', 'TRRUST'
    max_terms      : maximum rows shown (top by summed −log10 P across gene sets)
    """
    hm = ora_all[
        (ora_all['database_short'] == database_short) &
        (ora_all['adj_pval'] < ORA_PADJ_THRESH)
    ].copy()

    if hm.empty:
        print(f"  [skip] {database_short} heatmap: no significant terms")
        return

    hm['neg_log10_padj'] = -np.log10(hm['adj_pval'].clip(lower=1e-300))

    col_order = list(PATHWAY_LABEL_LONG.keys())
    pivot = hm.pivot_table(
        index='term',
        columns='gene_set_label',
        values='neg_log10_padj',
        aggfunc='max',
        fill_value=0,
    ).reindex(
        columns=[c for c in col_order if c in hm['gene_set_label'].unique()],
        fill_value=0,
    )

    # Sort by total signal, cap rows
    pivot = pivot.loc[pivot.sum(axis=1).sort_values(ascending=False).index]
    pivot = pivot.head(max_terms)

    def _clean_term(t: str) -> str:
        if database_short == 'Hallmarks':
            return t.replace('HALLMARK_', '').replace('_', ' ').title()
        return t[:60] + '…' if len(t) > 60 else t

    pivot.index = [_clean_term(t) for t in pivot.index]
    col_labels  = [PATHWAY_LABEL_LONG.get(c, c).split(' (')[0]
                   for c in pivot.columns]

    n_rows = len(pivot)
    n_cols = len(pivot.columns)
    fig_h  = max(4, n_rows * 0.30 + 2)
    fig_w  = max(6, n_cols * 1.8 + 2)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)
    sns.heatmap(
        pivot, ax=ax,
        cmap='YlOrRd',
        linewidths=0.3,
        linecolor='#e0e0e0',
        xticklabels=col_labels,
        yticklabels=True,
        cbar_kws={'label': '−log₁₀(adj. P)', 'shrink': 0.6},
        vmin=0,
    )
    db_label = _DB_TITLES.get(database_short, database_short)
    ax.set_title(
        f'{db_label} — ORA across all gene sets\n'
        f'(top {n_rows} significant terms, adj. P < {ORA_PADJ_THRESH})',
        fontsize=10, fontweight='bold',
    )
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.tick_params(axis='x', rotation=40, labelsize=10)
    ax.tick_params(axis='y', labelsize=10)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontweight('bold')

    savefig(fig, out_path)
