"""
de_plot_utils.py — shared plotting helpers for the Wilms Tumor DE analysis.

QC figure functions
-------------------
plot_library_sizes — horizontal bar chart of per-sample library sizes
plot_count_density — KDE density of log2-normalised counts across samples
plot_pca — PCA scatter coloured by sample_type and histology
plot_sample_correlation — Pearson correlation heatmap (all genes)

DE figure functions
-------------------
plot_volcano — volcano plot for any two-group DE result
plot_ma — MA plot for any two-group DE result
plot_heatmap — z-score heatmap for top N DEGs
plot_gene_bar  — horizontal bar chart for a custom gene panel (e.g. Wilms genes)

All functions:
  - accept a `palette` dict mapping direction/group → hex color
  - save PDF + PNG via savefig()
  - are designed to be called from any analysis script
"""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from sklearn.decomposition import PCA
import os
from utils.data_utils import savefig
from configuration import GDC_SAMPLE_TYPE_COLORS, PCA_GROUP_COLORS, PRIMARY_HISTOLOGY_COLORS

def plot_library_sizes(
    counts: pd.DataFrame,
    annot: pd.DataFrame,
    out_path: str,
    palette: dict | None = None,
    sample_col: str = 'sample_id',
    type_col: str = 'sample_type',
) -> None:
    """
    Horizontal bar chart of per-sample library sizes (millions of reads).

    Parameters
    ----------
    counts     : raw count matrix (genes × samples)
    annot      : sample annotation DataFrame
    out_path   : save path without extension
    palette    : sample_type → colour dict; defaults to GDC_SAMPLE_TYPE_COLORS
    sample_col : column in annot holding sample identifiers
    type_col   : column in annot holding sample-type labels (used for colouring)
    """
    palette = palette or GDC_SAMPLE_TYPE_COLORS

    lib = (counts.sum(axis=0) / 1e6).rename('lib_size')
    df  = lib.to_frame().join(annot.set_index(sample_col)[[type_col]])
    df  = df.sort_values('lib_size')
    colors = [palette.get(t, '#999') for t in df[type_col]]

    fig, ax = plt.subplots(figsize=(5, 8))
    ax.barh(range(len(df)), df['lib_size'], color=colors, linewidth=0)
    ax.set_yticks([])
    ax.set_xlabel('Mapped reads (millions)')
    ax.set_title('Library Sizes')
    handles = [mpatches.Patch(color=v, label=k)
               for k, v in palette.items() if k in df[type_col].values]
    ax.legend(handles=handles, fontsize=8, loc='lower right')
    fig.tight_layout()
    savefig(fig, out_path)


def plot_count_density(
    log2_matrix: pd.DataFrame,
    annot: pd.DataFrame,
    out_path: str,
    palette: dict | None = None,
    sample_col: str = 'sample_id',
    type_col: str = 'sample_type',
) -> None:
    """
    KDE density of log2-normalised counts, one line per sample.

    Normal samples are drawn with solid thick lines; tumours as thin semi-transparent.
    Samples are rendered in four layered passes so Normal and rarer types sit on top.

    Parameters
    ----------
    log2_matrix : log2-normalised expression matrix (genes × samples)
    annot       : sample annotation DataFrame
    out_path    : save path without extension
    palette     : sample_type → colour dict; defaults to GDC_SAMPLE_TYPE_COLORS
    sample_col  : column in annot holding sample identifiers
    type_col    : column in annot holding sample-type labels
    """
    palette = palette or GDC_SAMPLE_TYPE_COLORS
    type_map = annot.set_index(sample_col)[type_col]
    normal_key = 'Solid Tissue Normal'

    fig, ax = plt.subplots(figsize=(7, 4))
    # Draw in layered passes: primary tumors first (background), then special groups on top
    for pass_label, lw, alpha in [
        ('Primary Tumor',       0.8, 0.50),
        ('Recurrent Tumor',     1.0, 0.70),
        ('Metastatic',          1.0, 0.70),
        ('Solid Tissue Normal', 1.1, 0.70),
    ]:
        for sid in log2_matrix.columns:
            stype = type_map.get(sid, 'Unknown')
            if stype != pass_label:
                continue
            color = palette.get(stype, '#999')
            log2_matrix[sid].plot.kde(
                ax=ax, color=color, alpha=alpha, linewidth=lw)

    handles = [mpatches.Patch(color=palette.get(k, '#999'), label=k)
               for k in [normal_key, 'Primary Tumor', 'Recurrent Tumor', 'Metastatic'] if k in type_map.values]
    ax.legend(handles=handles, fontsize=8)
    ax.set_xlabel('log2(normalized counts + 1)')
    ax.set_ylabel('Density')
    ax.set_title('Count Density')
    fig.tight_layout()
    savefig(fig, out_path)


def _make_pca_group(row: pd.Series,
                    type_col: str = 'sample_type',
                    hist_col: str = 'histology') -> str:
    """
    Map sample_type + histology columns to a combined PCA group label.

    Parameters
    ----------
    row      : single annotation row (pd.Series)
    type_col : column name for sample type
    hist_col : column name for histology

    Returns
    -------
    One of: 'Normal', 'Primary FHWT', 'Primary DAWT', 'Recurrent', 'Metastatic'
    """
    t = row[type_col]
    if t == 'Solid Tissue Normal':
        return 'Normal'
    if t == 'Primary Tumor':
        return f"Primary {row[hist_col]}"
    if t == 'Recurrent Tumor':
        return 'Recurrent'
    if t == 'Metastatic':
        return 'Metastatic'
    return t


def plot_pca(
    log2_matrix: pd.DataFrame,
    annot: pd.DataFrame,
    out_path: str,
    sample_col: str = 'sample_id',
    type_col: str = 'sample_type',
    hist_col: str = 'histology',
    qc_flag_col: str | None = 'qc_flag',
    palette: dict | None = None,
    fit_exclude_types: list[str] | None = None,
) -> None:
    """
    Single-panel PCA coloured by combined group (Normal / Primary FHWT /
    Primary DAWT / Recurrent / Metastatic).

    PCA is fitted on all samples EXCEPT those whose sample_type is in
    `fit_exclude_types` (default: ['Recurrent Tumor', 'Metastatic']).
    All samples are then projected and plotted.
    QC-flagged samples receive an orange ring overlay.

    Parameters
    ----------
    log2_matrix       : log2-normalised expression matrix (genes × samples), all samples
    annot             : annotation DataFrame including all samples
    out_path          : save path without extension
    sample_col        : column in annot holding sample identifiers
    type_col          : column in annot holding sample-type labels
    hist_col          : column in annot holding histology labels
    qc_flag_col       : column name for QC flag (bool); None to skip overlay
    palette           : group → colour dict; defaults to PCA_GROUP_COLORS
    fit_exclude_types : sample types excluded from PCA fitting
    """
    palette = palette or PCA_GROUP_COLORS
    fit_exclude_types = fit_exclude_types or ['Recurrent Tumor']

    col_list = list(log2_matrix.columns)
    ann = annot.set_index(sample_col).reindex(col_list)
    ann['_group'] = ann.apply(
        lambda r: _make_pca_group(r, type_col, hist_col), axis=1)

    # Fit PCA on primary + normal samples only
    fit_mask = ~ann[type_col].isin(fit_exclude_types)
    fit_ids  = ann[fit_mask].index.tolist()
    pca_obj  = PCA(n_components=2, random_state=42)
    pca_obj.fit(log2_matrix[fit_ids].T.values)
    pct = pca_obj.explained_variance_ratio_ * 100

    # Project ALL samples
    coords = pca_obj.transform(log2_matrix.T.values)

    # Group order (rarer groups drawn last → on top)
    group_order  = ['Primary FHWT', 'Primary DAWT', 'Normal', 'Recurrent', 'Metastatic']
    group_counts = ann['_group'].value_counts()

    fig, ax = plt.subplots(figsize=(9, 8))

    for grp in group_order:
        mask = ann['_group'] == grp
        if not mask.any():
            continue
        idx   = ann[mask].index.tolist()
        cidx  = [col_list.index(s) for s in idx]
        color = palette.get(grp, '#999')
        is_special = grp in ('Recurrent', 'Metastatic')
        n = group_counts.get(grp, 0)
        ax.scatter(
            coords[cidx, 0], coords[cidx, 1],
            color=color,
            s=130 if is_special else 100,
            alpha=1.0 if is_special else 0.78,
            edgecolors='white' if not is_special else 'k',
            linewidths=0.5 if not is_special else 0.8,
            zorder=3 if is_special else 2,
            label=f'{grp} (n={n})',
        )

    # QC-flagged ring overlay
    if qc_flag_col and qc_flag_col in ann.columns:
        flagged_ids = ann[ann[qc_flag_col].astype(bool)].index.tolist()
        if flagged_ids:
            fidx = [col_list.index(s) for s in flagged_ids if s in col_list]
            ax.scatter(coords[fidx, 0], coords[fidx, 1],
                       s=300, facecolors='none',
                       edgecolors='#FF7043', linewidths=1.8,
                       zorder=4, label=f'QC flagged (n={len(fidx)})')

    ax.set_xlabel(f'PC1 ({pct[0]:.1f}%)', fontsize=15, labelpad=8, fontweight='bold')
    ax.set_ylabel(f'PC2 ({pct[1]:.1f}%)', fontsize=15, labelpad=8, fontweight='bold')
    ax.tick_params(labelsize=11)
    n_genes = log2_matrix.shape[0]
    ax.set_title(
        f'PCA — {n_genes:,} expressed genes\n'
        f'(fit: Primary + Normal; projected: all {log2_matrix.shape[1]} samples)',
        fontsize=13, pad=10, fontweight='bold')
    ax.legend(fontsize=11, framealpha=0.9, markerscale=0.8, prop={'weight': 'bold'})
    ax.axhline(0, color='lightgrey', lw=0.5, zorder=0)
    ax.axvline(0, color='lightgrey', lw=0.5, zorder=0)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    fig.tight_layout()
    savefig(fig, out_path)


def plot_sample_correlation(
    log2_matrix: pd.DataFrame,
    annot: pd.DataFrame,
    out_path: str,
    outlier_thresh: float = 0.90,
    sample_col: str = 'sample_id',
    type_col: str = 'sample_type',
    hist_col: str = 'histology',
) -> list[str]:
    """
    Pearson correlation heatmap (all genes), samples sorted Normal | DAWT | FHWT.

    Group boundaries are marked with white separator lines in data coordinates.
    Within-tumour outlier detection flags samples whose mean within-tumour
    Pearson r falls below `outlier_thresh`.

    Parameters
    ----------
    log2_matrix     : log2-normalised expression matrix (genes × samples)
    annot           : sample annotation DataFrame
    out_path        : save path without extension
    outlier_thresh  : mean within-tumour Pearson r below which a sample is flagged
    sample_col      : column in annot holding sample identifiers
    type_col        : column in annot holding sample-type labels
    hist_col        : column in annot holding histology labels

    Returns
    -------
    List of sample_id strings for QC-flagged outliers (may be empty).
    """
    cor_mat = log2_matrix.corr(method='pearson')
    ann = annot.set_index(sample_col).reindex(cor_mat.columns)

    # ── Within-tumor outlier check ────────────────────────────────────────────
    tumor_ids = ann[ann[type_col] == 'Primary Tumor'].index.tolist()
    tumor_cor = cor_mat.loc[tumor_ids, tumor_ids]
    mean_tc   = tumor_cor.mean(axis=1)
    flagged   = mean_tc[mean_tc < outlier_thresh].sort_values().index.tolist()
    if flagged:
        print(f"  [warn] Within-tumor mean Pearson r < {outlier_thresh}:")
        for sid in flagged:
            print(f"    {sid}  r={mean_tc[sid]:.4f}  {ann.loc[sid, hist_col]}")
    else:
        print(f"  Correlation QC passed "
              f"(within-tumor r: {mean_tc.min():.3f}–{mean_tc.max():.3f})")

    # ── Sort: Normal first, then DAWT, then FHWT ─────────────────────────────
    ann_sorted = ann.sort_values([type_col, hist_col], ascending=[False, True])
    order      = ann_sorted.index.tolist()
    cor_sort   = cor_mat.loc[order, order]
    n          = len(order)

    # ── Group boundary positions (in heatmap data coordinates, 0…n) ──────────
    n_normal = int((ann_sorted[type_col] == 'Solid Tissue Normal').sum())
    n_dawt   = int(((ann_sorted[type_col] == 'Primary Tumor') &
                    (ann_sorted[hist_col] == 'DAWT')).sum())
    separators = [pos for pos in [n_normal, n_normal + n_dawt] if 0 < pos < n]

    # ── Heatmap ───────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 8))
    sns.heatmap(cor_sort, ax=ax, cmap='RdBu_r', vmin=0.85, vmax=1.0,
                xticklabels=False, yticklabels=False,
                cbar_kws={'label': 'Pearson r'})

    # Remove all spines — they render as gray border lines around the heatmap
    for spine in ax.spines.values():
        spine.set_visible(False)

    # White separator lines between Normal | DAWT | FHWT blocks
    for pos in separators:
        ax.axvline(pos, color='white', lw=2, zorder=5)
        ax.axhline(pos, color='white', lw=2, zorder=5)

    # ── Block labels on the diagonal of each group ────────────────────────────
    # Each label sits at the center of its correlation sub-block
    n_fhwt = n - n_normal - n_dawt
    block_defs = [
        ('Normal',  GDC_SAMPLE_TYPE_COLORS['Solid Tissue Normal'], 0,        n_normal,              n_normal),
        ('DAWT',    PRIMARY_HISTOLOGY_COLORS['DAWT'],              n_normal,  n_normal + n_dawt,     n_dawt),
        ('FHWT',   PRIMARY_HISTOLOGY_COLORS['FHWT'],              n_normal + n_dawt, n,              n_fhwt),
    ]
    for lbl, color, start, end, count in block_defs:
        if count < 1:
            continue
        mid = (start + end) / 2
        ax.text(mid, mid, f'{lbl}\n(n={count})',
                ha='center', va='center', fontsize=9, fontweight='bold',
                color='white', alpha=0.85, zorder=6)

    # ── Title ─────────────────────────────────────────────────────────────────
    ax.set_title('Sample Pearson Correlation — all genes\n'
                 '(sorted: Normal | DAWT | FHWT)')

    fig.tight_layout()
    savefig(fig, out_path)
    return flagged


def plot_volcano(
    res: pd.DataFrame,
    out_path: str,
    title: str = 'Volcano Plot',
    palette: dict | None = None,
    lfc_thresh: float = 1.0,
    padj_thresh: float = 0.05,
    top_label: int = 20,
    gene_col: str = 'gene_name',
) -> None:
    """
    Volcano plot from a DE results DataFrame.

    Parameters
    ----------
    res        : DataFrame with columns [gene_col, log2FoldChange, padj, direction]
    out_path   : save path without extension (PDF + PNG saved automatically)
    palette    : dict mapping direction → colour; defaults to UP=red, DOWN=blue, NS=grey
    lfc_thresh : |LFC| cutoff line (dashed)
    padj_thresh: adjusted p-value cutoff line (dashed)
    top_label  : number of top genes to annotate (split top/bottom LFC)
    gene_col   : column name for gene labels
    """
    palette = palette or {'UP': '#E63946', 'DOWN': '#0077A8', 'NS': '#cccccc'}
    res = res.copy()
    res['_neglog10p'] = -np.log10(res['padj'].clip(lower=1e-300))

    sig = res[res['direction'] != 'NS']
    label_genes = set(
        pd.concat([sig.nlargest(top_label, 'log2FoldChange'),
                   sig.nsmallest(top_label // 2, 'log2FoldChange')])[gene_col]
    )

    fig, ax = plt.subplots(figsize=(8, 7))
    for direction, grp in res.groupby('direction'):
        is_ns = direction == 'NS'
        ax.scatter(
            grp['log2FoldChange'], grp['_neglog10p'],
            c=palette.get(direction, '#999'),
            s=8 if is_ns else 12,
            alpha=0.35 if is_ns else 0.75,
            label=f'{direction} (n={len(grp)})' if not is_ns else 'NS',
            rasterized=True,
        )

    for _, row in res[res[gene_col].isin(label_genes)].iterrows():
        ax.annotate(row[gene_col],
                    xy=(row['log2FoldChange'], row['_neglog10p']),
                    fontsize=6, alpha=0.85,
                    xytext=(3, 3), textcoords='offset points')

    ax.axvline( lfc_thresh,  color='grey', lw=0.8, ls='--')
    ax.axvline(-lfc_thresh,  color='grey', lw=0.8, ls='--')
    ax.axhline(-np.log10(padj_thresh), color='grey', lw=0.8, ls='--')
    ax.set_xlabel('log2 Fold Change (shrunken)')
    ax.set_ylabel('-log10(adjusted p-value)')
    ax.set_title(title)
    ax.legend(fontsize=8, markerscale=2)
    fig.tight_layout()
    savefig(fig, out_path)


def plot_ma(
    res: pd.DataFrame,
    out_path: str,
    title: str = 'MA Plot',
    palette: dict | None = None,
    lfc_thresh: float = 1.0,
) -> None:
    """
    MA plot (mean expression vs log2FC) from a DE results DataFrame.

    Parameters
    ----------
    res        : DataFrame with columns [log2FoldChange, baseMean, direction]
    out_path   : save path without extension
    palette    : dict mapping direction → colour
    lfc_thresh : |LFC| cutoff lines (dashed)
    """
    palette = palette or {'UP': '#E63946', 'DOWN': '#0077A8', 'NS': '#cccccc'}

    fig, ax = plt.subplots(figsize=(7, 5))
    for direction, grp in res.groupby('direction'):
        is_ns = direction == 'NS'
        ax.scatter(
            np.log10(grp['baseMean'] + 1), grp['log2FoldChange'],
            c=palette.get(direction, '#999'),
            s=6,
            alpha=0.25 if is_ns else 0.75,
            label=f'{direction} (n={len(grp)})' if not is_ns else 'NS',
            rasterized=True,
        )

    ax.axhline( lfc_thresh,  color='grey', lw=0.8, ls='--')
    ax.axhline(-lfc_thresh,  color='grey', lw=0.8, ls='--')
    ax.axhline(0,             color='black', lw=0.6)
    ax.set_xlabel('log10(mean expression + 1)')
    ax.set_ylabel('log2 Fold Change')
    ax.set_title(title)
    ax.legend(fontsize=8, markerscale=2)
    fig.tight_layout()
    savefig(fig, out_path)


def plot_heatmap(
    res: pd.DataFrame,
    log2_matrix: pd.DataFrame,
    annot: pd.DataFrame,
    out_path: str,
    title: str = 'Top DEGs',
    group_col: str = 'histology',
    palette: dict | None = None,
    top_n: int = 50,
    gene_col: str = 'gene_name',
) -> None:
    """
    z-score heatmap of top N DEGs (by |LFC|), samples sorted by group_col.

    Parameters
    ----------
    res         : DE results DataFrame with [gene_col, log2FoldChange, direction]
    log2_matrix : log2-normalized expression matrix (genes × samples)
    annot       : sample annotation DataFrame with 'sample_id' and group_col columns
    out_path    : save path without extension
    group_col   : column in annot used to colour the sample bar and sort columns
    palette     : dict mapping group_col values → colour
    top_n       : total number of genes shown (split top/bottom by LFC)
    gene_col    : column name for gene identifiers
    """
    palette = palette or {}
    sig = res[res['direction'] != 'NS']
    if sig.empty:
        print(f'  [skip] No significant DEGs for heatmap: {title}')
        return

    top_genes = pd.concat([
        sig.nlargest(top_n, 'log2FoldChange')[gene_col],
        sig.nsmallest(top_n, 'log2FoldChange')[gene_col],
    ]).drop_duplicates().tolist()
    top_genes = [g for g in top_genes if g in log2_matrix.index][:top_n]
    if not top_genes:
        print(f'  [skip] Genes not found in log2 matrix: {title}')
        return

    mat   = log2_matrix.loc[top_genes]
    mat_z = (mat
             .subtract(mat.mean(axis=1), axis=0)
             .divide(mat.std(axis=1).replace(0, 1), axis=0))

    ann_idx = annot.set_index('sample_id')
    order   = (ann_idx.reindex(mat_z.columns)
                      .sort_values(group_col)
                      .index.tolist())
    mat_z = mat_z[[s for s in order if s in mat_z.columns]]

    sample_colors = [palette.get(ann_idx.loc[s, group_col], '#999')
                     for s in mat_z.columns]

    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(mat_z, ax=ax, cmap='RdBu_r', vmin=-3, vmax=3,
                xticklabels=False, yticklabels=True,
                cbar_kws={'label': 'z-score'})
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)

    # Sample color bar above heatmap
    for i, color in enumerate(sample_colors):
        ax.add_patch(plt.Rectangle(
            (i, -0.5), 1, 0.5, color=color, clip_on=False,
            transform=ax.transData))

    if palette:
        handles = [mpatches.Patch(color=c, label=g) for g, c in palette.items()]
        ax.legend(handles=handles, loc='upper right', fontsize=8,
                  bbox_to_anchor=(1.15, 1.1))

    ax.set_title(title)
    fig.tight_layout()
    savefig(fig, out_path)


def plot_gene_bar(
    res: pd.DataFrame,
    gene_list: list[str],
    out_path: str,
    title: str = 'Gene Panel',
    palette: dict | None = None,
    lfc_thresh: float = 1.0,
    padj_thresh: float = 0.05,
    gene_col: str = 'gene_name',
    x_label: str = 'log2 Fold Change',
) -> None:
    """
    Horizontal bar chart of log2FC for a custom gene list (e.g. Wilms panel).
    Bars coloured by direction; significant genes marked with *.

    Parameters
    ----------
    res         : DE results DataFrame with [gene_col, log2FoldChange, padj, direction]
    gene_list   : ordered list of gene names to display
    out_path    : save path without extension
    title       : figure title
    palette     : dict mapping direction → colour
    lfc_thresh  : |LFC| cutoff lines (dashed)
    padj_thresh : significance threshold for * annotation
    gene_col    : column name for gene identifiers in res
    x_label     : x-axis label (default 'log2 Fold Change')
    """
    palette = palette or {'UP': '#E63946', 'DOWN': '#0077A8', 'NS': '#aaaaaa'}
    wt = res[res[gene_col].isin(gene_list)].copy()
    if wt.empty:
        print(f'  [skip] No matching genes for bar chart: {title}')
        return

    wt = wt.sort_values('log2FoldChange')
    # Color by LFC sign, alpha by significance — all genes shown regardless of direction label
    up_color   = palette.get('UP',   '#E63946')
    down_color = palette.get('DOWN', '#0077A8')
    ns_color   = palette.get('NS',   '#aaaaaa')
    bar_colors = []
    bar_alphas = []
    for lfc, d, p in zip(wt['log2FoldChange'], wt['direction'], wt['padj']):
        is_sig = (p < padj_thresh) and (abs(lfc) > lfc_thresh)
        bar_colors.append(up_color if lfc >= 0 else down_color)
        bar_alphas.append(0.90 if is_sig else 0.35)
    sig_marker = ['*' if p < padj_thresh and abs(lfc) > lfc_thresh else ''
                  for lfc, p in zip(wt['log2FoldChange'], wt['padj'])]

    fig, ax = plt.subplots(figsize=(6, max(4, len(wt) * 0.45)))
    for i, (val, color, alpha) in enumerate(zip(wt['log2FoldChange'], bar_colors, bar_alphas)):
        ax.barh(i, val, color=color, alpha=alpha, edgecolor='white')
    ax.set_yticks(range(len(wt)))
    ax.set_yticklabels(
        [f'{g} {m}' for g, m in zip(wt[gene_col], sig_marker)], fontsize=9)
    ax.axvline(0,            color='black', lw=0.8)
    ax.axvline( lfc_thresh,  color='grey',  lw=0.7, ls='--')
    ax.axvline(-lfc_thresh,  color='grey',  lw=0.7, ls='--')
    ax.set_xlabel(x_label)
    ax.set_title(title + '\n(* padj < 0.05)')
    if palette:
        handles = [mpatches.Patch(color=c, label=g) for g, c in palette.items()
                   if g != 'NS']
        ax.legend(handles=handles, fontsize=8)
    fig.tight_layout()
    savefig(fig, out_path)
