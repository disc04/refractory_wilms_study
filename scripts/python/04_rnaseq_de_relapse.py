#!/usr/bin/env python3
"""
04_rnaseq_de_relapse.py — Wilms Tumor Analysis Project
=======================================================
Differential expression: Relapse vs Progression (primary tumors).

Design
------
- Primary tumors only; QC-flagged samples excluded
- Primary comparison: FHWT Relapse vs FHWT Progression — FEMALE-ONLY (n=36 vs n=4)
  All 4 FHWT progressors are female; 33/69 relapsers are male → perfect sex confound.
  Female-only removes the Y-chromosome gene artifact. Sex cannot be modelled as covariate.
  PyDESeq2: design = ~first_event (Progression as baseline, Relapse as test level)
- Supplementary: FHWT all-sex (for confound sensitivity check; not primary)
- Exploratory: DAWT Relapse vs DAWT Progression (n=15 vs n=3) — underpowered

LFC threshold: |log2FC| > 1.0, padj < 0.05

Inputs
------
data/processed/matrices/counts_filtered.csv
data/processed/tables/sample_annotation.csv

Outputs
-------
data/processed/tables/de_relapse_fhwt_full.csv
data/processed/tables/de_relapse_fhwt_sig.csv
data/processed/tables/de_relapse_dawt_full.csv (exploratory)
data/processed/tables/de_relapse_dawt_sig.csv (exploratory)
results/figures/de_rel_fhwt_01_volcano.pdf/.png
results/figures/de_rel_fhwt_02_ma.pdf/.png
results/figures/de_rel_fhwt_03_heatmap.pdf/.png
results/figures/de_rel_fhwt_04_wilms_genes.pdf/.png
results/figures/de_rel_dawt_01_volcano.pdf/.png  (exploratory)
"""

from __future__ import annotations
import os
import warnings
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from utils.de_utils import run_deseq2
from configuration import PATHS, WILMS_PANEL, RELAPSE_PALETTE
from utils.de_plot_utils import plot_volcano, plot_ma, plot_heatmap, plot_gene_bar

warnings.filterwarnings('ignore')

processed_tables_path, matrices_path, figures_path = (PATHS["processed_tables_path"],
                                                      PATHS["matrices_path"],
                                                      PATHS["figures_path"])

LFC_THRESH  = 1.0
PADJ_THRESH = 0.05
TOP_HEATMAP = 50

# Contrast: Relapse vs Progression (Progression=reference, PyDESeq2 alphabetical baseline)
# UP = higher in Relapse, DOWN = higher in Progression


def load_data(
    histology: str,
    exclude_qc_flagged: bool = True,
    female_only: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load counts and annotation filtered to one histology's primary tumours.

    Parameters
    ----------
    histology          : 'FHWT' or 'DAWT'
    exclude_qc_flagged : drop QC-flagged samples (default True)
    female_only        : restrict to female samples; recommended for FHWT because
                         all 4 progressors are female while 33/69 relapsers are male,
                         creating a perfect sex confound without this filter

    Returns
    -------
    counts_sub : filtered raw count matrix (genes × selected samples)
    pt         : sample annotation DataFrame for the selected samples
    """
    counts = pd.read_csv(os.path.join(matrices_path, 'counts_filtered_primary_samples.csv'),
                         index_col=0)
    annot  = pd.read_csv(os.path.join(processed_tables_path, 'cohort_primary_samples_metadata.csv'))

    pt = annot[
        (annot['sample_type'] == 'Primary Tumor') &
        (annot['histology'] == histology)
    ].copy()
    if exclude_qc_flagged:
        n_before = len(pt)
        pt = pt[~pt['qc_flag']]
        if n_before > len(pt):
            print(f"  QC-flagged excluded: {n_before - len(pt)} samples")

    # Keep only samples with valid First Event for this comparison
    pt = pt[pt['First Event'].isin(['Relapse', 'Progression'])].copy()

    if female_only:
        n_before = len(pt)
        pt = pt[pt['Gender'] == 'Female'].copy()
        print(f"  Female-only filter: {n_before} → {len(pt)} samples "
              f"(removed {n_before - len(pt)} males)")

    # Rename for formulaic compatibility (no spaces in column name)
    pt = pt.rename(columns={'First Event': 'first_event'})

    counts_sub = counts[pt['sample_id'].tolist()]
    fe_counts  = pt['first_event'].value_counts().to_dict()
    print(f"  {histology} — {fe_counts} (total n={len(pt)})")
    return counts_sub, pt


def print_summary(res: pd.DataFrame, label: str) -> None:
    """
    Print a concise DE result summary to stdout.

    Parameters
    ----------
    res   : full DE results DataFrame with 'direction' column
    label : human-readable label printed as a header line
    """
    sig = res[res['direction'] != 'NS']
    n_up   = (res['direction'] == 'UP').sum()
    n_down = (res['direction'] == 'DOWN').sum()
    print(f"\n  {label}")
    print(f"  Genes tested : {len(res):,}")
    print(f"  Significant  : {len(sig):,}  (UP={n_up}, DOWN={n_down})")
    if not sig.empty:
        print(f"  Top 5 UP:")
        print(res[res['direction'] == 'UP'].head(5)
              [['gene_name', 'log2FoldChange', 'padj']].to_string(index=False))
        print(f"  Top 5 DOWN:")
        print(res[res['direction'] == 'DOWN'].head(5)
              [['gene_name', 'log2FoldChange', 'padj']].to_string(index=False))


def main() -> None:  # noqa: D103
    print("\n" + "="*60)
    print("  04_rnaseq_de_relapse.py — Relapse vs Progression")
    print("="*60)

    log2 = pd.read_csv(os.path.join(matrices_path, 'counts_log2_primary_samples.csv'), index_col=0)

    # PRIMARY: FHWT Relapse vs Progression (female-only)
    # NOTE: all 4 FHWT progressors are female; 33/69 relapsers are male.
    # This creates a perfect sex confound (no male progressors → can't model sex
    # as covariate). Primary analysis restricted to females to remove the
    # Y-chromosome gene artifact. Full cohort run retained as supplementary.
    print("\n── FHWT: Relapse vs Progression (female-only, primary) ─────")
    counts_f, meta_f = load_data('FHWT', female_only=True)

    res_f = run_deseq2(counts_f, meta_f,
                       design='~first_event',
                       ref=('first_event', 'Progression'),
                       lfc_thresh=LFC_THRESH,
                       padj_thresh=PADJ_THRESH)
    print_summary(res_f, 'FHWT Relapse vs Progression [female-only] (UP=Relapse)')

    sig_f = res_f[res_f['direction'] != 'NS']
    res_f.to_csv(os.path.join(processed_tables_path, 'de_relapse_fhwt_full.csv'), index=False)
    sig_f.to_csv(os.path.join(processed_tables_path, 'de_relapse_fhwt_sig.csv'), index=False)
    print(f"\n  de_relapse_fhwt_full.csv : {len(res_f):,} genes")
    print(f"  de_relapse_fhwt_sig.csv  : {len(sig_f):,} genes")

    log2_f = log2[counts_f.columns]
    title_base = f'FHWT: Relapse vs Progression — {len(sig_f)} DEGs'

    plot_volcano(
        res_f,
        out_path=os.path.join(figures_path, 'de_rel_fhwt_01_volcano'),
        title=f'{title_base}\n(|LFC|>{LFC_THRESH}, padj<{PADJ_THRESH}; UP=Relapse)',
        palette=RELAPSE_PALETTE,
        lfc_thresh=LFC_THRESH,
        padj_thresh=PADJ_THRESH,
    )
    print("  Saved: de_rel_fhwt_01_volcano")

    plot_ma(
        res_f,
        out_path=os.path.join(figures_path, 'de_rel_fhwt_02_ma'),
        title='MA Plot — FHWT Relapse vs Progression',
        palette=RELAPSE_PALETTE,
        lfc_thresh=LFC_THRESH,
    )
    print("  Saved: de_rel_fhwt_02_ma")

    plot_heatmap(
        res_f, log2_f, meta_f,
        out_path=os.path.join(figures_path, 'de_rel_fhwt_03_heatmap'),
        title=f'Top {TOP_HEATMAP} DEGs — z-score (FHWT Relapse vs Progression)',
        group_col='first_event',
        palette={'Relapse': '#1f77b4', 'Progression': '#9467bd'},
        top_n=TOP_HEATMAP,
    )
    print("  Saved: de_rel_fhwt_03_heatmap")

    plot_gene_bar(
        res_f, WILMS_PANEL,
        out_path=os.path.join(figures_path, 'de_rel_fhwt_04_wilms_genes'),
        title='Wilms Panel — FHWT Progression vs Relapse',
        palette=RELAPSE_PALETTE,
        lfc_thresh=LFC_THRESH,
        padj_thresh=PADJ_THRESH,
        x_label='log2 Fold Change (Relapse vs Progression)',
    )
    print("  Saved: de_rel_fhwt_04_wilms_genes")

    # SUPPLEMENTARY: FHWT all-sex (for sex-confound sensitivity check)
    print("\n── FHWT: Relapse vs Progression (all-sex, supplementary) ───")
    print("  [!] Y-chromosome genes confound this; use female-only as primary")
    counts_fs, meta_fs = load_data('FHWT', female_only=False)
    res_fs = run_deseq2(counts_fs, meta_fs,
                        design='~first_event',
                        ref=('first_event', 'Progression'),
                        lfc_thresh=LFC_THRESH,
                        padj_thresh=PADJ_THRESH)
    sig_fs = res_fs[res_fs['direction'] != 'NS']
    res_fs.to_csv(os.path.join(processed_tables_path, 'de_relapse_fhwt_allsex_full.csv'), index=False)
    sig_fs.to_csv(os.path.join(processed_tables_path, 'de_relapse_fhwt_allsex_sig.csv'), index=False)
    print(f"  de_relapse_fhwt_allsex_full.csv : {len(res_fs):,} genes, {len(sig_fs)} sig")

    # EXPLORATORY: DAWT Relapse vs Progression
    print("\n── DAWT: Relapse vs Progression (exploratory, n_prog=3) ─────")
    print("  [!] UNDERPOWERED — treat as hypothesis-generating only")
    counts_d, meta_d = load_data('DAWT')

    res_d = run_deseq2(counts_d, meta_d,
                       design='~first_event',
                       ref=('first_event', 'Progression'),
                       lfc_thresh=LFC_THRESH,
                       padj_thresh=PADJ_THRESH)
    print_summary(res_d, 'DAWT Relapse vs Progression')

    sig_d = res_d[res_d['direction'] != 'NS']
    res_d.to_csv(os.path.join(processed_tables_path, 'de_relapse_dawt_full.csv'), index=False)
    sig_d.to_csv(os.path.join(processed_tables_path, 'de_relapse_dawt_sig.csv'), index=False)

    plot_volcano(
        res_d,
        out_path=os.path.join(figures_path, 'de_rel_dawt_01_volcano'),
        title=f'DAWT Relapse vs Progression [EXPLORATORY, n_prog=3]\n'
              f'{len(sig_d)} DEGs (|LFC|>{LFC_THRESH}, padj<{PADJ_THRESH})',
        palette=RELAPSE_PALETTE,
        lfc_thresh=LFC_THRESH,
        padj_thresh=PADJ_THRESH,
    )
    print("  Saved: de_rel_dawt_01_volcano")
    print("\n── Complete ─────────────────────────────────────────────────")



if __name__ == '__main__':
    main()
