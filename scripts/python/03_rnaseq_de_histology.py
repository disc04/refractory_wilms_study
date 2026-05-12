#!/usr/bin/env python3
"""
03_rnaseq_de_histology.py — Wilms Tumor Analysis Project
=========================================================
Differential expression: FHWT vs DAWT (primary tumors).
Standard two-group comparison, the baseline analysis for the cohort.

Design
------
- Primary tumors only (dropped recurrent and metastatic); QC-flagged samples excluded
- PyDESeq2: design = ~ histology (FHWT vs DAWT, DAWT as reference)
- LFC threshold: |log2FC| > 1.0, padj < 0.05

Inputs
------
data/processed/matrices/counts_filtered.csv
data/processed/tables/sample_annotation.csv

Outputs
-------
data/processed/tables/de_histology_full.csv
data/processed/tables/de_histology_sig.csv
results/figures/de_hist_01_volcano.pdf/.png
results/figures/de_hist_02_ma.pdf/.png
results/figures/de_hist_03_heatmap.pdf/.png
results/figures/de_hist_04_wilms_genes.pdf/.png
"""

from __future__ import annotations
import os
import warnings
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from utils.de_utils import run_deseq2
from configuration import PATHS, PRIMARY_HISTOLOGY_COLORS, WILMS_PANEL, DE_PALETTE
from utils.de_plot_utils import plot_volcano, plot_ma, plot_heatmap, plot_gene_bar

warnings.filterwarnings('ignore')
processed_tables_path, matrices_path, figures_path = (PATHS["processed_tables_path"],
                                                      PATHS["matrices_path"],
                                                      PATHS["figures_path"])

LFC_THRESH  = 1.0
PADJ_THRESH = 0.05
TOP_HEATMAP = 50

# FHWT is reference: UP = higher in DAWT (anaplastic program), DOWN = higher in FHWT
# This convention aligns Step 1 UP genes with Step 3 AT UP genes for direct overlap

# ── Data loading ──────────────────────────────────────────────────────────────

def load_data(exclude_qc_flagged: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load filtered counts and annotation for the DAWT vs FHWT comparison.

    Parameters
    ----------
    exclude_qc_flagged : drop samples flagged during QC (default True)

    Returns
    -------
    counts_pt : filtered raw count matrix (genes × primary tumour samples)
    pt        : sample annotation DataFrame for the selected primary tumour samples
    """
    counts = pd.read_csv(os.path.join(matrices_path, 'counts_filtered_primary_samples.csv'),
                         index_col=0)
    annot  = pd.read_csv(os.path.join(processed_tables_path, 'cohort_primary_samples_metadata.csv'))

    pt = annot[annot['sample_type'] == 'Primary Tumor'].copy()
    if exclude_qc_flagged:
        n_before = len(pt)
        pt = pt[~pt['qc_flag']]
        print(f"  QC-flagged excluded: {n_before - len(pt)} samples")

    counts_pt = counts[pt['sample_id'].tolist()]
    print(f"  Histology: {pt['histology'].value_counts().to_dict()}")
    return counts_pt, pt


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:  # noqa: D103
    print("\n" + "="*60)
    print("  03_rnaseq_de_histology.py — FHWT vs DAWT")
    print("="*60)

    # Load
    print("\n── Step 1: Load data ────────────────────────────────────────")
    counts, meta = load_data(exclude_qc_flagged=True)

    # DE
    print("\n── Step 2: PyDESeq2 (DAWT vs FHWT, FHWT reference) ─────────")
    res = run_deseq2(counts, meta,
                     design='~histology',
                     ref=('histology', 'FHWT'),
                     lfc_thresh=LFC_THRESH,
                     padj_thresh=PADJ_THRESH)

    sig = res[res['direction'] != 'NS']
    print(f"\n  Total genes tested : {len(res):,}")
    print(f"  Significant DEGs   : {len(sig):,}  "
          f"(UP={len(sig[sig['direction']=='UP'])}, "
          f"DOWN={len(sig[sig['direction']=='DOWN'])})")
    print(f"\n  Top 10 UP (DAWT > FHWT — anaplastic program):")
    print(res[res['direction'] == 'UP'].head(10)
          [['gene_name', 'log2FoldChange', 'padj']].to_string(index=False))
    print(f"\n  Top 10 DOWN (FHWT > DAWT):")
    print(res[res['direction'] == 'DOWN'].head(10)
          [['gene_name', 'log2FoldChange', 'padj']].to_string(index=False))

    # Save tables
    print("\n── Step 3: Save tables ──────────────────────────────────────")
    res.to_csv(os.path.join(processed_tables_path, 'de_histology_full.csv'), index=False)
    sig.to_csv(os.path.join(processed_tables_path, 'de_histology_sig.csv'), index=False)
    print(f"  de_histology_full.csv : {len(res):,} genes")
    print(f"  de_histology_sig.csv  : {len(sig):,} genes")

    # Figures
    print("\n── Step 4: Figures ──────────────────────────────────────────")
    log2 = pd.read_csv(os.path.join(matrices_path, 'counts_log2_primary_samples.csv'), index_col=0)
    log2 = log2[counts.columns]

    plot_volcano(
        res,
        out_path=os.path.join(figures_path, 'de_hist_01_volcano'),
        title=f'DAWT vs FHWT — {len(sig)} DEGs (|LFC|>{LFC_THRESH}, padj<{PADJ_THRESH})\n'
              f'(FHWT reference; UP = anaplastic program)',
        palette=DE_PALETTE,
        lfc_thresh=LFC_THRESH,
        padj_thresh=PADJ_THRESH,
    )
    print("  Saved: de_hist_01_volcano")

    plot_ma(
        res,
        out_path=os.path.join(figures_path, 'de_hist_02_ma'),
        title='MA Plot — DAWT vs FHWT (FHWT reference)',
        palette=DE_PALETTE,
        lfc_thresh=LFC_THRESH,
    )
    print("  Saved: de_hist_02_ma")

    plot_heatmap(
        res, log2, meta,
        out_path=os.path.join(figures_path, 'de_hist_03_heatmap'),
        title=f'Top {TOP_HEATMAP} DEGs — z-score (DAWT vs FHWT)',
        group_col='histology',
        palette=PRIMARY_HISTOLOGY_COLORS,
        top_n=TOP_HEATMAP,
    )
    print("  Saved: de_hist_03_heatmap")

    plot_gene_bar(
        res, WILMS_PANEL,
        out_path=os.path.join(figures_path, 'de_hist_04_wilms_genes'),
        title='Wilms Panel Genes — DAWT vs FHWT',
        palette=DE_PALETTE,
        lfc_thresh=LFC_THRESH,
        padj_thresh=PADJ_THRESH,
        x_label='log2 Fold Change (DAWT vs FHWT)',
    )
    print("  Saved: de_hist_04_wilms_genes")
    print("\n── Complete ─────────────────────────────────────────────────")


if __name__ == '__main__':
    main()
