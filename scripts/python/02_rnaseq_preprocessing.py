#!/usr/bin/env python3
"""
02_rnaseq_preprocessing.py — Wilms Tumor Analysis Project
==========================================================
QC, filtering, and normalization of TARGET-WT RNA-seq counts.
Produces a clean counts matrix and log2-normalised matrix
for all downstream differential expression analyses.

Inputs
------
data/processed/tables/gdc_counts_matrix.csv — long-format counts (gene × sample)
data/processed/tables/cohort_metadata.csv — sample annotation

Outputs
-------
data/processed/matrices/counts_filtered.csv   — filtered raw counts (genes × samples)
data/processed/matrices/counts_normalized.csv — DESeq2 size-factor normalized (125 samples)
data/processed/matrices/counts_log2.csv       — log2(norm + 1), for DE / heatmap (125 samples)
data/processed/matrices/counts_log2_viz.csv   — log2(norm + 1), all 136 samples for PCA viz
data/processed/tables/cohort_primary_samples_metadata.csv — final sample annotation (125 samples)
results/figures/qc_01_library_sizes.pdf/.png
results/figures/qc_02_count_density.pdf/.png
results/figures/qc_03_pca.pdf/.png
results/figures/qc_04_sample_correlation.pdf/.png
"""

from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from pydeseq2.preprocessing import deseq2_norm

from configuration import PATHS, EXCLUDE_USI
from utils.data_utils import load_gdc_counts_matrix
from utils.de_plot_utils import (plot_library_sizes, plot_count_density,
                                 plot_pca, plot_sample_correlation)

processed_tables_path, matrices_path, figures_path = (PATHS["processed_tables_path"],
PATHS["matrices_path"], PATHS["figures_path"])

QC_FLAG_SAMPLES = [
    'TARGET-50-PAKGED-01A', 'TARGET-50-CAAAAM-01A',
    'TARGET-50-PAJMVC-01A', 'TARGET-50-PAJNAV-01A',
]

MIN_COUNTS  = 10
MIN_SAMPLES = 10


# ── Gene filtering ─────────────────────────────────────────────────────────────

def filter_genes(matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Keep genes with ≥ MIN_COUNTS reads in ≥ MIN_SAMPLES samples.

    Parameters
    ----------
    matrix : raw count matrix (genes × samples)

    Returns
    -------
    Filtered matrix retaining only genes that pass the low-count threshold.
    """
    keep = (matrix >= MIN_COUNTS).sum(axis=1) >= MIN_SAMPLES
    filtered = matrix.loc[keep]
    print(f"  After low-count filter (≥{MIN_COUNTS} in ≥{MIN_SAMPLES} samples): "
          f"{filtered.shape[0]:,} genes")
    return filtered


# ── Normalization ──────────────────────────────────────────────────────────────

def normalize(matrix: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    DESeq2 median-of-ratios normalization.

    Parameters
    ----------
    matrix : raw count matrix (genes × samples) with integer-compatible values

    Returns
    -------
    norm : size-factor normalized count matrix (genes × samples)
    log2 : log2(norm + 1) matrix (genes × samples)
    """
    norm_T, sf = deseq2_norm(matrix.T.astype(int))
    norm = norm_T.T
    log2 = np.log2(norm + 1)
    print(f"  Size factors: {sf.min():.3f} – {sf.max():.3f}")
    return norm, log2


def normalize_extra(counts_extra: pd.DataFrame,
                    counts_ref: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Normalize additional samples (Recurrent / Metastatic) by projecting them
    onto the normalization scale of the primary cohort.

    Size factors are computed using the gene-wise log geometric means from
    the primary cohort raw counts as the reference, following the DESeq2
    median-of-ratios method. This ensures the extra samples are comparable
    to the primary analysis matrices without their expression profiles
    influencing the primary size factors.

    Parameters
    ----------
    counts_extra : raw counts for extra samples (genes × extra_samples)
    counts_ref   : raw counts for the primary cohort (genes × primary_samples)

    Returns
    -------
    norm  : size-factor normalized counts
    log2  : log2(norm + 1)
    """
    # Gene-wise log geometric mean from the primary cohort raw counts
    log_ref = np.log(counts_ref.where(counts_ref > 0, np.nan))
    log_geo_means = log_ref.mean(axis=1)                    # NaN for all-zero genes
    valid_genes = log_geo_means[np.isfinite(log_geo_means)].index

    # Per-sample size factor: exp(median log-ratio) across valid genes
    size_factors = {}
    for sid in counts_extra.columns:
        col = counts_extra.loc[valid_genes, sid]
        log_col = np.log(col.where(col > 0, np.nan))
        log_ratios = log_col - log_geo_means.loc[valid_genes]
        size_factors[sid] = float(np.exp(np.nanmedian(log_ratios)))

    sf = pd.Series(size_factors)
    print(f"  Extra-sample size factors ({len(sf)} samples): "
          f"{sf.min():.3f} – {sf.max():.3f}")

    norm = counts_extra.div(sf, axis=1)
    log2 = np.log2(norm + 1)
    return norm, log2


# ── Primary cohort builder ─────────────────────────────────────────────────────

def build_primary_samples_matrices(
    counts: pd.DataFrame,
    meta_all_samples: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Subset counts to the primary analysis cohort (Primary Tumor + Solid Tissue
    Normal, minus EXCLUDE_USI), apply gene filter and DESeq2 normalization.

    Parameters
    ----------
    counts           : full protein-coding count matrix (genes × all samples)
    meta_all_samples : annotation DataFrame for all samples (includes 'usi',
                       'sample_id', 'sample_type' columns)

    Returns
    -------
    annot          : sample annotation for the 125-sample primary cohort
    matrix_primary : filtered raw counts (genes × 125 samples)
    norm_primary   : size-factor normalized counts (genes × 125 samples)
    log2_primary   : log2(norm + 1) (genes × 125 samples)
    """
    cohort_primary = meta_all_samples[
        ~meta_all_samples['usi'].isin(EXCLUDE_USI) &
        meta_all_samples['sample_type'].isin(['Primary Tumor', 'Solid Tissue Normal'])
    ].copy()
    print(f"\n  Samples after exclusion:")
    print(cohort_primary['sample_type'].value_counts().to_string())

    # Align annotation to available matrix columns
    shared = cohort_primary['sample_id'].isin(counts.columns)
    annot = cohort_primary[shared].copy()
    matrix_primary = counts[annot['sample_id'].values]

    # Gene filter and normalization on primary cohort only
    matrix_primary = filter_genes(matrix_primary)
    norm_primary, log2_primary = normalize(matrix_primary)

    # QC flag
    annot['qc_flag'] = annot['sample_id'].isin(QC_FLAG_SAMPLES)
    return annot, matrix_primary, norm_primary, log2_primary


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    print("\n" + "="*60)
    print("  02_rnaseq_preprocessing.py")
    print("="*60)

    # ── Step 1: Load counts (single pass, all samples) ───────────────────────
    print("\n Load & pivot counts ──────────────────────────────")
    meta_all_samples = pd.read_csv(os.path.join(processed_tables_path, 'cohort_metadata.csv'))
    counts = load_gdc_counts_matrix(os.path.join(processed_tables_path, 'gdc_counts_matrix.csv'))

    # ── Step 2: Filter genes and normalize primary cohort ────────────────────
    print("\n Filter genes & normalize primary cohort ───────")
    cohort_primary, matrix_primary, norm_primary, log2_primary = \
        build_primary_samples_matrices(counts, meta_all_samples)

    # ── Step 3: Build viz matrix (all 136 samples) ───────────────────────────
    # Restrict full counts to the gene set established by the primary filter,
    # then normalize extra samples using the primary cohort as reference.
    print("\n Build full matrix: primary + anaplastic transformation + metastatic ─────────")
    counts_filtered_all = counts.reindex(matrix_primary.index)

    extra_ids = [s for s in counts_filtered_all.columns
                 if s not in cohort_primary['sample_id'].values]
    print(f"  Extra samples (Recurrent / Metastatic): {len(extra_ids)}")
    matrix_extra = counts_filtered_all[extra_ids]
    _, log2_extra = normalize_extra(matrix_extra, matrix_primary)

    # Viz matrix: primary (125) + extra → all 136
    log2_viz = pd.concat([log2_primary, log2_extra], axis=1)

    # Full annotation aligned to viz matrix columns (with qc_flag)
    meta_viz = meta_all_samples[
        meta_all_samples['sample_id'].isin(log2_viz.columns)].copy()
    meta_viz['qc_flag'] = meta_viz['sample_id'].isin(QC_FLAG_SAMPLES)

    # ── Step 4: Save matrices ─────────────────────────────────────────────────
    print("\n Saving matrices ────────────────────────────────────")
    matrix_primary.to_csv(os.path.join(matrices_path, 'counts_filtered_primary_samples.csv'))
    norm_primary.to_csv(os.path.join(matrices_path, 'counts_normalized_primary_samples.csv'))
    log2_primary.to_csv(os.path.join(matrices_path, 'counts_log2_primary_samples.csv'))
    log2_viz.to_csv(os.path.join(matrices_path, 'counts_log2_all.csv'))
    cohort_primary.to_csv(
        os.path.join(processed_tables_path, 'cohort_primary_samples_metadata.csv'),
        index=False)

    print(f"  counts_filtered  : {matrix_primary.shape[0]:,} × {matrix_primary.shape[1]}")
    print(f"  counts_log2_viz  : {log2_viz.shape[0]:,} × {log2_viz.shape[1]}")
    print(f"  sample_annotation: {cohort_primary.shape[0]} samples "
          f"({cohort_primary['qc_flag'].sum()} QC-flagged)")

    # ── Step 5: QC figures ────────────────────────────────────────────────────
    print("\n QC figures ──────────────────────────────────────")
    plot_library_sizes(
        counts_filtered_all, meta_viz,
        os.path.join(figures_path, 'qc_01_library_sizes'))
    print("  Saved: qc_01_library_sizes")

    plot_count_density(
        log2_viz, meta_viz,
        os.path.join(figures_path, 'qc_02_count_density'))
    print("  Saved: qc_02_count_density")

    # PCA: fit on Primary + Normal (125), project all 136, plot all
    plot_pca(
        log2_viz, meta_viz,
        os.path.join(figures_path, 'qc_03_pca'),
        fit_exclude_types=['Recurrent Tumor', 'Metastatic'],
    )
    print("  Saved: qc_03_pca")

    plot_sample_correlation(
        log2_primary, cohort_primary,
        os.path.join(figures_path, 'qc_04_sample_correlation'))
    print("  Saved: qc_04_sample_correlation")
    print("\n Complete ─────────────────────────────────────────────────")


if __name__ == '__main__':
    main()
