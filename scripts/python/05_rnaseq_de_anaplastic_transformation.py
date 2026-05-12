#!/usr/bin/env python3
"""
05_rnaseq_de_anaplastic_transformation.py — Wilms Tumor Analysis Project
=========================================================================
Paired longitudinal DE: Primary → Recurrent (anaplastic transformation).

Patients
--------
4 patients with FHWT primary tumors confirmed to have transformed to anaplastic
histology at relapse (pathology_subtype = "anaplastic Wilms"). Each contributes
one primary (01A) and one recurrent (02A) sample.

    TARGET-50-PAJPDC  Female  Stage II  (anaplastic Wilms primary)
    TARGET-50-PAJNGH  Male    Stage I   (anaplastic Wilms primary)
    TARGET-50-PALFME  Male    Stage II  (anaplastic Wilms primary)
    TARGET-50-PAJNTJ  Male    Stage II  (anaplastic Wilms primary)

Excluded: TARGET-50-PALJIP — pathology_subtype = "relapse Wilms" (not anaplastic Wilms);
primary tumor retained in main FHWT cohort.

Design
------
PyDESeq2: ~ patient + time_point (patient=blocking factor, n=4 pairs)
  - time_point: Primary (reference) → Recurrent (test)
  - Positive LFC = upregulated in anaplastic recurrence
  - n=8 samples total; residual df=3 (tight — interpret with caution, LFC shrinkage applied)

Inputs
------
data/processed/tables/cohort_metadata.csv
data/processed/tables/gdc_counts_matrix.csv (raw long-format, all samples)
data/processed/matrices/counts_filtered.csv (gene universe, 17,341 genes)
data/processed/tables/de_relapse_fhwt_sig.csv (Step 2 results for overlap)

Outputs
-------
data/processed/matrices/counts_at_raw.csv — AT raw count matrix (8 samples)
data/processed/tables/de_at_full.csv
data/processed/tables/de_at_sig.csv
data/processed/tables/de_at_vs_relapse_overlap.csv — Step 2 gene overlap
results/figures/de_at_01_volcano.pdf/.png
results/figures/de_at_02_ma.pdf/.png
results/figures/de_at_03_heatmap.pdf/.png
results/figures/de_at_04_wilms_genes.pdf/.png
results/figures/de_at_05_overlap.pdf/.png
"""

from __future__ import annotations
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from utils.data_utils import savefig, load_gdc_counts_matrix
from utils.de_utils import assign_direction
from utils.de_plot_utils import plot_volcano, plot_ma, plot_heatmap, plot_gene_bar
from configuration import (PATHS, AT_PALETTE, TIME_PALETTE, WILMS_PANEL,
                           ANAPLASTIC_TRANSFORMATION_PATIENTS)

warnings.filterwarnings('ignore')

processed_tables_path, matrices_path, figures_path = (PATHS["processed_tables_path"],
                                                      PATHS["matrices_path"],
                                                      PATHS["figures_path"])

AT_USI = ANAPLASTIC_TRANSFORMATION_PATIENTS
LFC_THRESH  = 1.0
PADJ_THRESH = 0.05
TOP_HEATMAP = 40


# ── 1. Build AT data ───────────────────────────────────────────────

def build_at_metadata() -> pd.DataFrame:
    """
    Load cohort metadata and subset to the four anaplastic transformation patients.

    Derives 'time_point' (Primary / Recurrent) from the sample suffix (-01A / -02A)
    and 'patient' from the USI.

    Returns
    -------
    meta_at : annotation DataFrame for AT samples with 'patient' and 'time_point' columns
    """
    meta_full = pd.read_csv(os.path.join(processed_tables_path, 'cohort_metadata.csv'))
    meta_at = meta_full[meta_full['usi'].isin(AT_USI)].copy()

    # time_point from sample suffix (-01A=Primary, -02A=Recurrent)
    meta_at['time_point'] = meta_at['sample_id'].str.extract(r'-(0[12]A)$')[0] \
        .map({'01A': 'Primary', '02A': 'Recurrent'})
    meta_at['patient'] = meta_at['usi'].str.replace('TARGET-50-', '', regex=False)

    print(meta_at[['sample_id', 'patient', 'time_point', 'Gender', 'Stage',
                   'pathology_subtype']].sort_values(['patient', 'time_point']) \
          .to_string(index=False)
          )
    return meta_at

def check_complete_pairs(meta_at: pd.DataFrame) -> None:
    """
    Warn if any patient is missing their paired Primary or Recurrent sample.

    Parameters
    ----------
    meta_at : AT annotation DataFrame with 'patient' and 'time_point' columns
    """
    paired = meta_at.groupby('patient')['time_point'].nunique()
    complete_pairs = paired[paired == 2].index.tolist()
    if len(complete_pairs) < len(AT_USI):
        missing = set(meta_at['patient'].unique()) - set(complete_pairs)
        print(f"  [warn] Incomplete pairs excluded: {missing}")
        meta_at = meta_at[meta_at['patient'].isin(complete_pairs)]

    print(f"\n  Paired samples: {len(complete_pairs)} patients × 2 = {len(meta_at)} samples")

def build_at_counts_matrix(
    meta_at: pd.DataFrame,
) -> tuple[pd.DataFrame, list[str]]:
    """
    Load raw counts for AT samples from the GDC long-format matrix.

    Aligns the result to the primary analysis gene universe (from
    counts_filtered_primary_samples.csv) and saves counts_at_raw.csv.

    Parameters
    ----------
    meta_at : AT annotation DataFrame with a 'sample_id' column

    Returns
    -------
    matrix        : raw count matrix (genes × AT samples)
    valid_samples : sample IDs successfully loaded (subset of meta_at['sample_id'])
    """
    expected_sids = meta_at['sample_id'].tolist()

    # Gene universe from primary analysis — ensures consistent gene-set across all DE scripts.
    # usecols=[0] reads only the index column (gene names) without loading the full 17k × 125 matrix.
    gene_universe = pd.read_csv(
        os.path.join(matrices_path, 'counts_filtered_primary_samples.csv'),
        index_col=0, usecols=[0]).index

    matrix = load_gdc_counts_matrix(os.path.join(processed_tables_path, 'gdc_counts_matrix.csv'),
                                    sample_ids=expected_sids,
                                    gene_universe=gene_universe)
    matrix.to_csv(os.path.join(matrices_path, 'counts_at_raw.csv'))
    print(f"  Saved: counts_at_raw.csv "
          f"({matrix.shape[0]:,} genes × {matrix.shape[1]} samples)")
    # Align column order and drop any samples missing from the loaded matrix
    valid_samples = [s for s in expected_sids if s in matrix.columns]
    matrix = matrix[valid_samples]

    return matrix, valid_samples

# ── 2. PyDESeq2 paired analysis ───────────────────────────────────────────────

def run_paired_deseq2(counts: pd.DataFrame,
                      meta: pd.DataFrame) -> pd.DataFrame:
    """
    Paired design: ~ patient + time_point
    Blocking by patient removes inter-individual variance.
    Reference: Primary; Test: Recurrent (positive LFC = up in anaplastic).
    """
    counts_T = counts.T.astype(int)
    meta_dds = meta.set_index('sample_id')[['patient', 'time_point']]\
                   .reindex(counts_T.index)

    # Relevel: Primary first → reference baseline for time_point
    meta_dds['time_point'] = pd.Categorical(meta_dds['time_point'],
                                            categories=['Primary', 'Recurrent'])

    dds = DeseqDataSet(
        counts=counts_T,
        metadata=meta_dds,
        design='~patient + time_point',
        quiet=True,
    )
    dds.deseq2()

    stat = DeseqStats(
        dds,
        contrast=['time_point', 'Recurrent', 'Primary'],
        alpha=PADJ_THRESH,
        quiet=True,
    )
    stat.summary()
    stat.lfc_shrink(coeff='time_point[T.Recurrent]')

    res = stat.results_df.reset_index().rename(columns={'index': 'gene_name'})
    res = assign_direction(res, lfc_thresh=LFC_THRESH, padj_thresh=PADJ_THRESH)
    res = res.sort_values(['padj', 'log2FoldChange'], ascending=[True, False])
    return res


# ── 3. DEG overlap with Step 2 (FHWT relapse) ────────────────────────────────

def analyse_overlap(res_at: pd.DataFrame,
                    sig_rel: pd.DataFrame) -> pd.DataFrame:
    """
    Compare AT DEGs with Step 2 FHWT relapse DEGs.

    Parameters
    ----------
    res_at  : full AT DE results (all genes, with 'direction' column)
    sig_rel : significant FHWT relapse DEGs (de_relapse_fhwt_sig.csv)

    Returns
    -------
    DataFrame with direction in each comparison, or empty DataFrame if
    either input is empty.
    """
    if sig_rel.empty:
        print("  [skip overlap] relapse sig table is empty")
        return pd.DataFrame()

    sig_at = res_at[res_at['direction'] != 'NS'][
        ['gene_name', 'log2FoldChange', 'padj', 'direction']].copy()
    sig_at.columns = ['gene_name', 'lfc_at', 'padj_at', 'dir_at']

    sig_rel = sig_rel[['gene_name', 'log2FoldChange', 'padj', 'direction']].copy()
    sig_rel.columns = ['gene_name', 'lfc_rel', 'padj_rel', 'dir_rel']

    overlap = sig_at.merge(sig_rel, on='gene_name')
    concordant = overlap[overlap['dir_at'] == overlap['dir_rel']]
    discordant = overlap[overlap['dir_at'] != overlap['dir_rel']]

    print(f"\n  DEG overlap: AT sig={len(sig_at)}, Relapse sig={len(sig_rel)}, "
          f"overlap={len(overlap)}")
    print(f"  Concordant (same direction): {len(concordant)}")
    print(f"  Discordant (opposite direction): {len(discordant)}")
    if not concordant.empty:
        print(f"  Concordant genes:")
        print(concordant[['gene_name', 'lfc_at', 'lfc_rel', 'dir_at']].to_string(index=False))

    overlap['concordance'] = np.where(
        overlap['dir_at'] == overlap['dir_rel'], 'concordant', 'discordant')
    return overlap


def plot_overlap_bar(overlap: pd.DataFrame, out_path: str,
                     sig_at: int, sig_rel: int) -> None:
    """Bar chart: AT-only / Relapse-only / concordant overlap / discordant overlap."""
    if overlap.empty:
        return

    concordant = overlap[overlap['concordance'] == 'concordant']
    discordant = overlap[overlap['concordance'] == 'discordant']

    cats   = ['AT only', 'Relapse only', 'Concordant', 'Discordant']
    counts = [
        sig_at  - len(overlap),
        sig_rel - len(overlap),
        len(concordant),
        len(discordant),
    ]
    colors = ['#d62728', '#1f77b4', '#2ca02c', '#ff7f0e']

    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(cats, counts, color=colors, edgecolor='white')
    ax.bar_label(bars, padding=3, fontsize=9)
    ax.set_ylabel('Number of DEGs')
    ax.set_title('DEG Overlap: Anaplastic Transformation vs FHWT Relapse\n'
                 '(Concordant = same direction in both comparisons)')
    ax.set_ylim(0, max(counts) * 1.2)
    fig.tight_layout()
    savefig(fig, out_path)

def plot_overlap_bar_charts(
    sig_at_up: pd.DataFrame,
    conc_up: pd.DataFrame,
    dawt_up: pd.DataFrame,
    disc: pd.DataFrame,
    out_path: str,
) -> None:
    """
    Bar chart comparing AT UP genes against de novo DAWT UP genes.

    Four bars: AT-only, DAWT-only, concordant (UP in both), discordant (AT-UP ∩ DAWT-DOWN).

    Parameters
    ----------
    sig_at_up : AT-significant UP genes DataFrame
    conc_up   : genes UP in both AT and DAWT (concordant)
    dawt_up   : de novo DAWT UP genes DataFrame
    disc      : genes UP in AT but DOWN in DAWT (discordant)
    out_path  : save path without extension
    """
    cats = ['AT only', 'DAWT only', 'AT∩DAWT (concordant)', 'AT∩DAWT (discordant)']
    counts = [
        len(sig_at_up) - len(conc_up) - len(disc),
        len(dawt_up) - len(conc_up),
        len(conc_up),
        len(disc),
    ]
    colors = ['#d62728', '#ab47bc', '#ff7f0e', '#aec7e8']

    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(cats, counts, color=colors, edgecolor='white')
    ax.bar_label(bars, padding=3, fontsize=9)
    ax.set_ylabel('Number of genes')
    ax.set_title('AT Recurrent↑ vs de novo DAWT↑')
    ax.set_ylim(0, max(counts) * 1.25)
    ax.tick_params(axis='x', labelsize=8)
    fig.tight_layout()
    savefig(fig, out_path)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:  # noqa: D103
    print("\n" + "="*60)
    print("  05_rnaseq_de_anaplastic_transformation.py")
    print("="*60)

    # ── Step 1: Build AT metadata ─────────────────────────────────────────────
    print("\n── Step 1: Build AT metadata ────────────────────────────────")
    meta_at = build_at_metadata()
    check_complete_pairs(meta_at)

    # ── Step 2: Load raw counts for AT samples ───────────────────────────────
    print("\n── Step 2: Load AT raw counts ───────────────────────────────")
    matrix, valid_samples = build_at_counts_matrix(meta_at)
    meta_at = meta_at[meta_at['sample_id'].isin(valid_samples)].copy()

    # ── Step 3: PyDESeq2 paired ───────────────────────────────────────────────
    print("\n── Step 3: PyDESeq2 paired (~patient + time_point) ─────────")
    print(f"  Design: ~patient + time_point  "
          f"(n={len(meta_at[meta_at['time_point']=='Primary'])} Primary, "
          f"{len(meta_at[meta_at['time_point']=='Recurrent'])} Recurrent)")
    res = run_paired_deseq2(matrix, meta_at)

    sig = res[res['direction'] != 'NS']
    n_up   = (res['direction'] == 'UP').sum()
    n_down = (res['direction'] == 'DOWN').sum()
    print(f"\n  Genes tested : {len(res):,}")
    print(f"  Significant  : {len(sig):,}  (UP={n_up} Recurrent↑, DOWN={n_down} Primary↑)")

    if n_up > 0:
        print(f"\n  Top 10 UP (anaplastic recurrence↑):")
        print(res[res['direction'] == 'UP'].head(10)
              [['gene_name', 'log2FoldChange', 'padj']].to_string(index=False))
    if n_down > 0:
        print(f"\n  Top 10 DOWN (original FHWT primary↑):")
        print(res[res['direction'] == 'DOWN'].head(10)
              [['gene_name', 'log2FoldChange', 'padj']].to_string(index=False))

    # ── Step 4: Save tables ───────────────────────────────────────────────────
    print("\n── Step 4: Save tables ──────────────────────────────────────")
    res.to_csv(os.path.join(processed_tables_path, 'de_at_full.csv'), index=False)
    sig.to_csv(os.path.join(processed_tables_path, 'de_at_sig.csv'), index=False)
    print(f"  de_at_full.csv : {len(res):,} genes")
    print(f"  de_at_sig.csv  : {len(sig):,} DEGs")

    # ── Step 5: Figures ───────────────────────────────────────────────────────
    print("\n── Step 5: Figures ──────────────────────────────────────────")
    log2_viz = pd.read_csv(os.path.join(matrices_path, 'counts_log2_all.csv'),
                           index_col=0)
    at_cols = [s for s in meta_at['sample_id'] if s in log2_viz.columns]
    if len(at_cols) == len(meta_at):
        # Preferred: DESeq2-normalized log2 values computed in script 02
        log2_at = log2_viz[at_cols]
        print(f"  log2 matrix: using counts_log2_all.csv ({len(at_cols)} samples)")
    else:
        missing = set(meta_at['sample_id']) - set(log2_viz.columns)
        print(f"  [warn] {len(missing)} AT samples not in counts_log2_all.csv — "
              f"falling back to log2(raw+1). Re-run script 02 to fix.")
        log2_at = np.log2(matrix + 1)

    plot_volcano(
        res,
        out_path=os.path.join(figures_path, 'de_at_01_volcano'),
        title=f'Anaplastic Transformation: Recurrent vs Primary — {len(sig)} DEGs\n'
              f'(n=4 pairs, paired ~patient + time_point; UP=Recurrent)',
        palette=AT_PALETTE,
        lfc_thresh=LFC_THRESH,
        padj_thresh=PADJ_THRESH,
    )
    print("  Saved: de_at_01_volcano")

    plot_ma(
        res,
        out_path=os.path.join(figures_path, 'de_at_02_ma'),
        title='MA Plot — Anaplastic Transformation (Recurrent vs Primary)',
        palette=AT_PALETTE,
        lfc_thresh=LFC_THRESH,
    )
    print("  Saved: de_at_02_ma")

    plot_heatmap(
        res, log2_at, meta_at,
        out_path=os.path.join(figures_path, 'de_at_03_heatmap'),
        title=f'Top {TOP_HEATMAP} DEGs — Anaplastic Transformation\n'
              f'(Recurrent↑ = red, Primary↑ = cyan)',
        group_col='time_point',
        palette=TIME_PALETTE,
        top_n=TOP_HEATMAP,
    )
    print("  Saved: de_at_03_heatmap")

    plot_gene_bar(
        res, WILMS_PANEL,
        out_path=os.path.join(figures_path, 'de_at_04_wilms_genes'),
        title='Wilms Panel — Anaplastic Transformation',
        palette=AT_PALETTE,
        lfc_thresh=LFC_THRESH,
        padj_thresh=PADJ_THRESH,
        x_label='log2 Fold Change (Recurrent vs Primary)',
    )
    print("  Saved: de_at_04_wilms_genes")

    # ── Step 6: Overlap with FHWT relapse (Step 2) ───────────────────────────
    print("\n── Step 6: Overlap with FHWT relapse DEGs (Step 2) ─────────")

    relapse_path = os.path.join(processed_tables_path, 'de_relapse_fhwt_sig.csv')
    if not os.path.exists(relapse_path):
        print("  [skip overlap] de_relapse_fhwt_sig.csv not found — run 04 first")
        sig_rel_df = pd.DataFrame()
    else:
        sig_rel_df = pd.read_csv(relapse_path)

    overlap = analyse_overlap(res, sig_rel=sig_rel_df)
    if not overlap.empty:
        overlap.to_csv(os.path.join(processed_tables_path, 'de_at_vs_relapse_overlap.csv'),
                       index=False)
        plot_overlap_bar(
            overlap,
            out_path=os.path.join(figures_path, 'de_at_05_overlap'),
            sig_at=len(sig),
            sig_rel=len(sig_rel_df),
        )
        print("  Saved: de_at_05_overlap")

    # ── Step 7: Overlap with de novo DAWT identity (Step 1) ──────────────────
    # Key question: do AT recurrent tumors converge on the de novo DAWT
    # transcriptional program? Compares AT UP genes (anaplastic recurrence)
    # with DAWT UP genes (de novo anaplastic identity, FHWT-reference convention).

    print("\n── Step 7: Overlap with de novo DAWT identity (Step 1) ──────")
    dawt_path = os.path.join(processed_tables_path, 'de_histology_sig.csv')

    if not os.path.exists(dawt_path):
        print("  [skip] de_histology_sig.csv not found — run 03 first")
    else:
        sig_at_up   = sig[sig['direction'] == 'UP'][['gene_name', 'log2FoldChange', 'padj']]
        sig_at_down = sig[sig['direction'] == 'DOWN'][['gene_name', 'log2FoldChange', 'padj']]

        dawt_sig    = pd.read_csv(dawt_path)
        dawt_up     = dawt_sig[dawt_sig['direction'] == 'UP'][['gene_name', 'log2FoldChange', 'padj']]
        dawt_down   = dawt_sig[dawt_sig['direction'] == 'DOWN'][['gene_name', 'log2FoldChange', 'padj']]

        # Concordant: UP in both (acquired anaplastic program matches de novo DAWT)
        conc_up   = sig_at_up.merge(dawt_up.rename(
            columns={'log2FoldChange': 'lfc_dawt', 'padj': 'padj_dawt'}),
            on='gene_name')
        # Discordant: UP in AT but DOWN in DAWT (AT-specific, not DAWT-like)
        disc = sig_at_up.merge(dawt_down.rename(
            columns={'log2FoldChange': 'lfc_dawt', 'padj': 'padj_dawt'}), on='gene_name')

        print(f"\n  AT UP genes : {len(sig_at_up)}")
        print(f"  DAWT UP genes : {len(dawt_up)}")
        print(f"  Concordant (UP in both — acquired = de novo DAWT): {len(conc_up)}")
        print(f"  Discordant (AT-UP but DAWT-DOWN) : {len(disc)}")

        if not conc_up.empty:
            print(f"\n  Concordant genes (AT Recurrent↑ ∩ de novo DAWT↑):")
            print(conc_up[['gene_name', 'log2FoldChange', 'lfc_dawt']].to_string(index=False))

        # Save
        conc_up.to_csv(os.path.join(processed_tables_path, 'de_at_vs_dawt_concordant.csv'), index=False)
        plot_overlap_bar_charts(
            sig_at_up, conc_up, dawt_up, disc,
            out_path=os.path.join(figures_path, 'de_at_06_dawt_overlap'),
        )
        print("  Saved: de_at_06_dawt_overlap")

    print("\n── Complete ─────────────────────────────────────────────────")


if __name__ == '__main__':
    main()
