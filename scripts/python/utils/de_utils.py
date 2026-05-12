"""
de_utils.py — shared PyDESeq2 helpers for differential expression scripts.

Functions
---------
assign_direction  — label genes as UP / DOWN / NS based on LFC and adjusted p-value
run_deseq2        — run PyDESeq2 Wald test with LFC shrinkage (two-group design)
"""

from __future__ import annotations
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

DEFAULT_LFC_THRESH  = 1.0
DEFAULT_PADJ_THRESH = 0.05


def assign_direction(res: pd.DataFrame,
                     lfc_thresh: float  = DEFAULT_LFC_THRESH,
                     padj_thresh: float = DEFAULT_PADJ_THRESH) -> pd.DataFrame:
    """
    Label genes as UP / DOWN / NS and return the DataFrame.

    Parameters
    ----------
    res         : results DataFrame with 'log2FoldChange' and 'padj' columns
    lfc_thresh  : |log2 fold change| threshold (default 1.0)
    padj_thresh : adjusted p-value threshold (default 0.05)

    Returns
    -------
    Copy of res with 'direction' column added / overwritten.
    """
    res = res.copy()
    res['direction'] = 'NS'
    res.loc[(res['padj'] < padj_thresh) & (res['log2FoldChange'] >  lfc_thresh),
            'direction'] = 'UP'
    res.loc[(res['padj'] < padj_thresh) & (res['log2FoldChange'] < -lfc_thresh),
            'direction'] = 'DOWN'
    return res


def run_deseq2(counts: pd.DataFrame,
               meta: pd.DataFrame,
               design: str,
               ref: tuple[str, str],
               lfc_thresh: float  = DEFAULT_LFC_THRESH,
               padj_thresh: float = DEFAULT_PADJ_THRESH) -> pd.DataFrame:
    """
    Run PyDESeq2 Wald test with LFC shrinkage (two-group design).

    Parameters
    ----------
    counts      : genes × samples raw count matrix
    meta        : sample metadata with 'sample_id' column
    design      : model formula string, e.g. '~histology'
    ref         : (factor_column, reference_level) tuple
    lfc_thresh  : |log2 fold change| threshold for direction labeling
    padj_thresh : adjusted p-value threshold for direction labeling

    Returns
    -------
    DataFrame with gene_name, log2FoldChange, padj, baseMean, direction columns,
    sorted by padj ascending then log2FoldChange descending.

    Notes
    -----
    Reference leveling is done via pd.Categorical (putting ref_level first in
    categories) rather than DeseqDataSet's ref_level kwarg, which is unreliable
    across pydeseq2 versions.
    """
    factor, ref_level = ref
    counts_T = counts.T.astype(int)
    meta_dds = meta.set_index('sample_id')[[factor]].reindex(counts_T.index)

    # Relevel: reference level first so pydeseq2 treats it as baseline
    other_levels = [l for l in meta_dds[factor].unique() if l != ref_level]
    meta_dds[factor] = pd.Categorical(meta_dds[factor],
                                      categories=[ref_level] + sorted(other_levels))

    dds = DeseqDataSet(
        counts=counts_T,
        metadata=meta_dds,
        design=design,
        quiet=True,
    )
    dds.deseq2()

    test_level = [l for l in meta_dds[factor].unique() if l != ref_level][0]

    stat = DeseqStats(
        dds,
        contrast=[factor, test_level, ref_level],
        alpha=padj_thresh,
        quiet=True,
    )
    stat.summary()
    stat.lfc_shrink(coeff=f'{factor}[T.{test_level}]')

    res = stat.results_df.reset_index().rename(columns={'index': 'gene_name'})
    res = assign_direction(res, lfc_thresh=lfc_thresh, padj_thresh=padj_thresh)
    res = res.sort_values(['padj', 'log2FoldChange'], ascending=[True, False])
    return res
