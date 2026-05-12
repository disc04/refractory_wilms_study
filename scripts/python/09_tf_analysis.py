#!/usr/bin/env python3
"""
09_tf_analysis.py — Wilms Tumor Analysis Project
=================================================
TF regulatory network analysis for AT-specific UP genes.

Pipeline
--------
1. Load AT-specific UP gene set (165 genes)
2. Load TRRUST v2 TF–target database
   - Primary:  download trrust_rawdata.human.tsv from grnpedia.org (cached)
   - Fallback: reconstruct from GSEA Step 3 leading-edge results (no internet needed)
3. Cross-reference AT-specific genes against all TRRUST TF regulons
4. Focus analysis on RELA, NFKB1, JUN, STAT3 — build Venn / UpSet structure
5. Produce figures:
   enrich_11_tf_regulon_venn.png — 4-set UpSet plot (RELA / NFKB1 / JUN / STAT3)
   enrich_12_tf_regulon_heatmap.png — gene × TF membership heatmap for significant TFs
6. Save tables:
   tf_nfkb_direct_targets.csv — full regulon membership with DE stats
   tf_all_significant.csv — all significant TF–AT overlaps (Fisher exact p < 0.05)

Dependencies: pandas, numpy, matplotlib, seaborn, scipy, requests
Internet access: optional (TRRUST download); analysis runs offline from cached GSEA results.
"""

import os
import warnings
import requests
import pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from configuration import PATHS, CACHE_PATHS, FDR_THRESH, FOCUS_TFS, TRRUST_URL, MIN_OVERLAP
from utils.tf_plot_utils import fig_upset_nfkb, fig_heatmap
warnings.filterwarnings("ignore")

processed_tables_path, enrich_tables_path, figures_path = (PATHS["processed_tables_path"],
                                                           PATHS["enrich_tables_path"],
                                                           PATHS["figures_path"])
CACHE_TRRUST = CACHE_PATHS["CACHE_TRRUST"]

# 1. GENE LISTS
# ------------------------------------------------------------------------------
def load_at_specific_genes() -> tuple[list[str], pd.DataFrame]:
    """AT-specific UP gene list + DE statistics."""
    at_sig    = pd.read_csv(os.path.join(processed_tables_path, "de_at_sig.csv"))
    concordant = pd.read_csv(os.path.join(processed_tables_path, "de_at_vs_dawt_concordant.csv"))
    at_up     = set(at_sig[at_sig["direction"] == "UP"]["gene_name"])
    dawt_shared = set(concordant["gene_name"])
    specific  = sorted(at_up - dawt_shared)
    de_stats  = (at_sig[at_sig["direction"] == "UP"]
                 [["gene_name", "log2FoldChange", "padj"]]
                 .rename(columns={"gene_name": "gene"})
                 .set_index("gene"))
    print(f"  AT-specific UP genes: {len(specific)}")
    return specific, de_stats

# 2. TRRUST DATABASE
# ------------------------------------------------------------------------------
def load_trrust() -> dict[str, set]:
    """
    Load TRRUST v2 human TF–target pairs.
    Returns {TF_symbol: {target_gene, ...}}.

    Strategy:
      1. Use cached TSV if present.
      2. Try to download from grnpedia.org.
      3. Fall back to reconstructing from GSEA Step 3 leading-edge results.
    """
    if Path(CACHE_TRRUST).exists():
        print("  [cache] Loading TRRUST from cache")
        df = pd.read_csv(CACHE_TRRUST, sep="\t", header=None,
                         names=["tf", "target", "regulation", "pmids"])
    else:
        print(f"  Attempting TRRUST download from {TRRUST_URL}…")
        try:
            r = requests.get(TRRUST_URL, timeout=30)
            r.raise_for_status()
            with open(CACHE_TRRUST, "w") as f:
                f.write(r.text)
            df = pd.read_csv(CACHE_TRRUST, sep="\t", header=None,
                             names=["tf", "target", "regulation", "pmids"])
            print(f"  Downloaded {len(df)} TF–target pairs")
        except Exception as e:
            print(f"  Download failed ({e}) — reconstructing from GSEA leading-edge results")
            df = _reconstruct_trrust_from_gsea()

    regulons = {}
    for _, row in df.iterrows():
        tf = str(row["tf"]).strip()
        target = str(row["target"]).strip()
        regulons.setdefault(tf, set()).add(target)

    print(f"  TRRUST: {len(regulons)} TFs, "
          f"{sum(len(v) for v in regulons.values())} TF–target pairs")
    return regulons


def _reconstruct_trrust_from_gsea() -> pd.DataFrame:
    """
    Fallback: reconstruct TF→target mappings from GSEA step3 leading-edge genes.
    These are the AT-programme-relevant subset of each TRRUST gene set and are
    sufficient for the focused RELA/NFKB1/JUN/STAT3 analysis.
    """
    gsea = pd.read_csv(os.path.join(enrich_tables_path, "gsea_step3_TRRUST.csv"))
    rows = []
    for _, row in gsea.iterrows():
        # term format: "RELA human" / "JUN mouse"
        parts = str(row["term"]).split()
        tf = parts[0]
        leads = str(row.get("lead_genes", "")).split(";")
        for target in leads:
            target = target.strip()
            if target:
                rows.append({"tf": tf, "target": target,
                             "regulation": "Unknown", "pmids": ""})
    df = pd.DataFrame(rows).drop_duplicates()
    print(f"  Reconstructed {len(df)} TF–target pairs from GSEA leading edges")
    return df

# 3. OVERLAP ANALYSIS
# ------------------------------------------------------------------------------
def compute_overlaps(at_genes: list[str],
                     regulons: dict[str, set],
                     background_n: int = 20000) -> pd.DataFrame:
    """
    Fisher exact test for enrichment of each TF regulon in AT-specific genes.
    Returns sorted DataFrame with overlap stats.
    """
    at_set = set(at_genes)
    rows = []
    for tf, targets in regulons.items():
        overlap = at_set & targets
        if len(overlap) < MIN_OVERLAP:
            continue
        # 2×2 contingency: in AT + in regulon / not in AT + in regulon / etc.
        a = len(overlap)                             # AT ∩ regulon
        b = len(at_set) - a                          # AT \ regulon
        c = len(targets) - a                         # regulon \ AT
        d = background_n - a - b - c                 # neither
        _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({
            "tf": tf,
            "regulon_size": len(targets),
            "overlap_n": a,
            "overlap_pct_AT": round(100 * a / len(at_set), 1),
            "overlap_pct_regulon": round(100 * a / len(targets), 1),
            "pval": pval,
            "overlap_genes": ";".join(sorted(overlap)),
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    _, fdr, _, _ = multipletests(df["pval"], method="fdr_bh")
    df["fdr"] = fdr
    return df.sort_values("fdr").reset_index(drop=True)


def build_regulon_membership(at_genes: list[str],
                              regulons: dict[str, set],
                              focus_tfs: list[str]) -> pd.DataFrame:
    """
    Build a gene × TF binary membership matrix for focus TFs.
    Adds DE stats and community annotation.
    """
    at_set = set(at_genes)
    rows = []
    for gene in at_genes:
        row = {"gene": gene}
        for tf in focus_tfs:
            row[tf] = int(gene in regulons.get(tf, set()))
        rows.append(row)
    return pd.DataFrame(rows)


# MAIN
# ------------------------------------------------------------------------------
def main():
    print("\n" + "=" * 70)
    print("09_tf_analysis.py — TF regulon analysis for AT-specific UP genes")
    print("=" * 70)

    # ── 1. Gene list ──────────────────────────────────────────────────────────
    print("\n── Step 1: Load AT-specific UP genes ──")
    at_genes, de_stats = load_at_specific_genes()

    # ── 2. TRRUST database ────────────────────────────────────────────────────
    print("\n── Step 2: Load TRRUST v2 ──")
    regulons = load_trrust()

    # ── 3. Overlap analysis ───────────────────────────────────────────────────
    print("\n── Step 3: Fisher exact enrichment for all TF regulons ──")
    overlap_df = compute_overlaps(at_genes, regulons)
    sig = overlap_df[overlap_df["fdr"] < FDR_THRESH]
    print(f"  Significant TF regulons (FDR < {FDR_THRESH}): {len(sig)}")
    print(f"\n  Top 15:")
    print(sig.head(15)[["tf", "regulon_size", "overlap_n",
                         "overlap_pct_AT", "pval", "fdr"]].to_string(index=False))

    # ── 4. Focus analysis: RELA / NFKB1 / JUN / STAT3 ────────────────────────
    print("\n── Step 4: Focus analysis — RELA / NFKB1 / JUN / STAT3 ──")
    membership = build_regulon_membership(at_genes, regulons, FOCUS_TFS)

    # Annotate with DE stats and community
    comm_path = os.path.join(enrich_tables_path, "ppi_communities.csv")
    if Path(comm_path).exists():
        comm = pd.read_csv(comm_path).set_index("gene")
        comm_col = "community_label" if "community_label" in comm.columns else "community"
        membership["community"] = membership["gene"].map(
            lambda g: str(comm.loc[g, comm_col]) if g in comm.index else "isolated"
        )
    membership["log2FC"] = membership["gene"].map(
        lambda g: round(de_stats.loc[g, "log2FoldChange"], 3) if g in de_stats.index else None
    )
    membership["padj"] = membership["gene"].map(
        lambda g: float(f"{de_stats.loc[g, 'padj']:.3e}") if g in de_stats.index else None
    )

    # Summary
    nfkb = set(membership[(membership["RELA"] == 1) | (membership["NFKB1"] == 1)]["gene"])
    jun_only = set(membership[(membership["JUN"] == 1) &
                               (membership["RELA"] == 0) &
                               (membership["NFKB1"] == 0)]["gene"])
    all_four = set(membership[(membership[FOCUS_TFS].sum(axis=1) == len(FOCUS_TFS))]["gene"])

    print(f"\n  NF-κB direct targets (RELA|NFKB1): {len(nfkb)}/{len(at_genes)} "
          f"({100 * len(nfkb) / len(at_genes):.0f}%)")
    print(f"  JUN-only (AP-1, not NF-κB):          {len(jun_only)} genes: {sorted(jun_only)}")
    print(f"  All four TFs:                         {len(all_four)} genes: {sorted(all_four)}")

    # Save focus table
    focus_out = os.path.join(processed_tables_path, "tf_nfkb_direct_targets.csv")
    in_any = membership[FOCUS_TFS].sum(axis=1) > 0
    (membership[in_any]
     .sort_values(["RELA", "NFKB1", "JUN", "STAT3", "log2FC"],
                  ascending=[False, False, False, False, False])
     .to_csv(focus_out, index=False))
    print(f"\n  Saved tf_nfkb_direct_targets.csv ({in_any.sum()} genes)")

    # Save all-TF overlap table
    all_tf_out = os.path.join(processed_tables_path, "tf_all_significant.csv")
    overlap_df.to_csv(all_tf_out, index=False)
    print(f"  Saved tf_all_significant.csv ({len(overlap_df)} TFs tested, "
          f"{len(sig)} significant)")

    # ── 5. Figures ────────────────────────────────────────────────────────────
    print("\n── Step 5: Generating figures ──")
    fig_upset_nfkb(membership, de_stats, FOCUS_TFS, out_path=figures_path)
    fig_heatmap(membership, de_stats, out_path=figures_path, overlap_df = overlap_df)

    print("\n✓ Script complete.")
    print("  Figures: enrich_11_tf_regulon_upset, enrich_12_tf_regulon_heatmap")
    print("  Tables:  tf_nfkb_direct_targets.csv, tf_all_significant.csv")
    print("=" * 70)


if __name__ == "__main__":
    main()
