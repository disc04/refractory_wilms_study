#!/usr/bin/env python3
"""
08_druggability.py — Wilms Tumor Analysis Project
======================================================
Druggability screening for AT-specific genes.

Pipeline
--------
1. Query ChEMBL for top hub genes from PPI network → map to approved drugs / clinical candidates
2. Produce figures:
   enrich_10_druggability.png — druggability barplot for top hub genes
6. Save tables:
   ppi_druggability.csv — ChEMBL hits for hub genes
   supplementary_ppi_hub_genes.csv — curated supplementary table

Requires internet access (ChEMBL API).
Results cached to CSV; subsequent runs load from cache.
Dependencies: requests, pandas, matplotlib, seaborn, numpy
"""

import time
import os
import requests
import warnings
import pandas as pd
from pathlib import Path

from configuration import PATHS, CHEMBL_API, CACHE_PATHS, TOP_HUBS

from utils.ppi_plot_utils import fig_druggability

warnings.filterwarnings("ignore")

processed_tables_path, figures_path = PATHS["processed_tables_path"], PATHS["figures_path"]
CACHE_CHEMBL = CACHE_PATHS["CACHE_CHEMBL"]

# 1. CHEMBL DRUGGABILITY QUERY
# ------------------------------------------------------------------------------
def _chembl_target_query(gene: str) -> list[dict]:
    """Search ChEMBL targets by gene symbol; return target metadata."""
    url = f"{CHEMBL_API}/target/search.json"
    params = {
        "q": gene,
        "organism": "Homo sapiens",
        "target_type": "SINGLE PROTEIN",
        "limit": 5,
    }
    try:
        r = requests.get(url, params=params, timeout=30)
        if r.status_code != 200:
            return []
        hits = r.json().get("targets", [])
        return hits
    except Exception:
        return []


def _chembl_mechanisms(chembl_id: str) -> list[dict]:
    """Fetch mechanism-of-action entries for a given ChEMBL target ID."""
    url = f"{CHEMBL_API}/mechanism.json"
    params = {"target_chembl_id": chembl_id, "limit": 10}
    try:
        r = requests.get(url, params=params, timeout=30)
        return r.json().get("mechanisms", []) if r.status_code == 200 else []
    except Exception:
        return []


def _chembl_approved_drugs(chembl_id: str) -> list[dict]:
    """Fetch approved drugs acting on a given ChEMBL target."""
    url = f"{CHEMBL_API}/drug_indication.json"
    params = {"target_chembl_id": chembl_id, "max_phase": 4, "limit": 10}
    try:
        r = requests.get(url, params=params, timeout=30)
        return r.json().get("drug_indications", []) if r.status_code == 200 else []
    except Exception:
        return []


def _chembl_activity_count(chembl_id: str) -> int:
    """Count bioactivity records for a ChEMBL target (proxy for druggability)."""
    url = f"{CHEMBL_API}/activity.json"
    params = {"target_chembl_id": chembl_id, "limit": 1}
    try:
        r = requests.get(url, params=params, timeout=30)
        return r.json().get("page_meta", {}).get("total_count", 0) if r.status_code == 200 else 0
    except Exception:
        return 0


def query_chembl_druggability(hub_genes: list[str]) -> pd.DataFrame:
    """For each hub gene, get ChEMBL target metadata and druggability score."""
    if Path(CACHE_CHEMBL).exists():
        print("  [cache] Loading ChEMBL hits from cache")
        return pd.read_csv(CACHE_CHEMBL)

    print(f"  Querying ChEMBL for {len(hub_genes)} hub genes…")
    rows = []
    for i, gene in enumerate(hub_genes):
        print(f"    [{i+1}/{len(hub_genes)}] {gene}", end=" ")
        targets = _chembl_target_query(gene)
        if not targets:
            print("— no ChEMBL target")
            rows.append({"gene": gene, "chembl_id": None, "target_name": None,
                          "target_class": None, "activity_count": 0,
                          "max_clinical_phase": 0, "approved_drugs": None,
                          "drug_mechanisms": None, "druggability_tier": "Unknown"})
            time.sleep(0.3)
            continue

        # Take the top hit (first result)
        t = targets[0]
        chembl_id = t.get("target_chembl_id")
        target_name = t.get("pref_name")
        target_class = t.get("target_type")

        # Activity count
        n_activities = _chembl_activity_count(chembl_id)

        # Mechanism of action entries
        mechs = _chembl_mechanisms(chembl_id)
        mech_drugs = [m.get("molecule_chembl_id", "") for m in mechs]
        max_phase = max([m.get("max_phase", 0) or 0 for m in mechs], default=0)

        # Approved drug names (phase 4)
        drug_names = []
        for m in mechs:
            if (m.get("max_phase") or 0) >= 4:
                drug_names.append(m.get("mechanism_of_action", ""))

        druggability = (
            "Approved target" if max_phase >= 4
            else "Clinical candidate" if max_phase >= 1
            else "Bioactive compounds" if n_activities > 100
            else "Poorly drugged"
        )

        print(f"— {chembl_id} | phase={max_phase} | activities={n_activities}")
        rows.append({
            "gene": gene,
            "chembl_id": chembl_id,
            "target_name": target_name,
            "target_class": target_class,
            "activity_count": n_activities,
            "max_clinical_phase": max_phase,
            "approved_drugs": "; ".join(drug_names) if drug_names else None,
            "drug_mechanisms": "; ".join(mech_drugs[:5]) if mech_drugs else None,
            "druggability_tier": druggability,
        })
        time.sleep(0.5)   # polite API rate

    df = pd.DataFrame(rows)
    df.to_csv(CACHE_CHEMBL, index=False)
    return df

# 2. SUPPLEMENTARY TABLE
# ------------------------------------------------------------------------------
def save_supplementary_table(centrality: pd.DataFrame,
                              chembl: pd.DataFrame,
                              at_sig: pd.DataFrame) -> None:
    """Merged supplementary table: centrality + druggability + DE stats."""
    de_cols = at_sig[at_sig["direction"] == "UP"][
        ["gene_name", "log2FoldChange", "padj"]
    ].rename(columns={"gene_name": "gene", "log2FoldChange": "AT_log2FC",
                       "padj": "AT_padj"})

    supp = (centrality[centrality["degree"] > 0]
            .merge(chembl, on="gene", how="left")
            .merge(de_cols, on="gene", how="left"))

    out_cols = [
        "gene", "degree", "composite_hub_score", "hub_tier",
        "degree_centrality", "betweenness_centrality", "eigenvector_centrality",
        "AT_log2FC", "AT_padj",
        "chembl_id", "target_name", "druggability_tier",
        "max_clinical_phase", "activity_count", "approved_drugs",
    ]
    out_cols = [c for c in out_cols if c in supp.columns]
    supp = supp[out_cols].sort_values("composite_hub_score", ascending=False)

    path = os.path.join(processed_tables_path, "supplementary_ppi_hub_genes.csv")
    supp.to_csv(path, index=False)
    print(f"  Saved supplementary_ppi_hub_genes.csv ({len(supp)} genes)")
    chembl.to_csv(os.path.join(processed_tables_path, "ppi_druggability.csv"), index=False)

# MAIN
# ------------------------------------------------------------------------------
def main():
    print("\n" + "="*70)
    print("08_druggability.py — druggability screen for AT-specific genes from PPI network")
    print("="*70)

    # ── 1. Load hub gene list ──────────────────────────────────────────────────────────
    print("\n── Load centrality scores ──")
    centrality = pd.read_csv(os.path.join(processed_tables_path, "ppi_centrality.csv"))
    top10 = centrality[centrality["degree"] > 0].head(10)
    print("\n  Top 10 hub genes:")
    print(top10[["gene", "degree", "composite_hub_score", "hub_tier"]].to_string(index=False))

    # ── 2. ChEMBL druggability ────────────────────────────────────────────────
    print(f"\n── ChEMBL druggability (top {TOP_HUBS} hubs) ──")
    hub_genes = centrality[centrality["degree"] > 0].head(TOP_HUBS)["gene"].tolist()
    chembl = query_chembl_druggability(hub_genes)

    # Summary
    print("\n  Druggability summary:")
    print(chembl["druggability_tier"].value_counts().to_string())
    approved = chembl[chembl["druggability_tier"] == "Approved target"]
    if len(approved):
        print(f"\n  Approved targets ({len(approved)}):")
        print(approved[["gene", "target_name", "approved_drugs"]].to_string(index=False))

    # ── 3. Figure ────────────────────────────────────────────────────────────
    print("\n── Generating figure ──")
    fig_druggability(centrality, chembl)

    # ── 4. Supplementary table ────────────────────────────────────────────────
    print("\n── Saving supplementary_ppi_hub_genes table ──")
    at_sig = pd.read_csv(os.path.join(processed_tables_path, "de_at_sig.csv"))
    save_supplementary_table(centrality, chembl, at_sig)

    print("\n✓ Script complete. Figures: enrich_10  |  Tables: ppi_druggability.csv")
    print("="*70)


if __name__ == "__main__":
    main()
