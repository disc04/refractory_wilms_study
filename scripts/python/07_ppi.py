#!/usr/bin/env python3
"""
07_ppi.py — Wilms Tumor Analysis Project
======================================================
PPI network analysis for AT-specific genes.

Pipeline
--------
1. Build PPI network via STRING API (confidence >= 700) for AT-specific UP genes (n=163)
2. Compute network centrality metrics (degree, betweenness, eigenvector)
3. Derive a composite hub score; rank all genes
4. Query ChEMBL for top hub genes → map to approved drugs / clinical candidates
5. Produce figures:
   enrich_08_ppi_network.png — network graph (node size = degree, color = hub tier)
   enrich_09_ppi_centrality.png — top 25 genes by composite centrality
6. Save tables:
   ppi_centrality.csv — full centrality table for all network genes

Requires internet access (STRING API).
Results cached to CSV; subsequent runs load from cache.
Dependencies: networkx, requests, pandas, matplotlib, seaborn, numpy
"""

import os
import requests
import warnings
import pandas as pd
import networkx as nx
from pathlib import Path
from configuration import (PATHS, STRING_API, CACHE_PATHS, COMMUNITY_LABELS, TAXON_HUMAN,
                           STRING_SCORE, HUB_TIER_HIGH, HUB_TIER_MED)

from utils.ppi_plot_utils import fig_ppi_network, fig_centrality_barplot, fig_ppi_communities

warnings.filterwarnings("ignore")

processed_tables_path, figures_path = PATHS["processed_tables_path"], PATHS["figures_path"]
CACHE_PPI, CACHE_CENTRALITY, CACHE_COMMUNITIES = (CACHE_PATHS["CACHE_PPI"], CACHE_PATHS["CACHE_CENTRALITY"],
                                                  CACHE_PATHS["CACHE_COMMUNITIES"])




# 1. DE GENE LISTS
# ------------------------------------------------------------------------------
def load_at_specific_genes() -> list[str]:
    """AT recurrent UP minus 7 DAWT-concordant genes = AT-specific program."""
    at_sig = pd.read_csv(os.path.join(processed_tables_path, "de_at_sig.csv"))
    concordant = pd.read_csv(os.path.join(processed_tables_path, "de_at_vs_dawt_concordant.csv"))
    at_up = set(at_sig[at_sig["direction"] == "UP"]["gene_name"])
    dawt_shared = set(concordant["gene_name"])
    specific = sorted(at_up - dawt_shared)
    print(f"  AT-specific UP genes: {len(specific)}")
    return specific


# 2. STRING PPI NETWORK
# ------------------------------------------------------------------------------
def _string_map_identifiers(genes: list[str]) -> dict[str, str]:
    """Map gene symbols to STRING identifiers (returns {symbol: string_id})."""
    url = f"{STRING_API}/get_string_ids"
    params = {
        "identifiers": "\r".join(genes),
        "species": TAXON_HUMAN,
        "limit": 1,
        "echo_query": 1,
        "caller_identity": "wilms_ppi_analysis",
    }
    r = requests.post(url, data=params, timeout=60)
    r.raise_for_status()
    mapping = {}
    for entry in r.json():
        query = entry.get("queryItem", "")
        sid = entry.get("stringId", "")
        pref = entry.get("preferredName", query)
        if sid:
            mapping[pref] = sid
    return mapping


def _string_network(genes: list[str]) -> pd.DataFrame:
    """Fetch interaction network for gene list from STRING."""
    url = f"{STRING_API}/network"
    params = {
        "identifiers": "\r".join(genes),
        "species": TAXON_HUMAN,
        "required_score": STRING_SCORE,
        "caller_identity": "wilms_ppi_analysis",
        "add_nodes": 0,          # do not add new nodes
        "network_type": "functional",
    }
    r = requests.post(url, data=params, timeout=120)
    r.raise_for_status()
    data = r.json()
    if not data:
        return pd.DataFrame(columns=["gene_a", "gene_b", "score"])
    edges = pd.DataFrame(data)[["preferredName_A", "preferredName_B", "score"]]
    edges.columns = ["gene_a", "gene_b", "score"]
    return edges


def build_ppi_network(genes: list[str]) -> tuple[nx.Graph, pd.DataFrame]:
    """Build STRING PPI network; cache edges."""
    if Path(CACHE_PPI).exists():
        print("  [cache] Loading STRING edges from cache")
        edges = pd.read_csv(CACHE_PPI)
    else:
        print(f"  Querying STRING API (score >= {STRING_SCORE})…")
        edges = _string_network(genes)
        edges.to_csv(CACHE_PPI, index=False)
        print(f"  Retrieved {len(edges)} interactions")

    G = nx.Graph()
    G.add_nodes_from(genes)
    for _, row in edges.iterrows():
        if row.gene_a in genes and row.gene_b in genes:
            G.add_edge(row.gene_a, row.gene_b, weight=row.score / 1000)

    # Remove isolated nodes for network stats but keep for reference
    connected = [n for n in G.nodes if G.degree(n) > 0]
    print(f"  Network: {G.number_of_nodes()} nodes, "
          f"{G.number_of_edges()} edges, "
          f"{len(connected)} connected nodes")
    return G, edges

# 3. CENTRALITY ANALYSIS
# ------------------------------------------------------------------------------
def _apply_hub_tiers(df: pd.DataFrame) -> pd.DataFrame:
    """Assign hub tiers using quantiles computed on connected nodes only."""
    connected_scores = df.loc[df["degree"] > 0, "composite_hub_score"]
    q_high = connected_scores.quantile(HUB_TIER_HIGH)
    q_med  = connected_scores.quantile(HUB_TIER_MED)
    df["hub_tier"] = df["composite_hub_score"].apply(
        lambda s: "High" if s >= q_high else ("Medium" if s >= q_med else "Low")
    )
    return df


def compute_centrality(G: nx.Graph) -> pd.DataFrame:
    """Degree, betweenness, eigenvector centrality + composite hub score."""
    if Path(CACHE_CENTRALITY).exists():
        print("  [cache] Loading centrality from cache")
        df = pd.read_csv(CACHE_CENTRALITY)
        # Always reapply tiers — quantile fix may differ from cached version
        return _apply_hub_tiers(df)

    print("  Computing centrality metrics…")
    # Use largest connected component for eigenvector centrality
    lcc = G.subgraph(max(nx.connected_components(G), key=len)).copy()

    deg  = nx.degree_centrality(G)
    btw  = nx.betweenness_centrality(G, normalized=True, weight="weight")
    try:
        eig = nx.eigenvector_centrality(lcc, max_iter=1000, weight="weight")
    except nx.PowerIterationFailedConvergence:
        eig = {n: 0.0 for n in lcc.nodes}

    rows = []
    for gene in G.nodes:
        rows.append({
            "gene": gene,
            "degree": G.degree(gene),
            "degree_centrality": deg.get(gene, 0),
            "betweenness_centrality": btw.get(gene, 0),
            "eigenvector_centrality": eig.get(gene, 0),
            "in_lcc": gene in lcc.nodes,
        })
    df = pd.DataFrame(rows)

    # Normalize each metric to [0,1] and compute composite score
    for col in ["degree_centrality", "betweenness_centrality", "eigenvector_centrality"]:
        mx = df[col].max()
        df[f"{col}_norm"] = df[col] / mx if mx > 0 else 0

    df["composite_hub_score"] = (
        df["degree_centrality_norm"] * 0.4 +
        df["betweenness_centrality_norm"] * 0.4 +
        df["eigenvector_centrality_norm"] * 0.2
    )
    df = df.sort_values("composite_hub_score", ascending=False).reset_index(drop=True)

    df = _apply_hub_tiers(df)
    df.to_csv(CACHE_CENTRALITY, index=False)
    return df


# 4. COMMUNITY DETECTION
# ------------------------------------------------------------------------------

def detect_communities(G: nx.Graph) -> pd.DataFrame:
    """Louvain community detection; cache assignment table."""
    if Path(CACHE_COMMUNITIES).exists():
        print("  [cache] Loading community assignments")
        return pd.read_csv(CACHE_COMMUNITIES)

    print("  Running Louvain community detection…")
    import networkx.algorithms.community as nx_comm
    communities = nx_comm.louvain_communities(G, seed=42, weight="weight")
    communities = sorted(communities, key=len, reverse=True)

    mod = nx_comm.modularity(G, communities, weight="weight")
    print(f"  {len(communities)} communities  |  modularity Q = {mod:.3f}")
    for i, c in enumerate(communities):
        label = COMMUNITY_LABELS.get(i + 1, f"Module {i+1}")
        print(f"    C{i+1} ({label}, n={len(c)}): {', '.join(sorted(c)[:8])}…")

    rows = []
    for i, c in enumerate(communities):
        for g in c:
            rows.append({"gene": g, "community": i + 1,
                         "community_label": COMMUNITY_LABELS.get(i + 1, f"Module {i+1}")})
    df = pd.DataFrame(rows)
    df.to_csv(CACHE_COMMUNITIES, index=False)
    return df


# MAIN
# ------------------------------------------------------------------------------
def main():
    print("\n" + "="*70)
    print("07_ppi.py — AT-specific gene PPI network + druggability")
    print("="*70)

    # ── 1. Gene list ──────────────────────────────────────────────────────────
    print("\n── Step 1: Load AT-specific UP genes ──")
    at_genes = load_at_specific_genes()

    # ── 2. STRING PPI network ─────────────────────────────────────────────────
    print("\n── Step 2: Build STRING PPI network ──")
    G, edges = build_ppi_network(at_genes)

    # ── 3. Centrality ─────────────────────────────────────────────────────────
    print("\n── Step 3: Centrality analysis ──")
    centrality = compute_centrality(G)
    centrality.to_csv(os.path.join(processed_tables_path, "ppi_centrality.csv"), index=False)

    top10 = centrality[centrality["degree"] > 0].head(10)
    print("\n  Top 10 hub genes:")
    print(top10[["gene", "degree", "composite_hub_score", "hub_tier"]].to_string(index=False))

    # ── 5. Community detection ────────────────────────────────────────────────
    print("\n── Step 4: Community detection ──")
    connected_G = G.subgraph([n for n in G.nodes if G.degree(n) > 0]).copy()
    communities_df = detect_communities(connected_G)

    # ── 6. Figures ────────────────────────────────────────────────────────────
    print("\n── Step 5: Generating figures ──")
    fig_ppi_network(G, centrality)
    fig_ppi_communities(connected_G, centrality, communities_df)
    fig_centrality_barplot(centrality)

    print("\n✓ Script complete. Figures: enrich_08, 09 | Tables: ppi_centrality.csv")
    print("="*70)


if __name__ == "__main__":
    main()
