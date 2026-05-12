"""
ppi_plot_utils.py — shared plotting helpers for the Wilms Tumor PPI network analysis.

"""
import os
import warnings
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from configuration import (PATHS, COMMUNITY_LABELS, STRING_SCORE, TOP_HUBS,
                                 COMMUNITY_COLORS, TIER_COLORS, DRUG_COLORS)
from utils.data_utils import savefig
warnings.filterwarnings("ignore")

figures_path = PATHS["figures_path"]


def fig_ppi_network(G: nx.Graph, centrality: pd.DataFrame) -> None:
    """Network graph: node size ∝ hub score, colour = hub tier, labels for top 20 hubs."""
    cent_map = centrality.set_index("gene")
    connected = [n for n in G.nodes if G.degree(n) > 0]
    H = G.subgraph(connected).copy()

    if len(H.nodes) == 0:
        print("  No connected nodes — skipping network figure")
        return

    nodes_list = list(H.nodes)

    # Node attributes
    tiers  = [cent_map.loc[n, "hub_tier"] if n in cent_map.index else "Low"
               for n in nodes_list]
    scores = np.array([float(cent_map.loc[n, "composite_hub_score"])
                        if n in cent_map.index else 0.0
                        for n in nodes_list])
    node_colors = [TIER_COLORS.get(t, "#aec7e8") for t in tiers]

    # Fix: normalize sizes to [60, 500] — no more stop-sign nodes
    s_min, s_max = scores.min(), scores.max()
    norm = (scores - s_min) / (s_max - s_min) if s_max > s_min else np.zeros_like(scores)
    node_sizes = 60 + norm * 440

    # Top 20 labels
    top_labels = set(centrality[centrality["degree"] > 0].head(20)["gene"])
    labels = {n: n if n in top_labels else "" for n in nodes_list}

    # Fix: spring layout with reproducible seed — Kamada-Kawai clusters poorly here
    pos = nx.spring_layout(H, k=2.5 / np.sqrt(len(H.nodes)), seed=42, iterations=150)

    fig, ax = plt.subplots(figsize=(16, 14))

    # Fix: draw ALL edges in one call — per-edge loop breaks matplotlib autoscaling
    nx.draw_networkx_edges(H, pos, alpha=0.35, edge_color="#888888", width=0.7, ax=ax)

    nx.draw_networkx_nodes(H, pos, nodelist=nodes_list, node_color=node_colors,
                           node_size=node_sizes, alpha=0.92, ax=ax)

    # Shift label positions above each node.
    # node_size is in pt²; radius in pt = sqrt(size/π).
    # 1 pt ≈ 1/72 inch; with a 14-inch-tall axis spanning ~2.4 data units,
    # 1 pt ≈ 0.0048 data units → add one radius + small gap.
    pts_per_data = (max(pos[n][1] for n in nodes_list) -
                    min(pos[n][1] for n in nodes_list)) / (14 * 72) + 1e-6
    size_map = dict(zip(nodes_list, node_sizes))
    label_pos = {
        n: (pos[n][0],
            pos[n][1] + np.sqrt(size_map[n] / np.pi) * pts_per_data + 0.025)
        for n in nodes_list
    }
    nx.draw_networkx_labels(H, label_pos, labels=labels, font_size=18,
                            font_weight="bold", ax=ax)

    # Fix: set axis limits explicitly with padding
    xs = [pos[n][0] for n in nodes_list]
    ys = [pos[n][1] for n in nodes_list]
    pad_x = (max(xs) - min(xs)) * 0.12 + 0.05
    pad_y = (max(ys) - min(ys)) * 0.12 + 0.05
    ax.set_xlim(min(xs) - pad_x, max(xs) + pad_x)
    ax.set_ylim(min(ys) - pad_y, max(ys) + pad_y)

    # Size legend (small/medium/large dots)
    size_legend_handles = []
    for label, frac in [("Low hub", 0.0), ("Medium hub", 0.5), ("High hub", 1.0)]:
        sz = np.sqrt((60 + frac * 440) / np.pi)
        h = mpatches.Circle((0, 0), sz / 60, color="#999999", label=label)
        size_legend_handles.append(h)

    tier_patches = [mpatches.Patch(color=c, label=t) for t, c in TIER_COLORS.items()]
    legend1 = ax.legend(handles=tier_patches, title="Hub tier", loc="upper left",
                        framealpha=0.85, fontsize=18, title_fontsize=18) #, prop={'weight': 'bold'})
    ax.add_artist(legend1)

    n_edges = H.number_of_edges()
    ax.set_title(
        f"AT-specific UP genes PPI Network"
        f"(n={len(connected)} nodes, {n_edges} edges, STRING score ≥ {STRING_SCORE})",
        fontsize=22, pad=10, weight="bold")
    ax.axis("off")
    plt.tight_layout()
    savefig(fig, path=os.path.join(figures_path, "enrich_08_ppi_network.png"))
    print("  Saved enrich_08_ppi_network")


def fig_ppi_communities(G: nx.Graph, centrality: pd.DataFrame,
                        communities_df: pd.DataFrame) -> None:
    """Network with Louvain community layout: community centroids on a circle,
    members arranged within each community. SPP1 bottleneck highlighted."""
    import networkx.algorithms.community as nx_comm
    from scipy.spatial import ConvexHull

    cent_map = centrality.set_index("gene")
    comm_map = communities_df.set_index("gene")
    nodes_list = list(G.nodes)

    # Node sizes
    scores = np.array([float(cent_map.loc[n, "composite_hub_score"])
                        if n in cent_map.index else 0.0 for n in nodes_list])
    s_min, s_max = scores.min(), scores.max()
    norm_s = (scores - s_min) / (s_max - s_min) if s_max > s_min else np.zeros_like(scores)
    node_sizes = 80 + norm_s * 500

    # ── Community-aware layout ────────────────────────────────────────────────
    major_ids = list(range(1, len(COMMUNITY_LABELS) + 1))
    outer_r = 2.8
    centroids = {}
    for i, cid in enumerate(major_ids):
        angle = np.pi / 2 + 2 * np.pi * i / len(major_ids)
        centroids[cid] = np.array([outer_r * np.cos(angle), outer_r * np.sin(angle)])

    pos = {}
    for cid in major_ids:
        members = [n for n in nodes_list
                   if n in comm_map.index and int(comm_map.loc[n, "community"]) == cid]
        if not members:
            continue
        subG = G.subgraph(members)
        inner_r = 0.55 + 0.06 * len(members)
        if len(members) == 1:
            sub_pos = {members[0]: np.zeros(2)}
        elif subG.number_of_edges() > 0:
            sub_pos = nx.spring_layout(subG, seed=42 + cid, k=0.7, scale=inner_r)
        else:
            angs = np.linspace(0, 2 * np.pi, len(members), endpoint=False)
            sub_pos = {n: np.array([inner_r * 0.6 * np.cos(a), inner_r * 0.6 * np.sin(a)])
                       for n, a in zip(members, angs)}
        cx, cy = centroids[cid]
        for n, p in sub_pos.items():
            pos[n] = np.array([cx + p[0], cy + p[1]])

    # Minor communities at periphery
    minor_ids = [c for c in sorted(communities_df["community"].unique()) if c > len(COMMUNITY_LABELS)]
    base_angle = np.pi / 2 + np.pi / len(major_ids)
    for idx, cid in enumerate(minor_ids):
        members = [n for n in nodes_list
                   if n in comm_map.index and int(comm_map.loc[n, "community"]) == cid]
        mid_r = outer_r * 1.22
        angle = base_angle + idx * (2 * np.pi / max(len(minor_ids), 6))
        cx = mid_r * np.cos(angle); cy = mid_r * np.sin(angle)
        for j, n in enumerate(members):
            a2 = 2 * np.pi * j / max(len(members), 1)
            pos[n] = np.array([cx + 0.18 * np.cos(a2), cy + 0.18 * np.sin(a2)])

    rng = np.random.default_rng(42)
    for n in nodes_list:
        if n not in pos:
            pos[n] = rng.uniform(-0.3, 0.3, size=2)

    # ── Figure ───────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(18, 16))

    # Hull fills
    for cid in major_ids:
        members = [n for n in nodes_list
                   if n in comm_map.index and int(comm_map.loc[n, "community"]) == cid]
        if len(members) < 3:
            continue
        pts = np.array([[pos[n][0], pos[n][1]] for n in members])
        try:
            hull = ConvexHull(pts)
            hp = pts[hull.vertices]; hp = np.vstack([hp, hp[0]])
            centroid = pts.mean(axis=0)
            hp_exp = centroid + (hp - centroid) * 1.18
            ax.fill(hp_exp[:, 0], hp_exp[:, 1], alpha=0.12,
                    color=COMMUNITY_COLORS[cid])
            ax.plot(hp_exp[:, 0], hp_exp[:, 1], alpha=0.38,
                    color=COMMUNITY_COLORS[cid], lw=1.6, linestyle="--")
        except Exception:
            pass

    # Community label above each hull
    comm_sizes = {cid: len([n for n in nodes_list if n in comm_map.index
                             and int(comm_map.loc[n, "community"]) == cid])
                  for cid in major_ids}
    for cid in major_ids:
        cx, cy = centroids[cid]
        short_label = COMMUNITY_LABELS[cid].split(" / ")[0]
        ax.text(cx, cy + 0.62 + 0.065 * comm_sizes.get(cid, 8),
                short_label, ha="center", va="bottom", fontsize=20,   # fontsize=9,
                fontweight="bold", color=COMMUNITY_COLORS[cid], alpha=0.9)

    # Edges
    nx.draw_networkx_edges(G, pos, alpha=0.45, edge_color="#aaaaaa", width=0.6, ax=ax)

    # Nodes by community (important ones drawn last)
    for cid in reversed(major_ids):
        members = [n for n in nodes_list
                   if n in comm_map.index and int(comm_map.loc[n, "community"]) == cid]
        if not members:
            continue
        szs = [node_sizes[nodes_list.index(n)] for n in members]
        nx.draw_networkx_nodes(G, pos, nodelist=members,
                               node_color=COMMUNITY_COLORS[cid],
                               node_size=szs, alpha=0.90, ax=ax)

    small_nodes = [n for n in nodes_list
                   if n not in comm_map.index or int(comm_map.loc[n, "community"]) > len(COMMUNITY_LABELS)]
    if small_nodes:
        szs = [node_sizes[nodes_list.index(n)] for n in small_nodes]
        nx.draw_networkx_nodes(G, pos, nodelist=small_nodes,
                               node_color="#cccccc", node_size=szs, alpha=0.75, ax=ax)

    # Labels — top 4 per community + SPP1
    top_genes = set(["SPP1"])
    for cid in major_ids:
        members_in_comm = communities_df[communities_df["community"] == cid]["gene"].tolist()
        top4 = (centrality[centrality["gene"].isin(members_in_comm)]
                .sort_values("composite_hub_score", ascending=False)
                .head(4)["gene"].tolist())
        top_genes.update(top4)
    labels = {n: n if n in top_genes else "" for n in nodes_list}

    # Shift labels above each node (same approach as fig_ppi_network)
    ys_comm = [pos[n][1] for n in nodes_list]
    pts_per_data_c = (max(ys_comm) - min(ys_comm)) / (16 * 72) + 1e-6
    size_map_c = dict(zip(nodes_list, node_sizes))
    label_pos_c = {
        n: (pos[n][0],
            pos[n][1] + np.sqrt(size_map_c[n] / np.pi) * pts_per_data_c + 0.1)
        for n in nodes_list
    }
    nx.draw_networkx_labels(G, label_pos_c, labels=labels, font_size=18,
                            font_weight="bold", ax=ax)

    # SPP1 bottleneck annotation
    if "SPP1" in pos:
        sx, sy = pos["SPP1"]
        ring = plt.Circle((sx, sy), 0.13, fill=False, edgecolor="black", lw=2.5)
        ax.add_patch(ring)
        ax.annotate("bottleneck\n(highest betweenness)", xy=(sx, sy),
                    xytext=(sx + 0.50, sy + 0.40), fontsize=14, style="italic",
                    color="black", weight="bold",
                    arrowprops=dict(arrowstyle="->", color="black", lw=1.2))

    # Axis limits
    xs = [pos[n][0] for n in nodes_list]; ys = [pos[n][1] for n in nodes_list]
    px = (max(xs) - min(xs)) * 0.10 + 0.30
    py = (max(ys) - min(ys)) * 0.10 + 0.30
    ax.set_xlim(min(xs) - px, max(xs) + px)
    ax.set_ylim(min(ys) - py, max(ys) + py)

    # Legend
    comm_patches = [mpatches.Patch(color=COMMUNITY_COLORS[i],
                                   label=f"C{i}: {COMMUNITY_LABELS[i]}")
                    for i in sorted(COMMUNITY_LABELS)]
    ax.legend(handles=comm_patches, title="Network module", loc="upper right",
              framealpha=0.88, fontsize=16, title_fontsize=18)

    # Modularity
    all_c = [set(communities_df[communities_df["community"] == cid]["gene"])
             for cid in sorted(communities_df["community"].unique())]
    try:
        mod_q = nx_comm.modularity(G, [c & set(G.nodes) for c in all_c], weight="weight")
        mod_str = f"  |  modularity Q = {mod_q:.3f}"
    except Exception:
        mod_str = ""

    ax.set_title(
        f"AT-specific UP gene PPI network — Louvain community structure\n"
        f"(n={G.number_of_nodes()} nodes, {G.number_of_edges()} edges, "
        f"STRING ≥ {STRING_SCORE}{mod_str})",
        fontsize=26, pad=14, weight="bold") # fontsize=12, pad=14)
    ax.axis("off")
    plt.tight_layout()
    savefig(fig, path=os.path.join(figures_path, "enrich_08b_ppi_communities.png"))
    print("  Saved enrich_08b_ppi_communities")


def fig_centrality_barplot(centrality: pd.DataFrame) -> None:
    """4-panel horizontal barplot: degree, betweenness, eigenvector, composite score.
    Bars coloured by composite hub score gradient (white→#d62728) — tier labels as text."""
    top = centrality[centrality["degree"] > 0].head(25).copy()
    top = top.sort_values("composite_hub_score")   # ascending so best gene is at top

    # Color gradient mapped to composite score (light→dark red)
    cmap = plt.cm.get_cmap("YlOrRd")
    score_norm = (top["composite_hub_score"] - top["composite_hub_score"].min()) / (
        top["composite_hub_score"].max() - top["composite_hub_score"].min() + 1e-9
    )
    bar_colors = [cmap(0.25 + v * 0.70) for v in score_norm.values]

    genes = top["gene"].values
    metrics = [
        ("degree",                "Degree (interactions)",  "#4e79a7"),
        ("betweenness_centrality","Betweenness centrality",  "#f28e2b"),
        ("eigenvector_centrality","Eigenvector centrality",  "#59a14f"),
        ("composite_hub_score",   "Composite hub score",    None),     # gradient
    ]

    fig, axes = plt.subplots(1, 4, figsize=(18, 9), sharey=True)

    for ax, (col, xlabel, color) in zip(axes, metrics):
        vals = top[col].values
        colors = bar_colors if color is None else [color] * len(genes)
        ax.barh(genes, vals, color=colors, edgecolor="white", linewidth=0.4, height=0.72)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.tick_params(axis="y", labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        # Annotate top 5 values
        for i, (g, v) in enumerate(zip(genes, vals)):
            if i >= len(genes) - 5:
                ax.text(v * 1.02, i, f"{v:.2f}" if col != "degree" else f"{int(v)}",
                        va="center", fontsize=7, color="#333333")

    axes[0].set_ylabel("Gene", fontsize=10)

    # Hub tier tick annotation on y-axis (rightmost panel)
    tier_col_map = {"High": "#d62728", "Medium": "#ff7f0e", "Low": "#aaaaaa"}
    for i, (_, row) in enumerate(top.iterrows()):
        t = row["hub_tier"]
        axes[0].get_yticklabels()[i].set_color(tier_col_map.get(t, "black"))

    # Colorbar for composite score panel
    sm = plt.cm.ScalarMappable(cmap=cmap,
                                norm=plt.Normalize(top["composite_hub_score"].min(),
                                                   top["composite_hub_score"].max()))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes[3], shrink=0.4, pad=0.02, aspect=20)
    cbar.set_label("Composite hub score", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    # Tier legend
    patches = [mpatches.Patch(color=c, label=t) for t, c in tier_col_map.items()]
    axes[0].legend(handles=patches, title="Hub tier (gene colour)",
                   loc="lower right", fontsize=7.5, framealpha=0.85)

    fig.suptitle("Top 25 AT-specific hub genes — centrality metrics", fontsize=12, y=1.01)
    plt.tight_layout()
    savefig(fig, path=os.path.join(figures_path, "enrich_09_ppi_centrality.png"))
    print("  Saved enrich_09_ppi_centrality")


def fig_druggability(centrality: pd.DataFrame, chembl: pd.DataFrame) -> None:
    """Stacked bar showing druggability tier for hub genes, sorted by hub score."""
    merged = centrality[centrality["degree"] > 0].head(TOP_HUBS).merge(
        chembl, on="gene", how="left"
    )
    merged["druggability_tier"] = merged["druggability_tier"].fillna("Unknown")
    merged = merged.sort_values("composite_hub_score", ascending=True)

    tier_order = ["Approved target", "Clinical candidate",
                  "Bioactive compounds", "Poorly drugged", "Unknown"]

    fig, ax = plt.subplots(figsize=(10, 11))
    y_pos = np.arange(len(merged))

    for tier in tier_order:
        mask = merged["druggability_tier"] == tier
        ax.barh(
            y_pos[mask.values],
            merged.loc[mask, "composite_hub_score"].values,
            color=DRUG_COLORS[tier],
            label=tier,
            edgecolor="white",
            linewidth=0.4,
            height=0.75,
        )

    # Annotate with activity count
    for i, (_, row) in enumerate(merged.iterrows()):
        n = int(row.get("activity_count", 0) or 0)
        phase = row.get("max_clinical_phase", 0) or 0
        label = f"  Ph{int(phase)}" if phase > 0 else (f"  n={n}" if n > 0 else "")
        ax.text(row["composite_hub_score"] + 0.002, i, label,
                va="center", ha="left", fontsize=7.5, color="#444444")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(merged["gene"].values, fontsize=9)
    ax.set_xlabel("Composite hub score", fontsize=10)
    ax.set_title(f"AT-specific hub gene druggability (top {TOP_HUBS} by hub score)",
                 fontsize=11, pad=10)
    ax.legend(title="Druggability tier", loc="lower right",
              fontsize=8, framealpha=0.8)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    savefig(fig, path=os.path.join(figures_path, "enrich_10_druggability.png"))
    print("  Saved enrich_10_druggability")

