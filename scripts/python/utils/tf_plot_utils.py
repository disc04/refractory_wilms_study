import os
import warnings
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from configuration import FOCUS_TF_COLORS, FDR_THRESH
from utils.data_utils import savefig
warnings.filterwarnings("ignore")

def fig_upset_nfkb(membership: pd.DataFrame,
                   de_stats: pd.DataFrame,
                   focus_tfs: list[str],
                   out_path: str,) -> None:
    """
    UpSet-style bar chart showing gene counts per TF combination.
    Each bar is annotated with representative gene names.
    """
    # Build combination labels
    combos = {}
    for _, row in membership.iterrows():
        key = tuple(int(row[tf]) for tf in focus_tfs)
        if any(key):
            combos.setdefault(key, []).append(row["gene"])

    # Sort by descending count
    sorted_combos = sorted(combos.items(), key=lambda x: len(x[1]), reverse=True)

    fig, ax = plt.subplots(figsize=(14, 6))
    bar_colors = []
    labels = []
    counts = []

    for combo_key, genes in sorted_combos:
        active = [focus_tfs[i] for i, v in enumerate(combo_key) if v]
        label = "+".join(active)
        # Color by dominant TF (first active)
        color = FOCUS_TF_COLORS.get(active[0], "#888888") if active else "#888888"
        labels.append(label)
        counts.append(len(genes))
        bar_colors.append(color)

    bars = ax.bar(range(len(labels)), counts, color=bar_colors,
                  edgecolor="white", linewidth=0.5, width=0.65)

    # Annotate top-LFC genes per bar
    for i, (combo_key, genes) in enumerate(sorted_combos):
        in_de = [g for g in genes if g in de_stats.index]
        if in_de:
            top = sorted(in_de,
                         key=lambda g: de_stats.loc[g, "log2FoldChange"],
                         reverse=True)[:3]
            ax.text(i, counts[i] + 0.15, "\n".join(top),
                    ha="center", va="bottom", fontsize=6.5, color="#333333")

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=9)
    ax.set_ylabel("Number of AT-specific genes", fontsize=10)
    ax.set_title("Direct TF regulon membership of AT-specific UP genes\n"
                 "(RELA / NFKB1 / JUN / STAT3  ×  TRRUST v2)",
                 fontsize=11, pad=10)
    ax.spines[["top", "right"]].set_visible(False)

    # Legend
    patches = [mpatches.Patch(color=c, label=tf)
               for tf, c in FOCUS_TF_COLORS.items()]
    ax.legend(handles=patches, title="TF (dominant)", loc="upper right",
              fontsize=8, framealpha=0.85)

    plt.tight_layout()
    savefig(fig, os.path.join(out_path, "enrich_11_tf_regulon_upset.png"))
    print("  Saved enrich_11_tf_regulon_upset")


def fig_heatmap(membership: pd.DataFrame,
                de_stats: pd.DataFrame,
                overlap_df: pd.DataFrame,
                out_path: str,
                top_n_tfs: int = 20) -> None:
    """
    Gene × TF heatmap for the top significant TFs.
    Rows = AT-specific genes present in at least one top-TF regulon.
    Columns = top significant TFs by FDR.
    Rows ordered by decreasing sum of TF memberships, then by LFC.
    """
    # Top TFs
    top_tfs = overlap_df[overlap_df["fdr"] < FDR_THRESH].head(top_n_tfs)["tf"].tolist()
    if not top_tfs:
        print("  No significant TFs for heatmap — skipping")
        return

    # Build matrix: rows = genes in any top TF regulon, cols = top TFs
    regulon_dict = {}
    for _, row in overlap_df[overlap_df["tf"].isin(top_tfs)].iterrows():
        regulon_dict[row["tf"]] = set(row["overlap_genes"].split(";"))

    all_genes = sorted(set.union(*regulon_dict.values()))
    mat = pd.DataFrame(0, index=all_genes, columns=top_tfs)
    for tf, gset in regulon_dict.items():
        for g in gset:
            if g in mat.index:
                mat.loc[g, tf] = 1

    # Order rows: total membership desc, then LFC desc
    mat["_total"] = mat.sum(axis=1)
    mat["_lfc"] = [de_stats.loc[g, "log2FoldChange"] if g in de_stats.index else 0
                   for g in mat.index]
    mat = mat.sort_values(["_total", "_lfc"], ascending=[False, False])
    mat = mat.drop(columns=["_total", "_lfc"])

    # LFC annotation
    lfc_vals = [de_stats.loc[g, "log2FoldChange"] if g in de_stats.index else 0
                for g in mat.index]

    fig, axes = plt.subplots(1, 2, figsize=(max(10, len(top_tfs) * 0.55), max(8, len(all_genes) * 0.32)),
                             gridspec_kw={"width_ratios": [1, 6]})

    # LFC barplot
    ax_lfc = axes[0]
    colors = ["#d62728" if v > 0 else "#1f77b4" for v in lfc_vals]
    ax_lfc.barh(range(len(lfc_vals)), lfc_vals, color=colors,
                edgecolor="white", linewidth=0.3, height=0.8)
    ax_lfc.set_yticks(range(len(mat.index)))
    ax_lfc.set_yticklabels(mat.index, fontsize=7.5)
    ax_lfc.set_xlabel("log₂FC (AT recurrent vs primary)", fontsize=8)
    ax_lfc.axvline(0, color="black", linewidth=0.7)
    ax_lfc.spines[["top", "right"]].set_visible(False)
    ax_lfc.invert_yaxis()

    # TF membership heatmap
    ax_heat = axes[1]
    sns.heatmap(mat, ax=ax_heat, cmap="Reds", vmin=0, vmax=1,
                linewidths=0.3, linecolor="#dddddd",
                cbar=False, yticklabels=False, xticklabels=True)
    ax_heat.set_xticklabels(ax_heat.get_xticklabels(),
                             rotation=45, ha="right", fontsize=8)
    ax_heat.set_ylabel("")

    # Highlight FOCUS_TFS column labels
    for tick in ax_heat.get_xticklabels():
        if tick.get_text() in FOCUS_TF_COLORS:
            tick.set_color(FOCUS_TF_COLORS[tick.get_text()])
            tick.set_fontweight("bold")

    fig.suptitle(f"AT-specific genes in top {top_n_tfs} significant TF regulons (TRRUST v2)",
                 fontsize=11, y=1.01)
    plt.tight_layout()
    savefig(fig, os.path.join(out_path, "enrich_12_tf_regulon_heatmap.png"))
    print("  Saved enrich_12_tf_regulon_heatmap")

