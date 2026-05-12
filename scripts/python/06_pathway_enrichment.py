#!/usr/bin/env python3
"""
06_pathway_enrichment.py — Wilms Tumor Analysis Project
========================================================
Gene set enrichment analysis for all three DE comparisons.

Requires internet access (Enrichr API). Run on your local machine.
On first run all results are cached to CSV; subsequent runs load from cache.

───────────────────────────────────────────────────────────────────────────────
ORA — gseapy.enrichr(), hypergeometric test, BH correction
    6 gene sets × 4 databases
    Background: 17,341 expressed protein-coding genes (study universe)

    Gene sets
    ─────────
    1. at_up_specific   AT recurrent UP ∩ NOT de novo DAWT (163 genes)
    2. at_down          Genes suppressed during AT (32 genes)
    3. dawt_up          De novo anaplastic program (94 genes)
    4. fhwt_enriched    Favourable histology program (173 genes, DAWT DOWN)
    5. secondary_res    FHWT Relapse > Progression primary tumors (101 genes)
    6. primary_res      FHWT Progression > Relapse primary tumors (52 genes)

GSEA — gseapy.prerank(), Kolmogorov–Smirnov, gene permutation
    3 ranked gene lists × 4 databases
    Ranking metric: sign(LFC) × −log10(padj)   [NA padj → 0, sinks to middle]

    Ranked lists
    ────────────
    Step 1: DAWT vs FHWT   (positive = DAWT-enriched, de novo anaplastic)
    Step 2: Relapse vs Progression, FHWT female-only
            (positive = secondary resistance; negative = primary resistance)
    Step 3: AT Recurrent vs Primary, paired
            (positive = acquired anaplastic program)

Databases (Enrichr names)
    GO_Biological_Process_2023
    KEGG_2021_Human
    MSigDB_Hallmark_2020
    TRRUST_Transcription_Factors_2019

Significance thresholds
    ORA  : Adjusted P-value (BH) < 0.05
    GSEA : |NES| ≥ 1.5  AND  FDR q-value < 0.25

Inputs
──────
data/processed/tables/de_histology_full.csv
data/processed/tables/de_histology_sig.csv
data/processed/tables/de_relapse_fhwt_full.csv
data/processed/tables/de_relapse_fhwt_sig.csv
data/processed/tables/de_at_full.csv
data/processed/tables/de_at_sig.csv
data/processed/tables/de_at_vs_dawt_concordant.csv

Outputs
───────
data/processed/tables/enrichment/ora_{label}_{database}.csv
data/processed/tables/enrichment/gsea_step{n}_{database}.csv
data/processed/tables/enrichment/summary_ora_significant.csv
data/processed/tables/enrichment/summary_gsea_significant.csv
results/figures/enrich_01_ora_at_dotplot.pdf/.png
results/figures/enrich_02_ora_dawt_fhwt_dotplot.pdf/.png
results/figures/enrich_03_ora_resistance_dotplot.pdf/.png
results/figures/enrich_04_gsea_step1_nes.pdf/.png
results/figures/enrich_05_gsea_step2_nes.pdf/.png
results/figures/enrich_06_gsea_step3_nes.pdf/.png
results/figures/enrich_07_hallmarks_heatmap.pdf/.png
results/figures/enrich_08_kegg_heatmap.pdf/.png
"""

from __future__ import annotations
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import gseapy as gp

from configuration import PATHS
from utils.configuration_pathways import (PATHWAY_DATABASES, PATHWAY_DB_SHORT, PATHWAY_STEP_LABEL,
                                          ORA_PADJ_THRESH, GSEA_NES_THRESH, GSEA_FDR_THRESH,
                                          GSEA_PERMS, MIN_GENESET_SIZE, MAX_GENESET_SIZE)
from utils.enrichment_plot_utils import plot_ora_dotplot, plot_gsea_nes, plot_ora_heatmap

warnings.filterwarnings('ignore')

processed_tables_path, matrices_path, enrich_tables_path, figures_path = \
    (PATHS["processed_tables_path"],
     PATHS["matrices_path"],
     PATHS["enrich_tables_path"],
     PATHS["figures_path"])

# ── Connectivity check ────────────────────────────────────────────────────────

def check_enrichr_access() -> bool:
    """
    Probe the Enrichr API and return True if reachable.

    Prints a warning and returns False if the request times out or fails.
    """
    import urllib.request
    try:
        urllib.request.urlopen(
            'https://maayanlab.cloud/Enrichr/datasetStatistics', timeout=6)
        return True
    except Exception as e:
        print(f"\n  [!] Enrichr API not reachable: {e}")
        print("      Run this script on a machine with internet access.\n")
        return False


# ── Data loading ──────────────────────────────────────────────────────────────

def load_de_results() -> dict[str, pd.DataFrame]:
    """
    Load all three DE full result tables from disk.

    Returns
    -------
    dict with keys 'step1' (histology), 'step2' (relapse), 'step3' (AT),
    each mapping to the corresponding full DE results DataFrame.
    """
    return {
        'step1': pd.read_csv(os.path.join(processed_tables_path, 'de_histology_full.csv')),
        'step2': pd.read_csv(os.path.join(processed_tables_path, 'de_relapse_fhwt_full.csv')),
        'step3': pd.read_csv(os.path.join(processed_tables_path, 'de_at_full.csv')),
    }


def build_gene_sets(de: dict[str, pd.DataFrame]) -> dict[str, list[str]]:
    """
    Derive all 6 ORA gene sets from DE results.

    AT-specific UP excludes genes in the DAWT concordant set so the set captures
    what is unique to acquired anaplasia rather than the shared anaplastic module.

    Parameters
    ----------
    de : dict mapping step labels to full DE result DataFrames
         (must include 'step1', 'step2', 'step3' with a 'direction' column)

    Returns
    -------
    dict mapping gene_set_label → sorted list of gene names;
    keys: at_up_specific, at_down, dawt_up, fhwt_enriched, secondary_res, primary_res
    """
    # DAWT concordant genes (shared between AT UP and DAWT UP)
    dawt_conc_path = os.path.join(processed_tables_path, 'de_at_vs_dawt_concordant.csv')
    dawt_conc_genes = set()
    if os.path.exists(dawt_conc_path):
        dawt_conc_genes = set(
            pd.read_csv(dawt_conc_path)['gene_name'].dropna())

    sig1 = de['step1']
    sig2 = de['step2']
    sig3 = de['step3']

    at_up   = set(sig3.loc[sig3['direction'] == 'UP', 'gene_name'])
    at_down = set(sig3.loc[sig3['direction'] == 'DOWN', 'gene_name'])

    gene_sets = {
        'at_up_specific': sorted(at_up - dawt_conc_genes),
        'at_down':        sorted(at_down),
        'dawt_up':        sorted(sig1.loc[sig1['direction'] == 'UP',   'gene_name']),
        'fhwt_enriched':  sorted(sig1.loc[sig1['direction'] == 'DOWN', 'gene_name']),
        'secondary_res':  sorted(sig2.loc[sig2['direction'] == 'UP',   'gene_name']),
        'primary_res':    sorted(sig2.loc[sig2['direction'] == 'DOWN', 'gene_name']),
    }

    print("\n  Gene set sizes:")
    for label, genes in gene_sets.items():
        print(f"    {label:<20} n = {len(genes)}")
    return gene_sets


def build_ranked_lists(de: dict[str, pd.DataFrame]) -> dict[str, pd.Series]:
    """
    Build per-step ranked gene lists for GSEA prerank.

    Ranking metric: sign(LFC) × −log10(padj) + LFC × 1e-4 (tie-breaker).

    The primary component is a padj-based directional score.  The LFC × 1e-4
    tie-breaker separates NA-padj genes (DESeq2 outlier / low-count exclusions,
    which score 0 on the primary metric) without meaningfully shifting significant
    genes (max tie-breaker ≈ 0.002 vs primary > 1.3 at padj < 0.05).

    Parameters
    ----------
    de : dict mapping step labels to full DE result DataFrames
         (must include 'gene_name', 'log2FoldChange', 'padj' columns)

    Returns
    -------
    dict mapping step label → pd.Series (score, index=gene_name),
    sorted descending, duplicates removed (keep highest score).
    """
    ranked = {}
    for step, df in de.items():
        lfc  = df['log2FoldChange'].fillna(0)
        padj = df['padj'].fillna(1.0).clip(lower=1e-300)
        score = np.sign(lfc) * -np.log10(padj) + lfc * 1e-4
        s = pd.Series(score.values, index=df['gene_name'])
        s = s.sort_values(ascending=False)
        # Remove duplicates (keep first = highest score)
        s = s[~s.index.duplicated(keep='first')]
        ranked[step] = s
        print(f"    {step}: {len(s):,} genes ranked  "
              f"[max={s.max():.2f}, min={s.min():.2f}]")
    return ranked


def get_background(de: dict[str, pd.DataFrame]) -> list[str]:
    """
    Return the full gene universe for ORA (all expressed protein-coding genes).

    Uses the step1 gene list (17,341 genes) as the hypergeometric background.

    Parameters
    ----------
    de : dict mapping step labels to full DE result DataFrames

    Returns
    -------
    List of gene name strings (no NAs).
    """
    bg = de['step1']['gene_name'].dropna().tolist()
    print(f"  ORA background: {len(bg):,} genes")
    return bg


# ── ORA ───────────────────────────────────────────────────────────────────────

def _normalise_ora_columns(res: pd.DataFrame) -> pd.DataFrame:
    """
    Normalise ORA result column names to our internal convention regardless of
    gseapy version.  Works by lowercasing + collapsing spaces/hyphens first,
    then applying specific renames.
    """
    # Step 1: lowercase everything, replace spaces and hyphens with underscores
    res.columns = [c.strip().lower().replace(' ', '_').replace('-', '_')
                   for c in res.columns]
    # Step 2: map to our preferred names
    renames = {
        'adjusted_p_value': 'adj_pval',
        'p_value':          'pval',
        'gene_set':         'database',
        # 'term', 'overlap', 'genes', 'odds_ratio', 'combined_score' are already fine
    }
    return res.rename(columns=renames)


# All columns required by downstream figures; stale caches missing any of these
# are regenerated rather than silently producing broken plots.
_ORA_REQUIRED = {'gene_set_label', 'term', 'adj_pval', 'odds_ratio', 'database_short', 'genes'}


def run_ora_one(label: str, genes: list[str], background: list[str],
                databases: list[str]) -> pd.DataFrame:
    """
    Run ORA for one gene set against all databases.
    Returns combined DataFrame tagged with gene_set label.
    Uses CSV cache: skips API call if file is valid (has required columns).
    Invalid/stale caches (e.g. from a previous failed run) are regenerated.
    """
    cache_path = os.path.join(enrich_tables_path, f'ora_{label}.csv')
    if os.path.exists(cache_path):
        cached = pd.read_csv(cache_path)
        if _ORA_REQUIRED.issubset(set(cached.columns)):
            print(f"    [cache] {label}")
            return cached
        print(f"    [cache stale] {label} — missing columns, regenerating")
        os.remove(cache_path)

    if not genes:
        print(f"    [skip] {label}: empty gene set")
        return pd.DataFrame()

    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=databases,
        background=background,
        organism='human',
        outdir=None,
        verbose=False,
    )
    res = enr.results.copy()
    res['gene_set_label'] = label
    res = _normalise_ora_columns(res)
    res['database_short'] = res['database'].map(PATHWAY_DB_SHORT).fillna(res['database'])
    res.to_csv(cache_path, index=False)
    n_sig = (res['adj_pval'] < ORA_PADJ_THRESH).sum()
    print(f"    Saved: ora_{label}.csv  ({n_sig} significant)")
    return res


def run_all_ora(
    gene_sets: dict[str, list[str]],
    background: list[str],
) -> pd.DataFrame:
    """
    Run ORA for all 6 gene sets across all databases and concatenate results.

    Parameters
    ----------
    gene_sets  : dict mapping gene_set_label → list of gene names
    background : full gene universe for the hypergeometric test

    Returns
    -------
    Combined ORA results DataFrame with columns: gene_set_label, term, adj_pval,
    odds_ratio, genes, database, database_short.
    """
    print("\n── ORA: Running all gene sets ───────────────────────────────")
    frames = []
    for label, genes in gene_sets.items():
        df = run_ora_one(label, genes, background, PATHWAY_DATABASES)
        if not df.empty:
            # database_short added by run_ora_one; guard covers very old caches
            if 'database_short' not in df.columns:
                df['database_short'] = df['database'].map(PATHWAY_DB_SHORT).fillna(df['database'])
            frames.append(df)
    ora_all = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    sig = ora_all[ora_all['adj_pval'] < ORA_PADJ_THRESH]
    print(f"\n  Total ORA results  : {len(ora_all):,}")
    print(f"  Significant (FDR<{ORA_PADJ_THRESH}): {len(sig):,}")

    summary_path = os.path.join(enrich_tables_path, 'summary_ora_significant.csv')
    sig.sort_values(['gene_set_label', 'database', 'adj_pval'])\
       .to_csv(summary_path, index=False)
    print(f"  Saved: summary_ora_significant.csv")
    return ora_all


# ── GSEA ──────────────────────────────────────────────────────────────────────

def run_gsea_one(step: str, ranking: pd.Series, database: str) -> pd.DataFrame:
    """
    Run GSEA prerank for one step × one database.
    Returns results DataFrame tagged with step and database labels.
    Uses CSV cache.
    """
    db_short  = PATHWAY_DB_SHORT.get(database, database)
    cache_path = os.path.join(enrich_tables_path, f'gsea_{step}_{db_short}.csv')
    if os.path.exists(cache_path):
        print(f"    [cache] {step} × {db_short}")
        df = pd.read_csv(cache_path)
        return df

    pre = gp.prerank(
        rnk=ranking,
        gene_sets=database,
        outdir=None,
        min_size=MIN_GENESET_SIZE,
        max_size=MAX_GENESET_SIZE,
        permutation_num=GSEA_PERMS,
        seed=42,
        verbose=False,
    )
    # gseapy ≥1.0: .results is an internal dict; .res2d is the tidy DataFrame
    res = pre.res2d.reset_index().copy()
    # Normalise column names across gseapy versions
    col_map = {
        'Term':        'term',
        'ES':          'ES',
        'NES':         'NES',
        'NOM p-val':   'pval',
        'FDR q-val':   'fdr',
        'FWER p-val':  'fwer',
        'Tag %':       'tag_pct',
        'Gene %':      'gene_pct',
        'Signal %':    'signal_pct',
        'Lead_genes':  'lead_genes',
    }
    res = res.rename(columns={k: v for k, v in col_map.items() if k in res.columns})
    res['step']           = step
    res['database']       = database
    res['database_short'] = db_short
    res['step_label']     = PATHWAY_STEP_LABEL[step]

    res.to_csv(cache_path, index=False)
    sig = res[(res['NES'].abs() >= GSEA_NES_THRESH) & (res['fdr'] < GSEA_FDR_THRESH)]
    print(f"    {step} × {db_short}: {len(sig)} significant  (|NES|≥{GSEA_NES_THRESH}, FDR<{GSEA_FDR_THRESH})")
    return res


def run_all_gsea(ranked: dict[str, pd.Series]) -> pd.DataFrame:
    """
    Run GSEA prerank for all steps × all databases and concatenate results.

    Parameters
    ----------
    ranked : dict mapping step label → ranked pd.Series (score, index=gene_name)

    Returns
    -------
    Combined GSEA results DataFrame with columns: step, database_short, term,
    NES, fdr, ES, pval, and additional gseapy columns.
    """
    print("\n── GSEA prerank: Running all combinations ───────────────────")
    frames = []
    for step, ranking in ranked.items():
        for database in PATHWAY_DATABASES:
            frames.append(run_gsea_one(step, ranking, database))

    gsea_all = pd.concat([f for f in frames if not f.empty], ignore_index=True)
    sig = gsea_all[
        (gsea_all['NES'].abs() >= GSEA_NES_THRESH) &
        (gsea_all['fdr'] < GSEA_FDR_THRESH)
    ]
    summary_path = os.path.join(enrich_tables_path, 'summary_gsea_significant.csv')
    sig.sort_values(['step', 'database', 'NES'], ascending=[True, True, False])\
       .to_csv(summary_path, index=False)
    print(f"\n  Total GSEA results  : {len(gsea_all):,}")
    print(f"  Significant         : {len(sig):,}")
    print(f"  Saved: summary_gsea_significant.csv")
    return gsea_all


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:  # noqa: D103
    print("\n" + "=" * 60)
    print("  06_pathway_enrichment.py")
    print("=" * 60)

    # ── Connectivity check ────────────────────────────────────────────────────
    print("\n── Checking Enrichr API connectivity ────────────────────────")
    online = check_enrichr_access()
    if not online:
        # If all results cached, we can still run figures
        cached = [f for f in os.listdir(enrich_tables_path) if f.endswith('.csv')]
        if not cached:
            print("  No cached results found. Exiting.")
            return
        print(f"  {len(cached)} cached result files found — proceeding with figures only.")

    # ── Load DE results ───────────────────────────────────────────────────────
    print("\n── Load DE results ──────────────────────────────────────────")
    de = load_de_results()
    for step, df in de.items():
        print(f"  {step}: {len(df):,} genes  "
              f"(UP={( df['direction']=='UP').sum()}, "
              f"DOWN={(df['direction']=='DOWN').sum()})")

    # ── Build gene sets and ranking ───────────────────────────────────────────
    print("\n── Build gene sets (ORA) ────────────────────────────────────")
    gene_sets  = build_gene_sets(de)
    background = get_background(de)

    print("\n── Build ranked lists (GSEA) ────────────────────────────────")
    ranked = build_ranked_lists(de)

    # ── ORA ───────────────────────────────────────────────────────────────────
    if online or any(f.startswith('ora_') for f in os.listdir(enrich_tables_path)):
        ora_all = run_all_ora(gene_sets, background)
    else:
        ora_all = pd.DataFrame()

    # ── GSEA ──────────────────────────────────────────────────────────────────
    if online or any(f.startswith('gsea_') for f in os.listdir(enrich_tables_path)):
        gsea_all = run_all_gsea(ranked)
    else:
        gsea_all = pd.DataFrame()

    # ── Figures: ORA dot plots ─────────────────────────────────────────────
    if not ora_all.empty:
        print("\n── ORA figures ──────────────────────────────────────────────")
        print(f"  ora_all columns : {list(ora_all.columns)}")
        print(f"  ora_all shape   : {ora_all.shape}")

        # Panel 1: AT-specific program (AT UP and AT DOWN)
        plot_ora_dotplot(
            ora_all,
            gene_set_labels=['at_up_specific', 'at_down'],
            out_path=os.path.join(figures_path, 'enrich_01_ora_at_dotplot'),
            panel_title='ORA — Anaplastic Transformation gene sets\n'
                        '(AT-specific UP: acquired program; AT DOWN: lost identity)',
        )
        print("  Saved: enrich_01_ora_at_dotplot")

        # Panel 2: DAWT UP and FHWT-enriched (histological identity)
        plot_ora_dotplot(
            ora_all,
            gene_set_labels=['dawt_up', 'fhwt_enriched'],
            out_path=os.path.join(figures_path, 'enrich_02_ora_dawt_fhwt_dotplot'),
            panel_title='ORA — Histological identity gene sets\n'
                        '(DAWT UP: de novo anaplastic; FHWT-enriched: favourable histology)',
        )
        print("  Saved: enrich_02_ora_dawt_fhwt_dotplot")

        # Panel 3: Resistance gene sets
        plot_ora_dotplot(
            ora_all,
            gene_set_labels=['secondary_res', 'primary_res'],
            out_path=os.path.join(figures_path, 'enrich_03_ora_resistance_dotplot'),
            panel_title='ORA — Primary vs secondary resistance gene sets\n'
                        '(Secondary: Relapse UP; Primary: Progression UP)',
        )
        print("  Saved: enrich_03_ora_resistance_dotplot")

        # Summary heatmaps across all 6 gene sets — Hallmarks + KEGG
        # Hallmarks: best for oncogenic activation programs (UP-biased)
        # KEGG: broader pathway coverage, better for suppressed/DOWN sets
        for db_short, fig_idx in [('Hallmarks', '07'), ('KEGG', '08')]:
            fname = f'enrich_{fig_idx}_{db_short.lower()}_heatmap'
            plot_ora_heatmap(
                ora_all,
                out_path=os.path.join(figures_path, fname),
                database_short=db_short,
            )
            print(f"  Saved: {fname}")

    # ── Figures: GSEA NES bar charts ──────────────────────────────────────────
    if not gsea_all.empty:
        print("\n── GSEA figures ─────────────────────────────────────────────")
        for step, fig_idx in [('step1', '04'), ('step2', '05'), ('step3', '06')]:
            fname = f'enrich_{fig_idx}_gsea_{step}_nes'
            plot_gsea_nes(
                gsea_all, step,
                out_path=os.path.join(figures_path, fname),
            )
            print(f"  Saved: {fname}")

    # ── Print top results to console ──────────────────────────────────────────
    print("\n── Top ORA results (Hallmarks) ──────────────────────────────")
    if not ora_all.empty:
        top_hm = ora_all[
            (ora_all['database_short'] == 'Hallmarks') &
            (ora_all['adj_pval'] < ORA_PADJ_THRESH)
        ].sort_values('adj_pval')[
            ['gene_set_label', 'term', 'adj_pval', 'odds_ratio']
        ].head(30)
        if not top_hm.empty:
            print(top_hm.to_string(index=False))
        else:
            print("  No significant Hallmarks terms")

    print("\n── Top GSEA results (all databases) ─────────────────────────")
    if not gsea_all.empty:
        top_gsea = gsea_all[
            (gsea_all['NES'].abs() >= GSEA_NES_THRESH) &
            (gsea_all['fdr'] < GSEA_FDR_THRESH)
        ].sort_values('NES', ascending=False)[
            ['step', 'database_short', 'term', 'NES', 'fdr']
        ].head(30)
        if not top_gsea.empty:
            print(top_gsea.to_string(index=False))
        else:
            print("  No significant GSEA terms")

    print("  Outputs in data/processed/tables/enrichment/")
    print("\n── Complete ─────────────────────────────────────────────────")


if __name__ == '__main__':
    main()
