"""
configuration_pathways.py — constants for pathway enrichment analysis (scripts 06+).

Kept separate from configuration.py so that scripts 02–05 do not pay the
import cost for constants they never use, and dependencies are explicit.
"""

# ── ORA / GSEA thresholds ─────────────────────────────────────────────────────

ORA_PADJ_THRESH  = 0.05    # BH-adjusted p-value cutoff for ORA significance
GSEA_NES_THRESH  = 1.5     # minimum |NES| for GSEA significance
GSEA_FDR_THRESH  = 0.25    # FDR q-value cutoff for GSEA significance

# ── gseapy run parameters ─────────────────────────────────────────────────────

GSEA_PERMS       = 1000    # permutation count for prerank
TOP_N_PLOT       = 12      # max terms per panel in dot plots / NES bars
MIN_GENESET_SIZE = 10
MAX_GENESET_SIZE = 500

# ── Enrichr databases ─────────────────────────────────────────────────────────

PATHWAY_DATABASES = [
    'GO_Biological_Process_2023',
    'KEGG_2021_Human',
    'MSigDB_Hallmark_2020',
    'TRRUST_Transcription_Factors_2019',
]

PATHWAY_DB_SHORT = {
    'GO_Biological_Process_2023':      'GO:BP',
    'KEGG_2021_Human':                 'KEGG',
    'MSigDB_Hallmark_2020':            'Hallmarks',
    'TRRUST_Transcription_Factors_2019': 'TRRUST',
}

# ── Gene set labels ───────────────────────────────────────────────────────────

# Human-readable label for each ORA gene set (used as heatmap column headers)
PATHWAY_LABEL_LONG = {
    'at_up_specific': 'AT-specific UP (acquired anaplasia, n=163)',
    'at_down':        'AT DOWN (FHWT identity lost, n=32)',
    'dawt_up':        'De novo DAWT UP (anaplastic program, n=94)',
    'fhwt_enriched':  'FHWT-enriched (favourable histology, n=173)',
    'secondary_res':  'Secondary resistance UP (Relapse, n=101)',
    'primary_res':    'Primary resistance UP (Progression, n=52)',
}

# Human-readable label for each GSEA step (used in figure titles)
PATHWAY_STEP_LABEL = {
    'step1': 'DAWT vs FHWT (Step 1)',
    'step2': 'FHWT Relapse vs Progression (Step 2)',
    'step3': 'AT Recurrent vs Primary (Step 3)',
}

# ── Colour palette for gene set labels ───────────────────────────────────────

PATHWAY_PALETTE = {
    'at_up_specific': '#d62728',   # red — AT acquired program
    'at_down':        '#17becf',   # cyan — suppressed in AT
    'dawt_up':        '#9467bd',   # purple — de novo anaplastic
    'fhwt_enriched':  '#2ca02c',   # green — FHWT identity
    'secondary_res':  '#1f77b4',   # blue — secondary resistance
    'primary_res':    '#ff7f0e',   # orange — primary resistance
}
