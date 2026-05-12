"""Central config for the project. Handles:
   - Project-root-relative paths
   -- Environment variable / secrets loading from .env
   - Constants
   - Patients
   - Signatures
   - Labels and colors
 Usage:
    from configuration import PATHS, ENV
"""

from pathlib import Path
from dotenv import load_dotenv

# Project metadata
PROJECT_NAME = "Refractory Wilms Tumor Study"
VERSION      = "0.1.0"

# ─────────────────────────── paths ────────────────────────────────────────────
# Project root
# _CONF_DIR    = os.path.dirname(os.path.abspath(__file__))
# ROOT = os.path.dirname(os.path.dirname(os.path.dirname(_CONF_DIR)))
ROOT = Path(__file__).resolve().parent


# Paths
PATHS = {"root": ROOT,
         "clinical_data_path":  ROOT / "data" / "raw" / "GDCdata" / "TARGET-WT" / "Clinical" / "Clinical_Supplement",
         "rnaseq_data_path": ROOT / "data" / "raw" / "GDCdata" / "TARGET-WT"/
                                "Transcriptome_Profiling" / "Gene_Expression_Quantification",
        "matrices_path":  ROOT / "data" / "processed"/ "matrices",
        "processed_tables_path":  ROOT / "data" / "processed"/ "tables",
        "enrich_tables_path": ROOT / "data" / "processed"/ "tables" / "enrichment",
        "figures_path":  ROOT / "data" / "results" / "figures"}

CACHE_PATHS = {"CACHE_TRRUST": ROOT / "data" / "processed" / "tables" / "enrichment" / "trrust_rawdata_human.tsv",
                "CACHE_PPI": ROOT / "data" / "processed" / "tables" / "enrichment" / "ppi_string_edges.csv",
                "CACHE_CENTRALITY": ROOT / "data" / "processed" / "tables" / "enrichment" / "ppi_centrality.csv",
                "CACHE_CHEMBL": ROOT / "data" / "processed" / "tables" / "enrichment" / "ppi_chembl_hits.csv",
                "CACHE_COMMUNITIES": ROOT / "data" / "processed" / "tables" / "enrichment" / "ppi_communities.csv"
               }

# Create directories if they don't exist
for _dir in PATHS.values():
    if isinstance(_dir, Path) and _dir != ROOT:
        _dir.mkdir(parents=True, exist_ok=True)

# Load .env (secrets, API keys — never committed to git)
_env_path = ROOT / ".env"
load_dotenv(dotenv_path=_env_path, override=False)  # override=False: env vars set


# Environment variables / secrets
ENV = {
    # "API_KEY":    None  # os.getenv("API_KEY"),
    # "DB_URL":     None  # os.getenv("DB_URL"),
}

# ─────────────────────────── constants ────────────────────────────────────────
STRING_API = "https://string-db.org/api/json"
STRING_SCORE = 700  # high-confidence threshold (0–1000)
TAXON_HUMAN = 9606
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
TOP_HUBS = 30  # genes passed to ChEMBL
HUB_TIER_HIGH = 0.70  # composite score quantile for "high hub"
HUB_TIER_MED = 0.40
TRRUST_URL   = "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"
# Minimum overlap to report in all-TF analysis
MIN_OVERLAP  = 3
FDR_THRESH   = 0.05

# ─────────────────────────── patients ────────────────────────────────────────────

# excluded from histology differential expression analysis - anaplastic transformation and metastatic patients
EXCLUDE_USI = [
    'TARGET-50-PAJPDC', 'TARGET-50-PAJNGH', 'TARGET-50-PALFME',
    'TARGET-50-PAJNTJ', 'TARGET-50-PALJIP',   # anaplastic transformation
    'TARGET-50-PAJLUJ',                         # metastatic
]

ANAPLASTIC_TRANSFORMATION_PATIENTS = ['TARGET-50-PAJPDC',
                                      'TARGET-50-PAJNGH',
                                      'TARGET-50-PALFME',
                                      'TARGET-50-PAJNTJ',
                                      # TARGET-50-PALJIP excluded: pathology_subtype="relapse Wilms", not anaplastic Wilms
                                      ]

# ─────────────────────────── signatures ────────────────────────────────────────────

WILMS_PANEL = ['WT1', 'CTNNB1', 'WTX', 'TP53', 'SIX1', 'SIX2',
               'DROSHA', 'DICER1', 'IGF2', 'MYCN', 'DNMT3A']

KEY_DRIVERS = ['TP53', 'CTNNB1', 'WT1']

SIGNATURES: dict[str, list[str]] = {
    'P53_TARGETS':            ['CDKN1A','BAX','MDM2','GADD45A','BBC3','PMAIP1',
                               'TP53I3','SERPINB5','RRM2B','ZMAT3','SESN1',
                               'SESN2','TIGAR','FAS','TNFRSF10B'],
    'E2F_TARGETS':            ['MCM2','MCM3','MCM4','MCM5','MCM6','MCM7','PCNA',
                               'RFC4','CDC6','CDC25A','CCNE1','CCNE2','CDK1',
                               'CDK2','E2F1','E2F2','TOP2A','TYMS','RRM2','TK1',
                               'POLE','POLE2','UBE2C','DUT'],
    'MYC_TARGETS':            ['MYC','MYCN','MAX','HSPD1','NPM1','NOP56','NOP14',
                               'NOP10','PA2G4','LDHA','ENO1','PKM','SHMT2','PRDX3',
                               'TFDP1','ODC1','LDHB','EIF4E','EIF4G1','EIF4A1',
                               'PABPC1','TUFM'],
    'OXPHOS':                 ['ATP5F1A','ATP5F1B','ATP5MC1','ATP5MC3','COX4I1',
                               'COX5A','COX5B','COX6A1','COX6B1','NDUFA1','NDUFA5',
                               'NDUFA6','NDUFA8','NDUFA9','NDUFB7','NDUFB8','NDUFS1',
                               'NDUFS2','NDUFS7','SDHA','SDHB','UQCR10','UQCR11',
                               'UQCRC1','UQCRC2'],
    'EMT':                    ['SNAI1','SNAI2','ZEB1','ZEB2','TWIST1','VIM','FN1',
                               'CDH2','MMP2','MMP9','S100A4','TGFB1','TGFB2',
                               'TGFBR1','SERPINE1'],
    'EPITHELIAL':             ['CDH1','EPCAM','CLDN1','CLDN3','OCLN','TJP1','KRT8',
                               'KRT18','DSP','PKP3'],
    'BLASTEMAL_WT':           ['SIX1','SIX2','EYA1','LIN28B','DLK1','CITED1','MEOX1',
                               'PAX2','SALL1','MYCN','IGF2','WT1'],
    'IFNg_RESPONSE':          ['IFNGR1','IFNGR2','STAT1','IRF1','GBP1','GBP2','GBP4',
                               'GBP5','CXCL9','CXCL10','CXCL11','HLA-A','HLA-B',
                               'HLA-C','HLA-DRA','HLA-DRB1','TAP1','TAP2','PSMB8',
                               'PSMB9','MX1','OAS1','ISG15'],
    'CD8_TCELL':              ['CD8A','CD8B','GZMA','GZMB','GZMH','GZMK','PRF1',
                               'NKG7','CD3D','CD3E','CD3G','IFNG','TBX21'],
    'CHEMORESISTANCE':        ['ABCB1','ABCC1','ABCG2','GSTP1','HMOX1','NQO1','NRF2',
                               'TXN','TXNRD1','MGMT'],
    'P53_PATHWAY_HALL_subset':['CDKN1A','MDM2','BAX','GADD45A','PMAIP1','BBC3',
                               'TP53I3','RRM2B','ZMAT3','SESN1','SESN2','TIGAR',
                               'BTG2','CDKN1B','PML','SFN','TNFRSF10B'],
}

# TFs to highlight in focused analysis
FOCUS_TFS    = ["RELA", "NFKB1", "JUN", "STAT3"]

# Embedded fallback list of well-known X/Y/MT genes (used if no API is reachable)
EMBEDDED_XY_GENES = {
    # Y-chromosome
    'RPS4Y1','RPS4Y2','DDX3Y','USP9Y','EIF1AY','KDM5D','NLGN4Y','TXLNGY',
    'UTY','ZFY','PRKY','TBL1Y','TTTY14','TTTY15','TTTY10','SRY','TSPY1',
    'TSPY2','AMELY','TMSB4Y','PCDH11Y','BCORP1','ANOS2P','TGIF2LY','TBL1YP1',
    'AC011297.1','AC010889.1',
    # X-specific
    'XIST','TSIX','KDM6A','UBA1','ZFX','RPS4X','DDX3X','USP9X',
    # mitochondrial
    'MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6',
    'MT-CO1','MT-CO2','MT-CO3','MT-ATP6','MT-ATP8','MT-CYB',
    'MT-RNR1','MT-RNR2','MT-TF','MT-TV','MT-TL1','MT-TI','MT-TQ','MT-TM',
    'MT-TW','MT-TA','MT-TN','MT-TC','MT-TY','MT-TS1','MT-TD','MT-TK','MT-TG',
    'MT-TR','MT-TH','MT-TS2','MT-TL2','MT-TE','MT-TT','MT-TP',
}

# ─────────────────────────── labels and colors ────────────────────────────────────────────

GDC_SAMPLE_TYPE_COLORS = {
    'Primary Tumor':        '#FF7043',   # deep orange
    'Solid Tissue Normal':  '#0077A8',   # teal-blue
    'Recurrent Tumor':      '#AB47BC',   # violet
    'Metastatic':           '#26A69A',   # teal
}

PRIMARY_HISTOLOGY_COLORS = {
    'FHWT': '#2ca02c',
    'DAWT': '#d62728'
}

DE_PALETTE = {
    'UP':   PRIMARY_HISTOLOGY_COLORS['DAWT'],   # red  — DAWT/anaplastic
    'DOWN': PRIMARY_HISTOLOGY_COLORS['FHWT'],   # green — FHWT
    'NS':   '#cccccc',
}

AT_PALETTE = {
    'UP':  '#d62728',   # red — up in Recurrent/anaplastic
    'DOWN': '#17becf',  # cyan — down in Recurrent (up in Primary)
    'NS':  '#cccccc',
}
TIME_PALETTE = {
    'Primary':   '#2ca02c',   # green
    'Recurrent': '#d62728',   # red
}

# Combined group palette for PCA visualization (all 136 samples)
PCA_GROUP_COLORS = {
    'Normal':        '#1f77b4',   # blue
    'Primary FHWT':  '#2ca02c',   # green
    'Primary DAWT':  '#d62728',   # red
    'Recurrent':     '#ff7f0e',   # orange
    'Metastatic':    '#9467bd',   # purple
}

RELAPSE_PALETTE = {
    'UP':   '#1f77b4',   # blue — Relapse
    'DOWN': '#9467bd',   # purple — Progression
    'NS':   '#cccccc',
}

COMMUNITY_LABELS = {
    1: "NF-κB / Cytokine core",
    2: "ECM / Matrix remodeling",
    3: "EMT / Adhesion",
    4: "AP-1 / Stress TF",
    5: "Innate / Acute phase",
    6: "S100 / Annexin",
}

TIER_COLORS = {"High": "#d62728", "Medium": "#ff7f0e", "Low": "#aec7e8"}

COMMUNITY_COLORS = {
    1: "#e6293a",   # red — cytokine
    2: "#2ca02c",   # green — ECM
    3: "#1f77b4",   # blue — EMT
    4: "#9467bd",   # purple — TF
    5: "#17becf",   # teal — innate
    6: "#e377c2",   # pink — S100
}

DRUG_COLORS = {
    "Approved target":     "#2ca02c",
    "Clinical candidate":  "#1f77b4",
    "Bioactive compounds": "#ff7f0e",
    "Poorly drugged":      "#d3d3d3",
    "Unknown":             "#e8e8e8",
}

FOCUS_TF_COLORS = {"RELA": "#e6293a", "NFKB1": "#ff7f0e", "JUN": "#9467bd", "STAT3": "#1f77b4"}