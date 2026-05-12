# =============================================================================
# 00_setup_transcriptomics.R — Wilms Tumor Analysis Project
# Purpose : Install packages required for transcriptomics analysis.
#           Covers DE (DESeq2/edgeR/limma), pathway analysis, and
#           visualization workflows. TCGAbiolinks and GEOquery are
#           installed on demand by 01_data_download.R.
#           For mutation-only work use 00_setup_mutations.R instead.
# =============================================================================

# ── 0. Bootstrap ──────────────────────────────────────────────────────────────
if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

PROJECT_ROOT <- here::here()
LOGS         <- file.path(PROJECT_ROOT, "logs")
dir.create(LOGS, recursive = TRUE, showWarnings = FALSE)

message("Project root : ", PROJECT_ROOT)
message("Track        : transcriptomics")

# ── 1. Package lists ───────────────────────────────────────────────────────────
bioc_pkgs <- c(
  # RNA-seq / DE
  # Note: TCGAbiolinks and GEOquery are installed by 01_data_download.R
  "DESeq2",
  "apeglm",                  # LFC shrinkage for DESeq2 (lfcShrink type="apeglm")
  "edgeR",
  "limma",
  "tximport",
  # Genomic data structures
  "SummarizedExperiment",
  "GenomicRanges",
  "Biostrings",
  # Pathway analysis
  "clusterProfiler",
  "fgsea",
  "org.Hs.eg.db",
  "DOSE",
  "enrichplot",
  # Visualization
  "EnhancedVolcano",
  "ComplexHeatmap",
  "circlize",
  # Single-cell (future)
  "SingleCellExperiment"
)

cran_pkgs <- c(
  "here",
  "tidyverse",       # dplyr, tidyr, ggplot2, readr, purrr, stringr
  "janitor",         # clean column names
  "ggpubr",          # publication-ready ggplot2 helpers
  "patchwork",       # combine ggplot2 panels
  "ggridges",        # ridge plots (required by enrichplot::ridgeplot)
  "ggrepel",         # non-overlapping labels
  "viridis",
  "RColorBrewer",
  "pheatmap",
  "survival",
  "survminer",       # Kaplan-Meier plots
  "DT",              # interactive tables in Rmd
  "gt",              # publication-ready tables
  "renv",            # reproducibility
  "readxl",          # read Excel clinical files
  "msigdbr"          # MSigDB gene sets for fgsea / clusterProfiler
)

# ── 2. Install ─────────────────────────────────────────────────────────────────
message("\n── Installing CRAN packages ──")
installed_cran  <- rownames(installed.packages())
to_install_cran <- setdiff(cran_pkgs, installed_cran)
if (length(to_install_cran) > 0) {
  install.packages(to_install_cran, repos = "https://cloud.r-project.org")
} else {
  message("All CRAN packages already installed.")
}

message("\n── Installing Bioconductor packages ──")
BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE)

# ── 3. Verify ──────────────────────────────────────────────────────────────────
all_pkgs   <- c(bioc_pkgs, cran_pkgs)
pkg_status <- data.frame(
  package   = all_pkgs,
  installed = all_pkgs %in% rownames(installed.packages())
)
missing_pkgs <- pkg_status[!pkg_status$installed, "package"]

if (length(missing_pkgs) > 0) {
  warning("These packages failed to install:\n  ",
          paste(missing_pkgs, collapse = "\n  "))
} else {
  message("\n✓ All transcriptomics-track packages installed successfully!")
}
print(pkg_status)

# ── 4. Session info → log ──────────────────────────────────────────────────────
log_file <- file.path(LOGS, paste0("session_info_transcriptomics_", Sys.Date(), ".txt"))
writeLines(capture.output(sessionInfo()), log_file)
message("\nSession info saved to: ", log_file)

# ── 5. Check GDC connection ───────────────────────────────────────────────────
message("\n── Checking GDC API connection ──")
tryCatch({
  library(TCGAbiolinks)
  projects <- getGDCprojects()
  wt_proj  <- projects[projects$id == "TARGET-WT", ]
  if (nrow(wt_proj) > 0) {
    message("✓ GDC connection OK — TARGET-WT found: ", wt_proj$name)
  } else {
    message("GDC connected but TARGET-WT not found in project list.")
  }
}, error = function(e) {
  warning("GDC connection failed (check internet access): ", conditionMessage(e))
})

message("\n── Transcriptomics setup complete. Proceed to 01_data_download.R ──")
