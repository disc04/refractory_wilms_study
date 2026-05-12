# =============================================================================
# 00_setup_mutations.R — Wilms Tumor Analysis Project
# Purpose : Install packages required for mutation analysis only.
#           Covers scripts 02_mutations_maftools.R and 03_dnds_dndscv.R.
#           Use this when you only need
#           somatic mutation / dN/dS workflows.
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
message("Track        : mutations (maftools + dNdScv)")

# ── 1. Package lists ───────────────────────────────────────────────────────────
bioc_pkgs <- c(
  "maftools",
  "BSgenome.Hsapiens.UCSC.hg38"   # COSMIC SBS mutational signatures
)

cran_pkgs <- c(
  "here",
  "dplyr",
  "ggplot2",
  "ggrepel",     # dN/dS volcano labels
  "remotes"      # needed to install dndscv from GitHub
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

message("\n── Installing GitHub packages ──")
if (!requireNamespace("dndscv", quietly = TRUE))
  remotes::install_github("im3sanger/dndscv")

# ── 3. Verify ──────────────────────────────────────────────────────────────────
all_pkgs   <- c(bioc_pkgs, cran_pkgs, "dndscv")
pkg_status <- data.frame(
  package   = all_pkgs,
  installed = all_pkgs %in% rownames(installed.packages())
)
missing_pkgs <- pkg_status[!pkg_status$installed, "package"]

if (length(missing_pkgs) > 0) {
  warning("These packages failed to install:\n  ",
          paste(missing_pkgs, collapse = "\n  "))
} else {
  message("\n✓ All mutation-track packages installed successfully!")
}
print(pkg_status)

# ── 4. Session info → log ──────────────────────────────────────────────────────
log_file <- file.path(LOGS, paste0("session_info_mutations_", Sys.Date(), ".txt"))
writeLines(capture.output(sessionInfo()), log_file)
message("\nSession info saved to: ", log_file)

message("\n── Mutations setup complete. Run scripts/R/02_mutations_maftools.R ──")
