# =============================================================================
# 01_data_download.R — Wilms Tumor Analysis Project
# Purpose : Download all TARGET-WT data from GDC via TCGAbiolinks.
#           Downloads: RNA-seq (STAR counts), clinical data, somatic MAF.
#
# Output  : data/processed/rds/se_rnaseq.rds      — SummarizedExperiment
#           data/processed/rds/maf_raw.rds          — MAF data.frame
#           data/processed/tables/clinical.csv      — clinical metadata
#
# Runtime : ~20-40 min depending on bandwidth (RNA-seq ~3 GB total).
#           Checkpoint files are saved; re-run is safe (skips if exists).
# =============================================================================

# Install download-specific packages on first run
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
for (pkg in c("TCGAbiolinks", "GEOquery")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
}

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(here)
})

# =============================================================================
# Workaround: TCGAbiolinks bug with TARGET-WT
# GDCprepare internally calls TCGAquery_subtype() which tries to select a
# "disease_response" column that does not exist in TARGET-WT clinical data.
# The patch below adds a NA placeholder for this column so the merge succeeds.
# Ref: https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues
# =============================================================================
.patch_tcgabiolinks_target_wt <- function() {
  ns   <- getNamespace("TCGAbiolinks")
  fns  <- c("TCGAquery_subtype", "addClinicalData")   # either may host the bug

  for (fn_name in fns) {
    if (!exists(fn_name, envir = ns, inherits = FALSE)) next
    orig <- get(fn_name, envir = ns)

    patched <- local({
      orig_fn <- orig
      function(...) {
        result <- tryCatch(orig_fn(...), error = function(e) {
          if (grepl("disease_response", conditionMessage(e), fixed = TRUE)) {
            message("  [patch] Missing 'disease_response' column — returning NULL")
            return(NULL)
          }
          stop(e)
        })
        # Add missing column if data.frame was returned
        if (is.data.frame(result) &&
            !"disease_response" %in% names(result)) {
          result$disease_response <- NA_character_
        }
        result
      }
    })

    tryCatch({
      unlockBinding(fn_name, ns)
      assign(fn_name, patched, envir = ns)
      lockBinding(fn_name, ns)
      message("✓ Patched TCGAbiolinks::", fn_name,
              " for TARGET-WT compatibility")
    }, error = function(e) {
      message("  Could not patch ", fn_name, ": ", conditionMessage(e))
    })
  }
}

.patch_tcgabiolinks_target_wt()

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_RAW  <- here("data", "raw")
DATA_RDS  <- here("data", "processed", "rds")
DATA_TABS <- here("data", "processed", "tables")

dir.create(DATA_RAW,  recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_RDS,  recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_TABS, recursive = TRUE, showWarnings = FALSE)

# Set GDC download directory inside the project
GDCdata_dir <- file.path(DATA_RAW, "GDCdata")
dir.create(GDCdata_dir, recursive = TRUE, showWarnings = FALSE)

# ── Helper: checkpoint skip ───────────────────────────────────────────────────
checkpoint <- function(path, label) {
  if (file.exists(path)) {
    message("✓ Already exists, skipping: [", label, "] ", basename(path))
    return(TRUE)
  }
  return(FALSE)
}

# ── Helper: build SE directly from STAR count TSV files ──────────────────────
# Fallback when GDCprepare fails. Reads each downloaded TSV and assembles the
# SummarizedExperiment manually. Returns NULL if no files are found.
.build_se_from_star_files <- function(query, gdc_dir) {
  message("  Building SE manually from downloaded STAR count files...")

  manifest <- getResults(query)

  # GDC stores files at: <dir>/<project>/<case_id>/<uuid>/<filename>
  tsv_files <- list.files(gdc_dir,
    pattern    = "\\.rna_seq\\.augmented_star_gene_counts\\.tsv$",
    recursive  = TRUE,
    full.names = TRUE
  )
  if (length(tsv_files) == 0) {
    # Fall back to any TSV in the download dir
    tsv_files <- list.files(gdc_dir,
      pattern    = "\\.tsv$",
      recursive  = TRUE,
      full.names = TRUE
    )
  }
  if (length(tsv_files) == 0) {
    message("  No TSV files found in ", gdc_dir, " — cannot build SE manually.")
    return(NULL)
  }
  message("  Found ", length(tsv_files), " count files")

  # Read first file to get gene index and metadata columns
  first  <- read.delim(tsv_files[1], comment.char = "#",
                        stringsAsFactors = FALSE)
  gene_ids  <- first$gene_id
  gene_info <- first[, intersect(c("gene_id","gene_name","gene_type"),
                                  names(first)), drop = FALSE]
  count_col <- intersect(c("unstranded","expected_count"), names(first))[1]

  # Build count + TPM matrices
  n_genes   <- length(gene_ids)
  n_samples <- length(tsv_files)
  mat_counts <- matrix(NA_integer_,  nrow = n_genes, ncol = n_samples)
  mat_tpm    <- matrix(NA_real_,     nrow = n_genes, ncol = n_samples)
  sample_ids <- character(n_samples)
  file_ids   <- character(n_samples)

  # Determine which manifest column carries a unique sample-level barcode.
  # Priority: cases.samples.submitter_id > sample.submitter_id > file UUID
  sample_col <- intersect(
    c("cases.samples.submitter_id", "sample.submitter_id"),
    names(manifest)
  )[1]
  if (is.na(sample_col)) {
    message("  No sample-level barcode column found; will use file UUIDs as sample IDs")
  } else {
    message("  Using '", sample_col, "' as sample identifier")
  }

  for (i in seq_along(tsv_files)) {
    df <- read.delim(tsv_files[i], comment.char = "#",
                     stringsAsFactors = FALSE)
    mat_counts[, i] <- as.integer(df[[count_col]])
    tpm_col <- intersect(c("tpm_unstranded","tpm"), names(df))[1]
    if (!is.na(tpm_col)) mat_tpm[, i] <- as.numeric(df[[tpm_col]])

    # Primary match: UUID is the directory name one level above the file
    file_uuid <- basename(dirname(tsv_files[i]))
    idx        <- which(manifest$id == file_uuid)

    # Fallback: match by file basename
    if (length(idx) == 0) {
      bn  <- basename(tsv_files[i])
      idx <- which(manifest$file_name == bn)
      if (length(idx) > 0) file_uuid <- manifest$id[idx[1]]
    }

    if (length(idx) > 0) {
      sid <- if (!is.na(sample_col)) manifest[[sample_col]][idx[1]] else NA_character_
      # Use sample barcode if valid, else fall back to file UUID
      sample_ids[i] <- if (!is.na(sid) && nzchar(sid)) sid else file_uuid
      file_ids[i]   <- file_uuid
    } else {
      # Unknown file — use UUID derived from path
      sample_ids[i] <- file_uuid
      file_ids[i]   <- file_uuid
    }
  }

  # Guard: if sample_ids still has duplicates (e.g. two files per sample),
  # append a numeric suffix to guarantee uniqueness.
  if (anyDuplicated(sample_ids)) {
    message("  Duplicate sample IDs detected; appending file suffix to resolve")
    dup_mask <- duplicated(sample_ids) | duplicated(sample_ids, fromLast = TRUE)
    sample_ids[dup_mask] <- paste0(sample_ids[dup_mask], "_",
                                    substr(file_ids[dup_mask], 1, 8))
  }

  rownames(mat_counts) <- gene_ids
  colnames(mat_counts) <- sample_ids
  rownames(mat_tpm)    <- gene_ids
  colnames(mat_tpm)    <- sample_ids

  # Build colData — match rows by file UUID (always unique) then name by sample_ids
  col_data <- manifest[match(file_ids, manifest$id), , drop = FALSE]
  if (is.null(col_data) || nrow(col_data) == 0) {
    col_data <- data.frame(
      sample_id   = sample_ids,
      sample_type = "unknown",
      row.names   = sample_ids
    )
  } else {
    rownames(col_data) <- sample_ids
  }

  row_data <- S4Vectors::DataFrame(gene_info, row.names = gene_ids)

  se <- SummarizedExperiment(
    assays  = list(unstranded = mat_counts, tpm_unstrand = mat_tpm),
    rowData = row_data,
    colData = col_data
  )
  message("  Manual SE: ", nrow(se), " genes x ", ncol(se), " samples")
  se
}

# =============================================================================
# 1. RNA-seq (STAR counts) — TARGET-WT
# =============================================================================
se_path <- file.path(DATA_RDS, "se_rnaseq_raw.rds")

if (!checkpoint(se_path, "RNA-seq SE")) {

  message("\n── Step 1: Query RNA-seq data ──")
  query_rna <- GDCquery(
    project           = "TARGET-WT",
    data.category     = "Transcriptome Profiling",
    data.type         = "Gene Expression Quantification",
    workflow.type     = "STAR - Counts",
    access            = "open"
  )

  results_df <- getResults(query_rna)
  message("Files found: ", nrow(results_df))
  message("Samples:\n", paste(
    capture.output(table(results_df$sample_type)), collapse = "\n"))

  message("\n── Step 2: Download RNA-seq files ──")
  GDCdownload(
    query          = query_rna,
    method         = "api",
    files.per.chunk = 10,
    directory      = GDCdata_dir
  )

  message("\n── Step 3: Build SummarizedExperiment ──")
  se <- tryCatch(
    GDCprepare(query = query_rna, directory = GDCdata_dir),
    error = function(e) {
      msg <- conditionMessage(e)
      message("GDCprepare error: ", msg)
      if (grepl("disease_response|doesn't exist|don't exist", msg)) {
        message("Falling back to manual SE construction...")
        .build_se_from_star_files(query_rna, GDCdata_dir)
      } else {
        stop(e)
      }
    }
  )

  if (is.null(se)) stop("SE construction failed. Check GDCdata_dir: ", GDCdata_dir)

  message("SE dimensions: ", nrow(se), " genes x ", ncol(se), " samples")
  message("Assay names: ", paste(assayNames(se), collapse = ", "))
  message("colData fields: ", paste(names(colData(se)), collapse = ", "))

  saveRDS(se, se_path)
  message("✓ Saved: ", se_path)
}

# =============================================================================
# 2. Clinical data
# =============================================================================
clinical_path <- file.path(DATA_TABS, "clinical_raw.csv")

if (!checkpoint(clinical_path, "clinical")) {

  message("\n── Step 4: Download clinical data ──")
  query_clin <- GDCquery(
    project       = "TARGET-WT",
    data.category = "Clinical",
    data.type     = "Clinical Supplement",
    data.format   = "BCR Biotab"
  )
  GDCdownload(query_clin, directory = GDCdata_dir)
  clin_list <- GDCprepare(query_clin, directory = GDCdata_dir)

  # The prepare returns a list of data frames; use the patient table
  if (is.list(clin_list)) {
    # Flatten all available clinical tables
    for (nm in names(clin_list)) {
      out <- file.path(DATA_TABS, paste0("clinical_", nm, ".csv"))
      write.csv(clin_list[[nm]], out, row.names = FALSE)
      message("Saved: ", basename(out))
    }
    # Use patient table as primary clinical
    if ("clinical_patient" %in% names(clin_list)) {
      clinical <- clin_list[["clinical_patient"]]
    } else {
      clinical <- clin_list[[1]]
    }
  } else {
    clinical <- clin_list
  }

  write.csv(clinical, clinical_path, row.names = FALSE)
  message("✓ Saved clinical: ", clinical_path,
          " (", nrow(clinical), " rows)")
}

# =============================================================================
# 3. Somatic mutations (MAF) — open access masked somatic mutations
# =============================================================================
maf_path <- file.path(DATA_RDS, "maf_raw.rds")

if (!checkpoint(maf_path, "MAF")) {

  message("\n── Step 5: Download somatic mutation MAF ──")
  query_maf <- GDCquery(
    project       = "TARGET-WT",
    data.category = "Simple Nucleotide Variation",
    data.type     = "Masked Somatic Mutation",
    access        = "open"
  )

  files_maf <- getResults(query_maf)
  message("MAF files found: ", nrow(files_maf))

  GDCdownload(query_maf, directory = GDCdata_dir)
  maf_data <- GDCprepare(query_maf, directory = GDCdata_dir)

  saveRDS(maf_data, maf_path)
  message("✓ Saved MAF: ", maf_path,
          " (", nrow(maf_data), " mutations)")
}

# =============================================================================
# 4. Extract and save metadata from SE (for downstream use without SE load)
# =============================================================================
meta_path <- file.path(DATA_TABS, "sample_metadata.csv")

if (!checkpoint(meta_path, "metadata") && file.exists(se_path)) {
  se       <- readRDS(se_path)
  metadata <- as.data.frame(colData(se))
  write.csv(metadata, meta_path, row.names = TRUE)
  message("✓ Saved sample metadata: ", meta_path)
}

# =============================================================================
# Summary
# =============================================================================
message("\n", strrep("=", 60))
message("Download complete. Files:")
downloaded <- c(se_path, maf_path, clinical_path, meta_path)
for (f in downloaded) {
  exists <- file.exists(f)
  size   <- if (exists) paste0(round(file.size(f)/1e6, 1), " MB") else "missing"
  message(sprintf("  %s  %s  %s", if (exists) "\u2713" else "\u2717",
                  basename(f), size))
}
message(strrep("=", 60))
message("Next: run 02_preprocessing.R")
