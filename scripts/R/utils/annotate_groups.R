# =============================================================================
# utils/annotate_groups.R — Wilms Tumor Analysis Project
# Builds per-sample group annotation for the WES (MAF) cohort.
#
# Exports:
#   GROUP_COLORS             — named colour vector for consistent plotting
#   build_sample_annotation  — data.frame with group + reference columns
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(here)
})

GROUP_COLORS <- c(
  Normal      = "#1f77b4",
  PrimaryFHWT = "#2ca02c",
  PrimaryDAWT = "#d62728",
  Recurrent   = "#ff7f0e",
  Metastatic  = "#9467bd"
)

# ── Internal helpers ──────────────────────────────────────────────────────────

.barcode_to_usi <- function(b) {
  # TARGET-50-CAAAAA-01A-01D → TARGET-50-CAAAAA
  vapply(strsplit(b, "-"), function(x) paste(x[1:3], collapse = "-"),
         character(1))
}

.barcode_to_sample_id <- function(b) {
  # TARGET-50-CAAAAA-01A-01D → TARGET-50-CAAAAA-01A
  vapply(strsplit(b, "-"), function(x) paste(x[1:4], collapse = "-"),
         character(1))
}

# Infer GDC sample_type from barcode's 2-digit sample type code (4th element).
# Used as fallback when sample_metadata.csv does not cover a barcode.
.barcode_to_sample_type <- function(b) {
  code <- vapply(strsplit(b, "-"), function(x) substr(x[4], 1, 2), character(1))
  dplyr::recode(code,
    "01" = "Primary Tumor",
    "02" = "Recurrent Tumor",
    "06" = "Metastatic",
    "10" = "Blood Derived Normal",
    "11" = "Solid Tissue Normal",
    .default = NA_character_
  )
}

.find_clinical_file <- function(pattern) {
  base <- here("data", "raw", "GDCdata", "TARGET-WT",
               "Clinical", "Clinical_Supplement")
  hits <- Sys.glob(file.path(base, "*", pattern))
  if (length(hits) == 0) stop("Clinical file not found: ", pattern)
  hits[1]
}

# ── Main function ─────────────────────────────────────────────────────────────

#' Build per-sample group annotation for all barcodes present in the MAF.
#'
#' Three-source join (priority order in group assignment):
#'   1. Supplement  — Specimen Type contains "normal" → Normal
#'   2. Supplement  — Histological Subtype "relapse Wilms" → Recurrent
#'   3. GDC metadata / barcode code — sample_type → Normal / Recurrent / Metastatic
#'   4. Clinical files — Histologic Classification of Primary Tumor → PrimaryDAWT / FHWT
#'   5. Fallback → PrimaryUnk
#'
#' @param maf_barcodes  character vector — Tumor_Sample_Barcode values from MAF
#' @param metadata_path path to data/processed/tables/sample_metadata.csv
#' @return data.frame, one row per unique barcode, columns:
#'   Tumor_Sample_Barcode | sample_id | usi | sample_type |
#'   specimen_type | histological_subtype |   <- reference, not used for histology group
#'   primary_histology |                      <- DAWT | FHWT | NA  (from clinical only)
#'   group
build_sample_annotation <- function(maf_barcodes, metadata_path) {

  barcodes <- unique(maf_barcodes)

  ann <- data.frame(
    Tumor_Sample_Barcode = barcodes,
    usi                  = .barcode_to_usi(barcodes),
    sample_id            = .barcode_to_sample_id(barcodes),
    stringsAsFactors     = FALSE
  )

  # ── Source 1: GDC sample_metadata (sample_type) ───────────────────────────
  # sample_metadata.csv comes from colData(se_rnaseq) — covers RNA-seq samples.
  # WES-only samples fall back to barcode-inferred sample_type below.
  if (file.exists(metadata_path)) {
    meta <- read.csv(metadata_path, row.names = 1, check.names = FALSE,
                     stringsAsFactors = FALSE)
    meta$sample_id <- .barcode_to_sample_id(rownames(meta))

    keep_cols <- intersect(c("sample_id", "sample_type"), names(meta))
    meta_slim <- unique(meta[, keep_cols, drop = FALSE])

    ann <- merge(ann, meta_slim, by = "sample_id", all.x = TRUE)
  } else {
    message("  sample_metadata.csv not found — using barcode-inferred sample_type")
    ann$sample_type <- NA_character_
  }

  # Fallback: infer sample_type from barcode code for samples not in metadata
  missing_st <- is.na(ann$sample_type)
  if (any(missing_st)) {
    ann$sample_type[missing_st] <-
      .barcode_to_sample_type(ann$Tumor_Sample_Barcode[missing_st])
  }

  # ── Source 2: Supplement (Specimen Type + Histological Subtype) ───────────
  supp_path <- tryCatch(
    .find_clinical_file(paste0(
      "TARGET_WT_ClinicalData_Discovery_and_Validation_",
      "Percent_Tumor_Nuclei_and_Necrosis_Supplement_20230322.xlsx"
    )),
    error = function(e) { message("  Supplement not found: ", e$message); NULL }
  )

  if (!is.null(supp_path)) {
    supp <- as.data.frame(
      read_excel(supp_path, sheet = "WT Data", col_names = TRUE)
    )
    supp_slim <- supp[, c("TARGET - Sample ID", "Histological Subtype",
                           "Specimen Type")]
    names(supp_slim) <- c("sample_id", "histological_subtype", "specimen_type")
    supp_slim <- supp_slim[!is.na(supp_slim$sample_id), ]

    ann <- merge(ann,
                 supp_slim[, c("sample_id", "histological_subtype",
                                "specimen_type")],
                 by = "sample_id", all.x = TRUE)
  } else {
    ann$histological_subtype <- NA_character_
    ann$specimen_type        <- NA_character_
  }

  # ── Source 3: Discovery + Validation clinical (primary histology only) ─────
  # PrimaryDAWT / PrimaryFHWT are assigned ONLY from this source.
  clin_paths <- tryCatch(
    c(
      .find_clinical_file("TARGET_WT_ClinicalData_Discovery_20230322.xlsx"),
      .find_clinical_file("TARGET_WT_ClinicalData_Validation_20230322.xlsx")
    ),
    error = function(e) { message("  Clinical files not found: ", e$message); NULL }
  )

  if (!is.null(clin_paths)) {
    clin <- lapply(clin_paths, function(p) {
      df <- as.data.frame(
        read_excel(p, sheet = "Clinical Data", col_names = TRUE)
      )
      df[, c("TARGET USI", "Histologic Classification of Primary Tumor"),
         drop = FALSE]
    })
    clin <- do.call(rbind, clin)
    names(clin) <- c("usi", "primary_histology")
    clin <- clin[!duplicated(clin$usi) & !is.na(clin$usi), ]

    ann <- merge(ann, clin, by = "usi", all.x = TRUE)
  } else {
    ann$primary_histology <- NA_character_
  }

  # ── Group assignment — priority order ─────────────────────────────────────
  ann$group <- with(ann, mapply(
    function(specimen_type, histological_subtype, sample_type, primary_histology) {

      # 1. Supplement Specimen Type contains "normal"
      if (!is.na(specimen_type) &&
          grepl("normal", specimen_type, ignore.case = TRUE))
        return("Normal")

      # 2. Supplement Histological Subtype == "relapse Wilms"
      if (!is.na(histological_subtype) &&
          tolower(trimws(histological_subtype)) == "relapse wilms")
        return("Recurrent")

      # 3. GDC / barcode-inferred sample_type
      if (!is.na(sample_type)) {
        if (sample_type == "Solid Tissue Normal")    return("Normal")
        if (sample_type == "Blood Derived Normal")   return("Normal")
        if (sample_type == "Recurrent Tumor")        return("Recurrent")
        if (sample_type == "Metastatic")             return("Metastatic")
      }

      # 4. Primary histology — ONLY from clinical files
      if (!is.na(primary_histology)) {
        if (primary_histology == "DAWT") return("PrimaryDAWT")
        if (primary_histology == "FHWT") return("PrimaryFHWT")
      }

      "PrimaryUnk"
    },
    specimen_type, histological_subtype, sample_type, primary_histology,
    SIMPLIFY = TRUE, USE.NAMES = FALSE
  ))

  # Restore original barcode order and return
  ann[match(barcodes, ann$Tumor_Sample_Barcode), ]
}
