# =============================================================================
# 03_oncodrive.R — Wilms Tumor Analysis Project
# Purpose : Mutation hotspot / driver gene analysis using maftools::oncodrive.
#           Replaces dNdScv (incompatible with R-devel + low-TMB cohort).
#           oncodrive tests for clustering of mutations in functional hotspots,
#           which is more appropriate than dN/dS for small low-TMB cohorts.
#
# Input   : data/processed/rds/maf_object.rds     (from 02_mutations_maftools.R)
#           data/processed/tables/sample_annotation.csv
# Output  : data/processed/tables/oncodrive_all.csv
#           data/processed/tables/oncodrive_by_group.csv
#           data/processed/tables/oncodrive_wilms_panel.csv
#           results/figures/mut_09_oncodrive_all.pdf
#           results/figures/mut_10_oncodrive_wilms_panel.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(maftools)
  library(dplyr)
  library(ggplot2)
  library(here)
})

source(here("scripts", "R", "utils", "annotate_groups.R"))

DATA_RDS  <- here("data", "processed", "rds")
DATA_TABS <- here("data", "processed", "tables")
FIG_DIR   <- here("results", "figures")

MIN_SAMPLES_PER_GROUP <- 5

WILMS_GENES <- c("WT1", "CTNNB1", "WTX", "TP53", "SIX1", "SIX2",
                 "DROSHA", "DICER1", "IGF2", "MYCN", "DNMT3A")

# ── 1. Load inputs ────────────────────────────────────────────────────────────
message("── Loading MAF object and sample annotation ──")
maf_obj    <- readRDS(file.path(DATA_RDS, "maf_object.rds"))
sample_ann <- read.csv(file.path(DATA_TABS, "sample_annotation.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)

# Detect amino-acid change column
aa_col <- intersect(c("HGVSp_Short", "Protein_Change",
                       "amino_acid_change", "AAChange"),
                    names(maf_obj@data))[1]
if (is.na(aa_col)) {
  message("  No AA change column found — oncodrive hotspot detection may be sparse")
  aa_col <- NULL
} else {
  message("  Using AA column: ", aa_col)
}

# ── 2. Helper: run oncodrive safely ───────────────────────────────────────────
.run_oncodrive <- function(maf, label, min_mut = 3) {
  n_samp <- nrow(getSampleSummary(maf))
  message("  Running oncodrive: ", label, "  (n=", n_samp, " samples)")

  tryCatch({
    res <- oncodrive(maf = maf, AACol = aa_col,
                     minMut = min_mut, pvalMethod = "zscore")
    if (is.null(res) || nrow(res) == 0) {
      message("  No significant clusters found for '", label, "'")
      return(NULL)
    }
    res$label <- label
    res
  }, error = function(e) {
    message("  oncodrive failed for '", label, "': ", e$message)
    NULL
  })
}

# ── 3. Whole-cohort run ────────────────────────────────────────────────────────
message("\n── Oncodrive: whole cohort ──")
od_all <- .run_oncodrive(maf_obj, label = "all", min_mut = 2)

if (!is.null(od_all)) {
  write.csv(od_all, file.path(DATA_TABS, "oncodrive_all.csv"), row.names = FALSE)
  message("✓ oncodrive_all.csv  (", nrow(od_all), " genes)")

  message("\n  Wilms panel (whole cohort):")
  wilms_all <- od_all %>%
    filter(Hugo_Symbol %in% WILMS_GENES) %>%
    arrange(fdr)
  if (nrow(wilms_all) > 0)
    print(wilms_all[, intersect(c("Hugo_Symbol", "clusters", "muts",
                                   "pval", "fdr"), names(wilms_all))])
  else
    message("  No Wilms driver genes reached oncodrive minimum mutation threshold")
}

# ── 4. Per-group runs ──────────────────────────────────────────────────────────
message("\n── Oncodrive: per group ──")

group_counts <- sample_ann %>%
  filter(!is.na(group)) %>%
  count(group, name = "n_samples")
message("  Samples per group:")
print(group_counts)

groups_to_run <- group_counts %>%
  filter(n_samples >= MIN_SAMPLES_PER_GROUP, group != "Normal") %>%
  pull(group)

od_group_list <- lapply(groups_to_run, function(grp) {
  grp_barcodes <- sample_ann$Tumor_Sample_Barcode[sample_ann$group == grp]
  grp_maf      <- subsetMaf(maf = maf_obj, tsb = grp_barcodes)
  n_samp       <- length(grp_barcodes)

  if (n_samp < MIN_SAMPLES_PER_GROUP) {
    message("  Skipping '", grp, "': only ", n_samp, " samples")
    return(NULL)
  }

  # Lower minMut threshold for smaller groups
  min_mut <- if (n_samp < 10) 2 else 3
  res <- .run_oncodrive(grp_maf, label = grp, min_mut = min_mut)
  if (!is.null(res)) res$group <- grp
  res
})

od_by_group <- bind_rows(Filter(Negate(is.null), od_group_list))

if (nrow(od_by_group) > 0) {
  write.csv(od_by_group, file.path(DATA_TABS, "oncodrive_by_group.csv"),
            row.names = FALSE)
  message("✓ oncodrive_by_group.csv  (", nrow(od_by_group), " gene × group rows)")
} else {
  message("  No per-group results — oncodrive_by_group.csv not written")
}

# ── 5. Wilms driver panel summary ─────────────────────────────────────────────
message("\n── Wilms panel summary ──")

od_wilms <- bind_rows(
  if (!is.null(od_all)) od_all %>% mutate(group = "all") else NULL,
  if (nrow(od_by_group) > 0) od_by_group else NULL
) %>%
  filter(Hugo_Symbol %in% WILMS_GENES) %>%
  arrange(Hugo_Symbol, group)

if (nrow(od_wilms) > 0) {
  write.csv(od_wilms, file.path(DATA_TABS, "oncodrive_wilms_panel.csv"),
            row.names = FALSE)
  message("✓ oncodrive_wilms_panel.csv")
  print(od_wilms[, intersect(c("Hugo_Symbol", "group", "clusters",
                                "muts", "pval", "fdr"), names(od_wilms))])
} else {
  message("  No Wilms driver genes in oncodrive results")
}

# ── 6. Figures ─────────────────────────────────────────────────────────────────
message("\n── Oncodrive figures ──")

# Whole-cohort bubble plot
if (!is.null(od_all) && nrow(od_all) > 0) {
  tryCatch({
    pdf(file.path(FIG_DIR, "mut_09_oncodrive_all.pdf"), width = 10, height = 7)
    plotOncodrive(res = od_all, fdrCutOff = 0.1,
                  useFraction = TRUE, labelSize = 0.6)
    dev.off()
    message("Saved: mut_09_oncodrive_all.pdf")
  }, error = function(e) {
    message("  Oncodrive plot skipped: ", e$message)
  })
}

# Wilms panel: significance dot plot by group
if (nrow(od_wilms) > 0 && "group" %in% names(od_wilms) &&
    n_distinct(od_wilms$group) > 1) {
  tryCatch({
    plot_df <- od_wilms %>%
      filter(group != "all") %>%
      mutate(neg_logfdr = -log10(pmax(fdr, 1e-10)))

    p <- ggplot(plot_df,
                aes(x = reorder(Hugo_Symbol, -neg_logfdr),
                    y = neg_logfdr, colour = group, size = muts)) +
      geom_point(alpha = 0.8) +
      geom_hline(yintercept = -log10(0.1), linetype = "dashed",
                 colour = "grey50") +
      scale_colour_manual(values = GROUP_COLORS) +
      scale_size_continuous(range = c(2, 7), name = "# mutations") +
      labs(x = NULL, y = expression(-log[10](FDR)),
           title = "Wilms driver panel — oncodrive hotspot enrichment by group",
           colour = "Group") +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(FIG_DIR, "mut_10_oncodrive_wilms_panel.pdf"),
           p, width = 10, height = 5)
    message("Saved: mut_10_oncodrive_wilms_panel.pdf")
  }, error = function(e) {
    message("  Wilms panel plot skipped: ", e$message)
  })
}

message("\n── Driver analysis complete ──")
message("CSV outputs ready for Python transcriptomics integration:")
message("  ", file.path(DATA_TABS, "oncodrive_all.csv"))
message("  ", file.path(DATA_TABS, "oncodrive_by_group.csv"))
message("  ", file.path(DATA_TABS, "oncodrive_wilms_panel.csv"))
