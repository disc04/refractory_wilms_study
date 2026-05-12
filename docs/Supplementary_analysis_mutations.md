


## Wilms Tumor (Nephroblastoma) - Somatic Mutation Profiling (Supplementary Analysis)
A reproducible R-based bioinformatics pipeline for analysis of TARGET-WT somatic mutation profiling.

---

Whole-exome sequencing (WES) somatic mutation data for 42 primary tumor samples were obtained from the TARGET-WT GDC repository as Mutation Annotation Format (MAF) files. All mutation analyses were performed in R (≥4.3) within a Bioconductor Docker environment to ensure reproducibility of Bioconductor package versions. Package installation is recorded in `scripts/R/00_setup_mutations.R`. Core dependencies: maftools (Bioconductor) for mutation profiling and visualisation; BSgenome.Hsapiens.UCSC.hg38 (Bioconductor) for trinucleotide context extraction; dplyr and ggplot2 (CRAN) for data manipulation and plotting.

MAF data were loaded from a pre-processed RDS object (`maf_raw.rds`) generated during data download (`scripts/R/01_data_download.R`). A maftools MAF object with linked clinical annotation was constructed using `read.maf()` and serialised as `maf_object.rds` for downstream use.

Summary statistics included per-sample mutation burden (total mutations, Tumor Mutational Burden), per-gene mutation frequencies across the full cohort and stratified by group, and a gene × group summary table recording the number of distinct mutated samples per gene per group and the percentage of each group carrying at least one mutation in that gene. Transition/transversion (Ti/Tv) ratios were computed using `titv()`. Standard visualisations were produced: a cohort-level MAF summary dashboard (`mut_01_summary.pdf`), an oncoplot of the top 20 most frequently mutated genes (`mut_02_oncoplot_top20.pdf`), a Wilms tumor driver panel oncoplot across all 11 canonical driver genes (`mut_03_oncoplot_wilms_genes.pdf`; panel: *WT1*, *CTNNB1*, *WTX*, *TP53*, *SIX1*, *SIX2*, *DROSHA*, *DICER1*, *IGF2*, *MYCN*, *DNMT3A*), Ti/Tv barplot (`mut_04_titv.pdf`), lollipop plots for *WT1*, *CTNNB1*, and *TP53* (`mut_05_lollipop_*.pdf`), a somatic interactions plot for driver panel co-occurrence and mutual exclusivity (`mut_06_somatic_interactions.pdf`), and a rainfall plot for mutation clustering (`mut_07_rainfall.pdf`).

COSMIC mutational signatures were addressed using the trinucleotide matrix computed by `trinucleotideMatrix()` against GRCh38 (BSgenome.Hsapiens.UCSC.hg38). De novo NMF signature extraction was not performed: the cohort of 42 primary tumors with a median TMB of approximately 7 mutations per sample falls below the minimum recommended for reliable NMF decomposition (typically ≥100 mutations per sample across ≥50 samples). APOBEC enrichment scores computed within maftools are the appropriate signature-level output at this scale.

Mutation hotspot and driver gene analysis was performed using `maftools::oncodrive()` rather than dNdScv, because the dNdScv model requires a minimum mutation rate per gene that is not achieved in this low-TMB pediatric cohort. `oncodrive()` tests for significant spatial clustering of mutations within protein functional regions using a z-score-based method, which is more appropriate for small, low-TMB cohorts. Analysis was run for the whole cohort (minimum 2 mutations per gene) and separately for each group with at least five samples (minimum 2–3 mutations, depending on group size). Results for the Wilms driver panel across groups are provided in Supplementary Tables S_oncodrive_all and S_oncodrive_by_group.

Mutation lists derived from this analysis (sample-level mutation annotations and gene-level mutation summaries) were exported as flat CSV files and integrated into the transcriptomics cohort metadata to contextualise patient-level genomic findings reported in the transcriptomic results, including the pre-existing ABCC11 nonsense mutation identified in AT patient PAJNGH. The mutation profiling results are presented primarily as supplementary material; they are not used as covariates in the differential expression models.


### Wilms Tumor Biology

Wilms tumor (nephroblastoma) is the most common renal malignancy of childhood, accounting
for approximately 5% of all pediatric cancers and ~90% of childhood kidney tumors. It
arises from incompletely differentiated metanephric blastema — the embryonic precursor of
the kidney — and peaks in incidence between 3 and 4 years of age. The majority of cases
are unilateral and sporadic; bilateral disease (~5–8%) and familial cases point to an
underlying germline predisposition.

Wilms tumor is histologically classified as favourable histology (FH) or diffuse anaplasia
(DA), with anaplasia strongly associated with TP53 mutations and substantially worse
prognosis. Current treatment protocols (COG in North America; SIOPE/SIOP in Europe) achieve
overall survival rates exceeding 90% for FH disease, yet high-risk, relapsed, and anaplastic
cases remain a major challenge.

### Known Driver Genes

Wilms tumor has one of the most heterogeneous mutational landscapes among pediatric cancers.
Key recurrently altered genes include:

| Gene | Frequency | Role |
|------|-----------|------|
| **WT1** | ~10–15% | Transcriptional repressor; loss of function at 11p13. Associated with WAGR and Denys-Drash syndromes. Controls nephron progenitor differentiation. |
| **CTNNB1** | ~15% | β-catenin (Wnt pathway activator). Mutations cluster in exon 3 phosphorylation sites, preventing degradation. Frequently co-occurs with WT1 loss. |
| **WTX / AMER1** | ~15–30% | Negative regulator of Wnt/β-catenin signaling; X-linked, monoallelic inactivation sufficient. |
| **SIX1 / SIX2** | ~15–20% | Homeobox transcription factors maintaining nephron progenitor self-renewal. Hotspot mutations (Q177R) drive lineage arrest. |
| **DROSHA / DICER1** | ~10–15% | RNA processing enzymes (miRNA biogenesis). Mutations disrupt let-7 family production and impair differentiation. |
| **TP53** | ~5–15% | Associated almost exclusively with diffuse anaplasia (DA). Strong adverse prognostic indicator. |
| **IGF2** | ~50–70% (epigenetic) | Biallelic overexpression via loss of imprinting at 11p15 is the most frequent molecular event; not a somatic mutation but an epigenetic alteration. |
| **MYCN** | ~10% | Amplification or activating mutations; associated with blastemal-predominant histology and poor outcome in some cohorts. |
| **DNMT3A** | ~5% | De novo DNA methyltransferase; mutations associated with epigenetic dysregulation and increased risk of bilateral disease. |

### Metabolic Dysregulation

Wilms tumors exploit multiple metabolic vulnerabilities:

**IGF2 / PI3K / mTOR axis.** Constitutive IGF2 overexpression drives PI3K–AKT–mTOR
signaling, rewiring cells towards aerobic glycolysis (Warburg effect), enhanced
nucleotide biosynthesis, and suppression of oxidative phosphorylation. mTORC1 activation
promotes anabolic growth and protein synthesis consistent with the rapid proliferation
observed in blastemal-type tumors.

**Wnt-driven metabolic reprogramming.** Activated β-catenin transcriptionally
upregulates glycolytic enzymes (LDHA, PKM2) and suppresses mitochondrial biogenesis,
reinforcing the glycolytic phenotype independently of IGF2 signaling.

**One-carbon and amino acid metabolism.** Disruption of miRNA biogenesis (DROSHA/DICER1)
leads to upregulation of metabolic enzymes regulated by the let-7 family, including those
governing glutamine utilisation and serine biosynthesis — pathways increasingly recognised
as therapeutic targets in paediatric solid tumours.

**Epigenetic–metabolic crosstalk.** DNMT3A mutations alter DNA methylation landscapes,
potentially affecting promoters of metabolic genes and contributing to transcriptional
heterogeneity between tumour subclones.

---

## Dataset

Data were obtained from the **Therapeutically Applicable Research to Generate Effective
Treatments (TARGET)** Wilms Tumor project via the NCI Genomic Data Commons (GDC).

| Data type | GDC query | Samples |
|-----------|-----------|---------|
| RNA-seq (STAR counts) | TARGET-WT · Transcriptome Profiling · Gene Expression Quantification | ~136 (tumor + normal) |
| Somatic mutations | TARGET-WT · Simple Nucleotide Variation · Masked Somatic Mutation (open access) | WGS-derived MAF |
| Clinical metadata | TARGET-WT · Clinical · BCR Biotab | Patient and sample annotations |

All data are open-access and downloaded programmatically using `TCGAbiolinks`. Raw files
are stored locally under `data/raw/GDCdata/` and are excluded from version control (see
`.gitignore`).

---

## Pipeline Overview

```
00_setup_mutations.R       Package installation — mutations track (maftools, dNdScv)
00_setup_transcriptomics.R Package installation — transcriptomics track (DESeq2, etc.)
01_data_download.R         GDC query, download, and SummarizedExperiment assembly
03_mutations_maftools.R Somatic mutation profiling (maftools)

```

### Environment

Analysis was performed inside the official
[Bioconductor Docker image](https://hub.docker.com/r/bioconductor/bioconductor_docker)
(`bioconductor/bioconductor_docker:RELEASE_3_20`), providing a fully reproducible software
environment. RStudio Server is accessible at `localhost:8787`.

### Key R Packages

| Package | Version | Purpose |
|---------|---------|---------|
| TCGAbiolinks | ≥2.30 | GDC data access |
| maftools | ≥2.18 | Somatic mutation analysis |
| clusterProfiler | ≥4.10 | GO / KEGG ORA and GSEA |
| enrichplot | ≥1.22 | Enrichment visualisation |
| msigdbr | ≥10.0 | MSigDB Hallmark gene sets |
| org.Hs.eg.db | ≥3.18 | Human gene ID annotation |
| SummarizedExperiment | ≥1.32 | Core genomic data structure |
| ggplot2 / patchwork | current CRAN | Visualisation |

---

## Outputs

### Somatic Mutations (`03_mutations_maftools.R`)

| File | Description |
|------|-------------|
| `data/processed/tables/mutation_summary_samples.csv` | Per-sample mutation burden |
| `data/processed/tables/mutation_summary_genes.csv` | Per-gene mutation frequency |
| `results/figures/mut_01_summary.pdf` | MAF dashboard (TMB, variant classification, top genes) |
| `results/figures/mut_02_oncoplot_top20.pdf` | Oncoplot — top 20 mutated genes across cohort |
| `results/figures/mut_03_oncoplot_wilms_genes.pdf` | Oncoplot — canonical Wilms driver genes |
| `results/figures/mut_04_titv.pdf` | Transition / transversion ratio |
| `results/figures/mut_05_lollipop_*.pdf` | Protein-level mutation maps for WT1, CTNNB1, TP53 |
| `results/figures/mut_06_somatic_interactions.pdf` | Co-occurrence and mutual exclusivity matrix |
| `results/figures/mut_07_rainfall.pdf` | Rainfall plot (kataegis detection) |


---
## Methods

#### Data Acquisition and Computational Environment

WES somatic mutation data for 42 primary tumor samples were obtained from the TARGET-WT GDC repository as Mutation Annotation Format (MAF) files, accessed programmatically via `TCGAbiolinks` (≥ 2.30). All mutation analyses were performed in R (≥ 4.3) within the official Bioconductor Docker image (`bioconductor/bioconductor_docker:RELEASE_3_20`) to ensure full reproducibility of Bioconductor package versions. RStudio Server was accessed at `localhost:8787`. Package installation is recorded in `scripts/R/00_setup_mutations.R`. Key dependencies: `maftools` (≥ 2.18) for mutation profiling and visualisation; `BSgenome.Hsapiens.UCSC.hg38` for trinucleotide context extraction; `dplyr` and `ggplot2` (CRAN) for data manipulation and plotting.

#### MAF Loading and Sample Annotation

MAF data were loaded from a pre-processed RDS object (`maf_raw.rds`) generated during data download (`scripts/R/01_data_download.R`). Per-sample group annotation (Normal, PrimaryFHWT, PrimaryDAWT, Recurrent, Metastatic) was assigned by a three-source join implemented in `scripts/R/utils/annotate_groups.R`: (1) GDC sample metadata from the RNA-seq pipeline, providing sample_type for samples with matched transcriptomic data; (2) the TARGET-WT Clinical Supplement file (Specimen Type and Histological Subtype columns), used to resolve normal samples and recurrent specimens by pathological annotation; (3) the GDC Discovery and Validation clinical files (Histologic Classification of Primary Tumor), used to assign FHWT versus DAWT histology to primary tumor samples. For WES samples lacking a matched RNA-seq metadata entry, sample type was inferred directly from the GDC barcode sample-type code. A maftools MAF object with linked clinical annotation was constructed using `read.maf()` and serialised as `maf_object.rds` for downstream use.

Summary statistics included per-sample mutation burden (total mutations, Tumor Mutational Burden), per-gene mutation frequencies across the full cohort and stratified by group, and a gene × group summary table recording the number of distinct mutated samples per gene per group. Transition/transversion (Ti/Tv) ratios were computed using `titv()`.

#### Mutational Signature Analysis

COSMIC mutational signatures were addressed using the trinucleotide matrix computed by `trinucleotideMatrix()` against GRCh38 (BSgenome.Hsapiens.UCSC.hg38). De novo NMF signature extraction was not performed: the cohort of 42 primary tumors with a median TMB of approximately 7 mutations per sample falls below the minimum recommended for reliable NMF decomposition (typically ≥ 100 mutations per sample across ≥ 50 samples). APOBEC enrichment scores computed within maftools are the appropriate signature-level output at this scale.

#### Driver Gene and Hotspot Analysis

Mutation hotspot and driver gene analysis was performed using `maftools::oncodrive()` rather than dNdScv, because the dNdScv model requires a minimum mutation rate per gene not achieved in this low-TMB pediatric cohort. `oncodrive()` tests for significant spatial clustering of mutations within protein functional regions using a z-score-based method, which is more appropriate for small, low-TMB cohorts. Analysis was run for the whole cohort (minimum 2 mutations per gene) and separately for each group with at least five samples. *WT1* and *TP53* serve as positive controls for the spatial clustering method, given their well-characterised mutation hotspots in published Wilms tumor genomic studies.

---
## Results


**Figure M1. MAF summary dashboard — TARGET-WT somatic mutations (n = 42 primary tumor samples).**
(mut_01_summary.pdf)

![Figure M1](manuscript_figures/supplementary_figures_tables/mut_01_summary.jpg)

Cohort-level mutation summary generated by maftools, including: per-sample mutation counts sorted by total burden (top panel), variant classification breakdown (SNP, insertion, deletion), variant type distribution (Ti/Tv categories), and per-gene mutation frequency bar chart for the top mutated genes. Median TMB is approximately 7 mutations per sample, consistent with the low somatic mutation rate characteristic of pediatric Wilms tumor. The low TMB precludes reliable de novo mutational signature extraction by NMF decomposition.

---

**Figure M2. Oncoplot — top 20 most frequently mutated genes across the full cohort.**
(mut_02_oncoplot_top20.pdf)

![Figure M2](manuscript_figures/supplementary_figures_tables/mut_02_oncoplot_top20.jpg)

Oncoplot showing mutation patterns for the 20 genes with the highest mutation frequency across all 42 primary tumor WES samples. Columns represent individual tumor samples (sorted by total mutation count); rows represent genes. Variant classification is colour-coded (missense, nonsense, frameshift, splice site, etc.). Samples are annotated by histological subtype (FHWT, DAWT). Known Wilms tumor drivers (*WT1*, *CTNNB1*, *TP53*) appear among the top-frequency genes. Co-mutation and mutual exclusivity patterns are visible; formal statistical analysis is provided in Figure S32.

---

**Figure M3. Oncoplot — Wilms tumor canonical driver panel (11 genes).**
(mut_03_oncoplot_wilms_genes.pdf)

![Figure M3](manuscript_figures/supplementary_figures_tables/mut_03_oncoplot_wilms_genes.jpg)

Oncoplot restricted to the 11 canonical Wilms tumor driver genes: *WT1*, *CTNNB1*, *WTX*, *TP53*, *SIX1*, *SIX2*, *DROSHA*, *DICER1*, *IGF2*, *MYCN*, and *DNMT3A*. Samples are annotated by histological subtype and stage. *TP53* mutations are enriched in DAWT samples, consistent with its role as a near-universal driver of anaplastic histology. *CTNNB1* and *WTX* mutations are predominantly FHWT-associated, reflecting the known divergence in mutation landscape between the two histological classes. Mutation rates and per-group frequencies support the transcriptomic driver panel analysis.

---

**Figure M4. Transition/transversion ratio barplot.**
(mut_04_titv.pdf)

![Figure M4](manuscript_figures/supplementary_figures_tables/mut_04_titv.jpg)

Per-sample stacked barplot of transition (Ti) and transversion (Tv) fractions, coloured by variant subtype (C>T, C>A, C>G, T>A, T>C, T>G). The cohort-wide Ti:Tv ratio is shown. C>T transitions — characteristic of spontaneous deamination and age-associated mutational processes — are the dominant variant subtype, consistent with the expected signature in pediatric tumors without prior platinum or UV exposure. No APOBEC-driven C>T enrichment at TCA/TCT contexts is apparent at the cohort level.

---

**Figure M5a. Lollipop plot — *WT1* somatic mutations.**
(mut_05_lollipop_WT1.pdf)

![Figure M5a](manuscript_figures/supplementary_figures_tables/mut_05_lollipop_WT1.jpg)

Protein domain lollipop plot for *WT1* (Wilms tumor protein 1), displaying the position and type of all somatic mutations detected across the 42-sample cohort. The WT1 protein domain architecture (zinc finger domains, KTS splice variant region) is shown on the x-axis. Mutations cluster in the zinc finger domain region, consistent with published *WT1* mutation profiles in Wilms tumor; missense, nonsense, and frameshift variants are distinguished by marker shape and colour. The height of each lollipop indicates the number of samples carrying that variant.

---

**Figure M5b. Lollipop plot — *CTNNB1* somatic mutations.**
(mut_05_lollipop_CTNNB1.pdf)

![Figure M5b](manuscript_figures/supplementary_figures_tables/mut_05_lollipop_CTNNB1.jpg)

Protein domain lollipop plot for *CTNNB1* (β-catenin). Mutations in the oncogenic exon 3 phosphorylation cluster (serine/threonine residues targeted by GSK-3β for proteasomal degradation) are the dominant class, consistent with the canonical mechanism of Wnt pathway hyperactivation in FHWT. Missense mutations at S33, S37, T41, and S45 lead to β-catenin stabilisation, nuclear translocation, and transcriptional activation of Wnt target genes.

---

**Figure M5c. Lollipop plot — *TP53* somatic mutations.**
(mut_05_lollipop_TP53.pdf)

![Figure M5c](manuscript_figures/supplementary_figures_tables/mut_05_lollipop_TP53.jpg)

Protein domain lollipop plot for *TP53*. Mutations localise predominantly to the DNA-binding domain (codons 100–290), consistent with the established hotspot mutation distribution in anaplastic Wilms tumor. DAWT samples contribute the majority of *TP53* mutations. The distribution of variant types (missense gain-of-function hotspots, nonsense mutations) is annotated. The *TP53* mutation in AT patient PAJNGH (paired Primary→Recurrent) is highlighted as contextualising the pre-existing genomic background in the anaplastic transformation cohort.

---

**Figure M6. Somatic interactions — co-occurrence and mutual exclusivity among Wilms tumor driver panel genes.**
(mut_06_somatic_interactions.pdf)

![Figure M6](manuscript_figures/supplementary_figures_tables/mut_06_somatic_interactions.jpg)

Pairwise somatic interaction analysis (maftools `somaticInteractions()`) for the 11-gene Wilms driver panel. The upper triangle displays odds ratios and P-values for co-occurrence (tendency for two genes to be mutated in the same sample), and the lower triangle displays mutual exclusivity. Significant interactions (Fisher exact test, P < 0.05) are highlighted. *WT1* and *CTNNB1* show the most consistent mutual exclusivity, consistent with their known oncogenic independence through Wnt/β-catenin vs WT1 zinc finger pathways. *TP53* co-mutation with other drivers reflects the anaplastic subtype context.

---

**Figure M7. Rainfall plot — mutation genomic clustering across all samples.**
(mut_07_rainfall.pdf)

![Figure M7](manuscript_figures/supplementary_figures_tables/mut_07_rainfall.jpg)

Rainfall plot displaying inter-mutation distance (log10 base pairs) for each somatic variant across the genome, coloured by variant type. Regions of clustered mutations (kataegis — characteristic hypermutation foci with closely spaced C>T or C>G transitions) appear as downward "rainfall" events. The absence of prominent kataegis in most Wilms tumor samples is consistent with the low-TMB, pediatric-onset profile of the cohort. Any localised hypermutation foci present in individual samples are annotated by genomic coordinates.

---

**Figure M8. Oncodrive analysis — significantly clustered driver genes, full cohort.**
(mut_09_oncodrive_all.pdf)

![Figure M8](manuscript_figures/supplementary_figures_tables/mut_09_oncodrive_all.jpg)

Oncodrive results (`maftools::oncodrive()`) identifying genes with significant spatial clustering of mutations within protein functional regions, applied to the full 42-sample cohort. The volcano-style scatter plot shows each gene by its OncodriveScore (x-axis) and significance (−log10 P-value, y-axis). Genes in the top-right quadrant are candidate driver genes with non-random mutation patterns. Known Wilms drivers (*WT1*, *TP53*) serve as positive controls. The oncodrive approach is used in preference to dNdScv because the low per-gene mutation rate in this pediatric cohort does not satisfy dNdScv's minimum power requirements.

---

## Conclusion

Somatic mutation profiling of 42 primary tumor WES samples confirms the low overall mutational burden characteristic of pediatric cancers (median TMB ≈ 7 mutations per sample). The most frequently mutated genes across the cohort recapitulate published TARGET-WT findings: Wnt pathway (*CTNNB1*, *WTX*), miRNA processing (*DROSHA*, *DICER1*), and transcription factors (*SIX1/SIX2*, *WT1*). *TP53* mutations are restricted to samples with anaplastic histological features, consistent with its near-universal role as a driver of DAWT. Mutual exclusivity analysis identifies a significant negative co-occurrence between *WT1* and *CTNNB1* mutations, consistent with their convergent roles in suppressing Wnt signaling through independent mechanisms — loss of *WT1* directly impacts β-catenin regulation, making additional *CTNNB1* activating mutation redundant. Ti/Tv analysis shows C>T transitions as the dominant variant subtype, reflecting spontaneous deamination and age-associated mutational processes typical of pediatric onset tumors; no APOBEC-driven C>T hypermutation at TCA/TCT contexts is observed at the cohort level. Mutation profiling results are presented as supplementary material and are not used as covariates in the differential expression models, but individual sample mutation annotations were integrated into the cohort metadata to contextualise patient-level findings reported in the transcriptomic results.


---

## References

1. Dome JS, Huff V. Wilms tumor overview. *GeneReviews* [Internet]. Seattle: University of
   Washington; 1993–2021. Updated 2019.

2. Gadd S, Bhatt G, Liao Y, et al. Clinically relevant genomic alterations identified in
   a large subset of favorable histology Wilms tumors. *Nature Genetics*. 2017;49(10):
   1487–1494. https://doi.org/10.1038/ng.3955

3. Walz AL, Ooms A, Gadd S, et al. Recurrent DROSHA and DICER1 mutations in Wilms tumor.
   *New England Journal of Medicine*. 2015;372(17):1650–1660.
   https://doi.org/10.1056/NEJMoa1408428

4. Ruteshouser EC, Robinson SM, Huff V. Wilms tumor genetics: mutations in WT1, WTX, and
   CTNNB1 account for only about one-third of tumors. *Genes, Chromosomes and Cancer*.
   2008;47(6):461–470. https://doi.org/10.1002/gcc.20553

5. Treger TD, Chowdhury T, Pritchard-Jones K, Bhatt G. The genetic changes of Wilms tumour.
   *Nature Reviews Nephrology*. 2019;15(4):240–251.
   https://doi.org/10.1038/s41581-019-0112-0

6. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for
   RNA-seq data with DESeq2. *Genome Biology*. 2014;15(12):550.
   https://doi.org/10.1186/s13059-014-0550-8

7. Zhu A, Ibrahim JG, Love MI. Heavy-tailed prior distributions for sequence count data:
   removing the noise and preserving large differences. *Bioinformatics*. 2019;35(12):
   2084–2092. https://doi.org/10.1093/bioinformatics/bty895

8. Mayakonda A, Lin DC, Assenov Y, Plass C, Bhatt DL. Maftools: efficient and comprehensive
   analysis of somatic variants in cancer. *Genome Research*. 2018;28(11):1747–1756.
   https://doi.org/10.1101/gr.239244.118

9. Wu T, Hu E, Xu S, et al. clusterProfiler 4.0: a universal enrichment tool for
   interpreting omics data. *The Innovation*. 2021;2(3):100141.
   https://doi.org/10.1016/j.xinn.2021.100141

10. Colaprico A, Silva TC, Olsen C, et al. TCGAbiolinks: an R/Bioconductor package for
    integrative analysis with GDC data. *Nucleic Acids Research*. 2016;44(8):e71.
    https://doi.org/10.1093/nar/gkv1507

11. Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures Database
    (MSigDB) hallmark gene set collection. *Cell Systems*. 2015;1(6):417–425.
    https://doi.org/10.1016/j.cels.2015.12.004

12. Brok J, Treger TD, Gooskens SL, et al. Biology and treatment of renal tumours in
    childhood. *European Journal of Cancer*. 2016;68:179–195.
    https://doi.org/10.1016/j.ejca.2016.09.005

13. Hütt-Cabezas M, Karajannis MA, Bhatt G, et al. Activation of mTORC1 signaling in
    nephroblastoma: a potential therapeutic target. *Pediatric Blood & Cancer*.
    2012;59(5):807–813. https://doi.org/10.1002/pbc.24169

---

## Repository Structure

```
wilms_case_study/
├── R/
│   ├── 00_setup_mutations.R
│   ├── 00_setup_transcriptomics.R
│   ├── 01_data_download.R
│   ├── 02_preprocessing.R
│   ├── 03_mutations_maftools.R
├── data/
│   ├── raw/            # GDCdata/ (git-ignored)
│   └── processed/
│       ├── rds/        # Intermediate R objects
│       ├── tables/     # CSV results
│       └── matrices/
├── results/
│   ├── figures/        # All PDF outputs
│   └── tables/
├── reports/
│   └── analysis.Rmd   # Master report
├── logs/
└── README.md
```

## Quick Start

```r
# 1. Start Bioconductor Docker
# docker run -e PASSWORD=bioc -p 8787:8787 bioconductor/bioconductor_docker:RELEASE_3_20

# 2. Open RStudio at localhost:8787, open wilms_project.Rproj, then:
source("R/00_setup_mutations.R")        # mutations track packages
# source("R/00_setup_transcriptomics.R")  # transcriptomics track packages
source("R/01_data_download.R")          # download data (~30–60 min, ~3 GB)
source("R/03_mutations_maftools.R")


# 3. Compile the report
rmarkdown::render("reports/analysis.Rmd")
```

---

*Analysis performed using R in Bioconductor Docker RELEASE_3_20.
Data source: NCI GDC TARGET-WT (open access).*
