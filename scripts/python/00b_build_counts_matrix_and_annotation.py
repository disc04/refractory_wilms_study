#!/usr/bin/env python3
"""
00b_build_counts_matrix_and_annotation.py
===========================================
Run once to build the counts matrix and cohort metadata annotation table.
Builds cohort metadata annotation table.
Concatenates per-sample STAR `augmented_star_gene_counts.tsv` RNAseq files
(as downloaded from GDC for TARGET-WT) into a single genes × samples counts matrix.

Inputs
------
- data/raw/GDCdata/TARGET-WT/Transcriptome_Profiling/Gene_Expression_Quantification/<file_id>/<filename>.tsv
- data/processed/tables/sample_metadata.csv (GDC manifest with file_id <-> sample_id)
- data/raw/GDCdata/TARGET-WT/Clinical/ Clinical_Supplement/5e3502fb-7211-4716-8355-8eaaa90c0c43'/
                   'TARGET_WT_ClinicalData_Validation_20230322.xlsx'
- data/raw/GDCdata/TARGET-WT/Clinical/ Clinical_Supplement/'72bd5aba-4b34-4d74-861c-5586b8321ab7/
TARGET_WT_ClinicalData_Discovery_and_Validation_Percent_Tumor_Nuclei_and_Necrosis_Supplement_20230322.xlsx'


Outputs
-------
- data/processed/tables/gdc_counts_matrix.csv (genes × samples, unstranded counts, metadata)
- data/processed/tables/cohort_metadata.csv
"""

from __future__ import annotations
import os
import pandas as pd
from configuration import rnaseq_data_path, clinical_data_path, processed_tables_path
from utils.data_utils import extract_tables_from_folders, load_xlsx, load_metadata_csv

CLINICAL_COLUMNS = [
    'TARGET USI',
    'Gender',
    'First Event',
    'Age at Diagnosis in Days',
    'Vital Status',
    'Event Free Survival Time in Days',
    'Overall Survival Time in Days',
    'Stage',
    'Histologic Classification of Primary Tumor',
    'Reason for Death',
]
PATHOLOGY_COLUMNS = [
        'TARGET - Sample ID',
        'Histological Subtype',
        'Specimen Type',
        'Adequate Tissue (Y/N)',
        'Pathology: Diagnosis Confirmed'
    ]

# Merge Clinical Annotation ───────────────────────────────────────────────────────────
def build_cohort_annotation() -> pd.DataFrame:
    """
    Prepares sample metadata as follows:
    1. Sample_metadata.csv → keep 4 cols, rename 2 (n=136, RNA-seq cohort)
    2. Enrich with validation metadata → keep 11 cols, rename 1 → left-merge on usi
    3. Enrich with pathology metadata → keep selected cols, rename 4 → left-merge on sample_id
    """
    sample_metadata_path = os.path.join(processed_tables_path, 'sample_metadata.csv')
    raw_meta = load_metadata_csv(sample_metadata_path)
    meta = raw_meta[['cases.submitter_id', 'file_id', 'cases', 'sample_type']].copy()
    meta.index.name = 'sample_id'
    meta = meta.reset_index()
    meta = meta.rename(columns={'cases.submitter_id': 'usi'})

    clinical_path = os.path.join(clinical_data_path, 'd3002b03-71e8-4395-90ab-cd5c1020efe1',
                                 'TARGET_WT_ClinicalData_Discovery_20230322.xlsx')
    clinical_val = load_xlsx(clinical_path, sheet_hints=None)
    clinical_val = clinical_val[CLINICAL_COLUMNS].rename(columns={'TARGET USI': 'usi',
                                                                  'Histologic Classification of Primary Tumor': 'histology'})
    combined = meta.merge(clinical_val,  how='left', on='usi')


    pathology_path = os.path.join(clinical_data_path,'72bd5aba-4b34-4d74-861c-5586b8321ab7',
     'TARGET_WT_ClinicalData_Discovery_and_Validation_Percent_Tumor_Nuclei_and_Necrosis_Supplement_20230322.xlsx')
    pathology = load_xlsx(pathology_path, sheet_hints=None)
    pathology = pathology[PATHOLOGY_COLUMNS].rename(columns ={
        'TARGET - Sample ID':  'sample_id',
        'Histological Subtype': 'pathology_subtype',
        'Specimen Type':        'pathology_specimen_type'})
    pathology['Pathology: Diagnosis Confirmed'] = \
        pathology['Pathology: Diagnosis Confirmed'].str.capitalize()

    combined = combined.merge(pathology,  how='left', on='sample_id').drop_duplicates('sample_id')
    return combined

# Build Mutations Annotation ───────────────────────────────────────────────────────────
def build_mutations_annotation():
    """Builds mutations annotation table."""

    def aggregate_mutations(group):
        mutated_genes = sorted(list(set(group['Hugo_Symbol'].tolist())))
        mutations_comment = group.to_dict(orient='records')  # usi already excluded
        return pd.Series({
            'mutated_genes': mutated_genes,
            'mutations_comment': mutations_comment
        })

    raw = pd.read_csv(os.path.join(processed_tables_path, 'mutations_annotated.csv'))
    raw = raw[['usi', 'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
               'Strand', 'Variant_Classification', 'Variant_Type']]
    annotation = raw.groupby('usi').apply(aggregate_mutations, include_groups=False).reset_index()

    return annotation

# Concatenate count tables ───────────────────────────────────────────────────────────
def build_counts_matrix(coutns_path, metadata) -> None:
    """Concatenates RNAseq .tsv` files into a single genes × samples counts matrix."""
    counts = extract_tables_from_folders(coutns_path, skip_rows=1,
                                         use_columns=['gene_id', 'gene_name', 'gene_type', 'unstranded'])
    invalid_gene_ids = ['N_unmapped', 'N_multimapping', 'N_noFeature', 'N_ambiguous']
    counts = counts[~counts['gene_id'].isin(invalid_gene_ids)]

    counts = enrich_with_sample_id(counts, metadata[['sample_id', 'file_id']])
    print(f'Cohort counts matrix shape: {counts.shape}')
    counts.to_csv(os.path.join(processed_tables_path, 'gdc_counts_matrix.csv'), index=False)

def enrich_with_sample_id(counts: pd.DataFrame, cohort_metadata: pd.DataFrame) -> pd.DataFrame:
    """Adds sample_id column to counts matrix."""
    counts = counts.merge(cohort_metadata, how='left', on='file_id')
    return counts

if __name__ == '__main__':
    # Check n=136 samples in the Wilms RNA-seq cohort

    # Create cohort metadata annotation table
    cohort_metadata = build_cohort_annotation()
    # Enrich with mutations metadata
    mutations_annotation = build_mutations_annotation()
    cohort_metadata = cohort_metadata.merge(mutations_annotation, how='left', on='usi')
    cohort_metadata.to_csv(os.path.join(processed_tables_path, 'cohort_metadata.csv'), index=False)

    print(f'{cohort_metadata.shape[0]} samples found in cohort metadata')
    print(f'  -> {processed_tables_path}/cohort_metadata.csv')

    # Create cohort counts matrix
    folders = [f for f in os.listdir(rnaseq_data_path) if not f.startswith('.')]
    print(f'{len(folders)} files found in {rnaseq_data_path}')
    build_counts_matrix(coutns_path=rnaseq_data_path, metadata=cohort_metadata)

    print(f'  -> {processed_tables_path}/gdc_counts_matrix.csv')
    print("\n── Complete ─────────────────────────────────────────────────")