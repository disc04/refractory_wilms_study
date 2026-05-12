#!/usr/bin/env python3
"""
01_clinical_eda.py — Wilms Tumor Analysis Project
==================================================
Exploratory analysis of clinical metadata for the RNA-seq cohort (n=136).

Outputs
-------
results/figures/eda_01_cohort_overview.pdf
"""

from __future__ import annotations
import os
import pandas as pd
import matplotlib.pyplot as plt

from configuration import PATHS, PRIMARY_HISTOLOGY_COLORS, ANAPLASTIC_TRANSFORMATION_PATIENTS
from utils.data_utils import savefig

figures_path, processed_tables_path = PATHS["figures_path"], PATHS["processed_tables_path"]

def plot_cohort(df):
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    fig.suptitle("Cohort Clinical Summary (n=136, 124 primary tumor, 5 recurrent, 1 metastatic, 6 normal tissue)",
                 fontsize=12)

    # A — Histology × Stage (primary tumors only)
    ax = axes[0]
    pt = df[df['sample_type'] == 'Primary Tumor']
    ct = pd.crosstab(pt['Stage'], pt['histology'])
    ct.plot(kind='bar', ax=ax, color=[PRIMARY_HISTOLOGY_COLORS.get(c, '#999') for c in ct.columns],
            edgecolor='white', width=0.7)
    ax.set_xlabel("Stage")
    ax.set_ylabel("n")
    ax.set_title("Histology × Stage (Primary Tumors)")
    ax.tick_params(axis='x', rotation=0)
    ax.legend(title='')

    # B — Pathology subtype + confirmation
    ax = axes[1]
    sub = df[df['pathology_subtype'].notna()].copy()
    sub['confirmed'] = sub['Pathology: Diagnosis Confirmed'].str.lower().eq('yes')
    ct_c = (pd.crosstab(sub['pathology_subtype'], sub['confirmed'])
            .rename(columns={True: 'Confirmed', False: 'Not confirmed'}))
    ct_c[[c for c in ['Confirmed', 'Not confirmed'] if c in ct_c.columns]].plot(
        kind='bar', ax=ax, color=['#2196F3', '#BBDEFB'], edgecolor='white', width=0.6, stacked=True)
    ax.set_xlabel("")
    ax.set_ylabel("n")
    ax.set_title("Pathology Subtype + Confirmation")
    ax.tick_params(axis='x', rotation=15)
    ax.legend()

    plt.tight_layout()
    savefig(fig, os.path.join(figures_path, 'eda_01_cohort_overview.png'))
    print(f"Saved → eda_01_cohort_overview.png")
    print("\n── Complete ─────────────────────────────────────────────────")


def main():
    df = pd.read_csv(os.path.join(processed_tables_path, 'cohort_metadata.csv'))

    # Console summary
    print(f"  Wilms Tumor RNA-seq Cohort — Clinical EDA (n={len(df)})")
    print(f"{'-'*55}")
    print(df['sample_type'].value_counts().rename('n').to_string())
    print(df['histology'].value_counts().rename('n').to_string())
    print()

    print("Anaplastic transformation patients (FHWT primary → anaplastic at relapse):")
    at = df[df['usi'].isin(ANAPLASTIC_TRANSFORMATION_PATIENTS)]
    cols = ['usi', 'sample_type', 'histology',
            'pathology_subtype', 'Pathology: Diagnosis Confirmed']
    print(at[cols].sort_values(['usi', 'sample_type']).to_string(index=False))

    plot_cohort(df)




if __name__ == '__main__':
    main()
