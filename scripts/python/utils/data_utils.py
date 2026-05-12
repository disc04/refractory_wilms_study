from __future__ import annotations
import os
import glob
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from configuration import PATHS

rnaseq_data_path, clinical_data_path, processed_tables_path = (
    PATHS["rnaseq_data_path"], PATHS["clinical_data_path"], PATHS["processed_tables_path"])

_GDC_NA_VALUES = ['[Not Available]', '[Not Applicable]', '[Unknown]',
                  '[Discrepancy]', 'N/A', 'NA', '']

cohort_annotation_path = os.path.join(processed_tables_path, 'cohort_metadata_rnaseq.csv')

def load_xlsx(path: str, sheet_hints: list[str] | None = None) -> pd.DataFrame | None:
    """Load first matching sheet from an XLSX; drop all-empty rows/cols."""
    try:
        xl   = pd.ExcelFile(path, engine='openpyxl')
        sheets = xl.sheet_names
        target = sheets[0]
        if sheet_hints:
            for hint in sheet_hints:
                matches = [s for s in sheets if hint.lower() in s.lower()]
                if matches:
                    target = matches[0]
                    break
        df = xl.parse(target, header=0, na_values=_GDC_NA_VALUES,
                      keep_default_na=True)
        df = df.dropna(how='all').dropna(axis=1, how='all')
        # Strip whitespace from string columns
        for col in df.select_dtypes('object').columns:
            df[col] = df[col].str.strip() if hasattr(df[col], 'str') else df[col]
        return df
    except Exception as exc:
        print(f'  [error] Cannot read {os.path.basename(path)}: {exc}')
        return None


def extract_tables_from_folders(
    folders_path: str,
    skip_rows: int,
    use_columns: list[str],
) -> pd.DataFrame:
    """
    Concatenate per-sample tab-delimited files from a GDC folder tree into one DataFrame.

    Each immediate sub-directory of `folders_path` is treated as one sample.
    The directory name is added as a 'file_id' column (equivalent to sample_id).

    Parameters
    ----------
    folders_path : path to the parent directory containing one sub-folder per sample
    skip_rows    : number of header rows to skip in each file (passed to pd.read_csv)
    use_columns  : column names to retain; applied to all files

    Returns
    -------
    Concatenated DataFrame with one row per gene per sample and a 'file_id' column.

    Notes
    -----
    Assumes each sub-folder contains exactly one file and that folder names are unique.
    """
    frames = []
    folders = os.listdir(folders_path)
    folders = [f for f in folders if not f.startswith('.')]

    for folder in folders:
        filename = os.listdir(os.path.join(folders_path, folder))[0]
        file_path = os.path.join(folders_path, folder, filename)
        df = pd.read_csv(file_path, sep='\t', skiprows=skip_rows,
                         usecols=use_columns)
        df.columns = use_columns
        df['file_id'] = folder
        frames.append(df)
    concatenated = pd.concat(frames, axis=0)
    return concatenated

def load_annotated_count_matrix() -> pd.DataFrame:
    """
    Load counts_matrix.csv and merge with cohort annotation on 'file_id'.

    Returns
    -------
    Wide-format DataFrame with count columns and annotation columns joined.
    """
    counts = pd.read_csv(os.path.join(processed_tables_path, 'counts_matrix.csv'))
    annotation = pd.read_csv(cohort_annotation_path)
    counts_annotated = counts.merge(annotation, how='left', on='file_id')
    return counts_annotated


def savefig(fig: plt.Figure, path: str) -> None:
    """
    Save a matplotlib figure to disk at 150 dpi and close it.

    Parameters
    ----------
    fig  : figure to save
    path : output file path (including extension, e.g. 'results/figures/plot.png')
    """
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)


def _find_files(directory: str, pattern: str) -> list[str]:
    """Recursive glob under directory matching pattern."""
    return sorted(glob.glob(os.path.join(directory, '**', pattern), recursive=True))


def load_metadata_csv(path: str) -> pd.DataFrame | None:
    """
    Load the GDC metadata CSV produced by 01_data_download.R.

    Parameters
    ----------
    path : absolute path to the metadata CSV file

    Returns
    -------
    DataFrame with samples as rows, or None if the file is not found.
    """
    if not os.path.exists(path):
        print(f'  [missing] metadata file not found at:\n    {path}')
        return None
    df = pd.read_csv(path, index_col=0, low_memory=False,
                     na_values=_GDC_NA_VALUES, keep_default_na=True)
    print(f'  Loaded metadata file: {df.shape[0]} rows × {df.shape[1]} columns')
    return df

def _sep(title: str = '') -> None:
    """Print a horizontal separator line, optionally with a centred title."""
    print(f'\n{"─" * 62}')
    if title:
        print(f'  {title}')
        print(f'{"─" * 62}')


def load_gdc_counts_matrix(
    counts_path: str,
    sample_ids: list[str] | None = None,
    gene_universe: pd.Index | None = None,
) -> pd.DataFrame:
    """
    Chunked load of long-format GDC counts CSV → genes × samples raw count matrix.

    Filters to protein-coding genes during the read (memory-efficient).
    Optionally subsets samples and/or aligns to a pre-defined gene universe.

    Parameters
    ----------
    counts_path  : path to gdc_counts_matrix.csv (columns: sample_id, gene_name,
                   gene_type, unstranded)
    sample_ids   : if provided, only these samples are kept during load;
                   saves memory when loading a small subset (e.g. 8 AT samples)
    gene_universe: if provided, reindex result to this gene set after pivoting;
                   use to align a subset cohort to the primary analysis gene universe

    Returns
    -------
    DataFrame of shape (genes, samples) with raw integer-compatible counts
    """
    print("  Loading counts (chunked, protein-coding filter)…")
    sid_set = set(sample_ids) if sample_ids is not None else None
    chunks = []
    for chunk in pd.read_csv(
        counts_path,
        usecols=['sample_id', 'gene_name', 'gene_type', 'unstranded'],
        chunksize=500_000,
    ):
        pc = chunk[chunk['gene_type'] == 'protein_coding']
        if sid_set is not None:
            pc = pc[pc['sample_id'].isin(sid_set)]
        chunks.append(pc)

    pc = pd.concat(chunks, ignore_index=True).dropna(subset=['sample_id'])
    print(f"  Protein-coding genes: {pc['gene_name'].nunique():,} "
          f"across {pc['sample_id'].nunique()} samples")

    matrix = (
        pc.groupby(['gene_name', 'sample_id'])['unstranded']
        .sum()
        .unstack(level='sample_id')
    )

    if gene_universe is not None:
        matrix = matrix.reindex(gene_universe)
        n_missing = matrix.isna().all(axis=1).sum()
        if n_missing:
            print(f"  [warn] {n_missing} genes in universe have no counts in this subset")
        matrix = matrix.fillna(0)
        print(f"  After gene universe alignment: {matrix.shape[0]:,} genes")

    print(f"  Matrix: {matrix.shape[0]:,} genes × {matrix.shape[1]} samples")
    return matrix