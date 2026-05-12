#!/usr/bin/env python3
"""
00a_download_gdc_data.py
========================
Download TARGET-WT RNA-seq and clinical data from NCI GDC.

Run this once before 00b_build_counts_matrix_and_annotation.py.
Uses the GDC REST API (https://api.gdc.cancer.gov) via `requests` — no
GDC client install needed. All data downloaded here is open-access.

Downloads
---------
1. RNA-seq STAR count TSVs
       → data/raw/GDCdata/TARGET-WT/Transcriptome_Profiling/
         Gene_Expression_Quantification/<file_uuid>/<filename>.tsv

2. Clinical Supplement XLSX files
       → data/raw/GDCdata/TARGET-WT/Clinical/
         Clinical_Supplement/<file_uuid>/<filename>.xlsx

3. Sample metadata CSV
       → data/processed/tables/sample_metadata.csv
         Maps file_uuid → sample barcode, sample type, case ID.
         Required by 00b_build_counts_matrix_and_annotation.py.

Runtime : ~30-60 min depending on bandwidth (~3 GB RNA-seq total).
Resume  : safe to re-run — existing files are skipped automatically.
"""

from __future__ import annotations

import sys
import json
import time
import requests
import pandas as pd
from pathlib import Path
from configuration import ROOT, PATHS

# ── GDC API constants ─────────────────────────────────────────────────────────
GDC_API    = "https://api.gdc.cancer.gov"
FILES_EP   = f"{GDC_API}/files"
DATA_EP    = f"{GDC_API}/data"
PROJECT    = "TARGET-WT"

CHUNK_SIZE  = 1 << 20   # 1 MB streaming chunks
MAX_RETRIES = 3
RETRY_WAIT  = 10         # seconds between retries

RNASEQ_DIR = PATHS["rnaseq_data_path"]
CLINICAL_DIR = PATHS["clinical_data_path"]
META_OUT = ROOT / "data" / "processed" / "tables" / "sample_metadata.csv"

# ── GDC API helpers ───────────────────────────────────────────────────────────

def _query_files(filters: dict, fields: list[str], size: int = 2000) -> list[dict]:
    """
    Query the GDC /files endpoint and return all matching hits.

    Parameters
    ----------
    filters : GDC filter dict (will be JSON-encoded)
    fields  : list of dot-notation field names to return
    size    : max results (GDC default is 10; set high to get all at once)

    Returns
    -------
    List of hit dicts from data.hits.
    """
    params = {
        "filters": json.dumps(filters),
        "fields":  ",".join(fields),
        "format":  "JSON",
        "size":    str(size),
    }
    try:
        r = requests.get(FILES_EP, params=params, timeout=60)
        r.raise_for_status()
    except requests.RequestException as exc:
        sys.exit(f"\n[error] GDC query failed: {exc}")

    payload    = r.json()
    hits       = payload["data"]["hits"]
    pagination = payload["data"]["pagination"]
    total      = pagination.get("total", len(hits))
    if total > size:
        print(f"  [warn] {total} files found but only {size} fetched; "
              f"increase `size` if files are missing.")
    return hits


def _download_file(
    file_id:   str,
    dest_dir:  Path,
    file_name: str,
    file_size: int = 0,
) -> str:
    """
    Download a single GDC file by UUID into dest_dir/<file_name>.

    Parameters
    ----------
    file_id   : GDC file UUID
    dest_dir  : directory to place the downloaded file (created if absent)
    file_name : destination filename
    file_size : expected size in bytes (0 = skip size check)

    Returns
    -------
    "new"  — freshly downloaded
    "skip" — file already present and size matches
    "fail" — download failed after all retries
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    out_path = dest_dir / file_name

    # Skip if already present with correct size
    if out_path.exists():
        if file_size == 0 or out_path.stat().st_size == file_size:
            return "skip"

    url = f"{DATA_EP}/{file_id}"
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            with requests.get(url, stream=True, timeout=300) as r:
                r.raise_for_status()
                with open(out_path, "wb") as fh:
                    for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                        fh.write(chunk)
            return "new"
        except requests.RequestException as exc:
            if out_path.exists():
                out_path.unlink()           # remove partial download
            if attempt < MAX_RETRIES:
                print(f"    [retry {attempt}/{MAX_RETRIES}] {exc}")
                time.sleep(RETRY_WAIT * attempt)
            else:
                print(f"    [failed] {file_id} — {exc}")
                return "fail"

    return "fail"   # unreachable, satisfies type checker


# ── RNA-seq download ──────────────────────────────────────────────────────────

def _flatten_rnaseq_hit(hit: dict) -> dict:
    """
    Flatten a GDC /files hit into a flat metadata record.

    GDC returns nested JSON; one file → one case → one sample (for bulk RNA-seq).
    Falls back gracefully when fields are absent.
    """
    cases   = hit.get("cases") or [{}]
    case    = cases[0]
    samples = case.get("samples") or [{}]
    sample  = samples[0]
    return {
        # columns expected by 00b_build_counts_matrix_and_annotation.py
        "sample_id":          sample.get("submitter_id", hit["file_id"]),
        "file_id":            hit["file_id"],
        "cases.submitter_id": case.get("submitter_id", ""),
        "cases":              case.get("submitter_id", ""),   # kept by script 00
        "sample_type":        sample.get("sample_type", ""),
        # internal — not written to CSV
        "_file_name":         hit["file_name"],
        "_file_size":         hit.get("file_size", 0),
    }


def download_rnaseq() -> pd.DataFrame:
    """
    Query and download all TARGET-WT RNA-seq STAR count TSVs from GDC.

    Returns
    -------
    DataFrame with columns matching the sample_metadata.csv format expected
    by 00b_build_counts_matrix_and_annotation.py (indexed by sample_id).
    """
    _sep("RNA-seq: querying GDC")

    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {
                "field": "cases.project.project_id", "value": PROJECT}},
            {"op": "=", "content": {
                "field": "data_category", "value": "Transcriptome Profiling"}},
            {"op": "=", "content": {
                "field": "data_type", "value": "Gene Expression Quantification"}},
            {"op": "=", "content": {
                "field": "analysis.workflow_type", "value": "STAR - Counts"}},
            {"op": "=", "content": {
                "field": "access", "value": "open"}},
        ],
    }
    fields = [
        "file_id", "file_name", "file_size",
        "cases.submitter_id",
        "cases.samples.submitter_id",
        "cases.samples.sample_type",
    ]

    hits    = _query_files(filters, fields)
    records = [_flatten_rnaseq_hit(h) for h in hits]
    total   = len(records)
    print(f"  {total} files found")
    print(f"  Destination: {RNASEQ_DIR}\n")

    counts = {"new": 0, "skip": 0, "fail": 0}
    for i, rec in enumerate(records, 1):
        dest   = RNASEQ_DIR / rec["file_id"]
        status = _download_file(
            rec["file_id"], dest, rec["_file_name"], rec["_file_size"]
        )
        counts[status] += 1

        if status == "new":
            mb = rec["_file_size"] / 1e6 if rec["_file_size"] else 0
            print(f"  [{i:3d}/{total}] ✓  {rec['sample_id']:<30s}"
                  f"  {rec['sample_type']:<25s}  {mb:5.1f} MB")
        elif status == "fail":
            print(f"  [{i:3d}/{total}] ✗  {rec['sample_id']:<30s}  FAILED")
        elif i % 25 == 0 or i == total:
            print(f"  [{i:3d}/{total}]    {counts['skip']} files already present, skipping")

    _print_counts(counts)

    meta = (pd.DataFrame(records)
              .drop(columns=["_file_name", "_file_size"])
              .set_index("sample_id"))
    return meta


# ── Clinical download ─────────────────────────────────────────────────────────

def download_clinical() -> None:
    """
    Query and download all TARGET-WT Clinical Supplement XLSX files from GDC.

    The two XLSX files contain:
    - TARGET_WT_ClinicalData_Discovery_*.xlsx   (main clinical table)
    - TARGET_WT_*_Percent_Tumor_Nuclei_*.xlsx   (pathology supplement)

    Both are read by 00b_build_counts_matrix_and_annotation.py via load_xlsx().
    """
    _sep("Clinical: querying GDC")

    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {
                "field": "cases.project.project_id", "value": PROJECT}},
            {"op": "=", "content": {
                "field": "data_category", "value": "Clinical"}},
            {"op": "=", "content": {
                "field": "data_type", "value": "Clinical Supplement"}},
        ],
    }
    hits = _query_files(filters, fields=["file_id", "file_name", "file_size"])

    # Keep only XLSX files (GDC also returns BCR XML and other formats)
    xlsx_hits = [h for h in hits if h["file_name"].lower().endswith(".xlsx")]
    print(f"  {len(hits)} Clinical Supplement files found, "
          f"{len(xlsx_hits)} are XLSX")
    print(f"  Destination: {CLINICAL_DIR}\n")

    if not xlsx_hits:
        print("  [warn] No XLSX files found via query. "
              "Try downloading manually from https://portal.gdc.cancer.gov/")
        return

    counts = {"new": 0, "skip": 0, "fail": 0}
    for h in xlsx_hits:
        dest   = CLINICAL_DIR / h["file_id"]
        status = _download_file(
            h["file_id"], dest, h["file_name"], h.get("file_size", 0)
        )
        counts[status] += 1
        icon = {"new": "✓", "skip": "—", "fail": "✗"}[status]
        label = {"new": "downloaded", "skip": "already present", "fail": "FAILED"}[status]
        print(f"  {icon}  {label:<18s}  {h['file_name']}")

    _print_counts(counts)


# ── Sample metadata ───────────────────────────────────────────────────────────

def save_metadata(meta: pd.DataFrame) -> None:
    """
    Write sample_metadata.csv in the format expected by
    00b_build_counts_matrix_and_annotation.py.

    CSV layout: index = sample_id (sample barcode), then columns
    cases.submitter_id, file_id, cases, sample_type.
    Loaded downstream with pd.read_csv(path, index_col=0).
    """
    META_OUT.parent.mkdir(parents=True, exist_ok=True)
    meta.to_csv(META_OUT)
    print(f"\n  Saved: {META_OUT}")
    print(f"  ({len(meta)} samples)\n")

    # Quick breakdown by sample type
    if "sample_type" in meta.columns:
        counts = meta["sample_type"].value_counts()
        for stype, n in counts.items():
            print(f"    {n:3d}  {stype}")


# ── Utilities ─────────────────────────────────────────────────────────────────

def _sep(title: str = "") -> None:
    print(f"\n{'─' * 62}")
    if title:
        print(f"  {title}")
        print(f"{'─' * 62}")


def _print_counts(counts: dict) -> None:
    print(f"\n  → {counts['new']} downloaded  "
          f"{counts['skip']} skipped  "
          f"{counts['fail']} failed")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    print("=" * 62)
    print("  GDC TARGET-WT Data Download")
    print("=" * 62)
    print(f"  Project root : {ROOT}")
    print(f"  API endpoint : {GDC_API}")

    # 1. RNA-seq
    meta = download_rnaseq()

    # 2. Save metadata before clinical (in case clinical query fails)
    _sep("Sample metadata")
    save_metadata(meta)

    # 3. Clinical XLSX
    download_clinical()

    # Summary
    _sep()
    print("  Download complete.\n")
    print(f"  RNA-seq files : {RNASEQ_DIR}")
    print(f"  Clinical files: {CLINICAL_DIR}")
    print(f"  Metadata      : {META_OUT}")
    print("=" * 62)


if __name__ == "__main__":
    main()
