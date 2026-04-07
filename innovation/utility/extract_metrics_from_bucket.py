#!/usr/bin/env python3
"""
Extract WGS metrics from a GCS bucket into a single TSV file.

Reads bam-metrics summary TSV and GC bias summary (Picard format)
for every sample folder in the bucket.

Bucket structure expected:
    {bucket}/{sample_id}/bam-metrics/*.bam_metric.summary.tsv
    {bucket}/{sample_id}/gc-bias/*.gc_bias.summary_metrics

Usage:
    python extract_wgs_metrics.py gs://my-bucket
    python extract_wgs_metrics.py gs://my-bucket -o my_metrics.tsv
"""

from __future__ import annotations

import argparse
import re
import sys
from io import StringIO
from typing import Any, Dict, List, Optional

import pandas as pd
from gcsfs import GCSFileSystem


# ========= GCS I/O =========
def find_sample_folders(gcs: GCSFileSystem, bucket_path: str) -> List[str]:
    try:
        all_paths = gcs.ls(bucket_path, detail=False)
        return [p.rsplit("/", 1)[-1] for p in all_paths if gcs.isdir(p)]
    except Exception:
        return []


# ========= Parsing =========
def parse_bam_metrics_tsv(gcs: GCSFileSystem, file_path: str) -> Optional[pd.Series]:
    """Read a plain TSV with a header row and a single data row."""
    try:
        with gcs.open(file_path, "r") as f:
            df = pd.read_csv(f, sep="\t", nrows=1)
        if df.empty:
            return None
        df.columns = df.columns.astype(str).str.strip()
        s = df.iloc[0]
        s.index = s.index.astype(str).str.strip()
        return s
    except Exception:
        return None


def parse_picard_single_metrics(gcs: GCSFileSystem, file_path: str) -> Optional[pd.Series]:
    """Parse a Picard-style metrics file (## METRICS CLASS header)."""
    try:
        with gcs.open(file_path, "r") as f:
            content = f.read()
        m = content.find("## METRICS CLASS")
        if m == -1:
            return None
        block = content[m:]
        df = pd.read_csv(StringIO(block), sep="\t", comment="#", nrows=1)
        if df.empty:
            return None
        df.columns = df.columns.astype(str).str.strip()
        s = df.iloc[0]
        s.index = s.index.astype(str).str.strip()
        return s
    except Exception:
        return None


# ========= Role classification =========
def classify_role(sample_id: str) -> str:
    sid = re.sub(r"-(sbx|nseq6)$", "", str(sample_id).strip(), flags=re.I)
    if re.search(r"(?i)(?:^|[-_])ref(?:$|[-_])", sid):
        return "Normal"
    if re.match(r"(?i)^h\d{5,}", sid):
        return "Tumor"
    if re.search(r"(?i)T(?:\d+|I+)?$", sid):
        return "Tumor"
    if re.search(r"(?i)R(?:\d+|I+)?$", sid):
        return "Normal"
    return "Unknown"


# ========= Main extraction =========
def extract_metrics(bucket: str) -> pd.DataFrame:
    gcs = GCSFileSystem()
    if not gcs.exists(bucket):
        print(f"Error: bucket not found: {bucket}", file=sys.stderr)
        sys.exit(1)

    samples = sorted(find_sample_folders(gcs, bucket))
    if not samples:
        print(f"Error: no sample folders found in {bucket}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(samples)} samples in {bucket}", file=sys.stderr)

    rows: List[pd.Series] = []
    for i, sample_id in enumerate(samples, start=1):
        print(f"  [{i}/{len(samples)}] {sample_id}", file=sys.stderr, flush=True)
        base = f"{bucket.rstrip('/')}/{sample_id}"

        merged: Dict[str, Any] = {}

        # --- Bam metrics (plain TSV) ---
        bam_files = gcs.glob(f"{base}/bam-metrics/*.bam_metric.summary.tsv")
        if bam_files:
            bam_row = parse_bam_metrics_tsv(gcs, bam_files[0])
            if bam_row is not None:
                merged.update({str(k).strip(): v for k, v in bam_row.to_dict().items()})

        # --- GC bias summary metrics (Picard format) ---
        gc_files = gcs.glob(f"{base}/gc-bias/*.gc_bias.summary_metrics")
        if gc_files:
            gc_row = parse_picard_single_metrics(gcs, gc_files[0])
            if gc_row is not None:
                for col in ("AT_DROPOUT", "GC_DROPOUT"):
                    if col in gc_row.index:
                        merged[col] = gc_row[col]

        if merged:
            merged["SampleID"] = sample_id
            merged["SampleType"] = classify_role(sample_id)
            rows.append(pd.Series(merged))

    if not rows:
        print("Error: no metrics data could be loaded for any sample.", file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(rows)
    df.columns = df.columns.astype(str).str.strip()

    # Put SampleID and SampleType first, then all other columns alphabetically
    fixed = ["SampleID", "SampleType"]
    other_cols = sorted(c for c in df.columns if c not in fixed)
    df = df[fixed + other_cols].sort_values("SampleID").reset_index(drop=True)
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Extract WGS metrics from a GCS bucket into a TSV file."
    )
    parser.add_argument("bucket", help="GCS bucket path (e.g. gs://my-bucket)")
    parser.add_argument("-o", "--output", default=None,
                        help="Output TSV file path (default: wgs_metrics_<bucket_name>.tsv)")
    args = parser.parse_args()

    bucket = args.bucket.rstrip("/")
    if not bucket.startswith("gs://"):
        print("Error: bucket must start with gs://", file=sys.stderr)
        sys.exit(1)

    if args.output:
        output_path = args.output
    else:
        bucket_name = bucket.replace("gs://", "").replace("/", "_")
        output_path = f"wgs_metrics_{bucket_name}.tsv"

    df = extract_metrics(bucket)

    df.to_csv(output_path, sep="\t", index=False)
    print(f"Wrote {len(df)} samples to {output_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
