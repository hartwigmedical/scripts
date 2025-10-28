#!/usr/bin/env python3
"""
make_downsample_execs.py

Reads a TSV with:
    alignment_uri   src_cov   tgt_cov

Outputs separate YAML execution blocks like:
    name: SAMPLE0123_ds
    workflow: downsample-crams
    version: 1.0.0-beta1
    params:
      alignment_file_path: "gs://bucket/run42"
      alignment_file: "SAMPLE0123.cram"
      sample_id: "SAMPLE0123"
      src_cov: "108"
      tgt_cov: "83"
      seed: "42"
    ---
"""

import sys
import csv
import os
import yaml  # pip install pyyaml

WORKFLOW_NAME = "downsample-crams"
WORKFLOW_VERSION = "1.0.0"
DEFAULT_SEED = "42"


def split_alignment_uri(alignment_uri: str):
    """Return (dir_path, filename, sample_id)"""
    filename = os.path.basename(alignment_uri.strip())
    dir_path = alignment_uri.strip()[: -(len(filename))].rstrip("/")
    if dir_path == "":
        dir_path = "."
    sample_id = filename.split(".")[0]
    return dir_path, filename, sample_id


def load_rows(tsv_path: str):
    with open(tsv_path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"alignment_uri", "src_cov", "tgt_cov"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            sys.stderr.write(
                f"ERROR: Missing required TSV columns: {', '.join(sorted(missing))}\n"
            )
            sys.exit(1)
        for row in reader:
            if not row["alignment_uri"].strip():
                continue
            yield {
                "alignment_uri": row["alignment_uri"].strip(),
                "src_cov": row["src_cov"].strip(),
                "tgt_cov": row["tgt_cov"].strip(),
            }


def make_execution(row: dict) -> dict:
    dir_path, filename, sample_id = split_alignment_uri(row["alignment_uri"])
    exec_name = f"{sample_id}_ds"

    return {
        "name": exec_name,
        "workflow": WORKFLOW_NAME,
        "version": WORKFLOW_VERSION,
        "params": {
            "alignment_file_path": dir_path,
            "alignment_file": filename,
            "sample_id": sample_id,
            "src_cov": row["src_cov"],
            "tgt_cov": row["tgt_cov"],
            "seed": DEFAULT_SEED,
        },
    }


def main():
    if len(sys.argv) != 2:
        sys.stderr.write(f"Usage: {sys.argv[0]} input.tsv > executions.yaml\n")
        sys.exit(1)

    tsv_path = sys.argv[1]
    rows = list(load_rows(tsv_path))

    for i, row in enumerate(rows):
        exec_block = make_execution(row)
        yaml.safe_dump(exec_block, stream=sys.stdout, sort_keys=False)
        # add separator after each block
        sys.stdout.write("---\n")


if __name__ == "__main__":
    main()
