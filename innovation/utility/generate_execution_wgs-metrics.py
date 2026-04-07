#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_execution_wgs-metrics.py

Create wgs-metrics execution YAML entries from a TSV file with sampleId and gsUri columns.

Input TSV format (tab-separated, with header):
  sampleId	gsUri
  COLO829T	gs://bucket/path/to/file.cram

The gsUri is the full GCS path to the CRAM file. The cram_uri (folder) and
cram_file (basename) are derived from it automatically.

Usage:
  python generate_execution_wgs-metrics.py -i samples.tsv -o execution.yaml
  python generate_execution_wgs-metrics.py -i samples.tsv --workflow wgs-metrics-v2 --version 1.0.0
"""

import argparse
import csv
import sys
from typing import List


def build_yaml_doc(sample_id: str, cram_uri: str, cram_file: str,
                   workflow: str, version: str) -> str:
    return (
        f'name: "{sample_id}"\n'
        f'workflow: "{workflow}"\n'
        f'version: "{version}"\n'
        f'params:\n'
        f'  cram_uri: "{cram_uri}"\n'
        f'  sample_id: "{sample_id}"\n'
        f'  cram_file: "{cram_file}"\n'
    )


def parse_gs_uri(gs_uri: str):
    """Split a full GCS path into (folder, filename)."""
    gs_uri = gs_uri.strip()
    if not gs_uri.startswith("gs://"):
        raise ValueError(f"Expected gs:// URI, got: {gs_uri}")
    last_slash = gs_uri.rfind("/")
    folder = gs_uri[:last_slash + 1]
    filename = gs_uri[last_slash + 1:]
    return folder, filename


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Create wgs-metrics execution YAML from a TSV with sampleId and gsUri columns."
    )
    ap.add_argument("--input", "-i", required=True,
                    help="TSV file with sampleId and gsUri columns.")
    ap.add_argument("--output", "-o", default="wgs-metrics-execution.yaml",
                    help="Output YAML file (default: wgs-metrics-execution.yaml).")
    ap.add_argument("--workflow", default="wgs-metrics-v2",
                    help="Workflow name (default: wgs-metrics-v2).")
    ap.add_argument("--version", default="1.0.0",
                    help="Workflow version (default: 1.0.0).")
    args = ap.parse_args()

    yaml_docs: List[str] = []
    seen = set()

    with open(args.input, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "sampleId" not in reader.fieldnames or "gsUri" not in reader.fieldnames:
            print("ERROR: TSV must have 'sampleId' and 'gsUri' columns.", file=sys.stderr)
            sys.exit(1)

        for row in reader:
            sample_id = row["sampleId"].strip()
            gs_uri = row["gsUri"].strip()
            if not sample_id or not gs_uri:
                continue
            if sample_id in seen:
                print(f"WARNING: Duplicate sampleId '{sample_id}', skipping.", file=sys.stderr)
                continue
            seen.add(sample_id)

            cram_uri, cram_file = parse_gs_uri(gs_uri)
            yaml_docs.append(build_yaml_doc(sample_id, cram_uri, cram_file, args.workflow, args.version))

    if not yaml_docs:
        print("No samples found. Nothing written.", file=sys.stderr)
        sys.exit(2)

    with open(args.output, "w", encoding="utf-8") as out:
        out.write("---\n".join(yaml_docs))

    print(f"Wrote {len(yaml_docs)} execution entries to {args.output}")


if __name__ == "__main__":
    main()
