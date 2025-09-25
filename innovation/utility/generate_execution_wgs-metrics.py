#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_wgs_metrics_execution.py

Create 'wgs-metrics-execution.yaml' entries for each sample found under
the *immediate subfolders* of each GCS folder listed in an input file.

Input file format: one GCS folder per line, e.g.:
  gs://my-bucket/parent-folder/

Behavior:
- For each line (parent folder), list its immediate subfolders.
- In each subfolder, find *.cram or *.bam files (non-recursive).
- Prefer .cram over .bam when both exist for the same sample id.
- Infer sample_id from filename: basename -> remove extension -> split("_")[0]
- Produce one YAML document per sample:
    name: "<sample_id>"
    workflow: "wgs-metrics"
    version: "0.1.0"
    params:
      bam_uri: "gs://bucket/path/to/subfolder/"
      sample_id: "<sample_id>"
      bam_file: "<filename>"

Auth:
  gcloud auth application-default login
Requires:
  pip install gcsfs
"""

import argparse
import sys
import os
import re
from typing import Dict, List, Tuple, Optional

try:
    import gcsfs
except ImportError:
    print("ERROR: gcsfs is required. Install with: pip install gcsfs", file=sys.stderr)
    sys.exit(1)


def normalize_gs_folder(gs_uri: str) -> str:
    """Ensure a GCS folder path like 'gs://bucket/prefix/' with a single trailing slash."""
    gs_uri = gs_uri.strip()
    if not gs_uri.startswith("gs://"):
        raise ValueError("GCS URI must start with gs:// -> got: {}".format(gs_uri))
    # Remove extra slashes at end, then add single slash
    return gs_uri.rstrip("/") + "/"


def path_join_gs(folder: str, name: str) -> str:
    """Join gs folder and a child name (no leading slash in name)."""
    if not folder.endswith("/"):
        folder += "/"
    return folder + name


def basename(path: str) -> str:
    """Return last path component from a gs:// path."""
    return path.rstrip("/").split("/")[-1]


def infer_sample_id(filename: str) -> str:
    """
    From filename (e.g., 'SAMPLE_XX.cram' or 'SAMPLE.cram'):
      - strip extension
      - split on '_' and take first token
    """
    base = re.sub(r"\.(cram|bam)$", "", filename, flags=re.IGNORECASE)
    return base.split("_")[0]


def list_immediate_children(fs: gcsfs.GCSFileSystem, folder: str) -> List[str]:
    """
    List immediate children (both dirs and files) of a gs:// folder.
    Returns paths with 'gs://' prefix where possible.
    """
    # gcsfs can list 'gs://bucket/prefix'; entries are typically like 'bucket/prefix/child'
    entries = fs.ls(folder)
    out = []
    for e in entries:
        # Normalize to 'gs://'
        if e.startswith("gs://"):
            out.append(e)
        else:
            out.append("gs://" + e)
    return out


def is_dir(fs: gcsfs.GCSFileSystem, path: str) -> bool:
    try:
        return fs.isdir(path)
    except FileNotFoundError:
        return False


def build_yaml_doc(name: str, workflow: str, version: str,
                   bam_uri: str, sample_id: str, bam_file: str) -> str:
    # Ensure bam_uri ends with a slash in YAML
    bam_uri = normalize_gs_folder(bam_uri)
    return (
        'name: "{name}"\n'
        'workflow: "{workflow}"\n'
        'version: "{version}"\n'
        'params:\n'
        '  bam_uri: "{bam_uri}"\n'
        '  sample_id: "{sample_id}"\n'
        '  bam_file: "{bam_file}"\n'
        '---\n'
    ).format(
        name=name,
        workflow=workflow,
        version=version,
        bam_uri=bam_uri,
        sample_id=sample_id,
        bam_file=bam_file,
    )


def find_candidate_files(fs: gcsfs.GCSFileSystem, subfolder: str) -> List[str]:
    """List .cram/.bam files directly under subfolder (non-recursive)."""
    children = list_immediate_children(fs, subfolder)
    files = []
    for p in children:
        if not p.endswith("/"):  # files won't end with slash
            if re.search(r"\.(cram|bam)$", p, flags=re.IGNORECASE):
                files.append(p)
    return files


def choose_by_preference(files: List[str]) -> Dict[str, Tuple[str, str]]:
    """
    From a list of gs://.../(file), group by sample_id and prefer .cram over .bam.
    Returns dict: sample_id -> (bam_uri, bam_file)
    """
    best: Dict[str, Tuple[str, str]] = {}
    for fullpath in files:
        file_name = basename(fullpath)
        sample_id = infer_sample_id(file_name)
        # bam_uri is folder containing the file
        bam_uri = fullpath[: -(len(file_name))]
        prefer = file_name.lower().endswith(".cram")
        if sample_id not in best:
            best[sample_id] = (bam_uri, file_name)
        else:
            # Prefer .cram if available
            already = best[sample_id][1].lower()
            if prefer and already.endswith(".bam"):
                best[sample_id] = (bam_uri, file_name)
    return best


def main() -> None:
    ap = argparse.ArgumentParser(description="Create wgs-metrics-execution.yaml from GCS folders.")
    ap.add_argument("--input", "-i", required=True,
                    help="Text file with one GCS folder per line (e.g., gs://bucket/path/).")
    ap.add_argument("--output", "-o", default="wgs-metrics-execution.yaml",
                    help="Output YAML file (default: wgs-metrics-execution.yaml).")
    ap.add_argument("--workflow", default="wgs-metrics", help="Workflow name (default: wgs-metrics).")
    ap.add_argument("--version", default="0.1.0", help="Workflow version (default: 0.1.0).")
    args = ap.parse_args()

    fs = gcsfs.GCSFileSystem()

    # Read input folders
    with open(args.input, "r", encoding="utf-8") as fh:
        folders = [line.strip() for line in fh if line.strip()]

    yaml_docs: List[str] = []
    total_samples = 0

    for parent in folders:
        parent = normalize_gs_folder(parent)
        try:
            children = list_immediate_children(fs, parent)
        except FileNotFoundError:
            print(f"WARNING: Not found or no access: {parent}", file=sys.stderr)
            continue

        # Only subfolders
        subfolders = [c for c in children if c.endswith("/") and is_dir(fs, c)]
        if not subfolders:
            print(f"INFO: No subfolders under {parent}", file=sys.stderr)

        for sub in subfolders:
            candidates = find_candidate_files(fs, sub)
            if not candidates:
                # Silent skip if a subfolder has no .bam/.cram
                continue

            # Prefer .cram over .bam per sample_id
            chosen = choose_by_preference(candidates)

            for sample_id, (bam_uri, bam_file) in chosen.items():
                doc = build_yaml_doc(
                    name=sample_id,
                    workflow=args.workflow,
                    version=args.version,
                    bam_uri=bam_uri,
                    sample_id=sample_id,
                    bam_file=bam_file,
                )
                yaml_docs.append(doc)
                total_samples += 1

    if not yaml_docs:
        print("No samples found. Nothing written.", file=sys.stderr)
        sys.exit(2)

    with open(args.output, "w", encoding="utf-8") as out:
        out.writelines(yaml_docs)

    print(f"Wrote {len(yaml_docs)} documents to {args.output}")


if __name__ == "__main__":
    main()
