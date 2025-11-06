#!/usr/bin/env python3
import argparse
from pathlib import PurePosixPath
import sys
import yaml

def require_gs(uri: str) -> None:
    if not uri.startswith("gs://"):
        raise ValueError(f"Path does not start with gs:// -> {uri}")

def split_gs(uri: str):
    """
    Split a gs:// URI into (bucket, path_parts:list[str]).
    Example: gs://bucket/a/b/file.cram -> ("bucket", ["a","b","file.cram"])
    """
    require_gs(uri)
    rest = uri[len("gs://"):]
    if "/" in rest:
        bucket, tail = rest.split("/", 1)
        parts = [p for p in tail.split("/") if p]
    else:
        bucket, parts = rest, []
    return bucket, parts

def gs_parent_dir(uri: str) -> str:
    """
    Return parent directory (no trailing slash).
    gs://bucket/a/b/file.cram -> gs://bucket/a/b
    gs://bucket/file.cram -> gs://bucket
    """
    bucket, parts = split_gs(uri)
    if not parts:
        raise ValueError(f"No object path in URI (no file): {uri}")
    parent_parts = parts[:-1]
    if parent_parts:
        return f"gs://{bucket}/" + "/".join(parent_parts)
    else:
        return f"gs://{bucket}"

def gs_filename(uri: str) -> str:
    bucket, parts = split_gs(uri)
    if not parts:
        raise ValueError(f"No object filename in URI: {uri}")
    return parts[-1]

def sample_id_from_filename(filename: str) -> str:
    """
    Derive sample_id from filename (without extension).
    Rule: take the stem, then take everything up to first underscore (if present).
    Works for 'SAMPLE123.cram', 'SAMPLE123-ref.cram', 'SAMPLE123_xyz.cram' -> SAMPLE123 or SAMPLE123-ref.
    """
    p = PurePosixPath(filename)
    stem = p.stem
    return stem.split("_")[0]

def make_entry(uri: str, workflow: str, version: str) -> dict:
    input_file = gs_filename(uri)
    input_uri = gs_parent_dir(uri)  # NO trailing slash
    sample_id = sample_id_from_filename(input_file)
    return {
        "name": sample_id,
        "workflow": workflow,
        "version": version,
        "params": {
            "input_uri": input_uri,
            "sample_id": sample_id,
            "input_file": input_file,
        },
    }

def main():
    ap = argparse.ArgumentParser(description="Generate execution YAML entries (separated by ---) from GCS CRAM/BAM paths.")
    ap.add_argument("input_txt", help="Text file with one gs:// path per line (to a .cram or .bam).")
    ap.add_argument("--workflow", default="flagstat", help="Workflow name (default: flagstat)")
    ap.add_argument("--version", default="0.0.1", help="Workflow version (default: 0.0.1)")
    ap.add_argument("--output", default="executions.yaml", help="Output YAML file (default: executions.yaml)")
    args = ap.parse_args()

    uris = []
    with open(args.input_txt) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            require_gs(s)
            fn = gs_filename(s).lower()
            if not (fn.endswith(".cram") or fn.endswith(".bam")):
                print(f"[WARN] Skipping non-CRAM/BAM: {s}", file=sys.stderr)
                continue
            uris.append(s)

    if not uris:
        raise SystemExit("No valid gs:// CRAM/BAM paths found.")

    entries = [make_entry(u, args.workflow, args.version) for u in uris]

    # Write as separate YAML documents separated by '---', no list dashes.
    with open(args.output, "w") as out:
        for i, doc in enumerate(entries):
            if i > 0:
                out.write("---\n")
            yaml.dump(doc, out, sort_keys=False)

    print(f"Wrote {len(entries)} YAML entries to {args.output}")

if __name__ == "__main__":
    main()
