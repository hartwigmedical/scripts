#!/usr/bin/env python3
"""
Create SUCCESS files in non-empty subfolders of a GCS bucket.
- Skips bucket root and immediate children by default (min_depth=2).

Usage:
  python create_success_files.py gs://your-bucket            # min_depth=2
  python create_success_files.py gs://your-bucket 3          # only depth>=3
"""

import sys
from typing import List, Tuple
import gcsfs

def norm_bucket(uri: str) -> str:
    return uri.replace("gs://", "").strip("/")

def rel_depth(dirpath: str, bucket: str) -> int:
    # e.g. dirpath="bucket/a/b/c" -> rel_parts=["a","b","c"] -> depth=3
    rel = dirpath[len(bucket):].strip("/")
    return 0 if not rel else len(rel.split("/"))

def create_success_markers(bucket_uri: str, min_depth: int = 2) -> None:
    bucket = norm_bucket(bucket_uri)
    fs = gcsfs.GCSFileSystem()

    created = 0
    scanned = 0

    for dirpath, _dirnames, filenames in fs.walk(bucket):
        depth = rel_depth(dirpath, bucket)

        # Enforce depth threshold (skip root and optionally first-level dirs)
        if depth < min_depth:
            continue

        scanned += 1

        # Non-empty means: at least one file other than an existing SUCCESS
        has_real_files = any(name.rsplit("/", 1)[-1] != "SUCCESS" and name != '' for name in filenames)
        if not has_real_files:
            continue

        # Only create if this folder doesn't already contain SUCCESS
        names_here = [f.rsplit("/", 1)[-1] for f in filenames]
        if "SUCCESS" in names_here:
            continue

        success_path = f"{dirpath}/SUCCESS"
        with fs.open(success_path, "w") as f:
            f.write("SUCCESS\n")
        created += 1
        print(f"Created {success_path}")

    print(f"\nScanned {scanned} directories (depth≥{min_depth}). Created {created} SUCCESS files.")

if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        print("Usage: python create_success_files.py gs://your-bucket [min_depth]", file=sys.stderr)
        sys.exit(1)
    bucket_uri = sys.argv[1]
    min_depth = int(sys.argv[2]) if len(sys.argv) == 3 else 2
    create_success_markers(bucket_uri, min_depth)