#!/usr/bin/env python3
"""
build_execution_yaml.py
Create an execution.yaml for every sample that has an …/aligner/ folder.

USAGE
  python build_execution_yaml.py gs://input-bucket gs://resources [-o execution.yaml]

NEEDS
  • Google Cloud SDK ≥ 473 (for `gcloud storage`)
  • PyYAML  (pip install pyyaml)
"""
import argparse, pathlib, re, subprocess, sys, yaml


def list_aligner_dirs(bucket: str) -> list[str]:
    """
    Return every unique gs://…/aligner/ prefix inside *bucket*.
    Uses `**` wildcard support in gcloud storage ls (recursive). :contentReference[oaicite:0]{index=0}
    """
    pattern = f"{bucket.rstrip('/')}/**/aligner/**"
    cmd = [
        "gcloud",
        "storage",
        "ls",
        "--recursive",
        pattern,
    ]  # example ls syntax :contentReference[oaicite:1]{index=1}

    try:
        output = subprocess.check_output(cmd, text=True)
    except subprocess.CalledProcessError as e:
        sys.exit(e.stderr or str(e))

    # Collapse object URLs to their parent …/aligner/ “directory”.
    dirs = {
        m.group(1)
        for line in output.splitlines()
        if (m := re.match(r"^(gs://.*/aligner/)", line))
    }
    return sorted(dirs)


def make_docs(aligner_dirs: list[str], resources_uri: str):
    """Yield one execution-yaml doc per sample."""
    for d in aligner_dirs:
        sample_id = pathlib.PurePosixPath(d.replace("gs://", "")).parts[-2]
        yield {
            "name": sample_id,
            "workflow": "rerun-sage",
            "version": "0.0.1",
            "params": {
                "run_uri": d,
                "resources_uri": resources_uri,
                "sample_id": sample_id,
            },
        }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input_bucket", help="gs://input-bucket")
    p.add_argument("resources_bucket", help="gs://resources")
    p.add_argument("-o", "--output", default="execution.yaml")
    args = p.parse_args()

    aligner_dirs = list_aligner_dirs(args.input_bucket)
    if not aligner_dirs:
        sys.exit("nothing with an aligner/ folder found")

    with open(args.output, "w") as fh:
        yaml.safe_dump_all(
            make_docs(aligner_dirs, args.resources_bucket),
            fh,
            explicit_start=True,
            sort_keys=False,
        )

    print(f"wrote {len(aligner_dirs)} entries → {args.output}")


if __name__ == "__main__":
    main()