#!/usr/bin/env python3
# build_yaml_from_batch.py – fastq tumor-only YAML builder for pipeline v6
# Python 3.7.3

import argparse, json, os, subprocess, sys, yaml
from typing import List, Dict

TARGET_BUCKET = "gs://targeted-pipeline-output-prod-1/"
FASTQ_BUCKET  = "gs://oncopanel-prod-fastq-backup/"
OPTIONS = {  # unchanged
    "gcp": {"project": "hmf-crunch", "region": "europe-west4"},
    "serviceAccount": {
        "kubernetesServiceAccount": "hmf-crunch-sa",
        "gcpEmailAddress": "hmf-crunch@hmf-crunch.iam.gserviceaccount.com",
    },
    "cluster": "research-cluster-prod-1",
    "cmek": ("projects/hmf-database/locations/europe-west4/keyRings/hmf-database/"
             "cryptoKeys/hmf-database-20191001"),
    "image": ("europe-west4-docker.pkg.dev/hmf-build/hmf-docker/"
              "pipeline5:6.0.8-beta.5"),
    "outputBucket": "panel-v6-tmp-diag",
    "batch": {"size": "8"},
    "argumentOverrides": {
        "image_project": "hmf-pipeline-prod-e45b00f2",
        "image_name": "pipeline5-6-0-202504281206-202504281437-private",
        "ref_genome_version": 38,
        "use_target_regions": True,
        "cost_center_label": "innovation",
    },
}

class Dumper(yaml.Dumper):
    def increase_indent(self, flow=False, indentless=False):
        return super().increase_indent(flow, False)

def sh(cmd: List[str]) -> str:
    return subprocess.check_output(cmd, text=True)

# ───── GCS helpers ────────────────────────────────────────────────
def root_dirs(bucket: str) -> List[str]:
    """List ONLY top-level folders in a bucket (no recursion)."""
    return [ln.strip() for ln in sh(["gcloud", "storage", "ls", bucket]).splitlines()
            if ln.strip().endswith("/")]

def gcs_ls_recursive(prefix: str) -> List[str]:
    return [ln.strip() for ln in sh(["gcloud", "storage", "ls", "--recursive", prefix]).splitlines()
            if ln.strip() and not ln.startswith("TOTAL:")]

# ───── core logic ────────────────────────────────────────────────
def extract_barcodes(batch: str) -> List[str]:
    """grep top-level dirs by batch date, take third field as barcode."""
    dirs = [d for d in root_dirs(TARGET_BUCKET)
            if os.path.basename(d.rstrip("/")).startswith(batch)]
    bars = [os.path.basename(d.rstrip("/")).split("_")[2]
            for d in dirs
            if len(os.path.basename(d.rstrip("/")).split("_")) >= 3]
    return sorted(set(bars))

def fastqs_for_barcode(bar: str) -> List[str]:
    return [f for f in gcs_ls_recursive(f"{FASTQ_BUCKET}{bar}/") if f.endswith(".fastq.gz")]

def pair_fastqs(fqs: List[str]) -> List[Dict[str, str]]:
    fqs = sorted(fqs)
    return [{"read1": fqs[i], "read2": fqs[i + 1]}
            for i in range(0, len(fqs) - 1, 2)
            if "_R1_" in fqs[i] and "_R2_" in fqs[i + 1]]

def reporting_id(bar: str) -> str:
    raw = sh(["api", "-j", "samples", f"barcode={bar}"])
    return json.loads(raw)[0]["reporting_id"]

# ───── main ───────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-b", "--batch", required=True, help="batch date, e.g. 250606")
    ap.add_argument("-o", "--output", required=True, help="output YAML file")
    args = ap.parse_args()

    bars = extract_barcodes(args.batch)
    if not bars:
        sys.exit(f"No barcodes found for batch {args.batch}")

    samples = []
    for bar in bars:
        fq = fastqs_for_barcode(bar)
        if len(fq) > 8:
            print(f"WARNING: {bar} has {len(fq)} fastqs (expected 8)", file=sys.stderr)
        rid = reporting_id(bar)
        samples.append({"name": rid, "tumors": [{"name": rid, "fastq": pair_fastqs(fq)}]})

    with open(args.output, "w") as out:
        yaml.dump({"options": OPTIONS, "samples": samples},
                  out, Dumper=Dumper, sort_keys=False, default_flow_style=False)
    print("YAML written →", args.output)

if __name__ == "__main__":
    main()