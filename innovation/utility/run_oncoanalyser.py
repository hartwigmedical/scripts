#!/usr/bin/env python3
"""
run_oncoanalyser.py

Launch one right-sized GCP VM per sample, run nf-core/oncoanalyser on it with the
local Nextflow executor, sync results to an output bucket, and self-delete the VM.

Model
-----
- You provide a *manifest* CSV: a standard oncoanalyser samplesheet plus optional helper
  columns (`coverage`, `vm_type`, `sequencing_type`) that are stripped before the real
  samplesheet is written. `coverage` sizes the WORK DISK only; the machine is fixed by
  workload shape (see MACHINE SIZING). `vm_type` (or --machine-type) pins a specific machine.
- Rows are grouped by `group_id`; each group_id == one analysis == one VM.
- For each group the script writes samplesheet.csv (helper cols removed) + resources.config to
  gs://<staging-bucket>/<staging-prefix>/<run-id>/<group_id>/, renders a startup script, and
  creates a VM (Ubuntu 22.04) that runs the whole thing and then deletes itself.
- Status is reported by marker objects in the output bucket:
       gs://<output-bucket>/<runs-prefix>/<run-id>/<group_id>/_SUCCESS   or   _FAILED
  plus startup.log and nextflow.log for debugging.

Requirements (run side)
-----------------------
- google-cloud-cli (`gcloud`) installed and authenticated; Compute Engine API enabled.
- A service account the VMs run as, with at minimum:
    roles/storage.objectAdmin   on the input + output + staging buckets
    roles/compute.instanceAdmin.v1  (so the VM can delete itself)
- Reference data: either pre-stage it once and pass --reference-config gs://.../refdata.config
  (recommended), or omit it and let oncoanalyser auto-stage every run (slower, re-downloads).

Manifest example (tumor-only WGS from BAM)
------------------------------------------
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath,coverage
FFPE1,FFPE1,FFPE1-T,tumor,dna,bam,gs://my-in/FFPE1-T.bam,34
FFPE2,FFPE2,FFPE2-T,tumor,dna,bam,gs://my-in/FFPE2-T.bam,22

FASTQ / RNA inputs
------------------
FASTQ rows need an `info` column (library_id + lane, ';'-separated) and the paired R1;R2 paths
in a single `filepath` field, also ';'-separated. RNA is just sequence_type=rna (tumor-only is
fine) and runs under the default --mode wgts. Mixed DNA+RNA share a group_id. Only gzip'd,
non-interleaved paired-end FASTQ is supported. Index files are auto-fetched for BAM/CRAM only.

group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath,coverage
P1,P1,P1-T,tumor,dna,fastq,library_id:S1;lane:001,gs://my-in/P1-T_L001_R1.fastq.gz;gs://my-in/P1-T_L001_R2.fastq.gz,34
P1,P1,P1-RNA,tumor,rna,fastq,library_id:S1;lane:001,gs://my-in/P1-RNA_L001_R1.fastq.gz;gs://my-in/P1-RNA_L001_R2.fastq.gz,
"""

import argparse
import csv
import datetime as dt
import json
import os
import re
import secrets
import subprocess
import sys
import tempfile
import time
import urllib.request
from collections import OrderedDict

# --------------------------------------------------------------------------------------
# MACHINE SIZING -- by workload shape, not by coverage (coverage only picks a DISK size).
#
# oncoanalyser's nf-core labels cap any single task at 12 cpu / 72 GB (`process_high`), and
# Nextflow's local executor admits a task only if cpus AND memory fit, so a 72 GB ask lets the
# box run floor(RAM/72) heavy tasks at once. CPU utilisation therefore peaks at ~50-60% on any
# standard (4 GB/vCPU) shape, no matter how many cores you buy: memory, not cores, was the
# binding constraint. Measured over 141 runs, coverage-normalised: n2d-standard-32 -> -48 gained
# ~26% (1 -> 2 heavy slots) while -48 -> -64 gained nothing and cost 33% more.
#
# The fix is write_resources_config(): ask for the memory each tool actually uses (from
# execution_trace peak_rss) and give critical-path tools the thread counts HMF's own pipeline5
# gives them. Once the asks are honest the box becomes CPU-bound -- the correct state -- so one
# machine serves every from-BAM/CRAM run and highmem shapes are pure waste.
#
#   BAM/CRAM, DNA only: 48 cores. The post-REDUX burst needs ~830 core-min of side work while
#     ESVEE (the long pole) runs; at 48 cores that lands almost exactly in ESVEE's shadow.
#   BAM/CRAM, WITH RNA: 64 cores. The whole from-BAM RNA arm is one process, ISOFOX, whose only
#     data dependency is the RNA BAM -- so it is ready at t=0, concurrent with REDUX. It still
#     lost that race because the REDUX pair claimed every core, then landed in the post-REDUX
#     burst and overflowed ESVEE's shadow. ISOFOX_RESERVE_CPUS is held back from the REDUX split
#     and the box grows by exactly that much, so REDUX keeps its measured-optimal 36/12 and
#     ISOFOX runs alongside it. (The "48 -> 64 bought nothing" result above was measured on
#     DNA-only runs, where the box had no admissible work left; it does not transfer here.)
#   FASTQ: 80 cores. Alignment is many independent BWAMEM2_ALIGN shards that genuinely scale out.
#
# N4D = AMD EPYC Turin (5th gen) vs N2D's Milan (3rd gen): cheaper per hour AND faster per
# core. It requires Hyperdisk (no pd-*), which also lifts the disk ceiling from N2D+PD's
# ~1200 MiB/s to 2400 MiB/s.
# --------------------------------------------------------------------------------------
DEFAULT_MACHINE = ("n4d-standard-48", 48, 192)   # (machine_type, cpus, mem_gb) -- BAM/CRAM groups
RNA_MACHINE     = ("n4d-standard-64", 64, 256)   # BAM/CRAM + an rna row
FASTQ_MACHINE   = ("n4d-standard-80", 80, 320)   # groups with a fastq row

# Cores held back from the REDUX split for ISOFOX on RNA-bearing groups. Deliberately the same
# absolute count ISOFOX already gets on a 48-core box (plan_cpus frac(0.33) -> 16): this adds
# cores so ISOFOX no longer waits, it does not re-size ISOFOX itself. Absolute rather than a
# fraction so the reserve -- and therefore the REDUX split -- is stable across box sizes.
#
# 16 is NOT yet measured: ISOFOX is absent from PROCESS_CPUS_MEASURED, so its thread count has
# never been checked against real %cpu draw (VIRUSBREAKEND is the cautionary case, 4.3 cores
# from a 12-core grant). Take %cpu off the first WGTS execution_trace; if it draws ~6, drop this
# to 8 and move ISOFOX into PROCESS_CPUS_MEASURED.
ISOFOX_RESERVE_CPUS = 16

# Persistent WORK disk (/mnt/work: refdata staging, localized inputs, the Nextflow work dir,
# Docker storage, outputs before upload). This IS genuinely coverage-driven -- it must hold the
# inputs plus every intermediate BAM. coverage <= max_cov picks the row, last row is catch-all.
DISK_TIERS = [
    # max_cov, disk_gb
    (10,      500),
    (50,      750),
    (95,     1200),
    (150,    2000),
    (10_000, 3000),
]
DEFAULT_DISK_INDEX = 1   # BAM/CRAM group with no coverage given

# FASTQ groups carry localized FASTQ + per-lane BAMs + merged BAM + REDUX BAM simultaneously.
# NOTE: disk is NOT the alignment bottleneck -- a 190x sample succeeded on 4 TB (H00005238) and
# failed on 6 TB (H00006429); that failure mode was RAM (see the BWAMEM2_ALIGN notes).
FASTQ_DISK_TIERS = [
    (40,     1500),
    (95,     2500),
    (10_000, 4000),
]
FASTQ_DEFAULT_DISK_INDEX = 0

# Boot disk size (GB). Ephemeral, and holds ONLY the OS + apt packages: all heavy state lives on
# the separate persistent WORK disk mounted at /mnt/work. That disk has auto-delete OFF so it
# survives a Spot preemption -- the relaunch re-attaches it and `nextflow -resume` continues
# from the last completed step instead of re-aligning from zero.
BOOT_DISK_GB = 100

# --------------------------------------------------------------------------------------
# Cost estimation. Prices are read live from the Cloud Billing catalog for the run's region, so
# they cannot go stale: cpus*<Core SKU> + mem_gb*<Ram SKU> for the VM, plus the disk capacity
# and (for Hyperdisk) provisioned IOPS/throughput SKUs. Spot uses the real
# 'Spot Preemptible ...' SKUs rather than a fixed discount factor.
#
# Still only a burn-rate estimate, NOT invoice-accurate: it ignores sustained-use discounts
# (which make the real bill LOWER), network egress and GCS storage. The source of truth is the
# billing export filtered by the labels the script attaches (pipeline=oncoanalyser,
# run=<run-id>, group=<group>). If the catalog is unreachable the estimate is simply omitted.
# --------------------------------------------------------------------------------------
COMPUTE_SERVICE = "services/6F81-5844-456A"      # Compute Engine, in the Cloud Billing catalog
CURRENCY = "EUR"
HOURS_PER_MONTH = 730
# Hyperdisk Balanced bundles 3000 IOPS + 140 MiB/s; only the excess over that is metered.
HYPERDISK_FREE_IOPS = 3000
HYPERDISK_FREE_THROUGHPUT = 140
# disk type -> (capacity, IOPS, throughput) SKU description prefixes. pd-* bills capacity only;
# its performance is a fixed function of size, so there is nothing else to price.
DISK_SKU_PREFIX = {
    "hyperdisk-balanced": ("Hyperdisk Balanced Capacity", "Hyperdisk Balanced IOPS",
                           "Hyperdisk Balanced Throughput"),
    "hyperdisk-extreme":  ("Hyperdisk Extreme Capacity", "Hyperdisk Extreme IOPS", None),
    "pd-balanced": ("Balanced PD Capacity", None, None),
    "pd-ssd":      ("SSD backed PD Capacity", None, None),
    "pd-standard": ("Storage PD Capacity", None, None),
}

_SKU_CACHE = {}
_SHAPE_CACHE = {}


def region_of(zone):
    """europe-west4-a -> europe-west4"""
    return zone.rsplit("-", 1)[0]


def billing_skus(region):
    """Every Compute Engine SKU that applies to `region`, from the Cloud Billing catalog.

    Cached per region. Returns [] on any failure (no auth, API disabled, no network): a cost
    estimate is a nice-to-have and must never block a launch."""
    if region in _SKU_CACHE:
        return _SKU_CACHE[region]
    _SKU_CACHE[region] = []
    tok = subprocess.run(["gcloud", "auth", "print-access-token"], capture_output=True, text=True)
    if tok.returncode != 0:
        print("  (warning) no gcloud access token; cost estimates disabled")
        return []
    skus, page = [], ""
    while True:
        url = (f"https://cloudbilling.googleapis.com/v1/{COMPUTE_SERVICE}/skus"
               f"?currencyCode={CURRENCY}&pageSize=5000" + (f"&pageToken={page}" if page else ""))
        req = urllib.request.Request(url, headers={"Authorization": f"Bearer {tok.stdout.strip()}"})
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                data = json.load(resp)
        except Exception as exc:                    # any failure -> no estimate, never an abort
            print(f"  (warning) could not read the Cloud Billing catalog ({exc}); "
                  f"cost estimates disabled")
            return []
        skus += [s for s in data.get("skus", []) if region in (s.get("serviceRegions") or [])]
        page = data.get("nextPageToken") or ""
        if not page:
            break
    _SKU_CACHE[region] = skus
    return skus


def sku_price(sku):
    """The SKU's unit price (last tier), in CURRENCY per its usage unit."""
    tiers = ((sku.get("pricingInfo") or [{}])[0].get("pricingExpression") or {}).get("tieredRates")
    if not tiers:
        return None
    unit = tiers[-1].get("unitPrice") or {}
    return int(unit.get("units") or 0) + int(unit.get("nanos") or 0) / 1e9


def machine_hourly(machine_type, cpus, mem_gb, region, spot):
    """EUR/hour for the VM itself: cpus*<Core SKU> + mem_gb*<Ram SKU>. None if not priceable.

    Matches only the PREDEFINED instance SKUs of the machine family ('N4D Instance Core running
    in ...', 'N2D AMD Instance ...'), so the Custom / Sole Tenancy / Commitment / Confidential
    variants can never be picked up by accident. RAM is billed per GiB for most families and
    per GB for a few (C4D); that difference is well inside this estimate's error bar."""
    fam = machine_type.split("-", 1)[0].upper()
    spot_prefix = "Spot Preemptible " if spot else ""
    pat = re.compile(rf"^{spot_prefix}{fam}(?: [A-Z0-9]+)? Instance (Core|Ram) running in ")
    rates = {}
    for sku in billing_skus(region):
        cat = sku.get("category") or {}
        if cat.get("resourceFamily") != "Compute":
            continue
        if cat.get("usageType") != ("Preemptible" if spot else "OnDemand"):
            continue
        m = pat.match(sku.get("description") or "")
        if m:
            rates[m.group(1)] = sku_price(sku)
    if rates.get("Core") is None or rates.get("Ram") is None:
        return None
    return cpus * rates["Core"] + mem_gb * rates["Ram"]


def disk_hourly(disk_type, disk_gb, region, iops=0, throughput=0):
    """EUR/hour for one disk: capacity, plus metered Hyperdisk IOPS/throughput above baseline."""
    cap_sku, iops_sku, thr_sku = DISK_SKU_PREFIX.get(disk_type, (None, None, None))
    if cap_sku is None:
        return None
    prices = {}
    for sku in billing_skus(region):
        if (sku.get("category") or {}).get("resourceFamily") != "Storage":
            continue
        desc = sku.get("description") or ""
        for key, want in (("cap", cap_sku), ("iops", iops_sku), ("thr", thr_sku)):
            if want and desc.startswith(want + " in "):   # excludes Regional/HA/Confidential/pool
                prices[key] = sku_price(sku)
    if prices.get("cap") is None:
        return None
    per_month = prices["cap"] * disk_gb
    if prices.get("iops"):
        per_month += max(0, (iops or 0) - HYPERDISK_FREE_IOPS) * prices["iops"]
    if prices.get("thr"):
        per_month += max(0, (throughput or 0) - HYPERDISK_FREE_THROUGHPUT) * prices["thr"]
    return per_month / HOURS_PER_MONTH


def estimate_hourly(machine_type, cpus, mem_gb, disk_gb, disk_type, region, spot,
                    iops=0, throughput=0, boot_gb=0):
    """Estimated all-in EUR/hour for one VM (machine + work disk + boot disk), None if the
    catalog yielded no price. The boot disk is billed at baseline performance only -- the
    provisioned IOPS/throughput go to the work disk, where all the I/O happens."""
    machine = machine_hourly(machine_type, cpus, mem_gb, region, spot)
    work = disk_hourly(disk_type, disk_gb, region, iops, throughput)
    boot = disk_hourly(disk_type, boot_gb, region)
    if machine is None or work is None or boot is None:
        return None
    return round(machine + work + boot, 4)


def machine_shape(machine_type, zone, project, fallback=None):
    """(cpus, mem_gb) for a machine type, read from the Compute API so no local table can drift.

    `fallback` is used when the API cannot be reached -- the built-in defaults know their own
    shape. Without one (a user-pinned --machine-type / vm_type) an unresolvable machine is
    fatal, because the generated resources.config has to match the real box."""
    key = (machine_type, zone)
    if key not in _SHAPE_CACHE:
        res = subprocess.run(
            ["gcloud", "compute", "machine-types", "describe", machine_type,
             "--zone", zone, "--project", project, "--format=value(guestCpus,memoryMb)"],
            capture_output=True, text=True)
        parts = res.stdout.split() if res.returncode == 0 else []
        if len(parts) == 2 and all(p.isdigit() for p in parts):
            _SHAPE_CACHE[key] = (int(parts[0]), int(parts[1]) // 1024)
        elif fallback:
            print(f"  (warning) could not read the shape of {machine_type} in {zone}; "
                  f"assuming {fallback[0]} vCPU / {fallback[1]} GB")
            _SHAPE_CACHE[key] = tuple(fallback)
        else:
            sys.exit(f"machine type '{machine_type}' could not be resolved in zone {zone}: "
                     f"{(res.stderr or '').strip()}")
    return _SHAPE_CACHE[key]

# Columns that exist only to drive sizing/params and must NOT go into the real samplesheet.
HELPER_COLUMNS = {"coverage", "vm_type", "sequencing_type"}

# Accepted values for oncoanalyser's --sequencing_type
SEQUENCING_TYPES = {"ILLUMINA", "ULTIMA", "SBX"}

STARTUP_TEMPLATE = r"""#!/bin/bash
set -uo pipefail

GROUP_ID="%%GROUP_ID%%"
REVISION="%%REVISION%%"
GENOME="%%GENOME%%"
MODE="%%MODE%%"
SEQUENCING_TYPE="%%SEQUENCING_TYPE%%"
NXF_VER="%%NXF_VER%%"
RUN_ROOT="%%RUN_ROOT%%"
RUN_ID="%%RUN_ID%%"
STAGING_PREFIX="%%STAGING_PREFIX%%"
REFERENCE_CONFIG_GCS="%%REFERENCE_CONFIG_GCS%%"
PANEL_CONFIG_GCS="%%PANEL_CONFIG_GCS%%"
ZONE="%%ZONE%%"
KEEP_ON_FAILURE="%%KEEP_ON_FAILURE%%"
EXTRA_ARGS="%%EXTRA_ARGS%%"
LOCALIZE_INPUTS="%%LOCALIZE_INPUTS%%"
MACHINE_TYPE="%%MACHINE_TYPE%%"
DISK_GB="%%DISK_GB%%"
SPOT="%%SPOT%%"
HOURLY_RATE="%%HOURLY_RATE%%"
CURRENCY="%%CURRENCY%%"
COST_RECORD_URI="%%COST_RECORD_URI%%"

START_EPOCH=$(date +%s)
WORKDIR=/mnt/work
LOG=/var/log/oncoanalyser-startup.log
exec > >(tee -a "$LOG") 2>&1
echo "[$(date -u)] === oncoanalyser startup for ${GROUP_ID} (run ${RUN_ID}) ==="

# All run artifacts are namespaced by a unique run-id so separate runs never clobber each other.
RUN_PREFIX="${RUN_ROOT%/}/${GROUP_ID}"

report() {  # $1 = SUCCESS|FAILED
  echo "$1 $(date -u)" > /tmp/_STATUS
  gcloud storage cp /tmp/_STATUS "${RUN_PREFIX}/_${1}" || true
  gcloud storage cp "$LOG" "${RUN_PREFIX}/startup.log" || true
}

self_delete() {
  # Flip the work disk to auto-delete BEFORE deleting the instance, so a clean finish reclaims it
  # with the VM -- no orphan, no launcher-side race. A Spot PREEMPTION never runs this path: GCE
  # deletes the instance while auto-delete is still off, so the disk SURVIVES for the relaunch.
  gcloud compute instances set-disk-auto-delete "$(hostname)" --zone "$ZONE" \
    --device-name=work --auto-delete --quiet || true
  gcloud compute instances delete "$(hostname)" --zone "$ZONE" --quiet || true
}

# --- mount the persistent work disk at /mnt/work ------------------------------------
# Attached as device-name=work -> /dev/disk/by-id/google-work. Format ONLY on first use: a disk
# that already carries a filesystem is a preemption relaunch, and its work/ cache is exactly what
# lets -resume skip the completed steps. util-linux + e2fsprogs ship in the base image.
WORK_DEV=/dev/disk/by-id/google-work
for _ in $(seq 1 30); do [ -b "$WORK_DEV" ] && break; echo "[$(date -u)] waiting for work disk ${WORK_DEV}..."; sleep 2; done
if [ ! -b "$WORK_DEV" ]; then
  echo "[$(date -u)] ERROR: work disk ${WORK_DEV} never appeared"; report FAILED
  if [ "$KEEP_ON_FAILURE" != "true" ]; then self_delete; fi
  exit 1
fi
if ! blkid "$WORK_DEV" >/dev/null 2>&1; then
  echo "[$(date -u)] formatting fresh work disk ${WORK_DEV} (ext4)"
  mkfs.ext4 -q -F "$WORK_DEV"
fi
mkdir -p "$WORKDIR"
mount "$WORK_DEV" "$WORKDIR"
echo "[$(date -u)] mounted work disk: $(df -h "$WORKDIR" | tail -1)"
cd "$WORKDIR"

export DEBIAN_FRONTEND=noninteractive
apt-get update -y
apt-get install -y openjdk-17-jre-headless docker.io wget curl apt-transport-https ca-certificates gnupg

# Google Cloud CLI (Ubuntu base image does not ship it)
curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg \
  | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg
echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
  > /etc/apt/sources.list.d/google-cloud-sdk.list
apt-get update -y && apt-get install -y google-cloud-cli

# Put Docker's storage (image layers + the container /tmp scratch where sambamba sort spills)
# on the big persistent work disk, not the 100 GB boot disk. Bonus: images survive preemption,
# so relaunches skip re-pulling them. daemon.json must be in place before the daemon starts.
mkdir -p "$WORKDIR/docker" /etc/docker
cat > /etc/docker/daemon.json <<DOCKEREOF
{ "data-root": "$WORKDIR/docker" }
DOCKEREOF
systemctl restart docker

export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_VER="$NXF_VER"   # pin runtime to [25.04, 26.04): nf-schema needs >=25.04.0; 26.04+ defaults to the v2 config parser dev-3.0.0 can't parse
curl -s https://get.nextflow.io | bash
mv nextflow /usr/local/bin/

# Verify the version pin actually took effect before committing to a long, expensive run.
ACTUAL_VER=$(nextflow -version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
echo "[$(date -u)] Nextflow runtime: ${ACTUAL_VER} (pinned NXF_VER=${NXF_VER})"
if [ "${ACTUAL_VER}" != "${NXF_VER}" ]; then
  echo "[$(date -u)] ERROR: Nextflow ${ACTUAL_VER} is running but ${NXF_VER} was pinned. Aborting before wasting compute."
  report FAILED
  if [ "$KEEP_ON_FAILURE" != "true" ]; then self_delete; fi
  exit 1
fi

# Pin the exact pipeline revision (branch) up front so the run resolves to dev-3.0.0 (or whatever was set)
nextflow pull nf-core/oncoanalyser -revision "$REVISION"

# Pull generated inputs
gcloud storage cp "${STAGING_PREFIX%/}/samplesheet.csv" ./samplesheet.csv
gcloud storage cp "${STAGING_PREFIX%/}/resources.config" ./resources.config

REF_ARG=""
if [ -n "$REFERENCE_CONFIG_GCS" ]; then
  gcloud storage cp "$REFERENCE_CONFIG_GCS" ./reference_data.config
  REF_ARG="-config reference_data.config"
fi

# Targeted (panel) runs need a second config carrying params.panel_data_paths +
# params.ref_data_panel_data_path. Later -config wins on key collisions, so this is layered
# on top of the reference config (they define disjoint keys).
if [ -n "$PANEL_CONFIG_GCS" ]; then
  gcloud storage cp "$PANEL_CONFIG_GCS" ./panel_data.config
  REF_ARG="$REF_ARG -config panel_data.config"
fi

# Localize gs:// inputs to local disk with gcloud (parallel + resumable), then point the
# samplesheet at the local copies. This avoids Nextflow's single-stream foreign-file copier,
# which has no resume and stalls on large (100 GB+) CRAM/BAM inputs.
INPUT_SHEET=samplesheet.csv
if [ "$LOCALIZE_INPUTS" = "true" ]; then
  echo "[$(date -u)] localizing gs:// inputs to /mnt/work/inputs (robust parallel download)"
  mkdir -p /mnt/work/inputs
  python3 - <<'PYEOF'
import csv, os, subprocess, sys
INDIR = "/mnt/work/inputs"
os.makedirs(INDIR, exist_ok=True)

def remote_size(src):
    r = subprocess.run(["gcloud", "storage", "du", src], capture_output=True, text=True)
    try:
        return int(r.stdout.split()[0]) if r.returncode == 0 else None
    except (ValueError, IndexError):
        return None

def cp(src, required=True):
    dst = os.path.join(INDIR, src.rsplit("/", 1)[-1])
    # On a preemption relaunch the inputs already sit on the persistent work disk. Skip the
    # re-download only when the local size matches the remote object exactly -- a truncated
    # partial left by a preemption mid-localize won't match, so it gets re-fetched.
    if os.path.exists(dst):
        rs = remote_size(src)
        if rs is not None and os.path.getsize(dst) == rs:
            return dst
    err = None if required else subprocess.DEVNULL
    rc = subprocess.run(["gcloud", "storage", "cp", src, dst], stderr=err).returncode
    if rc != 0 and required:
        sys.exit("ERROR: failed to localize " + src)
    return dst if rc == 0 else None

def cp_tree(src, dst):
    # Directory inputs (amber_dir, purple_dir, ... -- see FileType.groovy) are whole prefixes, not
    # single objects, so `cp` cannot fetch them. rsync is idempotent, so a preemption relaunch
    # re-uses whatever already sits on the persistent work disk instead of re-downloading.
    os.makedirs(dst, exist_ok=True)
    rc = subprocess.run(["gcloud", "storage", "rsync", "-r", src.rstrip("/"), dst]).returncode
    if rc != 0:
        sys.exit("ERROR: failed to localize directory " + src)
    return dst

def localize_one(src, fetched):
    # Download one gs:// file (+ its co-located index for BAM/CRAM); return the local path.
    local = cp(src, required=True)
    low = src.lower()
    # also pull a co-located index so oncoanalyser finds it next to the local alignment file
    cands = []
    if low.endswith(".cram"):
        cands = [src + ".crai", src[:-5] + ".crai"]
    elif low.endswith(".bam"):
        cands = [src + ".bai", src[:-4] + ".bai"]
    for idx in cands:
        if idx not in fetched:
            cp(idx, required=False)   # ok if a given index naming form doesn't exist
            fetched.add(idx)
    return local

with open("samplesheet.csv") as fh:
    reader = csv.DictReader(fh); cols = reader.fieldnames; rows = list(reader)

fetched = set()
for row in rows:
    fp = (row.get("filepath") or "").strip()
    if not fp:
        continue
    # Existing-output directory rows resume a run without the Nextflow work/ cache. Local copies
    # are named <sample_id>.<filetype> because tool subdirectory basenames collide: sage/somatic
    # and pave/somatic would otherwise both land on inputs/somatic.
    ft_row = (row.get("filetype") or "").strip().lower()
    if ft_row.endswith("_dir"):
        if fp.startswith("gs://"):
            row["filepath"] = cp_tree(fp, os.path.join(INDIR, row.get("sample_id", "") + "." + ft_row))
        continue
    # FASTQ rows carry a ';'-separated R1;R2 pair in the single filepath field; BAM/CRAM is one
    # path. Localize each gs:// part independently and re-join with ';' so the pair is preserved.
    new_parts = []
    for part in (p.strip() for p in fp.split(";")):
        new_parts.append(localize_one(part, fetched) if part.startswith("gs://") else part)
    row["filepath"] = ";".join(new_parts)
    # For REDUX BAM/CRAM inputs, also localize the sibling REDUX TSVs so oncoanalyser can
    # infer them from the (now local) alignment file's directory and skip re-running REDUX.
    if (row.get("filetype") or "").strip().lower() in ("bam_redux", "cram_redux"):
        src = fp.split(";")[0].strip()
        if src.startswith("gs://"):
            base = src.rsplit(".", 1)[0]  # .../<sample>.redux.bam -> .../<sample>.redux
            for suffix, req in ((".bqr.tsv", True), (".jitter_params.tsv", True),
                                (".ms_table.tsv.gz", True), (".duplicate_freq.tsv", False)):
                cp(base + suffix, required=req)

with open("samplesheet.local.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)
print("wrote samplesheet.local.csv with local input paths")
PYEOF
  if [ ! -s samplesheet.local.csv ]; then
    echo "[$(date -u)] ERROR: input localization failed"; report FAILED
    if [ "$KEEP_ON_FAILURE" != "true" ]; then self_delete; fi
    exit 1
  fi
  INPUT_SHEET=samplesheet.local.csv
fi

echo "[$(date -u)] launching nextflow"
nextflow run nf-core/oncoanalyser \
  -revision "$REVISION" \
  -profile docker \
  $REF_ARG \
  -config resources.config \
  --mode "$MODE" \
  --genome "$GENOME" \
  --sequencing_type "$SEQUENCING_TYPE" \
  --input "$INPUT_SHEET" \
  --outdir ./output \
  $EXTRA_ARGS \
  -resume
RC=$?
echo "[$(date -u)] nextflow exited with code ${RC}"

# Always try to ship outputs + logs, even on failure
gcloud storage rsync -r ./output "${RUN_PREFIX}/output" || true
[ -f ./.nextflow.log ] && gcloud storage cp ./.nextflow.log "${RUN_PREFIX}/nextflow.log" || true

# --- cost record (self-computed estimate from wall time x injected hourly rate) ---
END_EPOCH=$(date +%s)
ELAPSED=$((END_EPOCH - START_EPOCH))
STATUS=$([ "$RC" -eq 0 ] && echo SUCCESS || echo FAILED)
if [ -n "$HOURLY_RATE" ]; then
  EST_COST=$(awk "BEGIN{printf \"%.2f\", ${HOURLY_RATE} * ${ELAPSED} / 3600}")
  HOURS=$(awk "BEGIN{printf \"%.2f\", ${ELAPSED} / 3600}")
  cat > /tmp/cost.json <<EOF
{"group":"${GROUP_ID}","machine_type":"${MACHINE_TYPE}","disk_gb":${DISK_GB},"spot":${SPOT},"sequencing_type":"${SEQUENCING_TYPE}","status":"${STATUS}","elapsed_seconds":${ELAPSED},"elapsed_hours":${HOURS},"hourly_rate":${HOURLY_RATE},"currency":"${CURRENCY}","est_cost":${EST_COST}}
EOF
  gcloud storage cp /tmp/cost.json "${COST_RECORD_URI}" || true
  echo "[$(date -u)] estimated cost: ${EST_COST} ${CURRENCY} (${HOURS} h x ${HOURLY_RATE} ${CURRENCY}/h)"
fi

if [ "$RC" -eq 0 ]; then
  report SUCCESS
  self_delete
else
  report FAILED
  if [ "$KEEP_ON_FAILURE" = "true" ]; then
    echo "[$(date -u)] failure: keeping VM for debugging (delete it manually when done)."
  else
    self_delete
  fi
fi
"""

# Runs during the ~30s grace when a Spot VM is preempted. It writes a _PREEMPTED marker so the
# launcher can tell a preemption apart from a real pipeline _FAILED and relaunch the group. The
# guard on the `preempted` metadata flag keeps it from firing on the normal post-run self-delete.
# Best-effort only: the grace is short and gcloud may not be installed if preemption hits early in
# boot -- the launcher's "VM vanished with no terminal marker" check is the backstop for that.
SHUTDOWN_TEMPLATE = r"""#!/bin/bash
RUN_PREFIX="%%RUN_PREFIX%%"
PREEMPTED=$(curl -s -H "Metadata-Flavor: Google" \
  "http://metadata.google.internal/computeMetadata/v1/instance/preempted" 2>/dev/null || echo "")
if [ "$PREEMPTED" = "TRUE" ]; then
  echo "PREEMPTED $(date -u)" > /tmp/_PREEMPTED
  gcloud storage cp /tmp/_PREEMPTED "${RUN_PREFIX%/}/_PREEMPTED" || true
  gcloud storage cp /var/log/oncoanalyser-startup.log "${RUN_PREFIX%/}/startup.preempted.log" || true
fi
"""


def sanitize(name: str) -> str:
    """RFC1035-safe instance name fragment."""
    s = re.sub(r"[^a-z0-9-]", "-", name.lower())
    return re.sub(r"-+", "-", s).strip("-")[:40] or "g"


RUN_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,49}$")


def validate_run_id(run_id):
    """Allow a human-friendly run-id, but keep it safe for GCS paths and derivable VM names."""
    if not RUN_ID_RE.match(run_id):
        sys.exit("--run-id must be 1-50 chars of letters, digits, '.', '_' or '-' and start with a "
                 "letter or digit (e.g. 'hcc-celllines-fastq-20260608')")
    return run_id


def group_is_fastq(group_rows):
    """True if any row is a FASTQ input, meaning the VM must align (and needs the fatter tiers)."""
    return any((r.get("filetype") or "").strip().lower() == "fastq" for r in group_rows)


def group_has_rna(group_rows):
    """True if the group carries an RNA sample, so ISOFOX runs alongside REDUX.

    This is the WGTS path: it earns a wider box plus an ISOFOX core reserve. DNA-only groups must
    stay on exactly the tuning they have today -- see the RNA note in the MACHINE SIZING block."""
    return any((r.get("sequence_type") or "dna").strip().lower() == "rna" for r in group_rows)


# Below this tumor:normal work ratio the default ~3:1 REDUX core split is wrong -- either the
# normal came in unusually deep, or it is the LARGER of the pair. Warn loudly and re-derive.
# In practice a T/N pair sits at 2.5-3.5 (HG008: 150x/54x tumor, size ratio 3.01, measured REDUX
# runtime ratio 3.17), so anything under 2 means the assumption behind the split has broken.
REDUX_RATIO_WARN_BELOW = 2.0
DNA_SAMPLE_TYPES = ("tumor", "normal", "donor")


def group_dna_inputs(group_rows):
    """{sample_type: {"paths": [gs://...], "filetypes": {...}}} for this group's DNA alignments.

    FASTQ rows carry ';'-joined R1;R2 in one field and one row per lane, so paths accumulate."""
    out = {}
    for r in group_rows:
        if (r.get("sequence_type") or "dna").strip().lower() != "dna":
            continue
        st = (r.get("sample_type") or "").strip().lower()
        if st not in DNA_SAMPLE_TYPES:
            continue
        ft = (r.get("filetype") or "").strip().lower()
        for p in (r.get("filepath") or "").split(";"):
            p = p.strip()
            if p.startswith("gs://"):
                e = out.setdefault(st, {"paths": [], "filetypes": set()})
                e["paths"].append(p)
                e["filetypes"].add(ft)
    return out


def gcs_object_sizes(paths, chunk=100):
    """{gs://path: bytes} via batched `gcloud storage ls -l`. Missing/unreadable paths are simply
    absent from the result -- callers must treat this as best-effort and fall back, never abort."""
    sizes = {}
    paths = sorted(set(paths))
    for i in range(0, len(paths), chunk):
        res = subprocess.run(["gcloud", "storage", "ls", "-l", *paths[i:i + chunk]],
                             capture_output=True, text=True)
        # Partial success is normal (one bad path -> ERROR on stderr, good lines still on stdout),
        # so parse whatever came back regardless of the exit code.
        for line in (res.stdout or "").splitlines():
            parts = line.split()
            if len(parts) >= 3 and parts[0].isdigit() and parts[-1].startswith("gs://"):
                sizes[parts[-1]] = int(parts[0])
    return sizes


def redux_work_ratio(group_rows, sizes):
    """(ratio, basis, detail): tumor REDUX work / the largest other DNA sample's work.

    Ratio only -- never absolute coverage. Compression per x varies far too much across
    platforms and BAM-vs-CRAM to infer depth from bytes, but a ratio between two files of the
    SAME filetype cancels that out, and it beat declared coverage as a predictor on HG008
    (size 3.01 vs coverage 2.78, against a measured REDUX runtime ratio of 3.17).

    Returns ratio=None when the group is tumor-only (nothing to balance against) or when
    neither basis is available.
    """
    inputs = group_dna_inputs(group_rows)
    others = [st for st in inputs if st != "tumor"]
    if "tumor" not in inputs:
        # RNA-only, or a normal/donor-only group. REDUX either does not run or has one task.
        return None, "single-sample", "no DNA tumor sample in this group"
    if not others:
        return None, "single-sample", "tumor-only: no second DNA sample to balance against"

    # 1. file size -- but only when both sides are the same filetype, since a CRAM is roughly
    #    2-2.5x smaller than the equivalent BAM and the ratio would be meaningless.
    t_ft = inputs["tumor"]["filetypes"]
    t_bytes = sum(sizes.get(p, 0) for p in inputs["tumor"]["paths"])
    if t_bytes:
        best = 0.0
        for st in others:
            if inputs[st]["filetypes"] != t_ft:
                best = 0.0
                break
            o_bytes = sum(sizes.get(p, 0) for p in inputs[st]["paths"])
            best = max(best, o_bytes)
        if best:
            return (t_bytes / best, "file size",
                    f"tumor {t_bytes / 2**30:.1f} GiB vs largest other {best / 2**30:.1f} GiB")

    # 2. declared coverage on both sides
    cov = {}
    for r in group_rows:
        st = (r.get("sample_type") or "").strip().lower()
        c = (r.get("coverage") or "").strip()
        if st in DNA_SAMPLE_TYPES and c:
            try:
                cov[st] = max(cov.get(st, 0.0), float(c))
            except ValueError:
                pass
    o_cov = max((cov.get(st, 0.0) for st in others), default=0.0)
    if cov.get("tumor") and o_cov:
        return (cov["tumor"] / o_cov, "coverage column",
                f"tumor {cov['tumor']:g}x vs largest other {o_cov:g}x")

    return None, None, "no file sizes and no coverage column for both samples"


def plan_redux_split(cpus, ratio, budget=None):
    """(tumor_cpus, other_cpus) so co-resident REDUX tasks finish as close together as possible.

    Cores are allocated proportional to work, which needs no special cases: it lands on exactly
    36/12 for the usual ~3:1 pair on a 48-core box, collapses to 24/24 when the pair is balanced,
    and inverts correctly when the normal is the larger sample. ratio=None (tumor-only) hands the
    whole box to the single REDUX task. A floor stops a wild ratio starving either side.

    `budget` is how many of the box's cores REDUX may claim; it defaults to all of them. On
    RNA-bearing groups the caller passes cpus - ISOFOX_RESERVE_CPUS so ISOFOX stays admissible
    while REDUX runs. Note the pair has always summed to exactly its budget, which is why
    --redux-split alone could never open that gap -- every T:N is a re-slice of the same total."""
    budget = cpus if budget is None else max(4, min(cpus, budget))
    if ratio is None:
        return budget, max(2, budget // 8)
    floor = max(2, budget // 8)
    t = int(round(budget * ratio / (ratio + 1.0)))
    t -= t % 2                                   # keep even core counts
    return max(floor, min(budget - floor, t)), budget - max(floor, min(budget - floor, t))


def parse_isofox_reserve(spec):
    """Resolve --isofox-reserve into a core count ('auto' -> ISOFOX_RESERVE_CPUS, '0' -> off)."""
    spec = (spec or "auto").strip().lower()
    if spec == "auto":
        return ISOFOX_RESERVE_CPUS
    try:
        n = int(spec)
    except ValueError:
        sys.exit(f"--isofox-reserve must be 'auto' or a whole number of cores; got {spec!r}")
    if n < 0:
        sys.exit(f"--isofox-reserve cannot be negative; got {spec!r}")
    return n


def parse_redux_split(spec, cpus, ratio, budget=None):
    """Resolve --redux-split into (tumor_cpus, other_cpus), within `budget` cores."""
    spec = (spec or "auto").strip().lower()
    if spec == "auto":
        return plan_redux_split(cpus, ratio, budget)
    if ":" not in spec:
        sys.exit(f"--redux-split must be 'auto' or 'T:N' (e.g. '3:1', '1:1'); got {spec!r}")
    try:
        a, b = (float(x) for x in spec.split(":", 1))
    except ValueError:
        sys.exit(f"--redux-split must be 'auto' or 'T:N' with numeric parts; got {spec!r}")
    if a <= 0 or b <= 0:
        sys.exit(f"--redux-split parts must both be > 0; got {spec!r}")
    return plan_redux_split(cpus, a / b)


def group_max_coverage(group_rows):
    """Highest `coverage` value declared in the group, or None if the column is absent/blank."""
    covs = []
    for r in group_rows:
        c = (r.get("coverage") or "").strip()
        if c:
            try:
                covs.append(float(c))
            except ValueError:
                pass
    return max(covs) if covs else None


def pick_tier(group_rows, zone, project, machine_override=""):
    """Return (machine_type, cpus, mem_gb, disk_gb) for one group.

    The MACHINE is fixed by workload shape (FASTQ groups align and scale out; BAM/CRAM groups
    do not); the one exception is an RNA arm, which adds a genuinely independent process
    (ISOFOX) that can run from t=0 and so has real cores to spend. Coverage drives only the
    DISK. An explicit `vm_type` column, or --machine-type, still wins -- its shape is read from
    the Compute API so the generated resources.config matches the real box.
    """
    fastq = group_is_fastq(group_rows)
    if fastq:
        default = FASTQ_MACHINE      # already wide enough to absorb the RNA arm
    elif group_has_rna(group_rows):
        default = RNA_MACHINE
    else:
        default = DEFAULT_MACHINE
    machine_type = default[0]
    cpus, mem_gb = machine_shape(machine_type, zone, project, fallback=default[1:])

    vm_type = machine_override.strip()
    for r in group_rows:                       # a per-row vm_type beats the global --machine-type
        vt = (r.get("vm_type") or "").strip()
        if vt:
            vm_type = vt
            break
    if vm_type:
        machine_type = vm_type
        cpus, mem_gb = machine_shape(vm_type, zone, project)

    tiers = FASTQ_DISK_TIERS if fastq else DISK_TIERS
    mx = group_max_coverage(group_rows)
    if mx is None:
        disk_gb = tiers[FASTQ_DEFAULT_DISK_INDEX if fastq else DEFAULT_DISK_INDEX][1]
    else:
        disk_gb = next((d for cov, d in tiers if mx <= cov), tiers[-1][1])
    return machine_type, cpus, mem_gb, disk_gb


def pick_sequencing_type(group_rows, default):
    """Per-group sequencing_type from the manifest column, else the global default."""
    for r in group_rows:
        st = (r.get("sequencing_type") or "").strip().upper()
        if st:
            if st not in SEQUENCING_TYPES:
                sys.exit(f"sequencing_type '{st}' must be one of {sorted(SEQUENCING_TYPES)}")
            return st
    return default


def read_manifest(path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []
        rows = list(reader)
    groups = OrderedDict()
    for r in rows:
        gid = (r.get("group_id") or "").strip()
        if not gid:
            sys.exit("Every row needs a non-empty group_id")
        groups.setdefault(gid, []).append(r)
    sheet_cols = [c for c in fieldnames if c not in HELPER_COLUMNS]
    return groups, sheet_cols


def write_samplesheet(group_rows, sheet_cols, out_path):
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sheet_cols, extrasaction="ignore")
        w.writeheader()
        for r in group_rows:
            w.writerow({c: r.get(c, "") for c in sheet_cols})


def parse_purple_purity_range(spec):
    """Parse a '--purple-purity-range MIN,MAX' value into a PURPLE ext.args string.

    Returns '' when unset, else '-min_purity MIN -max_purity MAX'. Aborts on a malformed
    or out-of-bounds range so a typo fails fast at launch rather than producing a bad fit."""
    spec = (spec or "").strip()
    if not spec:
        return ""
    parts = spec.split(",")
    if len(parts) != 2:
        sys.exit(f"--purple-purity-range must be 'MIN,MAX' (got {spec!r})")
    try:
        lo, hi = float(parts[0]), float(parts[1])
    except ValueError:
        sys.exit(f"--purple-purity-range values must be numbers (got {spec!r})")
    if not (0.0 <= lo <= hi <= 1.0):
        sys.exit(f"--purple-purity-range must satisfy 0 <= MIN <= MAX <= 1 (got {spec!r})")
    return f"-min_purity {lo:g} -max_purity {hi:g}"


def plan_cpus(cpus, isofox_cpus=None):
    """Per-process thread counts, scaled to the box.

    oncoanalyser's nf-core labels cap every task at 12 threads (`process_high`), far below what
    these tools use in HMF's own pipeline5 (HmfTool.java: REDUX 32, SAGE 32, ESVEE 32, AMBER 16,
    COBALT 16, TEAL 32), leaving ~60% of a 48-core box idle. These fractions re-spend that idle
    capacity on the critical path. The REDUX split is the biggest win: tumor and normal start
    together but a measured T/N run took 2h14m vs 42m, so the box sat at ~25% utilisation for 84
    minutes; splitting cores ~3:1 (the measured work ratio) roughly halves the stage.
    """
    def frac(f, lo=2):
        return max(lo, min(cpus, int(round(cpus * f / 2.0)) * 2))   # round to an even count
    plan = {
        "esvee":         frac(0.50),   # longest post-REDUX task -> the post-REDUX critical path
        "sage_somatic":  frac(0.50),
        "amber":         frac(0.33),
        "cobalt":        frac(0.33),
        "isofox":        frac(0.33),
        # BWAMEM2_ALIGN is many independent shards; 12 threads measured 12.33 cores of actual draw
        # (103% efficiency), so leave it and let concurrency (not thread count) fill the box.
        "bwamem2":       12,
        # small enough to slot into the cores left over while REDUX tumor is still running
        "light":         frac(0.17),
        "tiny":          frac(0.08),
    }
    # Tools that provably do NOT scale with the threads they are handed. Granting more than they
    # can use only reserves slots that block peers, so these are ABSOLUTE counts (clamped to fit a
    # small VM): they do not grow with the box, because the tool's parallelism doesn't.
    for k, v in PROCESS_CPUS_MEASURED.items():
        plan[k] = max(2, min(cpus, v))
    # On RNA groups ISOFOX is pinned to the reserve that was held back from the REDUX split, so
    # the two numbers cannot drift apart and leave ISOFOX unable to start (or over-committed).
    if isofox_cpus:
        plan["isofox"] = max(2, min(cpus, isofox_cpus))
    return plan


# Measured core draw, from the %cpu column of a real T/N WGS execution_trace. VIRUSBREAKEND is the
# worst offender: 4.3 cores from a 12-core grant, so the 16 the fraction table would hand it
# reserves a third of a 48-core box for something that cannot use it.
PROCESS_CPUS_MEASURED = {
    "virusbreakend": 6,    # measured 4.3 of 12 granted (427%)
    "sage_germline": 8,    # measured 6.1 of 12 (613%) -- and only runs ~30s
    "teal":          6,    # measured 4.4 of 6  (440%)
    "cider":         4,    # measured 2.2 of 6  (217%)
    "pave_somatic":  4,    # measured 3.1 of 6  (312%)
    "pave_germline": 4,    # measured 2.4 of 6  (239%)
    "lilac":         4,    # measured 2.3 of 6  (226%)
}


# Memory per process, in GB. Properties of the DATA, not the box, so they do not scale with the
# machine -- resourceLimits clamps them down on a small VM. Sized from measured execution_trace
# peak_rss on a T/N WGS run, cross-checked against HMF pipeline5, and split into two tiers because
# "observed peak" is only a safe basis for one of them.
#
# TIER 1, STRUCTURALLY BOUNDED: memory is set by the reference or by read depth, both visible in
# advance, so the observed peak plus headroom is safe.
#     process        peak_rss   HMF    here
#     REDUX          25.8 GB    64 GB   64   (scales with depth, which the manifest declares)
#     AMBER          16.8 GB    24 GB   32   (BAF at a fixed het-site list)
#     COBALT         12.5 GB    24 GB   24   (read-depth ratios over fixed 1kb bins)
#     SAGE_GERMLINE   6.8 GB    64 GB   24   (hotspots + coding only; TMB-independent)
#     PAVE_GERMLINE   7.9 GB    32 GB   24
#     BAMTOOLS        2.8 GB    24 GB   16   (streams the BAM)
#     TEAL_PREP       3.9 GB    32 GB   24   (telomeric reads only)
#     CIDER          13.8 GB    24 GB   24   (IG/TR loci only)
#     LILAC             ~8 GB   24 GB   24   (HLA region only)
#     BWAMEM2_ALIGN  35.5 GB    72 GB   48   (fixed -K batch + fixed sort buffer; see note below)
#
# TIER 2, BIOLOGY-SCALING: memory tracks somatic variant or SV burden -- a property of the TUMOUR,
# invisible in the manifest, and several times higher on a hypermutated or chromothriptic sample.
# So never size these below oncoanalyser's own stock label, however low the observed peak, and give
# them a working escalation path (see the oom_retry note in write_resources_config).
#     process        peak_rss   stock   here
#     ESVEE          86.4 GB    72 GB   96   (SV count; already the one genuine hog)
#     SAGE_SOMATIC   12.4 GB    72 GB   72   (variant count -- 12.4 GB here says nothing about TMB)
#     PAVE_SOMATIC   20.1 GB    36 GB   36   (annotates every somatic call)
#     PURPLE          5.5 GB    36 GB   96   (loads the whole somatic VCF; --purple-mem-gb)
#     LINX_SOMATIC     <1 GB    12 GB   16   (SV chains)
#
# NOTE: LINX's 16 GB is a JVM-only figure and does NOT transfer to LINX_VISUALISER, which forks
# task.cpus concurrent R children -- that one is sized separately, see linx_vis_mem below.
#
# BWAMEM2_ALIGN is the biggest single win, and the stock label is what blocks it. Measured across
# 24 shards of a T/N WGS run: peak_rss 35.0-35.9 GB (a 2.6% spread, while shard durations varied
# 1.8x) because memory is input-size independent by construction -- the module passes bwa-mem2
# `-K 100000000`, fixing the batch at 100M bases, and our sambamba `-m 16GiB` sort buffer is fixed
# too. At the stock 72 GB a 320 GB box admits only 4 shards (48 of 80 cores); at 48 GB it admits 6
# (74 cores, 435 MB/s of disk -- 54% of a 800 MiB/s work disk, so I/O is not the new limit). That
# turns a ~13h alignment phase into ~6.5h.
PROCESS_MEM_GB = {
    # tier 1 -- structurally bounded, sized from observation
    "redux": 64, "amber": 32, "cobalt": 24, "sage_germline": 24, "pave_germline": 24,
    "bamtools": 16, "teal": 24, "cider": 24, "lilac": 24,
    "virusbreakend": 48, "isofox": 48, "sage_append": 24, "bwamem2": 48,
    # tier 2 -- biology-scaling, never below the stock label
    "esvee": 96, "sage_somatic": 72, "pave_somatic": 36, "linx": 16,
}


def write_resources_config(cpus, mem_gb, out_path, align_max_forks=0, purple_ext_args="",
                           redux_version="", purple_mem_gb=96, retune=True, redux_split=None,
                           isofox_cpus=None):
    # Leave headroom below the box's RAM for the OS, the Nextflow JVM and the Docker daemon.
    # This is also the executor pool, so it is the real cap on how much work runs at once.
    limit_mem = max(mem_gb - max(12, mem_gb // 12), 8)
    n = plan_cpus(cpus, isofox_cpus=isofox_cpus)
    m = PROCESS_MEM_GB
    # (tumor_cpus, other_cpus) from the measured tumor:normal work ratio -- see plan_redux_split.
    redux_tumor, redux_other = redux_split or plan_redux_split(cpus, 3.0)
    # Cores left unclaimed by the REDUX pair. On RNA groups this is the ISOFOX reserve, and the
    # emitted comment records it so a resources.config can be read back without the launcher.
    redux_spare = cpus - (redux_tumor + redux_other)
    isofox_note = ""
    if isofox_cpus and redux_spare > 0:
        isofox_note = (
            f"\n    // {redux_spare} core(s) are deliberately left unclaimed by the pair above so "
            f"ISOFOX ({n['isofox']}\n    // cpu / {m['isofox']} GB) stays admissible and runs "
            f"CONCURRENTLY with REDUX. ISOFOX depends only on\n    // the RNA BAM, so it is ready "
            f"at t=0; without this gap it queued out the whole REDUX stage.\n")

    # OOM escalation for the biology-scaling (tier 2) processes. These hmftools set -Xmx to 95% of
    # task.memory, so the heap ceiling sits BELOW the container's cgroup limit: a hypermutated
    # sample dies of a clean Java OutOfMemoryError (exit 1) and the kernel never sends SIGKILL
    # (137). base.config retries only on 130..145 + 104 + 175..177, so exit 1 falls through to
    # 'finish' and kills the run -- without adding 1 here the `* task.attempt` escalation below is
    # unreachable for the exact failure it exists to handle. Scoped to tier 2 only, because exit 1
    # is also the generic "something went wrong" code and retrying it everywhere would double the
    # cost of every genuinely broken sample.
    oom_retry = ("        errorStrategy = { task.exitStatus in "
                 "((130..145) + 104 + (175..177) + 1) ? 'retry' : 'finish' }\n"
                 "        maxRetries    = 2\n")

    # BWAMEM2_ALIGN sizing goes into the withName block below (which also carries the ulimit and
    # sambamba sort tunings), so build it separately. It gets the same exit-1 ladder: sambamba
    # exits 1 on its own failures ('sambamba-sort: Unable to write to stream'), which base.config
    # does not retry, so a shard that needs more escalates to 96 GB instead of killing the phase.
    align_retune = ""
    if retune:
        align_retune = (f"        cpus   = {n['bwamem2']}\n"
                        f"        memory = {{ {m['bwamem2']}.GB * task.attempt }}\n"
                        + oom_retry)
    # Cap on concurrent BWAMEM2_ALIGN shards; 0 = no cap (the fast path for many-lane samples).
    maxforks_line = f"        maxForks = {align_max_forks}\n" if align_max_forks and align_max_forks > 0 else ""
    # PURPLE loads the full somatic VCF into the JVM heap and OOMs at the default ~30 GB label on
    # hypermutated tumor-only samples (e.g. HCC1954, no matched normal to filter germline). A
    # large but bounded slice -- still 2.4x HMF's production 40 GB -- since the tail now overlaps
    # other work. Raise with --purple-mem-gb.
    purple_mem = min(max(purple_mem_gb, 16), limit_mem)
    # Optional fit-grid override (e.g. force purity into [0.95,1.0] off a half-purity alias).
    # PURPLE's module appends ${task.ext.args} verbatim and the pipeline sets none, so this is
    # additive.
    purple_args_line = f"        ext.args = '{purple_ext_args}'\n" if purple_ext_args else ""
    # Pin the REDUX container to a specific hmftools-redux release, so a run can use a newer redux
    # (e.g. 2.0.1 for new SBX samples) without waiting for a pipeline bump. Default: the
    # pipeline's own pin. Same mechanism as the STAR_ALIGN pin below.
    redux_block = (
        "\n"
        "    // Override the REDUX container to a specific hmftools-redux release.\n"
        f"    withName: 'REDUX' {{\n"
        f"        container = 'biocontainers/hmftools-redux:{redux_version}--hdfd78af_0'\n"
        "    }\n"
    ) if redux_version else ""
    # ---------------------------------------------------------------------------------
    # The retune block: two rules, applied to every process that matters.
    #   cpus   -- a constant, sized by plan_cpus() so the box is actually filled.
    #   memory -- what the tool really uses (PROCESS_MEM_GB), times task.attempt so the
    #             retry-on-OOM escalation still works. resourceLimits clamps it.
    #
    # `meta` is in scope inside a resource closure (directives are LazyMaps resolved against the
    # task context, which carries the process's input variables), and REDUX's meta.sample_type is
    # 'tumor' / 'normal' / 'donor' -- that is what lets tumor and normal get different core counts.
    # ---------------------------------------------------------------------------------
    # LINX_VISUALISER is the one process here that is NOT JVM-only: it forks task.cpus concurrent
    # Rscript children (Gviz ideogram panels), each with a ~1 GB resident floor. Nextflow turns
    # `memory` into a HARD `docker run --memory` cgroup cap while `cpus` is only a soft
    # --cpu-shares, so sizing threads against a JVM-derived memory figure gets the R children
    # SIGKILLed: at PROCESS_MEM_GB['linx'] (16 GB, measured from LINX_SOMATIC) an 80-core box gives
    # 14 threads ~1 GB each, where 3 of 6 children died in hmftools-linx:2.3 and all survived at
    # 2 GB. A killed child prints no R error, so LINX surfaces only
    #   [FATAL] error executing R script -> 'plotting error' -> exit 1
    # which looks like a data bug and is not one. Budget 2 GB per child + 4 GB of JVM headroom; the
    # resulting large -Xmx is virtual, not resident (LINX_SOMATIC held 17 GB vmem against 1 GB RSS).
    linx_vis_mem = max(m['linx'], 4 + 2 * n['light'])
    retune_block = ""
    if retune:
        retune_block = f"""
    // ---- REDUX: split the box between tumor and normal ----------------------------
    // They start together but the tumor is normally ~3x the work, so equal threads leave
    // the machine idle for the tail of the tumor task. Cores are allocated proportional to
    // the measured tumor:normal work ratio so both land together; this run resolved to
    // {redux_tumor}:{redux_other}. 'donor' is treated like 'normal'.
{isofox_note}    withName: 'REDUX' {{
        cpus   = {{ meta.sample_type == 'tumor' ? {redux_tumor} : {redux_other} }}
        memory = {{ {m['redux']}.GB * task.attempt }}
    }}

    // ---- post-REDUX burst ----------------------------------------------------------
    // ESVEE is the long pole here (56m of a 72m burst in the measured run) and it started
    // 16m late purely because AMBER/COBALT were each holding a 72 GB reservation they did
    // not need. Right-sized asks let it start immediately; the extra threads shorten it.
    withName: '.*:ESVEE_CALLING:ESVEE' {{
        cpus   = {n['esvee']}
        memory = {{ {m['esvee']}.GB * task.attempt }}
{oom_retry}    }}
    withName: '.*:SAGE_CALLING:SAGE_SOMATIC' {{
        cpus   = {n['sage_somatic']}
        memory = {{ {m['sage_somatic']}.GB * task.attempt }}
{oom_retry}    }}
    withName: '.*:SAGE_CALLING:SAGE_GERMLINE' {{
        cpus   = {n['sage_germline']}
        memory = {{ {m['sage_germline']}.GB * task.attempt }}
    }}
    withName: 'AMBER' {{
        cpus   = {n['amber']}
        memory = {{ {m['amber']}.GB * task.attempt }}
    }}
    withName: 'COBALT' {{
        cpus   = {n['cobalt']}
        memory = {{ {m['cobalt']}.GB * task.attempt }}
    }}
    withName: 'VIRUSBREAKEND' {{
        cpus   = {n['virusbreakend']}
        memory = {{ {m['virusbreakend']}.GB * task.attempt }}
    }}
    withName: 'ISOFOX' {{
        cpus   = {n['isofox']}
        memory = {{ {m['isofox']}.GB * task.attempt }}
    }}

    // ---- light tasks -----------------------------------------------------------------
    // Deliberately small so they fit in the cores left over while REDUX tumor is still
    // running (BAMTOOLS on the normal starts there) rather than queueing behind it.
    withName: '.*:BAMTOOLS_METRICS:BAMTOOLS' {{
        cpus   = {n['light']}
        memory = {{ {m['bamtools']}.GB * task.attempt }}
    }}
    withName: 'TEAL.*' {{
        cpus   = {n['teal']}
        memory = {{ {m['teal']}.GB * task.attempt }}
    }}
    withName: 'CIDER' {{
        cpus   = {n['cider']}
        memory = {{ {m['cider']}.GB * task.attempt }}
    }}
    withName: 'LILAC' {{
        cpus   = {n['lilac']}
        memory = {{ {m['lilac']}.GB * task.attempt }}
    }}
    withName: '.*:PAVE_ANNOTATION:PAVE_SOMATIC' {{
        cpus   = {n['pave_somatic']}
        memory = {{ {m['pave_somatic']}.GB * task.attempt }}
{oom_retry}    }}
    withName: '.*:PAVE_ANNOTATION:PAVE_GERMLINE' {{
        cpus   = {n['pave_germline']}
        memory = {{ {m['pave_germline']}.GB * task.attempt }}
    }}
    withName: '.*:SAGE_APPEND:SAGE_APPEND_.*' {{
        cpus   = {n['light']}
        memory = {{ {m['sage_append']}.GB * task.attempt }}
    }}
    withName: '.*:LINX_ANNOTATION:LINX_.*' {{
        cpus   = {n['tiny']}
        memory = {{ {m['linx']}.GB * task.attempt }}
{oom_retry}    }}
    withName: '.*:LINX_PLOTTING:LINX_VISUALISER' {{
        cpus   = {n['light']}
        memory = {{ {linx_vis_mem}.GB * task.attempt }}
{oom_retry}    }}
"""

    with open(out_path, "w") as fh:
        fh.write(
            "process {\n"
            "    executor = 'local'\n"
            f"    resourceLimits = [ cpus: {cpus}, memory: {limit_mem}.GB, time: 72.h ]\n"
            "\n"
            f"{retune_block}"
            "    // The pipeline's STAR_ALIGN module pins STAR 2.7.3a, but the HMF\n"
            "    // star_index (star_index-gencode_38-2.7.11b) is built with 2.7.11b.\n"
            "    // The old binary aborts on the newer genomeParameters.txt\n"
            "    // ('unrecognized parameter name \"genomeType\"'), so pin the aligner\n"
            "    // container to a matching 2.7.11b build. No-op for DNA-only runs.\n"
            "    withName: '.*:STAR_ALIGN' {\n"
            "        container = 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_8'\n"
            "    }\n"
            "\n"
            "    // bwa-mem2 | sambamba view | sambamba sort. Two failure modes seen on deep\n"
            "    // tumors, both fixed here:\n"
            "    //  1. The sort spills the stream to hundreds of temp chunks and hits the\n"
            "    //     container's open-file limit ('sambamba-sort: ... Too many open files').\n"
            "    //     Fix: raise the open-file ceiling, and give sort a larger in-memory\n"
            "    //     buffer (ext.args3, -m) so it makes fewer, larger chunks.\n"
            "    //  2. 'sambamba-sort: Unable to write to stream' on a 190x 4-lane sample\n"
            "    //     (H00006429), on both 4 TB and 6 TB, while a 12-lane 190x sample\n"
            "    //     (H00005238) succeeded. Fixed by pinning sort --tmpdir to the task work\n"
            "    //     dir, which keeps temp chunks off RAM/tmpfs. --align-max-forks can also\n"
            "    //     cap concurrency, but the trace shows shard size does NOT drive memory\n"
            "    //     (peak_rss 35.0-35.9 GB over 24 shards whose durations varied 1.8x, since\n"
            "    //     bwa-mem2's -K 100000000 fixes the batch), so --tmpdir is the real fix.\n"
            "    withName: '.*:BWAMEM2_ALIGN' {\n"
            f"{maxforks_line}"
            f"{align_retune}"
            "        containerOptions = '--ulimit nofile=65535:65535'\n"
            "        ext.args3 = '-m 16GiB --tmpdir .'\n"
            "    }\n"
            "\n"
            "    // PURPLE OOMs loading the somatic VCF on hypermutated tumor-only samples at\n"
            "    // the default label memory. Give it most of the box (resourceLimits caps it\n"
            "    // on smaller tiers, so this is safe across all VM sizes).\n"
            f"    withName: '.*:PURPLE_CALLING:PURPLE' {{\n"
            f"        memory = {{ {purple_mem}.GB * task.attempt }}\n"
            f"{purple_args_line}"
            f"{oom_retry if retune else ''}"
            "    }\n"
            f"{redux_block}"
            "}\n"
            "\n"
            "// Pin the local executor's resource pool instead of letting Nextflow autodetect it.\n"
            "// It otherwise pools Runtime.availableProcessors() and the FULL physical RAM\n"
            "// (LocalPollingMonitor.configCpus/configMem), leaving no headroom for the OS, the\n"
            "// Nextflow JVM and the Docker daemon on a swapless VM. This pool -- not\n"
            "// resourceLimits -- is what decides how many tasks run at once.\n"
            "executor {\n"
            "    name   = 'local'\n"
            f"    cpus   = {cpus}\n"
            f"    memory = '{limit_mem} GB'\n"
            "}\n"
            "\n"
            "// oncoanalyser already enables trace/report/timeline into <outdir>/pipeline_info/.\n"
            "// Pin the field list so peak_rss and %cpu are always present -- they are what the\n"
            "// sizing above is derived from, and what any future retune has to be checked against.\n"
            "trace {\n"
            "    enabled   = true\n"
            "    overwrite = true\n"
            "    fields    = 'task_id,hash,name,status,exit,attempt,submit,start,complete,"
            "duration,realtime,%cpu,%mem,peak_rss,peak_vmem,rchar,wchar,cpus,memory'\n"
            "}\n"
        )


def run(cmd, dry_run=False):
    printable = " ".join(cmd)
    if dry_run:
        print("[dry-run]", printable)
        return ""
    print("[run]", printable)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit(f"command failed ({res.returncode}):\n{res.stderr}")
    return res.stdout


def run_soft(cmd):
    """Like run(), but never aborts the script -- used for best-effort housekeeping in the supervisor."""
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"  (warning) command failed ({res.returncode}): {' '.join(cmd)}")
    return res.returncode


def vm_instance_name(run_id, gid, attempt):
    """RFC1035-safe VM name; retries get a -rN suffix so a relaunch never name-clashes."""
    base = f"onco-{sanitize(run_id)[:24]}-{sanitize(gid)}"[:58]
    return base if attempt == 1 else f"{base}-r{attempt}"


def work_disk_name(run_id, gid):
    """Stable name for a group's persistent work disk -- identical across preemption relaunches,
    so the relaunch re-attaches the same disk (and its work/ cache) rather than making a new one."""
    return f"onco-{sanitize(run_id)[:24]}-{sanitize(gid)}-work"[:63]


def disk_exists(args, disk):
    """True only if the disk definitely exists (not on a transient API error)."""
    res = subprocess.run(
        ["gcloud", "compute", "disks", "describe", disk,
         "--zone", args.zone, "--project", args.project, "--format=value(name)"],
        capture_output=True, text=True)
    return res.returncode == 0


def sanitize_label(value):
    """GCP label value: lowercase, keep [a-z0-9_-] (underscores ARE valid here), <=63 chars."""
    return re.sub(r"[^a-z0-9_-]", "-", value.lower())[:63]


def build_labels(args, gid, run_id):
    """GCP labels for a VM (+ its boot disk). Values are sanitized to GCP's lowercase [a-z0-9_-] rules."""
    labels = {
        "pipeline": "oncoanalyser",
        "group": sanitize(gid),
        "run": sanitize(run_id),
    }
    if args.cost_center:
        labels["cost_center"] = sanitize_label(args.cost_center)
    return ",".join(f"{k}={v}" for k, v in labels.items())


def create_group_vm(args, run_id, gid, spec, attempt):
    """Create (or recreate) the VM for one group.

    Returns the instance name on success, or None if the create failed. A failed `gcloud instances
    create` must NEVER abort the batch: at 60 VMs a single transient API/capacity error mid-loop would
    otherwise strand every already-launched VM outside the --wait supervisor (no spot-relaunch, no cost
    report). The caller records the None and moves on.
    """
    instance = vm_instance_name(run_id, gid, attempt)

    # The persistent work disk (device-name=work, auto-delete=no) carries /mnt/work across
    # preemptions. Create it inline on first launch; attach the surviving one on a relaunch.
    # Deciding by existence (not by attempt #) also makes a manual same-run-id re-run reuse the
    # disk. In --dry-run we never hit the API, so just show the create form.
    wd = work_disk_name(run_id, gid)
    if (not args.dry_run) and disk_exists(args, wd):
        disk_arg = ["--disk", f"name={wd},device-name=work,auto-delete=no"]
    else:
        # Hyperdisk Balanced bills capacity, IOPS and throughput independently, so provision the
        # performance the pipeline actually needs rather than buying it via capacity. pd-* has no
        # such knobs -- its performance is a fixed function of size (pd-balanced 6 IOPS/GiB and
        # 0.28 MiB/s per GiB; pd-ssd 30 and 0.48), so the flags are omitted for those.
        create = (f"name={wd},size={spec['disk_gb']}GB,type={args.disk_type},"
                  f"device-name=work,auto-delete=no")
        if args.disk_type.startswith("hyperdisk"):
            create += f",provisioned-iops={args.disk_iops},provisioned-throughput={args.disk_throughput}"
        disk_arg = ["--create-disk", create]

    cmd = [
        "gcloud", "compute", "instances", "create", instance,
        "--project", args.project,
        "--zone", args.zone,
        "--machine-type", spec["machine_type"],
        "--image-family", args.image_family,
        "--image-project", args.image_project,
        "--boot-disk-size", f"{BOOT_DISK_GB}GB",
        "--boot-disk-type", args.disk_type,
        *disk_arg,
        "--scopes", "https://www.googleapis.com/auth/cloud-platform",
        "--metadata-from-file",
        f"startup-script={spec['startup_path']},shutdown-script={spec['shutdown_path']}",
        # labels must be lowercase [a-z0-9_-]; sanitize so a custom (e.g. mixed-case) run-id is valid
        "--labels", build_labels(args, gid, run_id),
    ]
    if args.service_account:
        cmd += ["--service-account", args.service_account]
    if args.spot:
        cmd += ["--provisioning-model", "SPOT", "--instance-termination-action", "DELETE"]

    if args.dry_run:
        print("[dry-run]", " ".join(cmd))
        return instance
    print("[run]", " ".join(cmd))
    # On a preemption relaunch the work disk may still read as "in use" for a few seconds while
    # the just-preempted VM finishes deleting (the disk detaches as part of that). Retry ONLY
    # that transient case so a fast relaunch doesn't wrongly trip _LAUNCH_FAILED; any other
    # create error fails immediately as before.
    for tries_left in range(5, -1, -1):
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode == 0:
            return instance
        err = (res.stderr or "").strip()
        if tries_left and any(s in err.lower() for s in ("already being used", "in use", "resourceinuse")):
            print(f"  (info) work disk for '{gid}' still attached to the preempted VM; "
                  f"retrying create in 15s ({tries_left} left)")
            time.sleep(15)
            continue
        # A non-zero exit does NOT always mean the VM wasn't created: `gcloud instances create` can
        # report a client-side/transient error AFTER the instance is already up server-side. Treating
        # that as a launch failure wrongly trips _LAUNCH_FAILED and drops the group from --wait
        # supervision even though the VM is running. Confirm with a describe before giving up.
        chk = subprocess.run(
            ["gcloud", "compute", "instances", "describe", instance,
             "--zone", args.zone, "--project", args.project, "--format=value(name)"],
            capture_output=True, text=True)
        if chk.returncode == 0:
            print(f"  (info) create reported an error for '{gid}' (exit {res.returncode}) but the VM "
                  f"exists server-side; treating as launched: {err}")
            return instance
        print(f"  (warning) VM create FAILED for group '{gid}' ({res.returncode}): {err}")
        return None


def instance_gone(args, instance):
    """True only if the instance definitely does not exist (not on a transient API error)."""
    res = subprocess.run(
        ["gcloud", "compute", "instances", "describe", instance,
         "--zone", args.zone, "--project", args.project, "--format=value(name)"],
        capture_output=True, text=True)
    if res.returncode == 0:
        return False
    return "not found" in (res.stderr or "").lower()


def delete_work_disk(args, run_id, gid):
    """Best-effort reclaim of a group's persistent work disk once the group is terminal.

    On SUCCESS / FAILED-without-keep the self-deleting VM already flips the disk to auto-delete
    and takes it down with the instance, so here we usually find it gone (or briefly still
    attached to a self-deleting VM) -- both are expected and stay silent. The case this call
    really covers is GAVE_UP_PREEMPTED / a failed relaunch, where no VM is left to clean it up.
    """
    wd = work_disk_name(run_id, gid)
    res = subprocess.run(["gcloud", "compute", "disks", "delete", wd,
                          "--zone", args.zone, "--project", args.project, "--quiet"],
                         capture_output=True, text=True)
    if res.returncode != 0:
        err = (res.stderr or "").lower()
        if not any(s in err for s in ("not found", "in use", "already being used", "resourceinuse")):
            print(f"  (warning) could not delete work disk {wd}: {(res.stderr or '').strip()}")


def run_name_prefix(run_id):
    """Common prefix of every VM / disk name derived from a run-id (see vm_instance_name)."""
    return f"onco-{sanitize(run_id)[:24]}-"


def gcloud_json(cmd):
    """Run a gcloud list/describe with --format=json; [] on any error (teardown must not abort)."""
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"  (warning) command failed ({res.returncode}): {' '.join(cmd)}\n"
              f"            {(res.stderr or '').strip()}")
        return []
    try:
        return json.loads(res.stdout or "[]")
    except json.JSONDecodeError:
        return []


def _zone_of(resource):
    """gcloud returns zone as a full URL; we only ever want the short name."""
    return (resource.get("zone") or "").rstrip("/").rsplit("/", 1)[-1]


def find_run_instances(project, run_id, groups=None):
    """Every VM of a run, across zones and retry attempts (-rN).

    Matched by label (run=<sanitized run-id>, set at create) OR by the derived name prefix, so a
    VM whose labels were stripped/edited is still found. `groups` limits it to those group ids.
    """
    prefix = run_name_prefix(run_id)
    filt = (f"(labels.pipeline=oncoanalyser AND labels.run={sanitize(run_id)}) "
            f"OR (name~^{prefix})")
    items = gcloud_json(["gcloud", "compute", "instances", "list", "--project", project,
                         "--filter", filt, "--format=json(name,zone,status,labels)"])
    out = []
    for it in items:
        name = it.get("name", "")
        gid = (it.get("labels") or {}).get("group") or name[len(prefix):].split("-r")[0]
        if groups and gid not in {sanitize(g) for g in groups}:
            continue
        out.append({"name": name, "zone": _zone_of(it), "status": it.get("status", ""), "group": gid})
    return sorted(out, key=lambda d: d["name"])


def disk_group_of_run(run_id, name):
    """The group id if `name` is a disk this run created, else None.

    Reconstructing the name from the derived group id (rather than trusting the prefix) keeps a
    run-id that shares its first 24 chars with another run from matching that run's disks.
    """
    prefix = run_name_prefix(run_id)
    if not name.startswith(prefix):
        return None
    rest = name[len(prefix):]
    if rest.endswith("-work"):
        gid = rest[:-len("-work")]
        return gid if work_disk_name(run_id, gid) == name else None
    # a boot disk (same name as its VM) -- normally auto-deleted with the instance, but reclaim
    # any that outlived it
    m = re.fullmatch(r"(.+?)(?:-r(\d+))?", rest)
    gid, attempt = m.group(1), int(m.group(2) or 1)
    return gid if vm_instance_name(run_id, gid, attempt) == name else None


def find_run_disks(project, run_id, groups=None):
    """Every disk of a run (persistent work disks + any orphaned boot disk), across zones."""
    items = gcloud_json(["gcloud", "compute", "disks", "list", "--project", project,
                         "--filter", f"name~^{run_name_prefix(run_id)}",
                         "--format=json(name,zone,sizeGb,status,users)"])
    out = []
    for it in items:
        name = it.get("name", "")
        gid = disk_group_of_run(run_id, name)
        if gid is None:
            continue
        if groups and gid not in {sanitize(g) for g in groups}:
            continue
        out.append({"name": name, "zone": _zone_of(it), "size_gb": it.get("sizeGb", "?"),
                    "group": gid, "in_use": bool(it.get("users"))})
    return sorted(out, key=lambda d: d["name"])


def gcs_rm_tree(uri):
    """Recursively delete a GCS prefix; a prefix that doesn't exist is fine, not an error."""
    res = subprocess.run(["gcloud", "storage", "rm", "--recursive", uri.rstrip("/") + "/"],
                         capture_output=True, text=True)
    if res.returncode == 0:
        print(f"  deleted {uri.rstrip('/')}/")
        return
    err = (res.stderr or "").strip()
    if "matched no objects" in err.lower() or "not found" in err.lower():
        print(f"  (nothing at {uri.rstrip('/')}/)")
        return
    print(f"  (warning) could not delete {uri.rstrip('/')}/: {err}")


def delete_run(args, run_id):
    """Tear down a run: its VMs, their persistent work disks and (optionally) its GCS folders.

    Deliberately discovery-based rather than manifest-based -- it finds whatever the run actually
    left in the project, including VMs from a manifest you no longer have, so an accidental launch
    can be cancelled with just its run-id.
    """
    groups = [g.strip() for g in args.delete_groups.split(",") if g.strip()]
    instances = find_run_instances(args.project, run_id, groups or None)
    disks = find_run_disks(args.project, run_id, groups or None)
    run_root = build_run_root(args.output_bucket, args.runs_prefix, run_id)
    staging_root = ""
    if args.staging_bucket:
        parts = [args.staging_bucket] + ([args.staging_prefix] if args.staging_prefix else []) + [run_id]
        staging_root = gcs_join(*parts)

    scope = f" (groups: {', '.join(groups)})" if groups else ""
    print(f"Teardown plan for run-id '{run_id}'{scope} in project {args.project}:\n")
    if instances:
        print("  VMs to delete:")
        for i in instances:
            print(f"    {i['name']:40s} {i['zone']:16s} {i['status']}")
    else:
        print("  VMs to delete: none found")
    if disks:
        print("  Disks to delete:")
        for d in disks:
            note = " (still attached; deleted after its VM)" if d["in_use"] else ""
            print(f"    {d['name']:40s} {d['zone']:16s} {d['size_gb']} GB{note}")
    else:
        print("  Disks to delete: none found")
    if args.delete_outputs:
        targets = [gcs_join(run_root, sanitize(g)) for g in groups] or [run_root]
        if staging_root:
            targets += [gcs_join(staging_root, g) for g in groups] or [staging_root]
        print("  GCS to delete (RECURSIVE, irreversible):")
        for t in targets:
            print(f"    {t}/")
    else:
        print(f"  GCS: kept ({run_root}/" + (f", {staging_root}/" if staging_root else "") +
              ") -- pass --delete-outputs to remove it too")

    if not (instances or disks or args.delete_outputs):
        print("\nNothing to do.")
        return
    if args.dry_run:
        print("\n[dry-run] no resources were touched.")
        return
    if not args.yes:
        if input(f"\nType the run-id to confirm deletion ('{run_id}'): ").strip() != run_id:
            print("Aborted; nothing was deleted.")
            return

    # A --wait supervisor in another terminal treats a vanished VM as a preemption and relaunches
    # it. Drop a terminal marker first so it stops instead (and warn, since a supervisor started
    # before this marker existed only recognises _SUCCESS / _FAILED).
    gids = {i["group"] for i in instances} | {d["group"] for d in disks}
    for gid in sorted(gids):
        write_marker(gcs_join(run_root, gid, "_CANCELLED"),
                     f"run {run_id} / group {gid} deleted via --delete")

    print("\nDeleting VMs...")
    by_zone = OrderedDict()
    for i in instances:
        by_zone.setdefault(i["zone"], []).append(i["name"])
    for zone, names in by_zone.items():
        run_soft(["gcloud", "compute", "instances", "delete", *names,
                  "--zone", zone, "--project", args.project, "--quiet"])
        print(f"  deleted in {zone}: {', '.join(names)}")

    print("Deleting disks...")
    for d in disks:
        res = subprocess.run(["gcloud", "compute", "disks", "delete", d["name"],
                              "--zone", d["zone"], "--project", args.project, "--quiet"],
                             capture_output=True, text=True)
        if res.returncode == 0:
            print(f"  deleted {d['name']} ({d['size_gb']} GB, {d['zone']})")
        elif "not found" in (res.stderr or "").lower():
            print(f"  (already gone) {d['name']}")  # boot disks go with their VM
        else:
            print(f"  (warning) could not delete {d['name']}: {(res.stderr or '').strip()}")

    if args.delete_outputs:
        print("Deleting GCS folders...")
        for t in ([gcs_join(run_root, sanitize(g)) for g in groups] or [run_root]):
            gcs_rm_tree(t)
        if staging_root:
            for t in ([gcs_join(staging_root, g) for g in groups] or [staging_root]):
                gcs_rm_tree(t)

    left_vms = find_run_instances(args.project, run_id, groups or None)
    left_disks = find_run_disks(args.project, run_id, groups or None)
    if left_vms or left_disks:
        print("\nStill present after teardown (a --wait supervisor may have relaunched the group; "
              "stop it with Ctrl-C and re-run --delete):")
        for i in left_vms:
            print(f"  VM   {i['name']} ({i['zone']}, {i['status']})")
        for d in left_disks:
            print(f"  disk {d['name']} ({d['zone']})")
    else:
        print(f"\nDone: no VMs or disks remain for run-id '{run_id}'{scope}. "
              f"Billing for it stops now.")
        if not args.delete_outputs:
            print(f"  GCS kept: {run_root}/" + (f" and {staging_root}/" if staging_root else ""))


def write_marker(uri, text):
    """Best-effort: upload a tiny status marker object to GCS."""
    fd, path = tempfile.mkstemp()
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(text + "\n")
        run_soft(["gcloud", "storage", "cp", path, uri])
    finally:
        os.unlink(path)


def gcs_join(*parts):
    return "/".join(p.strip("/") for p in parts) if not parts[0].endswith("//") \
        else parts[0] + "/".join(p.strip("/") for p in parts[1:])


# Terminal status markers an earlier attempt of a group may have left behind under
# gs://.../<run-id>/<group>/. They are status, never data.
STALE_MARKERS = ("_FAILED", "_SUCCESS", "_CANCELLED", "_PREEMPTED",
                 "_GAVE_UP_PREEMPTED", "_LAUNCH_FAILED")


def existing_output(run_root, gid):
    """What oncoanalyser results already sit under gs://.../<run-id>/<group>/output/ ?

    The startup script rsyncs ./output there only AFTER nextflow exits, so an empty/absent prefix
    means the run never got far enough to produce anything -- the only case where wiping a terminal
    marker and relaunching is safe. Anything present is real pipeline output that a from-scratch
    relaunch would rsync over (a _FAILED VM self-deletes its work disk, so the next attempt has no
    -resume cache and recomputes everything).

    Returns (state, entries):
        True,  [names]  output exists -- do NOT overwrite
        False, []       nothing there -- safe to clear markers and relaunch
        None,  []       the listing itself failed; unknown, so treat it as 'exists' (fail closed)
    """
    res = subprocess.run(["gcloud", "storage", "ls", f"{gcs_join(run_root, gid, 'output')}/"],
                         capture_output=True, text=True)
    if res.returncode == 0:
        entries = [ln.strip().rstrip("/").rsplit("/", 1)[-1]
                   for ln in res.stdout.splitlines() if ln.strip()]
        return (bool(entries), entries)
    err = (res.stderr or "").lower()
    if "matched no objects" in err or "not found" in err:
        return (False, [])
    return (None, [])


# Markers whose removal implies discarding a finished attempt's verdict. If output/ is already
# populated, clearing one of these and relaunching would overwrite real results, so it is refused.
PROTECTED_MARKERS = ("_FAILED", "_SUCCESS", "_GAVE_UP_PREEMPTED")


def clear_stale_markers(run_root, gid, dry_run=False, force=False):
    """Wipe a group's terminal markers just before (re)launching it, so the new VM starts clean.

    Relaunching into the same --run-id keeps the same output folder, and a marker from the previous
    attempt outlives it. Without this, --wait reads that stale marker on its first poll and calls the
    brand-new VM FAILED (or SUCCESS) while it is happily running -- and a leftover _PREEMPTED would
    make --spot --wait relaunch the group immediately, ending up with two VMs on one group. The
    watchdog reads the same markers, hence the false failure alerts.

    Safeguard: a _FAILED/_SUCCESS marker is only stale if the attempt that wrote it left nothing
    behind. If output/ already holds results, the group is BLOCKED instead of cleared -- relaunch it
    under a fresh --run-id, or pass --force-over-output to overwrite deliberately.

    Returns (cleared, ok, blocked): the marker names removed, False if any removal failed, and a
    reason string when the group must not be launched (nothing was deleted in that case).
    """
    base = gcs_join(run_root, gid)
    res = subprocess.run(["gcloud", "storage", "ls", f"{base}/"], capture_output=True, text=True)
    if res.returncode != 0:
        return [], True, ""  # folder does not exist yet: a first launch has nothing to clean
    listed = {line.strip() for line in res.stdout.splitlines()}
    cleared = [m for m in STALE_MARKERS if f"{base}/{m}" in listed]
    if not cleared:
        return [], True, ""

    protected = [m for m in cleared if m in PROTECTED_MARKERS]
    if protected and not force:
        has_output, entries = existing_output(run_root, gid)
        if has_output is None:
            return [], True, (f"could not list {base}/output/ to check for existing results; "
                              f"refusing to clear {', '.join(protected)} while that is unknown")
        if has_output:
            shown = ", ".join(sorted(entries)[:6]) + (", ..." if len(entries) > 6 else "")
            return [], True, (f"{base}/output/ already holds pipeline output ({shown}) alongside "
                              f"{', '.join(protected)}; a relaunch would overwrite it")

    uris = [f"{base}/{m}" for m in cleared]
    if dry_run:
        print("[dry-run]", " ".join(["gcloud", "storage", "rm"] + uris))
        return cleared, True, ""
    ok = run_soft(["gcloud", "storage", "rm"] + uris) == 0
    return cleared, ok, ""


def build_run_root(output_bucket, runs_prefix, run_id):
    """gs://<bucket>[/<runs_prefix>]/<run_id> — keeps run outputs out of the bucket root."""
    parts = [output_bucket]
    if runs_prefix:
        parts.append(runs_prefix)
    parts.append(run_id)
    return gcs_join(*parts)


def main():
    ap = argparse.ArgumentParser(description="Launch one oncoanalyser VM per sample.")
    ap.add_argument("--manifest", help="manifest CSV (samplesheet + optional coverage/vm_type). Required unless --report.")
    ap.add_argument("--project", required=True)
    ap.add_argument("--zone", default="", help="e.g. europe-west4-a. Required to launch; not needed "
                                               "for --report or --delete (which find resources in "
                                               "any zone).")
    ap.add_argument("--output-bucket", required=True, help="gs://bucket or gs://bucket/prefix")
    ap.add_argument("--staging-bucket", help="gs://bucket or gs://bucket/prefix for generated configs. Required unless --report.")
    ap.add_argument("--report", default="", metavar="RUN_ID",
                    help="don't launch; print the cost summary for a previous run-id and exit")
    ap.add_argument("--delete", default="", metavar="RUN_ID",
                    help="don't launch; tear down an existing run and exit. Deletes its VMs (all "
                         "groups, all retry attempts, any zone) and their persistent work disks, so "
                         "billing stops. Discovery is by GCP label/name, so no manifest is needed. "
                         "GCS outputs are KEPT unless --delete-outputs is given. Asks for "
                         "confirmation unless --yes; combine with --dry-run to preview.")
    ap.add_argument("--delete-groups", default="", metavar="G1,G2",
                    help="with --delete: limit the teardown to these group ids "
                         "(default: every group of the run)")
    ap.add_argument("--delete-outputs", action="store_true",
                    help="with --delete: ALSO recursively delete the run's GCS output folder "
                         "(gs://<bucket>/<runs-prefix>/<run-id>/) and, if --staging-bucket is given, "
                         "its staging folder. Irreversible.")
    ap.add_argument("--yes", "-y", action="store_true",
                    help="skip the interactive confirmation prompt of --delete")
    ap.add_argument("--runs-prefix", default="runs",
                    help="subfolder under the output bucket for run artifacts, keeping them out of "
                         "the bucket root (default: 'runs' -> gs://<bucket>/runs/<run-id>/...). "
                         "Pass '' to write at the bucket root.")
    ap.add_argument("--staging-prefix", default="staging",
                    help="subfolder under the staging bucket for the generated samplesheet/config, "
                         "keeping them out of the bucket root (default: 'staging' -> "
                         "gs://<staging-bucket>/staging/<run-id>/<group>/...). Pass '' to write at "
                         "the bucket root.")
    ap.add_argument("--run-id", default="",
                    help="logical name for this run's output folder (gs://<bucket>/<runs-prefix>/"
                         "<run-id>/...) and its --report key, e.g. 'hcc-celllines-fastq'. "
                         "Default: an auto timestamp+random id that can never collide. Allowed: "
                         "1-50 chars of letters/digits/.-_ . Reusing an existing run-id adds to / "
                         "overwrites that folder.")
    ap.add_argument("--force-over-output", action="store_true",
                    help="relaunch a group even though its folder already holds pipeline output "
                         f"alongside a terminal marker ({' '.join(PROTECTED_MARKERS)}). Without "
                         "this, such a group is BLOCKED and skipped so the existing results are "
                         "not overwritten; prefer a fresh --run-id.")
    ap.add_argument("--service-account", default="",
                    help="VM service account email. Omit to use the project's default "
                         "Compute Engine service account (which inherits project IAM).")
    ap.add_argument("--cost-center", default="innovation",
                    help="value for the 'cost_center' GCP label attached to every VM (+ its boot "
                         "disk), so spend is attributed in billing reports / BigQuery export. "
                         "Default: 'innovation'. Lowercased to satisfy GCP label rules. Pass '' to "
                         "omit the label.")
    ap.add_argument("--reference-config", default="", help="gs:// path to a staged reference_data.config (optional)")
    ap.add_argument("--panel-config", default="",
                    help="gs:// path to a staged panel config defining params.panel_data_paths and "
                         "params.ref_data_panel_data_path. Required for --mode targeted with a custom "
                         "panel; layered on top of --reference-config as a second -config.")
    ap.add_argument("--revision", default="dev-3.0.0-beta.15", help="oncoanalyser git branch/tag/commit (default: dev-3.0.0-beta.15)")
    ap.add_argument("--genome", default="GRCh38_hmf")
    ap.add_argument("--mode", default="wgts")
    ap.add_argument("--nextflow-version", default="25.04.8",
                    help="pinned NXF_VER on the VM. Must be in [25.04, 26.04): dev-3.0.0's bundled "
                         "nf-schema needs >=25.04.0, and >=26.04 defaults to the v2 config parser which "
                         "dev-3.0.0 fails to parse. 25.04.x or 25.10.x both work.")
    ap.add_argument("--sequencing-type", default="ILLUMINA", choices=sorted(SEQUENCING_TYPES),
                    help="default platform; overridden per-group by a 'sequencing_type' manifest column")
    ap.add_argument("--extra-args", default="", help="extra oncoanalyser args, e.g. '--processes_exclude virusinterpreter,orange'")
    ap.add_argument("--redux-version", default="", help="pin the hmftools-redux container to this release (e.g. '2.0.1'); default uses the pipeline's own pin")
    ap.add_argument("--purple-purity-range", default="", metavar="MIN,MAX",
                    help="constrain PURPLE's fit grid to this purity range, e.g. '0.95,1.0'. Injects "
                         "'-min_purity MIN -max_purity MAX' as ext.args on the PURPLE process. Use to "
                         "force a fit off a half-purity alias (e.g. a pure cell line stuck at 0.50). "
                         "Requires 0 <= MIN <= MAX <= 1. Omit to leave PURPLE's defaults (0.08-1.0).")
    ap.add_argument("--image-family", default="ubuntu-2204-lts")
    ap.add_argument("--image-project", default="ubuntu-os-cloud")
    ap.add_argument("--disk-type", default="hyperdisk-balanced",
                    help="boot + work disk type (default: hyperdisk-balanced). REQUIRED for N4D/C4D "
                         "machines, which do not support pd-* at all. Hyperdisk capacity is also "
                         "cheaper per GB than pd-balanced, with IOPS/throughput bought separately "
                         "(--disk-iops/--disk-throughput) instead of implied by size. Pass pd-balanced "
                         "or pd-ssd only when also pinning an N2D/E2 machine via --machine-type.")
    ap.add_argument("--disk-iops", type=int, default=20000,
                    help="provisioned IOPS on the work disk when --disk-type is hyperdisk-* "
                         "(default 20000; range 3000-160000). For reference the pd-balanced disks "
                         "this replaces gave only 6 IOPS/GiB, i.e. 12000 on a 2 TB disk, shared by "
                         "every tool doing indexed random reads over the BAM. Only the excess over "
                         "the free 3000 baseline is billed.")
    ap.add_argument("--disk-throughput", type=int, default=800, metavar="MIBPS",
                    help="provisioned throughput MiB/s on the work disk for hyperdisk-* "
                         "(default 800; range 140-2400, and capped by the machine: n4d-*-48 and "
                         "above allow 2400). Only the excess over the free 140 baseline is billed.")
    ap.add_argument("--machine-type", default="", metavar="NAME",
                    help="override the auto-selected machine for every group (default: "
                         f"{DEFAULT_MACHINE[0]} from BAM/CRAM, {FASTQ_MACHINE[0]} from FASTQ). Its "
                         "cpu/memory shape is read from the Compute API, so the generated "
                         "resources.config always matches the box. A per-row `vm_type` manifest "
                         "column still beats this.")
    ap.add_argument("--purple-mem-gb", type=int, default=96,
                    help="memory for PURPLE (default 96 GB). PURPLE loads the whole somatic VCF into "
                         "the JVM heap and OOMs on hypermutated tumor-only samples; raise this if you "
                         "hit that. Clamped to the box by resourceLimits.")
    ap.add_argument("--redux-split", default="auto", metavar="AUTO|T:N",
                    help="how to split cores between the tumor and normal REDUX tasks, which run "
                         "concurrently (default: auto). 'auto' allocates cores proportional to the "
                         "measured tumor:normal work ratio -- taken from input file size, falling "
                         "back to the `coverage` column -- which lands on the usual 3:1 for a "
                         "normal pair, collapses to 1:1 for a balanced pair, and inverts if the "
                         "normal is larger. A warning is printed whenever the ratio is under "
                         f"{REDUX_RATIO_WARN_BELOW:g}. Pass '1:1' or '3:1' to pin it.")
    ap.add_argument("--isofox-reserve", default="auto", metavar="AUTO|N",
                    help="cores held back from the REDUX split for ISOFOX on groups with an rna row "
                         f"(default: auto = {ISOFOX_RESERVE_CPUS}). ISOFOX depends only on the RNA "
                         "BAM, so it is ready before REDUX starts, but the REDUX pair otherwise "
                         "claims every core and ISOFOX queues out the whole stage. The reserve keeps "
                         "it admissible so the two run concurrently. '0' restores the old behaviour "
                         "(no reserve). No effect on DNA-only groups.")
    ap.add_argument("--no-retune", action="store_true",
                    help="skip the per-process cpu/memory overrides and fall back to oncoanalyser's "
                         "stock nf-core labels (12 cpu / 72 GB max). Use to A/B the retune, or if a "
                         "pipeline bump changes process names and a selector stops matching.")
    ap.add_argument("--spot", action="store_true", help="use Spot VMs (cheaper; can be preempted)")
    ap.add_argument("--max-preemption-retries", type=int, default=2,
                    help="for Spot VMs with --wait: how many times to relaunch a group whose VM is "
                         "preempted before giving up (default: 2 -> up to 3 attempts). A real pipeline "
                         "failure (_FAILED) is never retried.")
    ap.add_argument("--align-max-forks", type=int, default=0, metavar="N",
                    help="cap concurrent BWAMEM2_ALIGN shards (default 0 = no cap / full speed). "
                         "Set to 2 for few-lane deep samples (e.g. a 190x sample in <=4 lanes) where "
                         "3 concurrent sambamba sorts exhaust RAM and crash with 'Unable to write to stream'.")
    ap.add_argument("--keep-on-failure", action="store_true", help="do not delete VM if the run fails")
    ap.add_argument("--no-localize-inputs", action="store_true",
                    help="disable copying gs:// inputs to local disk first. Localizing (the default) "
                         "uses gcloud's parallel/resumable download and avoids Nextflow's single-stream "
                         "foreign-file copier stalling on large CRAM/BAM inputs.")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--wait", action="store_true", help="poll the output bucket until all samples finish")
    args = ap.parse_args()

    purple_ext_args = parse_purple_purity_range(args.purple_purity_range)
    isofox_reserve_arg = parse_isofox_reserve(args.isofox_reserve)

    if args.report:
        report_costs(build_run_root(args.output_bucket, args.runs_prefix, args.report), args.report)
        return
    if args.delete:
        delete_run(args, validate_run_id(args.delete))
        return
    if not args.manifest or not args.staging_bucket:
        sys.exit("--manifest and --staging-bucket are required unless --report or --delete is used")
    if not args.zone:
        sys.exit("--zone is required to launch a run (e.g. --zone europe-west4-a)")

    # Fail fast on machine/disk combinations GCP will reject, rather than after N VM creates.
    if args.machine_type:
        machine_shape(args.machine_type, args.zone, args.project)   # aborts if it doesn't exist
    hyperdisk_only = ("n4d-", "n4-", "c4d-", "c4-", "c3-")
    pinned = args.machine_type or DEFAULT_MACHINE[0]
    if args.disk_type.startswith("pd-") and pinned.startswith(hyperdisk_only):
        sys.exit(f"--disk-type {args.disk_type} cannot be used with {pinned}: the N4/N4D/C4/C4D "
                 f"families support Hyperdisk only. Use --disk-type hyperdisk-balanced, or pin an "
                 f"N2D machine with --machine-type.")
    if args.disk_type.startswith("hyperdisk"):
        if not 3000 <= args.disk_iops <= 160000:
            sys.exit("--disk-iops must be between 3000 and 160000 for hyperdisk-balanced")
        if not 140 <= args.disk_throughput <= 2400:
            sys.exit("--disk-throughput must be between 140 and 2400 MiB/s for hyperdisk-balanced")

    if args.run_id:
        run_id = validate_run_id(args.run_id)
    else:
        # Timestamp + random suffix => globally unique; this folder can never collide with a prior run.
        run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%d-%H%M%S") + "-" + secrets.token_hex(2)
    run_root = build_run_root(args.output_bucket, args.runs_prefix, run_id)

    # A custom run-id can name a folder that already exists; warn but don't block (re-runs are legit).
    if args.run_id and not args.dry_run and marker_exists(f"{run_root}/"):
        print(f"WARNING: {run_root}/ already exists -- its outputs may be overwritten or mixed "
              f"with this run. Use a different --run-id to keep runs separate.")
        print(f"  (each relaunched group's stale status markers ({' '.join(STALE_MARKERS)}) are "
              f"cleared just before its VM starts, so this attempt is judged on its own result.)")

    groups, sheet_cols = read_manifest(args.manifest)
    print(f"run-id={run_id}; outputs -> {run_root}/; {len(groups)} group(s): {', '.join(groups)}")

    # Size every DNA input up front in one batched pass, so the REDUX core split can be derived
    # from the real tumor:normal work ratio. Best-effort: if the objects are unreadable we fall
    # back to the coverage column, and failing that to the default 3:1 -- a launch never blocks
    # on this. Skipped entirely when the split is pinned, since then nothing consumes the sizes.
    input_sizes = {}
    if args.redux_split.strip().lower() == "auto" and not args.no_retune:
        all_paths = [p for rows in groups.values()
                     for e in group_dna_inputs(rows).values() for p in e["paths"]]
        if all_paths:
            input_sizes = gcs_object_sizes(all_paths)
            if not input_sizes:
                print("  (warning) could not read any input object sizes; REDUX split will fall "
                      "back to the coverage column or the 3:1 default")

    specs = OrderedDict()
    failed_launches = []
    blocked_launches = []
    tmpdir = tempfile.mkdtemp(prefix="oncoanalyser-")

    for gid, rows in groups.items():
        machine_type, cpus, mem_gb, disk_gb = pick_tier(rows, args.zone, args.project,
                                                        args.machine_type)

        # REDUX core split, from the real tumor:normal work ratio. Silent in the normal case;
        # loud when the assumption behind the default 3:1 has broken.
        ratio, basis, detail = redux_work_ratio(rows, input_sizes)
        # RNA groups hold back ISOFOX's cores so it can run alongside REDUX instead of behind it.
        # isofox_reserve stays 0 for DNA-only groups, which keeps their config byte-identical.
        isofox_reserve = isofox_reserve_arg if group_has_rna(rows) else 0
        isofox_reserve = min(isofox_reserve, max(0, cpus - 8))   # never starve REDUX itself
        redux_split = parse_redux_split(args.redux_split, cpus, ratio,
                                        budget=cpus - isofox_reserve)
        if isofox_reserve:
            print(f"  [{gid}] RNA group -> {machine_type}: REDUX {redux_split[0]}/{redux_split[1]} "
                  f"cpus + {isofox_reserve} reserved for ISOFOX, running concurrently")
        if args.redux_split.strip().lower() != "auto":
            print(f"  [{gid}] REDUX split pinned to {args.redux_split} "
                  f"-> tumor {redux_split[0]} / other {redux_split[1]} cpus")
        elif ratio is None and basis == "single-sample":
            print(f"  [{gid}] {detail}; REDUX gets all {redux_split[0]} cpus")
        elif ratio is None:
            print(f"  (warning) [{gid}] cannot determine the tumor:normal size ratio ({detail}); "
                  f"using the default {redux_split[0]}:{redux_split[1]} REDUX split. Add a "
                  f"`coverage` column, or pin it with --redux-split 1:1.")
        elif ratio < REDUX_RATIO_WARN_BELOW:
            print(f"  (warning) [{gid}] normal is unusually large relative to the tumor: "
                  f"ratio {ratio:.2f} by {basis} ({detail}). The default ~3:1 REDUX split assumes "
                  f">= {REDUX_RATIO_WARN_BELOW:g}, so cores were re-balanced to "
                  f"{redux_split[0]}:{redux_split[1]} to keep both tasks finishing together. "
                  f"Override with --redux-split 1:1 (or 3:1 to force the old behaviour).")
        seq_type = pick_sequencing_type(rows, args.sequencing_type)
        # keep generated configs out of the staging bucket root: gs://<bucket>/<staging-prefix>/<run-id>/<group>/
        staging_parts = [args.staging_bucket]
        if args.staging_prefix:
            staging_parts.append(args.staging_prefix)
        staging_parts += [run_id, gid]
        staging_prefix = gcs_join(*staging_parts)
        # disk cost spans both the persistent work disk (disk_gb, which carries the provisioned
        # Hyperdisk IOPS/throughput) and the fixed boot disk (baseline performance only)
        hourly = estimate_hourly(machine_type, cpus, mem_gb, disk_gb, args.disk_type,
                                 region_of(args.zone), args.spot, iops=args.disk_iops,
                                 throughput=args.disk_throughput, boot_gb=BOOT_DISK_GB)
        cost_record_uri = gcs_join(run_root, f"{sanitize(gid)}.cost.json")

        # generate + upload samplesheet and resources.config
        ss = os.path.join(tmpdir, f"{sanitize(gid)}.samplesheet.csv")
        rc = os.path.join(tmpdir, f"{sanitize(gid)}.resources.config")
        write_samplesheet(rows, sheet_cols, ss)
        write_resources_config(cpus, mem_gb, rc, align_max_forks=args.align_max_forks,
                               purple_ext_args=purple_ext_args, redux_version=args.redux_version,
                               purple_mem_gb=args.purple_mem_gb, retune=not args.no_retune,
                               redux_split=redux_split, isofox_cpus=isofox_reserve or None)
        run(["gcloud", "storage", "cp", ss, gcs_join(staging_prefix, "samplesheet.csv")], args.dry_run)
        run(["gcloud", "storage", "cp", rc, gcs_join(staging_prefix, "resources.config")], args.dry_run)

        # render startup + shutdown scripts to local files (reused verbatim on any relaunch)
        startup = (
            STARTUP_TEMPLATE
            .replace("%%GROUP_ID%%", gid)
            .replace("%%REVISION%%", args.revision)
            .replace("%%GENOME%%", args.genome)
            .replace("%%MODE%%", args.mode)
            .replace("%%SEQUENCING_TYPE%%", seq_type)
            .replace("%%NXF_VER%%", args.nextflow_version)
            .replace("%%RUN_ROOT%%", run_root)
            .replace("%%RUN_ID%%", run_id)
            .replace("%%STAGING_PREFIX%%", staging_prefix)
            .replace("%%REFERENCE_CONFIG_GCS%%", args.reference_config)
            .replace("%%PANEL_CONFIG_GCS%%", args.panel_config)
            .replace("%%ZONE%%", args.zone)
            .replace("%%KEEP_ON_FAILURE%%", "true" if args.keep_on_failure else "false")
            .replace("%%EXTRA_ARGS%%", args.extra_args)
            .replace("%%LOCALIZE_INPUTS%%", "false" if args.no_localize_inputs else "true")
            .replace("%%MACHINE_TYPE%%", machine_type)
            .replace("%%DISK_GB%%", str(disk_gb))
            .replace("%%SPOT%%", "true" if args.spot else "false")
            .replace("%%HOURLY_RATE%%", "" if hourly is None else str(hourly))
            .replace("%%CURRENCY%%", CURRENCY)
            .replace("%%COST_RECORD_URI%%", cost_record_uri)
        )
        startup_path = os.path.join(tmpdir, f"{sanitize(gid)}.startup.sh")
        with open(startup_path, "w") as fh:
            fh.write(startup)
        shutdown_path = os.path.join(tmpdir, f"{sanitize(gid)}.shutdown.sh")
        with open(shutdown_path, "w") as fh:
            fh.write(SHUTDOWN_TEMPLATE.replace("%%RUN_PREFIX%%", gcs_join(run_root, gid)))

        spec = {
            "machine_type": machine_type, "cpus": cpus, "mem_gb": mem_gb, "disk_gb": disk_gb,
            "seq_type": seq_type, "hourly": hourly,
            "startup_path": startup_path, "shutdown_path": shutdown_path,
            "attempts": 1, "instance": None,
        }
        # A relaunch into an existing run folder must not inherit the previous attempt's verdict:
        # clear the group's terminal markers immediately before the VM comes up.
        if args.run_id:
            cleared, markers_ok, blocked = clear_stale_markers(
                run_root, gid, args.dry_run, force=args.force_over_output)
            if blocked:
                blocked_launches.append(gid)
                print(f"  -> {gid}: BLOCKED, not launched -- {blocked}. Relaunch under a fresh "
                      f"--run-id, or pass --force-over-output to overwrite deliberately.")
                continue
            if cleared:
                print(f"  -> {gid}: cleared stale marker(s) from a previous attempt: "
                      f"{', '.join(cleared)}")
            if not markers_ok:
                print(f"  WARNING: {gid}: could not remove {gcs_join(run_root, gid)}/ marker(s) "
                      f"{', '.join(cleared)} -- --wait and the watchdog may read them as this "
                      f"attempt's result. Delete them by hand.")

        instance = create_group_vm(args, run_id, gid, spec, attempt=1)
        if instance is None:
            failed_launches.append(gid)
            write_marker(gcs_join(run_root, gid, "_LAUNCH_FAILED"),
                         "gcloud instances create failed at launch; VM was never started")
            print(f"  -> {gid}: LAUNCH FAILED (skipped; see warning above). Continuing with the rest.")
            continue
        spec["instance"] = instance
        specs[gid] = spec

        rate_str = "rate n/a (no catalog price)" if hourly is None else f"~{hourly} {CURRENCY}/h"
        print(f"  -> {gid}: {machine_type} ({cpus} vCPU / {mem_gb} GB / {disk_gb} GB disk), "
              f"sequencing_type={seq_type}, {rate_str}, as {spec['instance']}")

    if blocked_launches:
        print(f"\nWARNING: {len(blocked_launches)} group(s) were BLOCKED and never launched (their "
              f"run folder already holds output): {', '.join(blocked_launches)}")
        print("  (nothing was deleted for these; their existing results are intact.)")
    if failed_launches:
        print(f"\nWARNING: {len(failed_launches)} group(s) FAILED to launch and were skipped: "
              f"{', '.join(failed_launches)}")
        print("  (each has a _LAUNCH_FAILED marker under its output folder; the launched groups "
              "below are unaffected and are supervised normally.)")
    if not specs:
        print("\nNo VMs were launched (every group was blocked or failed to create). Nothing to "
              "supervise; exiting.")
        return

    print("\nLaunched:")
    for gid, spec in specs.items():
        print(f"  {gid:20s} {spec['machine_type']:16s} {spec['instance']}")
    rates = [s["hourly"] for s in specs.values() if s["hourly"] is not None]
    if rates:
        print(f"\nProjected burn rate while all VMs run: ~{round(sum(rates), 2)} {CURRENCY}/hour "
              f"(estimate; run id {run_id}). Final per-run cost via:  --report {run_id}")
    markers = "_SUCCESS | _FAILED" + (" | _PREEMPTED | _GAVE_UP_PREEMPTED" if args.spot else "")
    print(f"\nStatus markers ({markers}) will appear under:")
    for gid in specs:
        print(f"  {gcs_join(run_root, gid)}/")

    if args.spot and not args.wait and not args.dry_run:
        print("\nNOTE: --spot without --wait means NO auto-relaunch. A preempted group drops a "
              "_PREEMPTED marker but is not restarted, and its persistent work disk "
              "(onco-<run>-<group>-work) is left behind for you to delete or reuse. Add --wait "
              "to auto-relaunch on preemption (resuming from that disk) and reclaim it on completion.")

    if args.wait and not args.dry_run:
        wait_for_completion(args, run_id, run_root, specs)
        report_costs(run_root, run_id)


def report_costs(run_root, run_id):
    """Fetch and sum the per-group cost records written under the run root."""
    prefix = run_root
    listing = subprocess.run(["gcloud", "storage", "ls", f"{prefix}/"],
                             capture_output=True, text=True)
    if listing.returncode != 0 or not listing.stdout.strip():
        print(f"No cost records found at {prefix}/ "
              f"(run still in progress, or VMs were killed before finishing).")
        return
    uris = [u for u in listing.stdout.split() if u.endswith(".cost.json")]
    total = 0.0
    currency = CURRENCY
    rows = []
    for uri in uris:
        out = subprocess.run(["gcloud", "storage", "cat", uri], capture_output=True, text=True)
        if out.returncode != 0:
            continue
        try:
            rec = json.loads(out.stdout)
        except json.JSONDecodeError:
            continue
        rows.append(rec)
        total += float(rec.get("est_cost", 0) or 0)
        currency = rec.get("currency", currency)

    if not rows:
        print(f"No parseable cost records under {prefix}/")
        return
    print(f"\nEstimated cost for run {run_id}  (list-price estimate, excl. egress/storage):")
    print(f"  {'group':18s} {'machine':16s} {'status':8s} {'hours':>7s} {'est_cost':>10s}")
    for r in sorted(rows, key=lambda x: x.get("group", "")):
        print(f"  {r.get('group',''):18s} {r.get('machine_type',''):16s} "
              f"{r.get('status',''):8s} {float(r.get('elapsed_hours',0)):7.2f} "
              f"{float(r.get('est_cost',0)):10.2f}")
    print(f"  {'-'*61}")
    print(f"  TOTAL{'':56s}{total:10.2f} {currency}")
    print(f"\n  Authoritative cost (incl. egress/disk/discounts), once billing data lands:")
    print(f"  filter your billing report / BigQuery export by label  run={sanitize(run_id)}")
    print(f"  (or cost_center=<value> / pipeline=oncoanalyser for the whole programme)")


def marker_exists(uri):
    res = subprocess.run(["gcloud", "storage", "ls", uri],
                         capture_output=True, text=True)
    return res.returncode == 0


def wait_for_completion(args, run_id, run_root, specs, interval=30):
    """Poll markers until every group reaches a terminal state.

    On a Spot preemption -- signalled by a _PREEMPTED marker, or (backstop) the VM vanishing with no
    terminal marker -- the group is relaunched on a fresh VM up to --max-preemption-retries times.
    A real pipeline _FAILED is never retried. Markers left by an earlier attempt of the same
    --run-id were wiped at launch (clear_stale_markers), so anything seen here belongs to this
    attempt. Terminal markers under each gs://.../<group>/:
        _SUCCESS              finished (after archived _PREEMPTED.rN if it was retried)
        _FAILED               pipeline error
        _CANCELLED            torn down with --delete
        _GAVE_UP_PREEMPTED    preempted and out of retries
        _LAUNCH_FAILED        the VM could not be (re)created at all (gcloud create error)
    """
    print("\nPolling for completion (Ctrl-C to stop; VMs keep running regardless)...")
    max_retries = args.max_preemption_retries
    pending = set(specs)
    results = {}
    while pending:
        for gid in list(pending):
            base = gcs_join(run_root, gid)
            spec = specs[gid]
            if marker_exists(f"{base}/_SUCCESS"):
                results[gid] = "SUCCESS"; pending.discard(gid)
                delete_work_disk(args, run_id, gid)  # self-deleting VM usually took it; backstop
                print(f"  [OK]   {gid} (attempt {spec['attempts']})")
                continue
            if marker_exists(f"{base}/_CANCELLED"):
                # someone tore this group down with --delete: stop supervising it, and never
                # treat the now-missing VM as a preemption to relaunch
                results[gid] = "CANCELLED"; pending.discard(gid)
                print(f"  [CANCELLED] {gid} (deleted via --delete)")
                continue
            if marker_exists(f"{base}/_FAILED"):
                results[gid] = "FAILED"; pending.discard(gid)
                # keep the disk for post-mortem when --keep-on-failure kept the VM; else reclaim it
                if not args.keep_on_failure:
                    delete_work_disk(args, run_id, gid)
                print(f"  [FAIL] {gid}  (pipeline error; see {base}/nextflow.log)")
                continue

            # Preemption: explicit marker, or (backstop) the VM disappeared without any terminal marker.
            preempt_marker = args.spot and marker_exists(f"{base}/_PREEMPTED")
            vanished = args.spot and not preempt_marker and instance_gone(args, spec["instance"])
            if not (preempt_marker or vanished):
                continue

            if preempt_marker:  # archive so the next attempt starts clean and history stays visible
                run_soft(["gcloud", "storage", "mv",
                          f"{base}/_PREEMPTED", f"{base}/_PREEMPTED.r{spec['attempts']}"])
            if spec["attempts"] > max_retries:
                results[gid] = "GAVE_UP_PREEMPTED"; pending.discard(gid)
                delete_work_disk(args, run_id, gid)  # no VM left to reclaim the surviving disk
                write_marker(f"{base}/_GAVE_UP_PREEMPTED",
                             f"gave up after {spec['attempts']} attempt(s); preempted repeatedly")
                print(f"  [GAVE UP] {gid}  preempted {spec['attempts']}x (limit {max_retries} retries)")
                continue
            spec["attempts"] += 1
            print(f"  [PREEMPTED] {gid}  relaunching (attempt {spec['attempts']}/{max_retries + 1})")
            new_instance = create_group_vm(args, run_id, gid, spec, spec["attempts"])
            if new_instance is None:
                results[gid] = "LAUNCH_FAILED"; pending.discard(gid)
                delete_work_disk(args, run_id, gid)  # no VM left to reclaim the surviving disk
                write_marker(f"{base}/_LAUNCH_FAILED",
                             "relaunch after preemption failed to create the VM")
                print(f"  [LAUNCH FAIL] {gid}  relaunch create failed; giving up on this group")
                continue
            spec["instance"] = new_instance
        if pending:
            time.sleep(interval)
    print("\nAll done:", results)
    return results


if __name__ == "__main__":
    main()