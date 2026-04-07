#!/usr/bin/env python3
import argparse
import csv
import re
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import mysql.connector  # pip3 install mysql-connector-python
import yaml  # pip3 install pyyaml


class _LiteralStr(str):
    pass


yaml.add_representer(
    _LiteralStr,
    lambda dumper, data: dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|"),
)


GCP_SECRET_PROJECT = "hmf-secrets"
DB_SECRET_NAME = "mysql-patients-sql-prod-1-reader"
DB_NAME = "hmfpatients534"

WORKFLOW_NAME = "cfdna-report"
WORKFLOW_VERSION = "0.1.0"

MASTER_SAMPLE_ID_COL = "timepoint"
MASTER_DATE_COL = "IntakeDate"


@dataclass(frozen=True)
class Config:
    execution_name: str
    wgs_purple_dir: str
    cfdna_wisp_dir: str
    cfdna_master_file: Path
    output_file: Path


def parse_args() -> Config:
    ap = argparse.ArgumentParser(
        description="Create cfdna-report execution YAML."
    )
    ap.add_argument("--execution_name", "-n", required=True,
                    help="Execution name.")
    ap.add_argument("--wgs_purple_dir", "-w", required=True,
                    help="GCP URI to WGS purple directory.")
    ap.add_argument("--cfdna_wisp_dir", "-c", required=True,
                    help="GCP URI to directory with cfdna WISP summary files.")
    ap.add_argument("--cfdna_master_file", "-m", required=True, type=Path,
                    help="Path to file with sampling dates for cfDNA samples.")
    ap.add_argument("--output_file", "-o", required=True, type=Path,
                    help="Path to output file.")
    args = ap.parse_args()
    return Config(args.execution_name, args.wgs_purple_dir, args.cfdna_wisp_dir, args.cfdna_master_file, args.output_file)


def main(config: Config) -> None:
    wgs_sample_id = wgs_sample_id_from_purple_dir(config.wgs_purple_dir)

    wgs_date, tumor_type = load_wgs_clinical_info(wgs_sample_id)

    wisp_files = list_wisp_files(config.cfdna_wisp_dir)

    relevant_wisp_files = select_relevant_files(wisp_files, wgs_sample_id)
    if not relevant_wisp_files:
        sys.exit(f"ERROR: no WISP files found for WGS sample {wgs_sample_id!r} in {config.cfdna_wisp_dir}.")
    sample_id_to_date = load_master(config.cfdna_master_file)

    manifest = build_manifest(
        wgs_sample_id, config.wgs_purple_dir, wgs_date, relevant_wisp_files, sample_id_to_date
    )
    doc = make_execution(config, manifest, tumor_type)

    with open(config.output_file, "w") as out_f:
        yaml.dump(doc, out_f, sort_keys=False)


def make_execution(config: Config, manifest: str, tumor_type: str) -> dict:
    return {
        "name": config.execution_name,
        "workflow": WORKFLOW_NAME,
        "version": WORKFLOW_VERSION,
        "params": {
            "wgs_input_bucket_uri": bucket_root(config.wgs_purple_dir),
            "cfdna_input_bucket_uri": bucket_root(config.cfdna_wisp_dir),
            "tumor_type": tumor_type,
            "sample_manifest": _LiteralStr(manifest),
        },
    }


def build_manifest(
        wgs_sample_id: str,
        wgs_purple_dir: str,
        wgs_date: str,
        wisp_files: list[str],
        sample_id_to_date: dict[str, str],
) -> str:
    lines = ["SampleId,SampleType,SamplingDate,WispSummaryTsv,PurpleDirectory"]
    lines.append(f"{wgs_sample_id},WGS,{wgs_date},,{wgs_purple_dir}")
    for wisp_uri in wisp_files:
        sample_id = sample_id_from_wisp(wisp_uri)
        if sample_id not in sample_id_to_date:
            sys.exit(f"ERROR: no date found in master file for cfDNA sample {sample_id!r}.")
        date = sample_id_to_date[sample_id]
        lines.append(f"{sample_id},CFDNA,{date},{wisp_uri},")
    return "\n".join(lines) + "\n"


def bucket_root(gcs_uri: str) -> str:
    parts = gcs_uri.replace("gs://", "").split("/")
    return f"gs://{parts[0]}"


def wgs_sample_id_from_purple_dir(wgs_purple_dir: str) -> str:
    pattern = f"{wgs_purple_dir.rstrip('/')}/*.purple.purity.tsv"
    cmd = ["gcloud", "storage", "ls", pattern]
    try:
        output = subprocess.check_output(cmd, text=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        sys.exit(f"ERROR listing purple purity files: {e.stderr or e}")
    matches = [line.strip() for line in output.splitlines() if line.strip()]
    if len(matches) != 1:
        sys.exit(f"ERROR: expected exactly one *.purple.purity.tsv in {wgs_purple_dir}, found {len(matches)}.")
    filename = matches[0].rsplit("/", 1)[-1]
    return re.sub(r"\.purple\.purity\.tsv$", "", filename)


def gcp_secret(secret_name: str) -> str:
    cmd = ["gcloud", "secrets", "versions", "access", "latest",
           f"--secret={secret_name}", f"--project={GCP_SECRET_PROJECT}"]
    try:
        return subprocess.check_output(cmd, text=True, stderr=subprocess.PIPE).strip()
    except subprocess.CalledProcessError as e:
        sys.exit(f"ERROR fetching secret {secret_name}: {e.stderr or e}")


def load_wgs_clinical_info(wgs_sample_id: str) -> tuple[str, str]:
    secret = gcp_secret(DB_SECRET_NAME)
    creds = {k: v.strip('"') for token in secret.split()
             if "=" in token for k, v in [token.split("=", 1)]}
    cnx = mysql.connector.connect(
        user=creds["user"],
        password=creds["password"],
        host=creds["host"],
        port=int(creds["port"]),
        database=DB_NAME,
    )
    cursor = cnx.cursor()
    cursor.execute(
        "SELECT sampleArrivalDate, primaryTumorType FROM clinical WHERE sampleId = %s",
        (wgs_sample_id,),
    )
    rows = cursor.fetchall()
    cursor.close()
    cnx.close()
    if not rows:
        sys.exit(f"ERROR: no clinical record found for WGS sample {wgs_sample_id!r}.")
    if len(rows) > 1:
        sys.exit(f"ERROR: multiple clinical records found for WGS sample {wgs_sample_id!r}.")
    biopsy_date, tumor_type = rows[0]
    if biopsy_date is None:
        sys.exit(f"ERROR: sampleArrivalDate is NULL in clinical record for WGS sample {wgs_sample_id!r}.")
    if tumor_type is None:
        sys.exit(f"ERROR: primaryTumorType is NULL in clinical record for WGS sample {wgs_sample_id!r}.")
    return str(biopsy_date), tumor_type


def select_relevant_files(wisp_files: list[str], wgs_sample_id: str) -> list[str]:
    return [f for f in wisp_files if f.rsplit("/", 1)[-1].split("_", 1)[0] == wgs_sample_id[:12]]


def list_wisp_files(cfdna_wisp_dir: str) -> list[str]:
    pattern = f"{cfdna_wisp_dir.rstrip('/')}/*.wisp.summary.tsv"
    cmd = ["gcloud", "storage", "ls", pattern]
    try:
        output = subprocess.check_output(cmd, text=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        sys.exit(f"ERROR listing WISP files: {e.stderr or e}")
    return [line.strip() for line in output.splitlines() if line.strip()]


def sample_id_from_wisp(wisp_uri: str) -> str:
    filename = wisp_uri.rsplit("/", 1)[-1]
    name = re.sub(r"\.wisp\.summary\.tsv$", "", filename)
    if "_" in name:
        return name.split("_", 1)[1]
    return name


def load_master(master_file: Path) -> dict[str, str]:
    delimiter = "\t" if master_file.suffix.lower() in (".tsv", ".txt") else ","
    with open(master_file, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if reader.fieldnames is None:
            sys.exit(f"ERROR: {master_file} appears to be empty.")
        missing = {MASTER_SAMPLE_ID_COL, MASTER_DATE_COL} - set(reader.fieldnames)
        if missing:
            sys.exit(f"ERROR: master file missing columns: {', '.join(sorted(missing))}")
        return {
            re.sub(r"T(?:II)?(P\d)", r"\1", row[MASTER_SAMPLE_ID_COL].strip()): datetime.strptime(row[MASTER_DATE_COL].strip(), "%d-%m-%Y %H:%M:%S").strftime("%Y-%m-%d")
            for row in reader
            if row[MASTER_SAMPLE_ID_COL].strip()
        }


if __name__ == "__main__":
    main(parse_args())