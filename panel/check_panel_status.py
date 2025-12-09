#!/usr/bin/env python3
import argparse
import logging
import subprocess

from pathlib import Path
from typing import List, Optional

TARGETED_PIPELINE_DIRECTORY = "gs://targeted-pipeline-output-prod-1"


def main(batch: Optional[str], input_directory: Optional[Path]) -> None:
    run_buckets = run_bash_command(["gcloud", "storage",  "ls", TARGETED_PIPELINE_DIRECTORY]).split("\n")

    if batch is None:
        non_empty_run_buckets = [bucket for bucket in run_buckets if bucket]
        most_recent_run = non_empty_run_buckets[-1]
        logging.info(f"Most recent run found: {most_recent_run}")
        batch = most_recent_run.strip().split("/")[3].split("_")[0]
        logging.info(f"Inferred batch: {batch}")

    if input_directory is None:
        input_directory = Path.home() / f"tmp_diag_{batch}"
        logging.info(f"Inferred input directory: {input_directory}")

    relevant_buckets = [bucket for bucket in run_buckets if bucket.startswith(f"{TARGETED_PIPELINE_DIRECTORY}/{batch}")]
    relevant_count_found = len(relevant_buckets)
    logging.info(f"Found {relevant_count_found} relevant runs")

    relevant_orange_pdfs = [
        bucket for bucket in run_bash_command(["gcloud", "storage", "ls", f"{TARGETED_PIPELINE_DIRECTORY}/{batch}*/orange/*.pdf"]).split("\n")
        if bucket
    ]
    log_or_warn(f"Found {len(relevant_orange_pdfs)} relevant Orange PDFs", len(relevant_orange_pdfs) == relevant_count_found)

    if input_directory is None or not input_directory.exists():
        logging.info(
            f"If all expected runs are done, provide a QC check directory with the option '-i' or create one with something like: "
            f"qc_check_oncopanel_batch {batch} {Path.home()}/tmp_diag_{batch}"
        )
    else:
        remarks_yaml = input_directory / f"{batch}_remarks.yaml"
        if remarks_yaml.exists():
            logging.info(f"Found remarks YAML {remarks_yaml}")

            remarks_path = f"gs://oncoact-panel-remarks-output/remarks-{batch}/remarks/remarks.tsv"
            if run_bash_command(["gcloud", "storage",  "ls", remarks_path]):
                remarks_output = run_bash_command(["gcloud", "storage",  "cat", remarks_path])
                logging.info(f"Found remarks output:\n{remarks_output}")
            else:
                logging.warning(f"No remarks output for batch {batch}")
        else:
            logging.warning(f"No remarks YAML {remarks_yaml}")


def run_bash_command(command: List[str]) -> str:
    return subprocess.run(command, capture_output=True, text=True).stdout


def log_or_warn(message: str, should_log: bool) -> None:
    if should_log:
        logging.info(message)
    else:
        logging.warning(message)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="check_panel_status.py", description="Check status for batch of OncoPanel runs")
    parser.add_argument("--batch", "-b", type=str, required=False, help="Batch for runs")
    parser.add_argument("--input-directory", "-i", type=Path, required=False, help="Local working directory")
    return parser.parse_args()


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s - [%(levelname)-8s] - %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
    )
    args = parse_args()
    main(args.batch, args.input_directory)
