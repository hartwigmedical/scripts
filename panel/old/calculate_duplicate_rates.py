#!/usr/bin/env python
import argparse
import logging
import re
import sys
from dataclasses import dataclass
from typing import List, Optional

from google.cloud import storage

from gcp.base import GCPPath
from gcp.client import GCPClient

MAPPED_LINE_REGEX = re.compile("^\d+ \+ 0 mapped \(\d+\.\d+%:N/A\)$")

@dataclass(frozen=True)
class FlagstatSummary(object):
    total_alignment_count: int
    duplicate_read_count: int
    mapped_alignment_count: int
    paired_in_sequencing_read_count: int

    @classmethod
    def from_text(cls, text: str) -> "FlagstatSummary":
        split_text = text.split("\n")
        total_alignment_count = cls._get_total_alignment_count(split_text[0])
        duplicate_read_count = cls._get_duplicate_read_count(split_text[3])
        mapped_alignment_count = cls._get_mapped_alignment_count(split_text[4])
        paired_in_sequencing_read_count = cls._get_paired_in_sequencing_read_count(split_text[5])
        summary = FlagstatSummary(
            total_alignment_count, duplicate_read_count, mapped_alignment_count, paired_in_sequencing_read_count,
        )
        return summary

    @classmethod
    def _get_total_alignment_count(cls, line: str) -> int:
        if not line.endswith(" + 0 in total (QC-passed reads + QC-failed reads)"):
            raise SyntaxError(f"Unexpected format for total alignment count line: '{line}'")
        return int(line.split(" ")[0])

    @classmethod
    def _get_duplicate_read_count(cls, line: str) -> int:
        if not line.endswith(" + 0 duplicates"):
            raise SyntaxError(f"Unexpected format for duplicate read count line: '{line}'")
        return int(line.split(" ")[0])

    @classmethod
    def _get_mapped_alignment_count(cls, line: str) -> int:
        if not MAPPED_LINE_REGEX.fullmatch(line):
            raise SyntaxError(f"Unexpected format for mapped alignment line: '{line}'")
        return int(line.split(" ")[0])

    @classmethod
    def _get_paired_in_sequencing_read_count(cls, line: str) -> int:
        if not line.endswith(" + 0 paired in sequencing"):
            raise SyntaxError(f"Unexpected format for paired-in-sequencing read count line: '{line}'")
        return int(line.split(" ")[0])


def main(non_umi_flagstat_path: GCPPath, umi_flagstat_path: Optional[GCPPath]) -> None:
    logging.debug(f"Start calculating duplicate rates for {non_umi_flagstat_path} and {umi_flagstat_path}")
    gcp_client = GCPClient(storage.Client())
    non_umi_flagstat_summary = get_flagstat_summary(non_umi_flagstat_path, gcp_client)
    if umi_flagstat_path is not None:
        umi_flagstat_summary = get_flagstat_summary(umi_flagstat_path, gcp_client)
    else:
        umi_flagstat_summary = None

    # See DEV-2565 for details
    non_umi_duplicate_read_count = non_umi_flagstat_summary.duplicate_read_count
    total_read_count = non_umi_flagstat_summary.paired_in_sequencing_read_count
    total_alignment_count = non_umi_flagstat_summary.total_alignment_count
    mapped_alignment_count = non_umi_flagstat_summary.mapped_alignment_count

    # Each unmapped read has exactly one alignment
    unmapped_read_count = total_alignment_count - mapped_alignment_count

    non_umi_duplicate_rate = non_umi_duplicate_read_count / total_read_count
    unmapped_rate = unmapped_read_count / total_read_count

    if umi_flagstat_summary is not None:
        umi_read_count = umi_flagstat_summary.paired_in_sequencing_read_count
        # All reads are either unmapped, duplicate or in the UMI-deduplicated BAM
        umi_duplicate_rate = 1 - (umi_read_count + unmapped_read_count) / total_read_count
    else:
        umi_duplicate_rate = None

    if umi_duplicate_rate is not None:
        logging.debug(f"Non-UMI duplicate rate: {non_umi_duplicate_rate:.4f}")
        logging.debug(f"UMI duplicate rate: {umi_duplicate_rate:.4f}")
        logging.debug(f"Unmapped rate: {unmapped_rate:.4f}")
        print(f"{non_umi_flagstat_path}\t{non_umi_duplicate_rate}\t{umi_duplicate_rate}\t{unmapped_rate}")
    else:
        logging.debug(f"Non-UMI duplicate rate: {non_umi_duplicate_rate:.4f}")
        logging.debug(f"Unmapped rate: {unmapped_rate:.4f}")
        print(f"{non_umi_flagstat_path}\t{non_umi_duplicate_rate}\t{unmapped_rate}")
    logging.debug("Finished calculating_duplicate_rates")


def get_flagstat_summary(gcp_flagstat_path: GCPPath, gcp_client: GCPClient) -> FlagstatSummary:
    flagstat_text = gcp_client.get_text(gcp_flagstat_path)
    return FlagstatSummary.from_text(flagstat_text)


def parse_args(sys_args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate non-UMI duplicate rate, UMI duplicate rate, and unmapped rate from flagstat output.\n"
            "The argument and calculation for the UMI-deduplicated BAM are optional.\n"
            "\n"
            "The duplicate rates are calculated as 'count of duplicate reads'/'total read count'.\n"
            "The 'unmapped rate' is calculated as 'count of unmapped reads'/'total read count'.\n"
            "\n"
            "Assumes that the non-UMI deduplicated BAM has been deduplicated by sambamba markdup, \n"
            "and the UMI deduplicated BAM by UMICollapse with dropping all duplicate and unmapped reads.\n"
            "\n"
            "When the flagstat file for a UMI-deduplicated BAM is provided, the output order is as follows:\n"
            "'${Non-UMI flagstat path}\t${Non-UMI duplicate rate}\t${UMI duplicate rate}\t${Unmapped rate}'\n"
            "When no flagstat file for a UMI-deduplicated BAM is provided, the output order is:\n"
            "'${Non-UMI flagstat path}\t${Non-UMI duplicate rate}\t${Unmapped rate}'"
        )
    )
    parser.add_argument(
        '--non_umi_flagstat', '-n', type=GCPPath.from_string, required=True, help="GCP path to flagstat for non-UMI-deduplicated BAM.",
    )
    parser.add_argument(
        '--umi_flagstat', '-u', type=GCPPath.from_string, help="Optional GCP path to flagstat for UMI-deduplicated BAM.",
    )
    return parser.parse_args(sys_args)


if __name__ == '__main__':
    logging.basicConfig(
        format="%(asctime)s - [%(levelname)-8s] - %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
    )
    args = parse_args(sys.argv[1:])
    main(args.non_umi_flagstat, args.umi_flagstat)
