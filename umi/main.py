#!/usr/bin/env python
import argparse
import logging
import sys
from collections import defaultdict

from typing import List, NamedTuple, Set, Dict, Optional, DefaultDict, Tuple

import pysam as pysam

# set up logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.DEBUG,
    datefmt='%Y-%m-%d %H:%M:%S')


class SingleReadSummary(NamedTuple):
    query_name: str
    reference_name: str  # contig or chromosome
    reference_start: int  # position
    explicit_umi: Optional[str]  # from BX tag, if present
    mapping_quality: int
    is_duplicate: bool
    is_reverse: bool
    is_unmapped: bool


class ReadPosition(NamedTuple):
    reference_name: str  # contig or chromosome
    reference_start: int  # position
    is_reverse: bool
    is_unmapped: bool


class PairedReadSummary(NamedTuple):
    query_name: str
    umi: Optional[str]
    worst_mapping_quality: int
    is_duplicate: bool


class PairedReadTracker(object):
    def __init__(self) -> None:
        self.__read_positions_to_paired_read_summaries = defaultdict(list)

    def add(self, reads: Set[SingleReadSummary], warnings: DefaultDict[Tuple[str, ...], List[str]]) -> None:
        # Order of arguments does not matter
        # umis = {read.explicit_umi for read in reads if read.explicit_umi is not None}
        # if len(umis) == 0:
        #     umis = {read.query_name.split(":")[-1] for read in reads}

        umis = {read.query_name.split(":")[-1] for read in reads}
        if len(umis) == 1:
            umi = umis.pop()
        elif len(umis) == 0:
            raise ValueError(f"No umi could be extracted {reads}")
        else:
            raise ValueError(f"Reads disagree on umi: umis={umis}, reads={reads}")

        worst_mapping_quality = min(read.mapping_quality for read in reads)

        read_summaries = {
            PairedReadSummary(read.query_name, umi, worst_mapping_quality, read.is_duplicate) for read in reads
        }
        if len(read_summaries) == 1:
            read_summary = read_summaries.pop()
        else:
            read_is_duplicate = any(read.is_duplicate for read in reads)
            read_summaries = {
                PairedReadSummary(read.query_name, umi, worst_mapping_quality, read_is_duplicate) for read in reads
            }
            assert len(read_summaries) == 1, f"Read summaries are still not the same: {reads}"
            warnings[tuple(sorted(list({read.query_name for read in reads})))].append(f"Not all reads agree on whether they are duplicate: {reads}")
            read_summary = read_summaries.pop()

        read_positions = frozenset({
            ReadPosition(read.reference_name, read.reference_start, read.is_reverse, read.is_unmapped) for read in reads
        })

        self.__read_positions_to_paired_read_summaries[read_positions].append(read_summary)

    def get_reads_to_be_undupped(self, warnings: DefaultDict[Tuple[str, ...], List[str]]) -> Set[str]:
        reads_to_be_undupped = set()
        for read_summaries in self.__read_positions_to_paired_read_summaries.values():
            umi_to_summaries = defaultdict(list)
            for summary in read_summaries:
                umi_to_summaries[summary.umi].append(summary)

            for umi, summaries_grouped_by_umi in umi_to_summaries.items():
                non_duplicate_read_pair_summaries = [summary for summary in summaries_grouped_by_umi if not summary.is_duplicate]
                if not non_duplicate_read_pair_summaries:
                    # a read should be unduplicated
                    best_quality_pair = sorted(summaries_grouped_by_umi, key=lambda x: x.worst_mapping_quality, reverse=True)[0]
                    reads_to_be_undupped.add(best_quality_pair.query_name)
                else:
                    if len(non_duplicate_read_pair_summaries) != 1:
                        warning_msg = (
                            f"Multiple non duplicate reads for same umi at same locations.\n"
                            f"{non_duplicate_read_pair_summaries}"
                        )
                        warnings[tuple(sorted(list({summary.query_name for summary in non_duplicate_read_pair_summaries})))].append(warning_msg)

        return reads_to_be_undupped


def main(input_bam_path: str, output_bam_path: str) -> None:
    print("START UNDUP")
    warnings = defaultdict(list)

    bam_file = pysam.AlignmentFile(input_bam_path, "rb")

    query_name_to_single_read_summaries: DefaultDict[str] = defaultdict(list)

    tracker = PairedReadTracker()

    for read in bam_file:
        # print(read)
        # tags = dict(read.get_tags())
        # umi = tags['BX'] if 'BX' in tags.keys() else None

        # print(
        #     f"("
        #     f"read_id={read.query_name}, "
        #     f"flag={read.flag}, "
        #     f"seq_name={read.reference_name}, "
        #     f"position={read.reference_start}, "
        #     f"mapping_quality={read.mapping_quality}, "
        #     f"next_read_seq_name={read.next_reference_name}, "
        #     f"next_read_position={read.next_reference_start}, "
        #     f"tags={tags}, "
        #     f"duplicate={read.is_duplicate}, "
        #     f"paired={read.is_paired}, "
        #     f"reverse={read.is_reverse}, "
        #     f"umi={umi}, "
        #     f")"
        # )

        query_name_to_single_read_summaries[read.query_name].append(SingleReadSummary(
            read.query_name,
            read.reference_name,
            read.reference_start,
            read.get_tag("BX") if read.has_tag("BX") else None,
            read.mapping_quality,
            read.is_duplicate,
            read.is_reverse,
            read.is_unmapped,
        ))

    for grouped_read_summaries in query_name_to_single_read_summaries.values():
        if not reads_match(grouped_read_summaries):
            grouped_reads_string = "\n".join(str(summary) for summary in grouped_read_summaries)
            warning_msg = f"Reads don't match perfectly:\n{grouped_reads_string}"
            warnings[(grouped_read_summaries[0].query_name,)].append(warning_msg)

        tracker.add(grouped_read_summaries, warnings)

    reads_to_be_undupped = tracker.get_reads_to_be_undupped(warnings)
    # print(reads_to_be_undupped)

    bam_file = pysam.AlignmentFile(input_bam_path, "rb")
    with pysam.AlignmentFile(output_bam_path, "wb", header=bam_file.header) as outf:
        for read in bam_file:
            if read.query_name in reads_to_be_undupped:
                if read.flag % 2048 < 1024:
                    warning_msg = f"Impossible bit flag for read to be undupped: {read.flag}\n{read}"
                    warnings[(read.query_name,)].append(warning_msg)
                else:
                    read.flag = read.flag - 1024
            outf.write(read)

    warning_keys = sorted(list(warnings.keys()), key=lambda x: (len(x), x))
    print("\n\n".join([f"{key}:\n" + "\n".join(warnings[key]) for key in warning_keys]))

    if warnings:
        raise ValueError("There was at least one warning")

    print("DONE WITH UNDUP")


def reads_match(reads: Set[SingleReadSummary]) -> bool:
    return len({(read.query_name, read.is_duplicate) for read in reads}) == 1


def parse_args(sys_args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=("Check whether duplicate reads from sambamba markdup are duplicates according to umi-tools group, "
                     "and then unmark them as duplicates if they aren't.")
    )
    parser.add_argument('input_bam', type=str, help='Input bam. Needs to have gone through umi-tools group and sambamba markdup already.')
    parser.add_argument('output_bam', type=str, help='The output bam. If file exists, it is deleted first.')
    return parser.parse_args(sys_args)


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    main(args.input_bam, args.output_bam)
