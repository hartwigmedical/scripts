#!/usr/bin/env python
import argparse
import logging
import sys
from collections import defaultdict

from typing import List, NamedTuple, Set, Dict, Optional

import pysam as pysam

# set up logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.DEBUG,
    datefmt='%Y-%m-%d %H:%M:%S')


class ReadPosition(NamedTuple):
    reference_name: str  # contig, chromosome
    reference_start: int  # position
    mapped_to_reverse: bool
    is_unmapped: bool


class PairedReadSummary(NamedTuple):
    query_name: str
    umi: Optional[str]
    worst_mapping_quality: int
    is_duplicate: bool


class PairedReadTracker(object):
    def __init__(self) -> None:
        self.__read_positions_to_paired_read_summaries = defaultdict(list)

    def add(self, left_read: pysam.AlignedSegment, right_read: pysam.AlignedSegment) -> None:
        # Order of arguments does not matter
        reads = {left_read, right_read}

        umis = {read.get_tag("BX") for read in reads if read.has_tag("BX")}
        if len(umis) == 1:
            umi = umis.pop()
        elif len(umis) == 0:
            umi = left_read.query_name.split(":")[-1]
        else:
            raise ValueError(f"Pair has different umis: {umis}")

        worst_mapping_quality = min(read.mapping_quality for read in reads)

        read_summaries = {
            PairedReadSummary(read.query_name, umi, worst_mapping_quality, read.is_duplicate) for read in reads
        }
        assert len(read_summaries) == 1, f"Read summaries are not the same"
        read_summary = read_summaries.pop()

        read_positions = frozenset({
            ReadPosition(read.reference_name, read.reference_start, read.is_reverse, read.is_unmapped) for read in reads
        })
        assert len(read_positions) == 2, f"Read Positions are the same"

        self.__read_positions_to_paired_read_summaries[read_positions].append(read_summary)

    def get_reads_to_be_undupped(self) -> Set[str]:
        reads_to_be_undupped = set()
        for read_summaries in self.__read_positions_to_paired_read_summaries.values():
            umi_to_summaries = defaultdict(list)
            for summary in read_summaries:
                umi_to_summaries[summary.umi].append(summary)

            for umi, summaries_grouped_by_umi in umi_to_summaries.items():
                non_duplicate_read_pair_summaries = [summary for summary in summaries_grouped_by_umi if not summary.is_duplicate]
                if not non_duplicate_read_pair_summaries:
                    # a read should be unduplicated
                    best_quality_pair = sorted(summaries_grouped_by_umi, key=lambda x: x.mapping_quality, reverse=True)[0]
                    reads_to_be_undupped.add(best_quality_pair.query_name)
                else:
                    assert len(non_duplicate_read_pair_summaries) == 1, (
                        f"Multiple non duplicate reads for same umi at same locations.\n"
                        f"{non_duplicate_read_pair_summaries}"
                    )

        return reads_to_be_undupped


def main(input_bam_path: str, output_bam_path: str) -> None:
    # TODO: remove counter limit
    bam_file = pysam.AlignmentFile(input_bam_path, "rb")

    query_name_to_single_read = {}
    query_names_of_found_pairs = set()

    tracker = PairedReadTracker()

    counter = 0
    for read in bam_file:
        print(read)
        tags = dict(read.get_tags())
        umi = tags['BX'] if 'BX' in tags.keys() else None

        print(
            f"("
            f"read_id={read.query_name}, "
            f"flag={read.flag}, "
            f"seq_name={read.reference_name}, "
            f"position={read.reference_start}, "
            f"mapping_quality={read.mapping_quality}, "
            f"next_read_seq_name={read.next_reference_name}, "
            f"next_read_position={read.next_reference_start}, "
            f"tags={tags}, "
            f"duplicate={read.is_duplicate}, "
            f"paired={read.is_paired}, "
            f"reverse={read.is_reverse}, "
            f"umi={umi}, "
            f")"
        )

        assert read.query_name not in query_names_of_found_pairs, "Read pair already seen!"

        if read.query_name not in query_name_to_single_read.keys():
            query_name_to_single_read[read.query_name] = read
        else:
            other_half = query_name_to_single_read.pop(read.query_name)
            assert reads_match(read, other_half), f"Reads don't match perfectly:\n{read}\n{other_half}"
            query_names_of_found_pairs.add(read.query_name)

            tracker.add(read, other_half)

        counter += 1
        if counter > 1000:
            break

    print("hello")
    print(query_name_to_single_read)
    print(query_names_of_found_pairs)
    reads_to_be_undupped = tracker.get_reads_to_be_undupped()
    print(reads_to_be_undupped)

    bam_file = pysam.AlignmentFile(input_bam_path, "rb")
    with pysam.AlignmentFile(output_bam_path, "wb", header=bam_file.header) as outf:
        for read in bam_file:
            if read.query_name in reads_to_be_undupped:
                if read.flag % 2048 < 1024:
                    error_msg = f"Impossible bit flag for read to be undupped: {read.flag}\n{read}"
                    raise ValueError(error_msg)
                read.flag = read.flag - 1024
            outf.write(read)

    print("DONE")


def reads_match(left: pysam.AlignedSegment, right: pysam.AlignedSegment) -> bool:
    return (
            left.query_name == right.query_name and
            left.reference_name == right.next_reference_name and
            left.reference_start == right.next_reference_start and
            left.next_reference_name == right.reference_name and
            left.next_reference_start == right.reference_start and
            left.is_duplicate == right.is_duplicate
    )


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
