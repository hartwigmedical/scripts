#!/usr/bin/env python3
import gzip
import logging
import sys

from copy import copy
from enum import Enum, auto
from functools import total_ordering
from pathlib import Path
from typing import Dict, Set, Tuple

# Meant for SAGE 3.3

HEADER_LINES = [
    "##fileformat=VCFv4.2",
    '##FILTER=<ID=PASS,Description="All filters passed">',
    '##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description="Average calculated base quality">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
    '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allelic frequency calculated from read context counts as (Full + Partial + Core + Realigned + Alt) / Coverage">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=RABQ,Number=2,Type=Integer,Description="Raw allelic base quality">',
    '##FORMAT=<ID=RAD,Number=2,Type=Integer,Description="Raw allelic depth">',
    '##FORMAT=<ID=RC_CNT,Number=7,Type=Integer,Description="Read context counts [Full, Partial, Core, Realigned, Alt, Reference, Total]">',
    '##FORMAT=<ID=RC_IPC,Number=1,Type=Integer,Description="Read context improper pair count">',
    '##FORMAT=<ID=RC_JIT,Number=3,Type=Integer,Description="Read context jitter [Shortened, Lengthened, QualityPenalty]">',
    '##FORMAT=<ID=RC_QUAL,Number=7,Type=Integer,Description="Read context quality [Full, Partial, Core, Realigned, Alt, Reference, Total]">',
    '##FORMAT=<ID=RDP,Number=1,Type=Integer,Description="Raw read depth">',
    '##FORMAT=<ID=SB,Number=1,Type=Float,Description="Strand bias - percentage of first-in-pair reads">',
    '##FORMAT=<ID=UMI_CNT,Number=6,Type=Integer,Description="UMI type counts [RefNone,RefSingle,RefDualStrand,AltNone,AltSingle,AltDualStrand]">',
    "##sageVersion=3.3",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFAKE",
]
FORMAT = "GT:ABQ:AD:AF:DP:RABQ:RAD:RC_CNT:RC_IPC:RC_JIT:RC_QUAL:RDP:SB"


class ReferenceGenomeVersion(Enum):
    V37 = auto()
    V38 = auto()


@total_ordering
class Chromosome(Enum):
    CHR1 = auto()
    CHR2 = auto()
    CHR3 = auto()
    CHR4 = auto()
    CHR5 = auto()
    CHR6 = auto()
    CHR7 = auto()
    CHR8 = auto()
    CHR9 = auto()
    CHR10 = auto()
    CHR11 = auto()
    CHR12 = auto()
    CHR13 = auto()
    CHR14 = auto()
    CHR15 = auto()
    CHR16 = auto()
    CHR17 = auto()
    CHR18 = auto()
    CHR19 = auto()
    CHR20 = auto()
    CHR21 = auto()
    CHR22 = auto()
    CHRX = auto()
    CHRY = auto()
    CHRMT = auto()

    def __str__(self) -> str:
        return self.name[3:]

    def __repr__(self) -> str:
        return f"Chromosome.{self.name}"

    def __lt__(self, other):
        if self.__class__ is other.__class__:
            return self.value < other.value
        return NotImplemented

    @classmethod
    def from_string(cls, chromosome: str) -> "Chromosome":
        if chromosome[:3] != "chr":
            chromosome = f"chr{chromosome}"

        if chromosome.upper() == "CHRM":
            return Chromosome.CHRMT
        else:
            return Chromosome[chromosome.upper()]

    def to_string(self, reference_genome_version: ReferenceGenomeVersion) -> str:
        if reference_genome_version == ReferenceGenomeVersion.V37:
            return self.name[3:]
        elif reference_genome_version == ReferenceGenomeVersion.V38:
            return f"chr{self.name[3:]}"
        else:
            raise NotImplementedError(f"Unrecognized reference genome version: {reference_genome_version}")


def main() -> None:
    logging.info("Start creating SAGE append VCF")

    logging.info("Parse CLI arguments")
    input_bed, output_vcf, reference_genome_version = parse_args()

    logging.info("Perform sanity checks")
    if not input_bed.exists():
        raise ValueError(f"Input BED file does not exist: {input_bed}")
    if input_bed.name[-4:] != ".bed":
        logging.warning("Name of input BED file does not end in '.bed' file extension.")
    if output_vcf.name[-7:] != ".vcf.gz":
        raise ValueError("Name of output VCF file does not end in '.vcf.gz' file extension.")

    logging.info("Read input BED file")
    chromosome_to_positions = read_input_bed(input_bed)

    total_position_count = sum(len(positions) for positions in chromosome_to_positions.values())
    logging.info(f"Making VCF with {total_position_count} positions")

    logging.info("Write output file")
    write_output_vcf(chromosome_to_positions, output_vcf, reference_genome_version)

    logging.info("Finished creating SAGE append VCF")


def parse_args() -> Tuple[Path, Path, ReferenceGenomeVersion]:
    if len(sys.argv) != 4:
        error_msg = "\n".join(
            [
                "Incorrect number of arguments provided. Should be:",
                "   create_sage_append_vcf_from_bed.py ${input_bed} ${output_vcf} ${reference_genome_version}",
                "e.g.:",
                "   create_sage_append_vcf_from_bed.py input.bed output.vcf.gz 37",
                "   create_sage_append_vcf_from_bed.py path/input.bed other_path/output.vcf.gz 38",
            ]
        )
        raise ValueError(error_msg)
    input_bed_arg = sys.argv[1]
    output_vcf_arg = sys.argv[2]
    reference_genome_version_arg = sys.argv[3]

    input_bed = Path(input_bed_arg)
    output_vcf = Path(output_vcf_arg)

    if reference_genome_version_arg == "37":
        reference_genome_version = ReferenceGenomeVersion.V37
    elif reference_genome_version_arg == "38":
        reference_genome_version = ReferenceGenomeVersion.V38
    else:
        raise ValueError(f"Invalid reference genome version argument: {reference_genome_version_arg}")

    return input_bed, output_vcf, reference_genome_version


def read_input_bed(input_bed: Path) -> Dict[Chromosome, Set[int]]:
    chromosome_to_positions: Dict[Chromosome, Set[int]] = {}
    with open(input_bed) as in_f:
        for line in in_f.readlines():
            split_line = line.split("\t")
            chromosome = Chromosome.from_string(split_line[0])
            start_position = int(split_line[1])
            end_position = int(split_line[2])

            if chromosome not in chromosome_to_positions.keys():
                chromosome_to_positions[chromosome] = set()

            # BED file are 0-indexed and have inclusive start and exclusive end
            for position in range(start_position + 1, end_position + 1):
                chromosome_to_positions[chromosome].add(position)
    return chromosome_to_positions


def write_output_vcf(
        chromosome_to_positions: Dict[Chromosome, Set[int]],
        output_vcf: Path,
        reference_genome_version: ReferenceGenomeVersion
) -> None:
    text = get_output_vcf_text(chromosome_to_positions, reference_genome_version)
    with gzip.open(output_vcf, "wt") as out_f:
        out_f.write(text)


def get_output_vcf_text(
        chromosome_to_positions: Dict[Chromosome, Set[int]], reference_genome_version: ReferenceGenomeVersion
) -> str:
    lines = copy(HEADER_LINES)
    for chromosome in sorted(chromosome_to_positions.keys()):
        for position in sorted(chromosome_to_positions[chromosome]):
            line = get_output_line(chromosome, position, reference_genome_version)
            lines.append(line)
    text = "\n".join(lines)
    return text


def get_output_line(chromosome: Chromosome, position: int, reference_genome_version: ReferenceGenomeVersion) -> str:
    ref_base = "A"
    alt_base = "T"
    line = (
        f"{chromosome.to_string(reference_genome_version)}\t{position}\t.\t{ref_base}\t{alt_base}\t.\tPASS\t.\t{FORMAT}"
        f"\t./.:0:0,0:0.0:0:0,0:0,0:0,0,0,0,0,0,0:0:0,0,0:0,0,0,0,0,0,0:0:0.0"
    )
    return line


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s - [%(levelname)-8s] - %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
    )
    main()
