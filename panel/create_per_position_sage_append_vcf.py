#!/usr/bin/env python3
import argparse
import gzip
import logging
import sys

from copy import copy
from dataclasses import dataclass
from enum import Enum, auto
from functools import total_ordering
from pathlib import Path
from typing import Dict, Set, List

# Meant for SAGE 3.4

HEADER_LINES = [
    "##fileformat=VCFv4.2",
    '##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description="Average calculated base quality">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
    '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allelic frequency calculated from read context counts as (Full + Partial + Core + Realigned + Alt) / Coverage">',
    '##FORMAT=<ID=AMQ,Number=2,Type=Integer,Description="Average map quality count (all,alt)">',
    '##FORMAT=<ID=ANM,Number=2,Type=Float,Description="Average NM count (all,alt)">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=RABQ,Number=2,Type=Integer,Description="Raw allelic base quality">',
    '##FORMAT=<ID=RAD,Number=2,Type=Integer,Description="Raw allelic depth">',
    '##FORMAT=<ID=RC_CNT,Number=7,Type=Integer,Description="Read context counts [Full, Partial, Core, Realigned, Alt, Reference, Total]">',
    '##FORMAT=<ID=RC_IPC,Number=1,Type=Integer,Description="Read context improper pair count">',
    '##FORMAT=<ID=RC_JIT,Number=3,Type=Integer,Description="Read context jitter [Shortened, Lengthened, QualityPenalty]">',
    '##FORMAT=<ID=RC_QUAL,Number=7,Type=Integer,Description="Read context quality [Full, Partial, Core, Realigned, Alt, Reference, Total]">',
    '##FORMAT=<ID=RDP,Number=1,Type=Integer,Description="Raw read depth">',
    '##FORMAT=<ID=RSB,Number=1,Type=String,Description="Read strand bias - percentage of forward-orientation reads">',
    '##FORMAT=<ID=SB,Number=1,Type=String,Description="Fragment strand bias - percentage of forward-orientation fragments">',
    '##FORMAT=<ID=UMI_CNT,Number=6,Type=Integer,Description="UMI type counts [TotalNone,TotalSingle,TotalDualStrand,AltNone,AltSingle,AltDualStrand]">',
    "##sageVersion=3.4",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFAKE",
]
FORMAT = "GT:ABQ:AD:AF:AMQ:ANM:DP:RABQ:RAD:RC_CNT:RC_IPC:RC_JIT:RC_QUAL:RDP:RSB:SB"


class ReferenceGenomeVersion(Enum):
    V37 = auto()
    V38 = auto()

    @classmethod
    def from_string(cls, string_value: str) -> "ReferenceGenomeVersion":
        if string_value == "37":
            return ReferenceGenomeVersion.V37
        elif string_value == "38":
            return ReferenceGenomeVersion.V38
        else:
            raise ValueError(f"Invalid reference genome version argument: {string_value}")


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


@dataclass(frozen=True)
class Config:
    input: Path
    output: Path
    reference_genome_version: ReferenceGenomeVersion


def main() -> None:
    logging.info("Start creating SAGE append VCF")

    logging.info("Parse CLI arguments")
    config = parse_args(sys.argv[1:])

    logging.info("Perform sanity checks")
    if not config.input.exists():
        raise ValueError(f"Input file does not exist: {config.input}")
    if not ends_with(config.input.name, ".bed") and not ends_with(config.input.name, ".vcf.gz"):
        raise ValueError("Name of input file does not end in '.bed' or '.vcf.gz' file extension.")
    if not ends_with(config.output.name, ".vcf.gz"):
        raise ValueError("Name of output VCF file does not end in '.vcf.gz' file extension.")

    logging.info("Read input file")
    if ends_with(config.input.name, ".bed"):
        chromosome_to_positions = read_input_bed(config.input)
    else:
        chromosome_to_positions = read_input_vcf(config.input)

    total_position_count = sum(len(positions) for positions in chromosome_to_positions.values())
    logging.info(f"Making VCF with {total_position_count} positions")

    logging.info("Write output file")
    write_output_vcf(chromosome_to_positions, config.output, config.reference_genome_version)

    logging.info("Finished creating SAGE append VCF")


def parse_args(arguments: List[str]) -> Config:
    parser = argparse.ArgumentParser(description="Create VCF with SNV per position for SAGE append")
    parser.add_argument("--input", "-i", type=str, help="BED file or VCF file of relevant positions")
    parser.add_argument("--output", "-o", type=str, help="Path for output VCF.")
    parser.add_argument("--reference_genome_version", "-r", type=str, help="Reference genome version")
    args = parser.parse_args(arguments)
    return Config(Path(args.input), Path(args.output), ReferenceGenomeVersion.from_string(args.reference_genome_version))


def read_input_bed(input_bed: Path) -> Dict[Chromosome, Set[int]]:
    chromosome_to_positions: Dict[Chromosome, Set[int]] = {}
    with open(input_bed) as in_f:
        for line in in_f.readlines():
            split_line = line.replace("\n", "").split("\t")
            chromosome = Chromosome.from_string(split_line[0])
            start_position = int(split_line[1])
            end_position = int(split_line[2])

            if chromosome not in chromosome_to_positions.keys():
                chromosome_to_positions[chromosome] = set()

            # BED file are 0-indexed and have inclusive start and exclusive end
            for position in range(start_position + 1, end_position + 1):
                chromosome_to_positions[chromosome].add(position)
    return chromosome_to_positions


def read_input_vcf(input_vcf: Path) -> Dict[Chromosome, Set[int]]:
    chromosome_to_positions: Dict[Chromosome, Set[int]] = {}
    with gzip.open(input_vcf, "rt") as in_f:
        for line in in_f.readlines():
            if line[0] != "#":
                split_line = line.replace("\n", "").split("\t")
                chromosome = Chromosome.from_string(split_line[0])
                start_position = int(split_line[1])
                ref = split_line[3]

                if chromosome not in chromosome_to_positions.keys():
                    chromosome_to_positions[chromosome] = set()

                for position in range(start_position, start_position + len(ref)):
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


def ends_with(input_string: str, end_string: str) -> bool:
    return input_string[-len(end_string):] == end_string


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s - [%(levelname)-8s] - %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
    )
    main()
