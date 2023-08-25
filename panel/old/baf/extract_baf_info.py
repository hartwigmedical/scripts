#!/usr/bin/env python
import argparse
import gzip
import sys
from typing import List, Dict, Tuple


def main(amber_output: str, loci: str, output_file: str) -> None:
    position_to_allele_counts: Dict[Tuple[str, str], Tuple[int, int]] = {}
    with gzip.open(amber_output, "rt") as amber_output_f:
        for line in amber_output_f:
            if line[:3] == "chr":
                chrom, pos, _, _, _, _, _, _, _, data = line.split("\t")
                allele_counts = data.split(":")[1]
                ref_count, alt_count = allele_counts.split(",")
                ref_count = int(ref_count)
                alt_count = int(alt_count)
                position_to_allele_counts[(chrom, pos)] = (ref_count, alt_count)

    with open(loci, "r") as loci_f:
        with open(output_file, "w") as output_f:
            for line in loci_f:
                if line[:3] == "chr":
                    chrom, pos, _, _, _, _, _, _ = line.split("\t")
                    if (chrom, pos) in position_to_allele_counts.keys():
                        ref_count, alt_count = position_to_allele_counts[(chrom, pos)]
                    else:
                        ref_count = 0
                        alt_count = 0
                    new_line = f"{ref_count}\t{alt_count}\n"
                    output_f.write(new_line)


def parse_args(sys_args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract BAF info from AMBER output")
    parser.add_argument('amber_output', type=str, help='AMBER output vcf.')
    parser.add_argument('loci', type=str, help='AMBER config file with loci.')
    parser.add_argument('output_file', type=str, help='Output file.')
    return parser.parse_args(sys_args)


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    main(args.amber_output, args.loci, args.output_file)


