#!/usr/bin/env python
import argparse
import logging
import sys

from typing import List

# set up logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.DEBUG,
    datefmt='%Y-%m-%d %H:%M:%S')


def main(input_bam: str, output_bam: str) -> None:
    print("hello")


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
