import argparse
import sys
from decimal import Decimal, ROUND_DOWN
from typing import List, NamedTuple


BASES_PER_READ = 150
BASES_PER_GBASE = 10**9


class Config(NamedTuple):
    read_count: int
    excluded_base_percentage: Decimal


def main(config: Config) -> None:
    if config.read_count < 0:
        raise ValueError(f"Read count cannot be negative: {config.read_count}")

    if config.excluded_base_percentage <= Decimal(0):
        raise ValueError(f"Excluded base percentage cannot be negative: {config.excluded_base_percentage}")

    if config.excluded_base_percentage >= Decimal(1):
        raise ValueError(f"Excluded base percentage cannot be 1 or greater: {config.excluded_base_percentage}")

    yield_in_gbase = (1-config.excluded_base_percentage) * config.read_count * BASES_PER_READ // BASES_PER_GBASE

    if yield_in_gbase.as_tuple()[0] != 0:
        raise ValueError(f"Resulting yield is negative: {yield_in_gbase}")

    if yield_in_gbase.as_tuple()[2] != 0:
        raise ValueError(f"Resulting yield is not in the form of an integer: {yield_in_gbase}")

    print(yield_in_gbase)


def parse_args(sys_args: List[str]) -> Config:
    parser = argparse.ArgumentParser(
        prog="get_yield_from_read_count_and_excluded_base_percentage.py",
        description=(
            "Calculate yield in gbases from the read count and PCT_EXC_TOTAL as output by Picard CollectWgsMetrics. "
            "Formula is ((1-exc) * reads * BASES_PER_READ / BASES_PER_GBASE). The result is rounded down to an integer."
        ),
    )
    parser.add_argument("--read_count", "-c", type=int, required=True, help="Read count.")
    parser.add_argument(
        "--excluded_base_percentage",
        "-e",
        type=Decimal,
        required=True,
        help="The value Picard CollectWgsMetrics would output as PCT_EXC_TOTAL.",
    )

    args = parser.parse_args(sys_args)

    config = Config(args.read_count, args.excluded_base_percentage)
    return config


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
