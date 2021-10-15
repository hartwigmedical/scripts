import argparse
import sys
from decimal import Decimal, ROUND_DOWN, InvalidOperation
from typing import List, NamedTuple, Optional


BASES_IN_GBASE = 10**9


class Config(NamedTuple):
    base_count: int
    round_like: Optional[Decimal]


def main(config: Config) -> None:
    gbases = Decimal(config.base_count) / BASES_IN_GBASE
    if config.round_like is not None:
        round_like = Decimal(config.round_like)
        sign, digits, exponent = round_like.as_tuple()
        if exponent > 0:
            round_like = round_like.quantize(Decimal("1"), rounding=ROUND_DOWN)

        gbases = gbases.quantize(round_like, rounding=ROUND_DOWN)
    print(gbases.normalize())


def parse_args(sys_args: List[str]) -> Config:
    parser = argparse.ArgumentParser(
        prog="base_count_to_gbase.py",
        description=(
            "Transform string representing a number of bases to a string representing a number of gbase."
        ),
    )
    parser.add_argument("--base_count", "-i", type=int, required=True, help="Base count.")
    parser.add_argument(
        "--round_like",
        "-r",
        type=Decimal,
        help="Round base count in gbase down to at most the same number of decimals as this number",
    )

    args = parser.parse_args(sys_args)

    config = Config(args.base_count, args.round_like)
    return config


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
