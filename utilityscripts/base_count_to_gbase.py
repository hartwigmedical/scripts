import argparse
import sys
from decimal import Decimal, ROUND_DOWN
from typing import List, NamedTuple, Optional


BASES_IN_GBASE = 10**9


class Config(NamedTuple):
    base_count: int
    round_like: Optional[Decimal]


def main(config: Config) -> None:
    if config.base_count < 0:
        raise ValueError(f"Base count cannot be negative: {config.base_count}")

    if config.round_like is not None:
        round_like = add_significant_digits_if_less_precise_than_integer(config.round_like)

        # Round down to the same number of decimals as round_like.
        gbases = (Decimal(config.base_count) / BASES_IN_GBASE).quantize(round_like, rounding=ROUND_DOWN)
    else:
        gbases = Decimal(config.base_count) / BASES_IN_GBASE

    # Strip as many zeroes from the right as possible
    normalized_gbases = gbases.normalize()

    # If we stripped so many zeroes that scientific notation is necessary to represent the number (e.g. 1200 -> 1.2e3),
    # increase the precision to get the result back to an integer (so back to 1200).
    pretty_gbases = add_significant_digits_if_less_precise_than_integer(normalized_gbases)

    print(pretty_gbases)


def add_significant_digits_if_less_precise_than_integer(number: Decimal) -> Decimal:
    """
    If input does not have enough significant digits to be cleanly represented as an integer,
    add significant zeroes until it does.

    In other words, if the exponent of the input Decimal object is greater than 0,
    return a Decimal object that is equal to the input as a number but that has exponent 0.
    In all other cases, return the input Decimal object.

    1.2e3 -> 1200
    1314 -> 1314
    14284.13214 -> 14284.13214
    1e-9 -> 1e-9
    3e4 -> 30000
    """
    exponent = number.as_tuple()[2]
    if exponent > 0:
        # This doesn't actually do any rounding
        result = number.quantize(Decimal("1"), rounding=ROUND_DOWN)
    else:
        result = number

    if result != number:
        raise ValueError(f"Accidentally rounded decimal: input={number}, result={result}")

    return result


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
