#!/usr/bin/env python3
import argparse

# We only need the first 9 headers for the output
cutoff_i = 9

parser = argparse.ArgumentParser(description="Stages the shallow sample sheet so it can be parsed in email")
parser.add_argument("input_tsv")

def main():
    args = parser.parse_args()

    with open(args.input_tsv, 'r') as input:
        data = [l.strip().split('\t')[:cutoff_i] for l in input.readlines()]

    print(*['\t'.join(r) for r in data], sep='\n')



if __name__ == '__main__':
    main()