import argparse
import csv
import os
import sys

from typing import List, NamedTuple

import openpyxl


def main():
    config = parse_args()
    check_arg_input_file_exists(config.input_report)
    data = load_tsv(config.input_report)
    outpath = generate_outpath(config.input_report)
    generate_workbook(data, outpath)

def load_tsv(inpath: str) -> List[List[str]]:
    rows = []
    with open(inpath) as file:
         tsv_file = csv.reader(file, delimiter="\t")
         for row in tsv_file:
             rows.append(row)
    return rows

def generate_outpath(inpath: str) -> str:
   return f"{os.path.expanduser('~')}/{inpath.split('/')[-1].split('.')[0]}.xlsx"

def generate_workbook(data: List[List[str]], outpath: str):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Sheet1"
    for row in data:
        ws.append(row[0:9])
    wb.save(outpath)

def check_arg_input_file_exists(path_input_file: str):
    if not os.path.isfile(path_input_file):
        print(f"ERROR: File Not Found '{path_input_file}'")
        sys.exit(1)

def parse_args() -> NamedTuple:
    parser = argparse.ArgumentParser(prog="convert_shallow_tsv_to_xlsx.py", description="Convert ShallowSeq tsv to xlsx format")
    parser.add_argument("--input_report", "-i", type=str, required=True, help="path to input report")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    main()