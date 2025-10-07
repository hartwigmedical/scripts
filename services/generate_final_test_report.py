#!/usr/bin/env python3

import argparse
import datetime
import os
import shutil
import sys

from typing import Dict, List, NamedTuple



def main():
    config = parse_args()
    check_arg_input_file_exists(config.input_report)
    if config.output_dir != None:
        check_arg_output_dir_exists(config.output_dir)
    else:
        config.output_dir = set_arg_output_dir(config.input_report)
    report_info = prepare_report_info(config)
    print("[INFO]: Input report parsed.")
    generate_html_report(report_info)
    print("[INFO]: html report generated.")
    print_command_convert_html_to_pdf_report(report_info['output_file'])


def check_args(config: NamedTuple):
    check_arg_input_file_exists(config.input_report)
    check_arg_output_dir_exists(config.output_dir)

def check_arg_input_file_exists(path_input_file: str):
    if not os.path.isfile(path_input_file):
        print(f"[ERROR]: File Not Found '{path_input_file}'")
        sys.exit(1)

def check_arg_output_dir_exists(path_output_dir: str):
    if not os.path.isdir(path_output_dir) or not path_output_dir.endswith("/"):
        print(f"[ERROR]: Directory Not Found '{path_output_dir}', make sure it ends with '/'")
        sys.exit(1)

def set_arg_output_dir(path_input_file: str):
    return "/".join(os.path.abspath(path_input_file).split("/")[:-1]) + "/"


def prepare_report_info(config: NamedTuple) -> Dict:
    report_info = parse_input_report(config.input_report)
    report_info['user_name'] = config.name
    report_info['output_file'] = generate_output_filepath(report_info['output_file'], config.output_dir)
    report_info['general_info'] = prepare_general_information(report_info['general_info'])
    report_info['sample_info'] = prepare_sample_information(report_info['sample_info'])
    return report_info

def parse_input_report(input_report: str) -> Dict:
    with open(input_report) as report_file:
        lines = report_file.readlines()

    sample_section_start_i = next(i for i, line in enumerate(lines) if line.strip() == "# Sample table for report:")

    output_filename = lines[0].split('(')[1].split('.')[0]
    # there is a hardcoded empty line before the sample section split, which we don't want to include in the general_info
    # hence we do sample_section_start_i - 1
    # this is very hacky
    general_info = lines[1:sample_section_start_i - 1]
    sample_info = lines[sample_section_start_i + 1:]
    return {'output_file': output_filename,
            'general_info': general_info,
            'sample_info': sample_info}

def generate_output_filepath(output_report: str, output_dir: str) -> str:
    return f"{output_dir}{output_report}.html"

def prepare_general_information(general_info: List[str]) -> List[str]:
    return [line.strip() for line in general_info]

def prepare_sample_information(sample_info: List[str]) -> List[List[str]]:
    return [[index+1] + line.strip().split("\t") for index, line in enumerate(sample_info)]


def generate_html_report(report_info: Dict):
    copy_template_file(report_info['output_file'])
    place_general_field_info(report_info['output_file'], report_info['general_info'], report_info['user_name'])
    place_sample_field_info(report_info['output_file'], report_info['sample_info'])

def copy_template_file(output_file: str):
    shutil.copy("/data/resources/ops/hartwig_sequencing_services/final_test_report_template.html", output_file)

def place_general_field_info(output_file: str, general_info: List[str], user_name: str):
    text = open_file(output_file)
    order_fields = ("{field_client_name}", "{field_client_email}", "{field_data_email}", "{field_dvo}",
                    "{field_sample_count}", "{field_output_level}", "{field_yield_generated}", "{field_create_date}",
                    "{field_checked_by}")
    for field, value in zip(order_fields, general_info + [get_date_today(), user_name]):
        text = text.replace(field, value)
    close_file(output_file, text)

def place_sample_field_info(output_file: str, sample_info: List[List[str]]):
    samples_formatted = ""
    for sample in sample_info:
        samples_formatted += "\t\t\t\t<tr>\n"
        for value in sample:
            samples_formatted += f'\t\t\t\t\t<td class="tg-value">{value}</td>\n'
        samples_formatted += "\t\t\t\t</tr>\n"
    text = open_file(output_file)
    text = text.replace("{field_samples}", samples_formatted)
    close_file(output_file, text)


def print_command_convert_html_to_pdf_report(output_file: str):
    out_file = os.path.abspath(output_file)
    print("Convert html to pdf format")
    print(f" wkhtmltopdf --enable-local-file-access {out_file} {out_file.split('.')[0]}.pdf")
    print("")
    print("And download the generated pdf report (Execute locally!)")
    print(f" gcloud compute scp ops-vm-prod-2:~/{out_file.split('.')[0].split('/')[-1]}.pdf ~/")
    print("")


def parse_args() -> NamedTuple:
    parser = argparse.ArgumentParser(prog="generate_hss_final_report.py", description="Generates Final Test Report")
    parser.add_argument("--input_report", "-i", type=str, required=True, help="path to input report")
    parser.add_argument("--name", "-n", type=str, required=True, help="name of user generating report")
    parser.add_argument("--output_dir", "-o", type=str, help="path to output directory, output file will share name with input")
    args = parser.parse_args()
    return args

def open_file(filepath: str) -> str:
    with open(filepath, 'r') as f:
        return f.read()

def close_file(filepath: str, text: str):
    with open(filepath, 'w') as f:
        f.write(text)

def get_date_today() -> str:
    return datetime.datetime.today().strftime("%B %d, %Y")


if __name__ == "__main__":
    main()
