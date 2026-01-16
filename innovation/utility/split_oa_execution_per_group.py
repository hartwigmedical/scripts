#!/usr/bin/env python3

import argparse, sys, yaml
import logging
from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import List

class Dumper(yaml.Dumper):
    def increase_indent(self, flow=False, indentless=False):
        return super().increase_indent(flow, False)

def str_representer(dumper, data):
    if '\n' in data:
        return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')
    return dumper.represent_scalar('tag:yaml.org,2002:str', data)

Dumper.add_representer(str, str_representer)


@dataclass(frozen=True)
class Config:
    input_path: Path
    output_dir: Path


def main(config: Config) -> None:
    logging.info("Started splitting OA execution YAML")

    with open(config.input_path) as in_f:
        input_yaml = yaml.safe_load(in_f)

    samplesheet_text = input_yaml["params"]["samplesheet"]
    samplesheet_lines = samplesheet_text.split("\n")

    group_id_to_lines = defaultdict(list)
    for line in samplesheet_lines[1:]:
        group_id = line.split(",")[0]
        if group_id:
            group_id_to_lines[group_id].append(line)

    samplesheet_header = samplesheet_lines[0]

    for group_id, lines in group_id_to_lines.items():
        new_file_name = config.output_dir /f"{config.input_path.stem}.{group_id}{config.input_path.suffix}"
        new_yaml = deepcopy(input_yaml)
        new_yaml["name"] = f"{input_yaml['name']}-{group_id}"
        new_yaml["params"]["samplesheet"] = "\n".join([samplesheet_header, *lines]) + "\n"
        new_file_name.parent.mkdir(parents=True, exist_ok=True)

        with open(new_file_name, "w") as out_f:
            yaml.dump(
                new_yaml,
                out_f,
                Dumper=Dumper,
                sort_keys=False,
                default_flow_style=False,
            )

    logging.info("Finished splitting OA execution YAML")


def parse_args(sys_args: List[str]) -> Config:
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, type=Path, help="Execution YAML to split")
    ap.add_argument("-o", "--output", required=True, type=Path, help="Output directory")
    args = ap.parse_args(sys_args)

    return Config(args.input, args.output)


if __name__ == "__main__":
    config = parse_args(sys.argv[1:])
    main(config)