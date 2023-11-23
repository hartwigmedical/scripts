#!/usr/bin/env python3
import argparse
import gzip
import logging
import statistics
import subprocess
from dataclasses import dataclass
from enum import Enum, auto
from functools import total_ordering

from pathlib import Path
from typing import List, Dict, Set, Optional


@dataclass
class Config:
    input_vcf: str
    sample_name: str
    regions_path: str
    blacklist_path: Optional[str]
    working_directory: Path
    output_directory: str


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


@dataclass(frozen=True, order=True)
class Coordinate:
    chromosome: Chromosome
    position: int


@dataclass(frozen=True, order=True)
class Region:
    chromosome: Chromosome
    chromosome_string: str
    start: int  # inclusive
    end: int  # inclusive
    name: str

    def get_coordinates(self) -> Set[Coordinate]:
        return {Coordinate(self.chromosome, position) for position in range(self.start, self.end+1)}


@dataclass(frozen=True)
class CoverageSummary:
    mean_depth: float
    median_depth: float
    _min_depth_to_percentage: Dict[int, float]

    def get_percentage_for_min_depth(self, min_depth: int) -> float:
        return self._min_depth_to_percentage[min_depth]


DEPTHS_TO_CHECK = [
    1, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000
]
OUTPUT_DELIMITER = "\t"
REGION_DELIMITER = "_"
GCP_PATH_START = "gs://"
DEFAULT_WORKING_DIRECTORY = Path("/data/tmp")


def main(config: Config) -> None:
    logging.info("Start summarizing coverage")

    local_config = get_local_only_config(config)

    do_local_processing(local_config)

    if is_gcp_path(config.output_directory):
        copy_directory_gcp(local_config.output_directory, config.output_directory)

    logging.info("Finished summarizing coverage")


def do_local_processing(local_config: Config) -> None:
    logging.info("Load regions")
    regions = get_regions(local_config.regions_path)
    logging.info("Load depths from input VCF")
    coordinate_to_depth = get_depths(local_config.input_vcf, local_config.sample_name)
    create_output_files(
        coordinate_to_depth, regions, local_config.sample_name, local_config.output_directory, "all")
    if local_config.blacklist_path is not None:
        logging.info("Write blacklist filtered output files")
        blacklisted_region_ids = get_blacklisted_region_ids(local_config.blacklist_path)
        non_blacklisted_regions = [region for region in regions if region.name not in blacklisted_region_ids]
        create_output_files(
            coordinate_to_depth,
            non_blacklisted_regions,
            local_config.sample_name,
            local_config.output_directory,
            "non_blacklist"
        )


def get_local_only_config(config: Config) -> Config:
    if is_gcp_path(config.input_vcf):
        local_input_vcf_path = f"{config.working_directory}/{get_file_name(config.input_vcf)}"
        copy_file_gcp(config.input_vcf, local_input_vcf_path)
    else:
        local_input_vcf_path = config.input_vcf
    if is_gcp_path(config.regions_path):
        local_regions_path = f"{config.working_directory}/{get_file_name(config.regions_path)}"
        copy_file_gcp(config.regions_path, local_regions_path)
    else:
        local_regions_path = config.regions_path
    if is_gcp_path(config.blacklist_path):
        local_blacklist_path = f"{config.working_directory}/{get_file_name(config.blacklist_path)}"
        copy_file_gcp(config.blacklist_path, local_blacklist_path)
    else:
        local_blacklist_path = config.blacklist_path
    if is_gcp_path(config.output_directory):
        local_output_directory = config.working_directory / "output"
    else:
        local_output_directory = config.output_directory
    
    local_config = Config(
        local_input_vcf_path,
        config.sample_name,
        local_regions_path,
        local_blacklist_path,
        config.working_directory,
        local_output_directory,
    )
    return local_config


def get_file_name(path: str) -> str:
    return path.split("/")[-1]


def is_gcp_path(path: str) -> bool:
    return path[:len(GCP_PATH_START)] == GCP_PATH_START


def copy_file_gcp(source: str, target: str) -> None:
    logging.info(f"Copy {source} to {target}")
    if not is_gcp_path(target):
        Path(target).parent.mkdir(parents=True, exist_ok=True)
    run_bash_command(["gsutil", "-m", "cp", source, target])


def copy_directory_gcp(source: str, target: str) -> None:
    logging.info(f"Copy {source} to {target}")
    if not is_gcp_path(target):
        Path(target).mkdir(parents=True, exist_ok=True)
    run_bash_command(["gsutil", "-m", "rsync", "-r", source, target])


def get_blacklisted_region_ids(blacklist_path: str) -> Set[str]:
    region_ids = set()
    with open(blacklist_path, "r") as in_f:
        for line in in_f:
            region_id = line.rstrip("\n")
            if region_id:
                region_ids.add(region_id)
    return region_ids


def create_output_files(
        coordinate_to_depth: Dict[Coordinate, int], 
        regions: List[Region], 
        sample_name: str, 
        output_directory: str, 
        region_description: str,
) -> None:
    logging.info(f"Create output file with coverage summary of {region_description} regions")
    region_to_coverage_summary = {
        region: get_coverage_summary([region], coordinate_to_depth) for region in regions
    }
    create_per_region_output_file(
        regions, region_to_coverage_summary, f"{output_directory}/{sample_name}.{region_description}.region.coverage.tsv"
    )
    logging.info(f"Create output file with coverage summary of {region_description} regions with less than 95% 100x")
    create_per_region_output_file(
        [region for region in regions if region_to_coverage_summary[region].get_percentage_for_min_depth(100) < 0.95],
        region_to_coverage_summary,
        f"{output_directory}/{sample_name}.{region_description}.region.less_than_95_percent_below_100.coverage.tsv",
    )
    logging.info(f"Create output file with coverage summary of {region_description} regions with less than 95% 200x")
    create_per_region_output_file(
        [region for region in regions if region_to_coverage_summary[region].get_percentage_for_min_depth(200) < 0.95],
        region_to_coverage_summary,
        f"{output_directory}/{sample_name}.{region_description}.region.less_than_95_percent_below_200.coverage.tsv",
    )
    if all(REGION_DELIMITER in region.name for region in regions):
        logging.info(f"Create output file with coverage summary per type for {region_description} regions")
        create_per_type_output_file(
            regions, coordinate_to_depth, f"{output_directory}/{sample_name}.{region_description}.type.coverage.tsv"
        )


def create_per_region_output_file(
        regions: List[Region], 
        region_to_coverage_summary: Dict[Region, CoverageSummary], 
        output_file_path: str,
) -> None:
    per_region_output_lines = [
        OUTPUT_DELIMITER.join(["chromosome", "start", "end", "name", "meanDepth", "medianDepth"] + [f"PCT_{depth}X" for depth in DEPTHS_TO_CHECK])
    ]
    for index, region in enumerate(regions):
        line = OUTPUT_DELIMITER.join(
            [
                region.chromosome_string,
                str(region.start),
                str(region.end),
                region.name,
                get_coverage_summary_line(region_to_coverage_summary[region]),
            ]
        )
        per_region_output_lines.append(line)
        if index % 1000 == 0 and index > 0:
            logging.info(f"Handled {index} regions")

    Path(output_file_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_file_path, "w") as out_f:
        out_f.write("\n".join(per_region_output_lines) + "\n")


def create_per_type_output_file(
        regions: List[Region], coordinate_to_depth: Dict[Coordinate, int], output_file_path: str
) -> None:
    per_type_output_lines = [
        OUTPUT_DELIMITER.join(["type", "meanDepth", "medianDepth"] + [f"PCT_{depth}X" for depth in DEPTHS_TO_CHECK])
    ]
    region_to_region_type = {region: region.name.split(REGION_DELIMITER)[0] for region in regions}
    for region_type in sorted(set(region_to_region_type.values())):
        logging.info(f"Handling {region_type}")
        coverage_summary = get_coverage_summary(
            [region for region in regions if region_to_region_type[region] == region_type], coordinate_to_depth
        )
        per_type_output_lines.append(OUTPUT_DELIMITER.join([region_type, get_coverage_summary_line(coverage_summary)]))

    logging.info(f"Handling TOTAL")
    coverage_summary = get_coverage_summary(regions, coordinate_to_depth)
    per_type_output_lines.append(OUTPUT_DELIMITER.join(["TOTAL", get_coverage_summary_line(coverage_summary)]))
    
    logging.info("Writing output file")
    Path(output_file_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_file_path, "w") as out_f:
        out_f.write("\n".join(per_type_output_lines) + "\n")


def get_coverage_summary_line(coverage_summary: CoverageSummary) -> str:
    line = OUTPUT_DELIMITER.join(
        [f"{coverage_summary.mean_depth:.1f}", str(coverage_summary.median_depth)] 
            + [f"{coverage_summary.get_percentage_for_min_depth(min_depth):.3f}" for min_depth in DEPTHS_TO_CHECK]
    )
    return line


def get_coverage_summary(regions: List[Region], coordinate_to_depth: Dict[Coordinate, int]) -> CoverageSummary:
    relevant_coordinates = {coordinate for region in regions for coordinate in region.get_coordinates()}
    depths = [coordinate_to_depth[coordinate] for coordinate in relevant_coordinates]
    mean_depth = statistics.mean(depths)
    median_depth = statistics.median(depths)
    min_depth_to_percentage = get_min_depth_to_percentage(depths)
    coverage_summary = CoverageSummary(mean_depth, median_depth, min_depth_to_percentage)
    return coverage_summary


def get_min_depth_to_percentage(depths: List[int]) -> Dict[int, int]:
    min_depth_to_count = get_min_depth_to_count(depths)
    return {min_depth: min_depth_to_count[min_depth] / len(depths) for min_depth in DEPTHS_TO_CHECK}


def get_min_depth_to_count(depths: List[int]) -> Dict[int, int]:
    # For each depth, determine the highest coverage threshold it surpasses, if any.
    # For each coverage threshold, count the depths for which this occurs in 'min_depth_to_count_first_threshold'.
    # The 'min_depth_to_count' values can then efficiently be determined from the values in 'min_depth_to_count_first_threshold'
    min_depth_to_count_first_threshold = {min_depth: 0 for min_depth in DEPTHS_TO_CHECK}
    for depth in depths:
        for min_depth in reversed(DEPTHS_TO_CHECK):
            if depth >= min_depth:
                min_depth_to_count_first_threshold[min_depth] += 1
                break
    
    min_depth_to_count = {min_depth: 0 for min_depth in DEPTHS_TO_CHECK}
    for min_depth in DEPTHS_TO_CHECK:
        relevant_counts = [
            min_depth_to_count_first_threshold[source_min_depth]
            for source_min_depth in DEPTHS_TO_CHECK if source_min_depth >= min_depth
        ]
        min_depth_to_count[min_depth] = sum(relevant_counts)
    
    return min_depth_to_count


def get_regions(regions_path: str) -> List[Region]:
    if regions_path.endswith(".tsv"):
        with open(regions_path, "r") as in_f:
            return [get_region_from_tsv_line(line) for line in in_f]
    elif regions_path.endswith(".bed"):
        with open(regions_path, "r") as in_f:
            return [get_region_from_bed_line(line) for line in in_f]
    else:
        raise ValueError(f"File extension of regions file is not '.tsv' or '.bed': {regions_path}")


def get_region_from_bed_line(line: str) -> Region:
    split_line = line.rstrip("\n").split("\t")
    chromosome_string = split_line[0]
    chromosome = Chromosome.from_string(chromosome_string)
    start = int(split_line[1]) + 1
    end = int(split_line[2])
    name = split_line[3]
    return Region(chromosome, chromosome_string, start, end, name)


def get_region_from_tsv_line(line: str) -> Region:
    split_line = line.rstrip("\n").split("\t")
    chromosome_string = split_line[0]
    chromosome = Chromosome.from_string(chromosome_string)
    start = int(split_line[1])
    end = int(split_line[2])
    name = split_line[3]
    return Region(chromosome, chromosome_string, start, end, name)


def get_depths(input_vcf: str, sample_name: str) -> Dict[Coordinate, int]:
    coordinate_to_depth = {}
    sample_index = -1
    depth_index = -1
    with gzip.open(input_vcf, "rt") as in_f:
        for line in in_f:
            if line[0:2] != "##":
                split_line = line.rstrip("\n").split("\t")
                if line[0] == "#":
                    sample_index = split_line.index(sample_name)
                else:
                    chromosome = Chromosome.from_string(split_line[0])
                    position = int(split_line[1])
                    coordinate = Coordinate(chromosome, position)
                    if coordinate in coordinate_to_depth.keys():
                        raise ValueError(f"Coordinate in VCF multiple times: {coordinate.chromosome}:{coordinate.position}")

                    if depth_index == -1:
                        info_format = split_line[8]
                        depth_index = info_format.split(":").index("DP")

                    sample_info = split_line[sample_index]
                    depth = int(sample_info.split(":")[depth_index])
                    coordinate_to_depth[coordinate] = depth
    return coordinate_to_depth


def run_bash_command(command: List[str]) -> str:
    return subprocess.run(command, capture_output=True, text=True).stdout


def parse_args() -> Config:
    parser = argparse.ArgumentParser(prog="summarize_sage_vcf_coverage.py", description="Summarize coverage from SAGE append VCF")
    parser.add_argument("--input_vcf", "-i", type=str, required=True, help="Input VCF. Can be local or GCP bucket path.")
    parser.add_argument("--sample_name", "-n", type=str, required=True, help="Sample name in VCF (and output files).")
    parser.add_argument("--regions", "-r", type=str, required=True, help="TSV or BED file defining the regions to look at. Can be local or GCP bucket path.")
    parser.add_argument("--output_directory", "-o", type=str, required=True, help="Output directory. Can be local or GCP bucket path.")
    parser.add_argument("--blacklist", "-b", type=str, default=None, help="Optional file with blacklisted region IDs. Can be local or GCP bucket path.")
    parser.add_argument("--working_directory", "-w", type=Path, default=DEFAULT_WORKING_DIRECTORY, help=f"Local working directory, if needed. Default value is {DEFAULT_WORKING_DIRECTORY}")
    args = parser.parse_args()
    return Config(args.input_vcf, args.sample_name, args.regions, args.blacklist, args.working_directory, args.output_directory)


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s - [%(levelname)-8s] - %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
    )
    main(parse_args())
