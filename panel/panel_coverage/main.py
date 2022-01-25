import argparse
import logging
import sys
from pathlib import Path
from typing import List

from config import PanelFileConfig, AnalysisTypeConfig, ProgramConfig
from analysis import do_analysis
from gcp.base import GCPPath
from util import assert_file_exists


def main(program_config: ProgramConfig) -> None:
    # See gs://hmf-crunch-experiments/210518_david_FUNC-79_panel-v1-coverage-analysis/config/
    # in hmf-crunch for panel config files

    set_up_logging()

    panel_file_config = PanelFileConfig(
        all_genes_tsv=program_config.panel_config_dir / "all_genes.38.tsv",
        baf_sites_list=program_config.panel_config_dir / "baf_points_list_v1.tsv",
        fusion_sites_list=program_config.panel_config_dir / "fusion_intron_list_v1.tsv",
        gene_list=program_config.panel_config_dir / "gene_list_v2.txt",
        hotspot_list=program_config.panel_config_dir / "hotspot_list_grch38_v1.tsv",
        msi_sites_list=program_config.panel_config_dir / "msi_list_v1.tsv",
        pgx_sites_list=program_config.panel_config_dir / "pgx_target_list_v1.tsv",
        tert_site=program_config.panel_config_dir / "tert_promoter_list_v1.tsv",
    )
    panel_file_config.validate()
    assert_file_exists(program_config.samtools)

    analysis_type_config = AnalysisTypeConfig(
        baf=True,
        exome=True,
        fusion=True,
        hotspot=True,
        msi=True,
        pgx=True,
        tert=True,
    )

    do_analysis(program_config, panel_file_config, analysis_type_config)


def set_up_logging() -> None:
    logging.basicConfig(
        format="%(asctime)s - [%(levelname)-8s] - %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
    )


def parse_args(sys_args: List[str]) -> ProgramConfig:
    parser = argparse.ArgumentParser(
        prog="panel_coverage",
        description=(
            "Determine coverage performance of panel bams at interesting sites."
        ),
    )
    parser.add_argument("--panel_config_dir", "-p", type=Path, required=True, help="Dir with panel config files.")
    parser.add_argument("--output_dir", "-o", type=GCPPath.from_string, required=True, help="Output GCP dir.")
    parser.add_argument("--samtools", "-s", type=Path, required=True, help="Samtools version 1.13 or greater.")
    parser.add_argument(
        "--working_dir", "-w", type=Path, required=True, help="Working dir to store intermediate files."
    )
    parser.add_argument(
        "--min_coverage", "-c", type=int, required=True, action="append", help="Min coverage. Can be specified multiple times."
    )
    parser.add_argument(
        "--bam", "-b", type=GCPPath.from_string, required=True, action="append", help="GCP path to bam. Can be specified multiple times."
    )
    args = parser.parse_args(sys_args)

    sorted_min_coverages: List[int] = sorted(args.min_coverage)
    config = ProgramConfig(
        args.panel_config_dir,
        args.output_dir,
        args.samtools,
        args.working_dir,
        tuple(sorted_min_coverages),
        tuple(args.bam),
    )
    return config


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
