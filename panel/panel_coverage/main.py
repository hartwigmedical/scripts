import argparse
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple

from config import PanelFileConfig, AnalysisTypeConfig, AnalysisConfig, Config
from analysis import do_analysis
from util import assert_file_exists


def main(config: Config) -> None:
    # See gs://hmf-crunch-experiments/210518_david_FUNC-79_panel-v1-coverage-analysis/config/ in hmf-crunch for panel config files

    panel_file_config = PanelFileConfig(
        all_genes_tsv=config.panel_config_dir / "all_genes.38.tsv",
        baf_sites_list=config.panel_config_dir / "baf_points_list_v1.tsv",
        fusion_sites_list=config.panel_config_dir / "fusion_intron_list_v1.tsv",
        gene_list=config.panel_config_dir / "gene_list_v2.txt",
        hotspot_list=config.panel_config_dir / "hotspot_list_grch38_v1.tsv",
        msi_sites_list=config.panel_config_dir / "msi_list_v1.tsv",
        pgx_sites_list=config.panel_config_dir / "pgx_target_list_v1.tsv",
        tert_site=config.panel_config_dir / "tert_promoter_list_v1.tsv",
    )
    panel_file_config.validate()
    assert_file_exists(config.samtools)
    for bam in config.bams:
        assert_file_exists(bam)

    analysis_type_config = AnalysisTypeConfig(
        baf=True,
        exome=True,
        fusion=True,
        hotspot=True,
        msi=True,
        pgx=True,
        tert=True,
    )

    config.working_dir.mkdir(parents=True, exist_ok=True)

    sample_with_depth_file_list: List[Tuple[str, Path]] = []
    for bam in config.bams:
        depth_file = config.working_dir / f"{bam.name}.depth"
        if not depth_file.exists():
            create_depth_file(config.samtools, bam, depth_file)
        if not depth_file.exists():
            raise FileNotFoundError(f"Depth file creation failed: {depth_file}")
        sample_with_depth_file_list.append((bam.stem, depth_file))

    for min_coverage in config.min_coverages:
        analysis_output_dir = config.output_dir / str(min_coverage)
        analysis_config = AnalysisConfig(
            tuple(sample_with_depth_file_list),
            panel_file_config,
            analysis_type_config,
            min_coverage,
            analysis_output_dir,
        )
        do_analysis(analysis_config)


def create_depth_file(samtools: Path, bam: Path, depth_file: Path) -> None:
    cli_args = [samtools, "depth", "-s", bam]
    with open(depth_file, "w") as depth_f:
        subprocess.run(cli_args, stdout=depth_f)


def parse_args(sys_args: List[str]) -> Config:
    parser = argparse.ArgumentParser(
        prog="panel_coverage",
        description=(
            "Determine coverage performance of panel bams at interesting sites."
        ),
    )
    parser.add_argument("--panel_config_dir", "-p", type=Path, required=True, help="Dir with panel config files.")
    parser.add_argument("--output_dir", "-o", type=Path, required=True, help="Output dir.")
    parser.add_argument("--samtools", "-s", type=Path, required=True, help="Samtools version 1.13 or greater.")
    parser.add_argument(
        "--working_dir", "-w", type=Path, required=True, help="Working dir to store intermediate files."
    )
    parser.add_argument(
        "--min_coverage", "-c", type=int, required=True, help="Min coverage. Can be specified multiple times."
    )
    parser.add_argument(
        "--bam", "-b", type=Path, required=True, help="Path to bam. Can be specified multiple times."
    )
    args = parser.parse_args(sys_args)

    sorted_min_coverages: List[int] = sorted(args.min_coverage)
    config = Config(
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
