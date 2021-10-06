import argparse
import sys
from pathlib import Path
from typing import List, Tuple

from config import PanelFileConfig, AnalysisTypeConfig, AnalysisConfig, Config
from analysis import do_analysis


def main(config: Config) -> None:
    # See gs://hmf-crunch-experiments/210518_david_FUNC-79_panel-v1-coverage-analysis/config/ in hmf-crunch for panel config files
    # Depth files are output of samtools depth
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

    analysis_type_config = AnalysisTypeConfig(
        baf=True,
        exome=True,
        fusion=True,
        hotspot=True,
        msi=True,
        pgx=True,
        tert=True,
    )

    # sample_with_depth_file_list: List[Tuple[str, str]] = [
    #     ("COLO829vP01R", Path.home() / "depths/COLO829vP01R.bam.depth"),
    #     ("COLO829vP01T", Path.home() / "depths/COLO829vP01T.bam.depth"),
    #     ("GIAB12878v1T", Path.home() / "depths/GIAB12878v1T.bam.depth"),
    #     ("GIAB12878v2T", Path.home() / "depths/GIAB12878v2T.bam.depth"),
    #     ("GIAB12878v3T", Path.home() / "depths/GIAB12878v3T.bam.depth"),
    #     ("COLO100TCP", Path.home() / "depths/COLO100TCPT.bam.depth"),
    #     ("COLOT50NG", Path.home() / "depths/COLOT50NG.bam.depth"),
    # ]
    sample_with_depth_file_list: List[Tuple[str, Path]] = [
        ("FR12232711.non_umi_dedup", Path.home() / "depths/FR12232711.non_umi_dedup.bam.depth"),
        ("FR16673161.non_umi_dedup", Path.home() / "depths/FR16673161.non_umi_dedup.bam.depth"),
        ("FR30728561.non_umi_dedup", Path.home() / "depths/FR30728561.non_umi_dedup.bam.depth"),
        ("FR30728562.non_umi_dedup", Path.home() / "depths/FR30728562.non_umi_dedup.bam.depth"),
    ]

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


def parse_args(sys_args: List[str]) -> Config:
    parser = argparse.ArgumentParser(
        prog="panel_coverage",
        description=(
            "Determine coverage performance of panel bams at interesting sites."
        ),
    )
    parser.add_argument("--panel_config_dir", "-p", type=Path, required=True, help="Dir with panel config files.")
    parser.add_argument("--output_dir", "-o", type=Path, required=True, help="Output dir.")
    parser.add_argument(
        "--min_coverage", "-c", type=int, required=True, help="Min coverage. Can be specificied multiple times."
    )
    args = parser.parse_args(sys_args)

    config = Config(args.panel_config_dir, args.output_dir, tuple(args.min_coverage))
    return config


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
