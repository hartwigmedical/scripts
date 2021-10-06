from pathlib import Path
from typing import List, Tuple

from config import PanelFileConfig, AnalysisTodoConfig, Config
from analysis import do_analysis


def main() -> None:
    # See gs://hmf-crunch-experiments/210518_david_FUNC-79_panel-v1-coverage-analysis/config/ in hmf-crunch for panel config files
    # Depth files are output of samtools depth
    panel_file_config = PanelFileConfig(
        all_genes_tsv=Path.home() / "config/all_genes.38.tsv",
        baf_sites_list=Path.home() / "config/baf_points_list_v1.tsv",
        fusion_sites_list=Path.home() / "config/fusion_intron_list_v1.tsv",
        gene_list=Path.home() / "config/gene_list_v2.txt",
        hotspot_list=Path.home() / "config/hotspot_list_grch38_v1.tsv",
        msi_sites_list=Path.home() / "config/msi_list_v1.tsv",
        pgx_sites_list=Path.home() / "config/pgx_target_list_v1.tsv",
        tert_site=Path.home() / "config/tert_promoter_list_v1.tsv",
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
    output_dir_50 = Path.home() / "coverage_analysis_50/"
    output_dir_100 = Path.home() / "coverage_analysis_100/"
    output_dir_200 = Path.home() / "coverage_analysis_200/"
    # output_dir_500 = Path.home() / "coverage_analysis_500/"

    analysis_todo_config = AnalysisTodoConfig(
        baf=True,
        exome=True,
        fusion=True,
        hotspot=True,
        msi=True,
        pgx=True,
        tert=True,
    )
    # analysis_todo_config = AnalysisTodoConfig(
    #     baf=False,
    #     exome=True,
    #     fusion=False,
    #     hotspot=False,
    #     msi=False,
    #     pgx=False,
    #     tert=False,
    # )
    config_50 = Config(
        tuple(sample_with_depth_file_list), panel_file_config, analysis_todo_config, 50, output_dir_50
    )
    config_100 = Config(
        tuple(sample_with_depth_file_list), panel_file_config, analysis_todo_config, 100, output_dir_100
    )
    config_200 = Config(
        tuple(sample_with_depth_file_list), panel_file_config, analysis_todo_config, 200, output_dir_200
    )
    # config_500 = Config(
    #     tuple(sample_with_depth_file_list), panel_file_config, analysis_todo_config, 500, output_dir_500
    # )
    #
    do_analysis(config_50)
    do_analysis(config_100)
    do_analysis(config_200)
    # do_analysis(config_500)


if __name__ == "__main__":
    main()
