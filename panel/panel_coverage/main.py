from typing import List, Tuple

from config import PanelFileConfig, AnalysisTodoConfig, Config
from analysis import do_analysis


def main() -> None:
    # See gs://func-79-panel-v1-coverage-analysis/config in hmf-crunch for panel config files
    # Depth files are output of samtools depth
    panel_file_config = PanelFileConfig(
        all_genes_tsv="/home/david/config/all_genes.38.tsv",
        baf_sites_list="/home/david/config/baf_points_list_v1.tsv",
        fusion_sites_list="/home/david/config/fusion_intron_list_v1.tsv",
        gene_list="/home/david/config/gene_list_v2.txt",
        hotspot_list="/home/david/config/hotspot_list_grch38_v1.tsv",
        msi_sites_list="/home/david/config/msi_list_v1.tsv",
        pgx_sites_list="/home/david/config/pgx_target_list_v1.tsv",
        tert_site="/home/david/config/tert_promoter_list_v1.tsv",
    )

    sample_with_depth_file_list: List[Tuple[str, str]] = [
        ("COLOT20NGU_umi_collapse", "/home/david/depths/COLOT20NGU.umi_collapse.bam.depth"),
        ("COLOT35NGU_umi_collapse", "/home/david/depths/COLOT35NGU.umi_collapse.bam.depth"),
        ("COLOT50NGU_umi_collapse", "/home/david/depths/COLOT50NGU.umi_collapse.bam.depth"),
        ("COLOT65NGU_umi_collapse", "/home/david/depths/COLOT65NGU.umi_collapse.bam.depth"),
        ("COLOT80NGU_umi_collapse", "/home/david/depths/COLOT80NGU.umi_collapse.bam.depth"),
        ("COLO100NGU_umi_collapse", "/home/david/depths/COLO100NGU.umi_collapse.bam.depth"),
        ("COLO150NGU_umi_collapse", "/home/david/depths/COLO150NGU.umi_collapse.bam.depth"),
        ("COLO200NGU_umi_collapse", "/home/david/depths/COLO200NGU.umi_collapse.bam.depth"),
        ("COLOT50NGU_no_dedup", "/home/david/depths/COLOT50NGU.bam.depth"),
        ("COLOT50NGU_sambamba", "/home/david/depths/COLOT50NGU.sambamba.bam.depth"),
        ("COLOT50NGU_umi_tools", "/home/david/depths/COLOT50NGU.umi_tools.bam.depth"),
    ]
    output_dir_100 = "/home/david/coverage_analysis_100/"
    output_dir_500 = "/home/david/coverage_analysis_500/"

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
    config_100 = Config(
        tuple(sample_with_depth_file_list), panel_file_config, analysis_todo_config, 100, output_dir_100
    )
    config_500 = Config(
        tuple(sample_with_depth_file_list), panel_file_config, analysis_todo_config, 500, output_dir_500
    )

    do_analysis(config_100)
    do_analysis(config_500)


if __name__ == "__main__":
    main()
