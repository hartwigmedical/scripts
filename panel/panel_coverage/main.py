from typing import List, Tuple

from config import PanelFileConfig, AnalysisTodoConfig, Config
from analysis import do_analysis


def main() -> None:
    # See gs://func-79-panel-v1-coverage-analysis/config in hmf-crunch for panel config files
    # Depth files are output of samtools depth
    panel_file_config = PanelFileConfig(
        all_genes_tsv="/Users/davidkoetsier/Documents/all_genes.38.tsv",
        baf_sites_list="/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/baf_points_list_v1.tsv",
        fusion_sites_list="/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/fusion_intron_list_v1.tsv",
        gene_list="/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/gene_list_v2.txt",
        hotspot_list="/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/hotspot_list_grch38_v1.tsv",
        msi_sites_list="/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/msi_list_v1.tsv",
        pgx_sites_list="/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/pgx_target_list_v1.tsv",
        tert_site="/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/tert_promoter_list_v1.tsv",
    )
    min_coverage = 100
    # min_coverage = 500

    # sample_with_depth_file_list: List[Tuple[str, str]] = [
    #     ("COLO829vP01R", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO829vP01R.bam.depth"),
    #     ("COLO829vP01T", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO829vP01T.bam.depth"),
    #     ("GIAB12878v1T", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/GIAB12878v1T.bam.depth"),
    #     ("GIAB12878v2T", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/GIAB12878v2T.bam.depth"),
    #     ("GIAB12878v3T", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/GIAB12878v3T.bam.depth"),
    #     ("GIAB23485v1T", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/GIAB23485v1T.bam.depth"),
    #     ("TSO500", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/TSO500.bam.depth"),
    #     ("COLO829vP01R_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO829vP01R.bam.all.depth"),
    #     ("COLO829vP01T_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO829vP01T.bam.all.depth"),
    #     ("GIAB12878v1T_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/GIAB12878v1T.bam.all.depth"),
    #     ("GIAB12878v2T_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/GIAB12878v2T.bam.all.depth"),
    #     ("GIAB12878v3T_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/GIAB12878v3T.bam.all.depth"),
    #     ("GIAB23485v1T_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/GIAB23485v1T.bam.all.depth"),
    #     ("TSO500_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/TSO500.bam.all.depth"),
    # ]
    # output_dir = "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/coverage_analysis/run3"

    # sample_with_depth_file_list = [
    #     ("COLOT20NG", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLOT20NG.bam.depth"),
    #     ("COLOT35NG", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLOT35NG.bam.depth"),
    #     ("COLOT50NG", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLOT50NG.bam.depth"),
    #     ("COLOT65NG", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLOT65NG.bam.depth"),
    #     ("COLOT80NG", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLOT80NG.bam.depth"),
    #     ("COLO100NG", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO100NG.bam.depth"),
    #     ("COLO150NG", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO150NG.bam.depth"),
    #     ("COLO200NG", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO200NG.bam.depth"),
    # ]
    # output_dir = "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/coverage_analysis/run2"

    # sample_with_depth_file_list: List[Tuple[str, str]] = [
    #     ("COLO100TCP", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO100TCPT.bam.depth"),
    #     ("COLOT50NG", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLOT50NG.bam.depth"),
    #     ("NA12878v1", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/NA12878v1T.bam.depth"),
    #     ("COLO100TCP_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO100TCPT.bam.all.depth"),
    #     ("COLOT50NG_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLOT50NG.bam.all.depth"),
    #     ("NA12878v1_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/NA12878v1T.bam.all.depth"),
    # ]
    # output_dir = "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/coverage_analysis/run4"

    sample_with_depth_file_list: List[Tuple[str, str]] = [
        ("COLO829vP01R", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO829vP01R.bam.depth"),
        ("COLO829vP01T", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO829vP01T.bam.depth"),
        ("GIAB12878v1T", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/GIAB12878v1T.bam.depth"),
        ("GIAB12878v2T", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/GIAB12878v2T.bam.depth"),
        ("GIAB12878v3T", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/GIAB12878v3T.bam.depth"),
        ("GIAB23485v1T", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/GIAB23485v1T.bam.depth"),
        ("TSO500", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/TSO500.bam.depth"),
        ("COLO829vP01R_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO829vP01R.bam.all.depth"),
        ("COLO829vP01T_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO829vP01T.bam.all.depth"),
        ("GIAB12878v1T_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/GIAB12878v1T.bam.all.depth"),
        ("GIAB12878v2T_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/GIAB12878v2T.bam.all.depth"),
        ("GIAB12878v3T_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/GIAB12878v3T.bam.all.depth"),
        ("GIAB23485v1T_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/GIAB23485v1T.bam.all.depth"),
        ("TSO500_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/TSO500.bam.all.depth"),
        ("COLOT20NG", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLOT20NG.bam.depth"),
        ("COLOT35NG", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLOT35NG.bam.depth"),
        ("COLOT50NG", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLOT50NG.bam.depth"),
        ("COLOT65NG", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLOT65NG.bam.depth"),
        ("COLOT80NG", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLOT80NG.bam.depth"),
        ("COLO100NG", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO100NG.bam.depth"),
        ("COLO150NG", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO150NG.bam.depth"),
        ("COLO200NG", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO200NG.bam.depth"),
        ("COLO100TCP", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO100TCPT.bam.depth"),
        ("NA12878v1", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/NA12878v1T.bam.depth"),
        ("COLO100TCP_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLO100TCPT.bam.all.depth"),
        ("COLOT50NG_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/COLOT50NG.bam.all.depth"),
        ("NA12878v1_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/NA12878v1T.bam.all.depth"),
    ]
    output_dir = "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/coverage_analysis/run5"

    # sample_with_depth_file_list: List[Tuple[str, str]] = [
    #     ("GIAB12878v1T_no_dedup", "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/bams/GIAB12878v1T.bam.all.depth"),
    # ]
    # output_dir = "/Users/davidkoetsier/Documents/LowPurityPanel/CoverageAnalysis/coverage_analysis/run6"

    # analysis_todo_config = AnalysisTodoConfig(
    #     baf=True,
    #     exome=True,
    #     fusion=True,
    #     hotspot=True,
    #     msi=True,
    #     pgx=True,
    #     tert=True,
    # )
    analysis_todo_config = AnalysisTodoConfig(
        baf=False,
        exome=True,
        fusion=False,
        hotspot=False,
        msi=False,
        pgx=False,
        tert=False,
    )
    config = Config(
        tuple(sample_with_depth_file_list), panel_file_config, analysis_todo_config, min_coverage, output_dir
    )

    do_analysis(config)


if __name__ == "__main__":
    main()
