import subprocess
from pathlib import Path
from typing import Set, List, Tuple, Dict

from config import AnalysisTypeConfig, Panel, PanelFileConfig, ProgramConfig
from genome import Interval
from output_writer import OutputWriter
from panel_reader import PanelReader
from util import parallel_get_sample_to_coverage_info, CoverageInfo


def do_analysis(
        program_config: ProgramConfig,
        panel_file_config: PanelFileConfig,
        analysis_type_config: AnalysisTypeConfig,
) -> None:
    program_config.working_dir.mkdir(parents=True, exist_ok=True)

    sample_to_depth_file = Dict[str, Path]
    for bam in program_config.bams:
        depth_file = program_config.working_dir / f"{bam.name}.depth"
        if not depth_file.exists():
            create_depth_file(program_config.samtools, bam, depth_file)
        if not depth_file.exists():
            raise FileNotFoundError(f"Depth file creation failed: {depth_file}")
        sample_to_depth_file[bam.stem] = depth_file

    panel = PanelReader.get_panel(panel_file_config)
    coverage_intervals = get_coverage_intervals(analysis_type_config, panel)

    sample_to_coverage_info = parallel_get_sample_to_coverage_info(
        sample_to_depth_file, coverage_intervals, program_config.min_coverages)
    ordered_sample_with_coverage_info_list = [
        (bam.stem, sample_to_coverage_info[bam.stem]) for bam in program_config.bams
    ]
    for min_coverage in program_config.min_coverages:
        output_dir = program_config.output_dir / str(min_coverage)
        write_analyses(
            ordered_sample_with_coverage_info_list, analysis_type_config, min_coverage, panel, output_dir,
        )


def get_coverage_intervals(analysis_type_config: AnalysisTypeConfig, panel: Panel) -> Set[Interval]:
    coverage_intervals: Set[Interval] = set()
    if analysis_type_config.baf:
        baf_intervals = {baf_site.get_site_interval() for baf_site in panel.baf_sites}
        coverage_intervals = coverage_intervals.union(baf_intervals)
    if analysis_type_config.exome:
        exome_intervals = {exon.interval for exon in panel.exons}
        coverage_intervals = coverage_intervals.union(exome_intervals)
    if analysis_type_config.fusion:
        fusion_intervals = {fusion_site.interval for fusion_site in panel.fusion_sites}
        coverage_intervals = coverage_intervals.union(fusion_intervals)
    if analysis_type_config.hotspot:
        hotspot_intervals = {hotspot.get_interval() for hotspot in panel.hotspots}
        coverage_intervals = coverage_intervals.union(hotspot_intervals)
    if analysis_type_config.msi:
        msi_intervals = {msi_site.get_site_interval() for msi_site in panel.msi_sites}
        coverage_intervals = coverage_intervals.union(msi_intervals)
    if analysis_type_config.pgx:
        pgx_intervals = {pgx_site.interval for pgx_site in panel.pgx_sites}
        coverage_intervals = coverage_intervals.union(pgx_intervals)
    if analysis_type_config.tert:
        coverage_intervals.add(panel.tert_site)
    return coverage_intervals


def write_analyses(
        sample_with_coverage_info_list: List[Tuple[str, CoverageInfo]],
        analysis_type_config: AnalysisTypeConfig,
        min_coverage: int,
        panel: Panel,
        output_dir: Path,
) -> None:
    if analysis_type_config.baf:
        OutputWriter.write_baf_coverage_file(
            sample_with_coverage_info_list, panel.baf_sites, output_dir)
    if analysis_type_config.exome:
        OutputWriter.write_exome_count_min_coverage_file(
            sample_with_coverage_info_list, panel.exons, min_coverage, output_dir)
        OutputWriter.write_exome_cumulative_coverage_analysis(
            sample_with_coverage_info_list, panel.exons, output_dir)
    if analysis_type_config.fusion:
        OutputWriter.write_fusion_count_min_coverage_file(
            sample_with_coverage_info_list, panel.fusion_sites, min_coverage, output_dir)
        OutputWriter.write_fusion_cumulative_coverage_file(
            sample_with_coverage_info_list, panel.fusion_sites, output_dir)
    if analysis_type_config.hotspot:
        OutputWriter.write_hotspot_coverage_file(
            sample_with_coverage_info_list, panel.hotspots, output_dir)
    if analysis_type_config.msi:
        OutputWriter.write_msi_count_min_coverage_file(
            sample_with_coverage_info_list, panel.msi_sites, min_coverage, output_dir)
        OutputWriter.write_msi_cumulative_coverage_file(
            sample_with_coverage_info_list, panel.msi_sites, output_dir)
    if analysis_type_config.pgx:
        OutputWriter.write_pgx_count_min_coverage_file(
            sample_with_coverage_info_list, panel.pgx_sites, min_coverage, output_dir)
        OutputWriter.write_pgx_cumulative_coverage_file(
            sample_with_coverage_info_list, panel.pgx_sites, output_dir)
    if analysis_type_config.tert:
        OutputWriter.write_tert_count_min_coverage_file(
            sample_with_coverage_info_list, panel.tert_site, min_coverage, output_dir)
        OutputWriter.write_tert_cumulative_coverage_file(
            sample_with_coverage_info_list, panel.tert_site, output_dir)


def create_depth_file(samtools: Path, bam: Path, depth_file: Path) -> None:
    cli_args = [samtools, "depth", "-s", bam]
    with open(depth_file, "w") as depth_f:
        subprocess.run(cli_args, stdout=depth_f)
