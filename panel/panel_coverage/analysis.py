from pathlib import Path
from typing import Set, List, Tuple

from config import Config, AnalysisTodoConfig, Panel
from genome import Interval
from output_writer import OutputWriter
from panel_reader import PanelReader
from util import parallel_get_sample_to_coverage_info, CoverageInfo


def do_analysis(config: Config) -> None:
    panel = PanelReader.get_panel(config.panel_file_config)

    coverage_intervals = get_coverage_intervals(config.analysis_todo_config, panel)

    sample_to_coverage_info = parallel_get_sample_to_coverage_info(
        config.sample_with_depth_file_pairs, coverage_intervals, config.min_coverage)
    sample_with_coverage_info_list = [
        (sample, sample_to_coverage_info[sample]) for sample, _ in config.sample_with_depth_file_pairs
    ]

    write_analyses(
        sample_with_coverage_info_list, config.analysis_todo_config, config.min_coverage, panel, config.output_dir,
    )


def get_coverage_intervals(analysis_todo_config: AnalysisTodoConfig, panel: Panel) -> Set[Interval]:
    coverage_intervals: Set[Interval] = set()
    if analysis_todo_config.baf:
        baf_intervals = {baf_site.get_site_interval() for baf_site in panel.baf_sites}
        coverage_intervals = coverage_intervals.union(baf_intervals)
    if analysis_todo_config.exome:
        exome_intervals = {exon.interval for exon in panel.exons}
        coverage_intervals = coverage_intervals.union(exome_intervals)
    if analysis_todo_config.fusion:
        fusion_intervals = {fusion_site.interval for fusion_site in panel.fusion_sites}
        coverage_intervals = coverage_intervals.union(fusion_intervals)
    if analysis_todo_config.hotspot:
        hotspot_intervals = {hotspot.get_interval() for hotspot in panel.hotspots}
        coverage_intervals = coverage_intervals.union(hotspot_intervals)
    if analysis_todo_config.msi:
        msi_intervals = {msi_site.get_site_interval() for msi_site in panel.msi_sites}
        coverage_intervals = coverage_intervals.union(msi_intervals)
    if analysis_todo_config.pgx:
        pgx_intervals = {pgx_site.interval for pgx_site in panel.pgx_sites}
        coverage_intervals = coverage_intervals.union(pgx_intervals)
    if analysis_todo_config.tert:
        coverage_intervals.add(panel.tert_site)
    return coverage_intervals


def write_analyses(
        sample_with_coverage_info_list: List[Tuple[str, CoverageInfo]],
        analysis_todo_config: AnalysisTodoConfig,
        min_coverage: int,
        panel: Panel,
        output_dir: Path,
) -> None:
    if analysis_todo_config.baf:
        OutputWriter.write_baf_coverage_file(
            sample_with_coverage_info_list, panel.baf_sites, output_dir)
    if analysis_todo_config.exome:
        OutputWriter.write_exome_count_min_coverage_file(
            sample_with_coverage_info_list, panel.exons, min_coverage, output_dir)
        OutputWriter.write_exome_cumulative_coverage_analysis(
            sample_with_coverage_info_list, panel.exons, output_dir)
    if analysis_todo_config.fusion:
        OutputWriter.write_fusion_count_min_coverage_file(
            sample_with_coverage_info_list, panel.fusion_sites, min_coverage, output_dir)
        OutputWriter.write_fusion_cumulative_coverage_file(
            sample_with_coverage_info_list, panel.fusion_sites, output_dir)
    if analysis_todo_config.hotspot:
        OutputWriter.write_hotspot_coverage_file(
            sample_with_coverage_info_list, panel.hotspots, output_dir)
    if analysis_todo_config.msi:
        OutputWriter.write_msi_count_min_coverage_file(
            sample_with_coverage_info_list, panel.msi_sites, min_coverage, output_dir)
        OutputWriter.write_msi_cumulative_coverage_file(
            sample_with_coverage_info_list, panel.msi_sites, output_dir)
    if analysis_todo_config.pgx:
        OutputWriter.write_pgx_count_min_coverage_file(
            sample_with_coverage_info_list, panel.pgx_sites, min_coverage, output_dir)
        OutputWriter.write_pgx_cumulative_coverage_file(
            sample_with_coverage_info_list, panel.pgx_sites, output_dir)
    if analysis_todo_config.tert:
        OutputWriter.write_tert_count_min_coverage_file(
            sample_with_coverage_info_list, panel.tert_site, min_coverage, output_dir)
        OutputWriter.write_tert_cumulative_coverage_file(
            sample_with_coverage_info_list, panel.tert_site, output_dir)
