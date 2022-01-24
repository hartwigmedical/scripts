import concurrent.futures
import logging
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Set, List, Tuple

from google.cloud import storage

from config import AnalysisTypeConfig, Panel, PanelFileConfig, ProgramConfig
from gcp.base import GCPPath
from gcp.client import GCPClient
from genome import Interval
from output_writer import OutputWriter
from panel_reader import PanelReader
from coverage_info import CoverageInfo, get_coverage_info


def do_analysis(
        program_config: ProgramConfig,
        panel_file_config: PanelFileConfig,
        analysis_type_config: AnalysisTypeConfig,
) -> None:
    program_config.working_dir.mkdir(parents=True, exist_ok=True)

    with ProcessPoolExecutor() as executor:
        future_to_bam = {
            executor.submit(analyze_bam, bam, program_config, panel_file_config, analysis_type_config): bam
            for bam in program_config.bams
        }
        for future in concurrent.futures.as_completed(future_to_bam):
            bam = future_to_bam[future]
            try:
                future.result()
            except Exception as exc:
                logging.info(f"BAM {bam} generated an exception: {exc}")
            else:
                logging.info(f"BAM {bam} handled successfully.")


def analyze_bam(
        bam_path: GCPPath,
        program_config: ProgramConfig,
        panel_file_config: PanelFileConfig,
        analysis_type_config: AnalysisTypeConfig,
) -> None:
    gcp_client = GCPClient(storage.Client())

    bam_file_name = bam_path.relative_path.split("/")[-1]
    sample_name = bam_file_name.replace(".bam", "")
    sample_working_dir = program_config.working_dir / sample_name
    sample_working_dir.mkdir(parents=True, exist_ok=True)
    local_bam_path = sample_working_dir / bam_file_name

    logging.info(f"Downloading bam file for sample: {sample_name}")
    gcp_client.download_file(bam_path, local_bam_path)
    gcp_client.download_file(bam_path.append_suffix(".bai"), local_bam_path.with_suffix(".bai"))

    logging.info("Getting samtools depth file")
    depth_file = get_depth_file(local_bam_path, sample_working_dir, program_config.samtools)

    logging.info("Getting panel from file")
    panel = PanelReader.get_panel(panel_file_config)

    logging.info("Getting interesting coverage intervals")
    coverage_intervals = get_coverage_intervals(analysis_type_config, panel)

    logging.info("Getting coverages")
    coverage_info = get_coverage_info(depth_file, coverage_intervals, program_config.min_coverages)

    logging.info("Writing output files")
    local_output_dir = sample_working_dir / "output"
    local_output_dir.mkdir(parents=True, exist_ok=True)
    write_analyses(
        [(sample_name, coverage_info)],
        analysis_type_config,
        program_config.min_coverages,
        panel,
        local_output_dir,
    )

    logging.info("Uploading output files")
    for local_output_file in local_output_dir.iterdir():
        gcp_client.upload_file(
            local_output_file,
            program_config.output_dir.append_suffix(f"/{sample_name}/{local_output_file.name}"),
        )

    logging.info(f"Cleaning up data for sample: {sample_name}")
    shutil.rmtree(sample_working_dir)


def get_depth_file(bam: Path, working_dir: Path, samtools: Path) -> Path:
    depth_file = working_dir / f"{bam.name}.depth"
    if not depth_file.exists():
        create_depth_file(samtools, bam, depth_file)
    if not depth_file.exists():
        raise FileNotFoundError(f"Depth file creation failed: {depth_file}")
    return depth_file


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
        min_coverages: Tuple[int, ...],
        panel: Panel,
        output_dir: Path,
) -> None:
    if analysis_type_config.baf:
        OutputWriter.write_baf_coverage_file(
            sample_with_coverage_info_list, panel.baf_sites, output_dir)
    if analysis_type_config.exome:
        OutputWriter.write_exome_cumulative_coverage_analysis(
            sample_with_coverage_info_list, panel.exons, output_dir)
        for min_coverage in min_coverages:
            OutputWriter.write_exome_count_min_coverage_file(
                sample_with_coverage_info_list, panel.exons, min_coverage, output_dir)
    if analysis_type_config.fusion:
        OutputWriter.write_fusion_cumulative_coverage_file(
            sample_with_coverage_info_list, panel.fusion_sites, output_dir)
        for min_coverage in min_coverages:
            OutputWriter.write_fusion_count_min_coverage_file(
                sample_with_coverage_info_list, panel.fusion_sites, min_coverage, output_dir)
    if analysis_type_config.hotspot:
        OutputWriter.write_hotspot_coverage_file(
            sample_with_coverage_info_list, panel.hotspots, output_dir)
    if analysis_type_config.msi:
        OutputWriter.write_msi_cumulative_coverage_file(
            sample_with_coverage_info_list, panel.msi_sites, output_dir)
        for min_coverage in min_coverages:
            OutputWriter.write_msi_count_min_coverage_file(
                sample_with_coverage_info_list, panel.msi_sites, min_coverage, output_dir)
    if analysis_type_config.pgx:
        OutputWriter.write_pgx_cumulative_coverage_file(
            sample_with_coverage_info_list, panel.pgx_sites, output_dir)
        for min_coverage in min_coverages:
            OutputWriter.write_pgx_count_min_coverage_file(
                sample_with_coverage_info_list, panel.pgx_sites, min_coverage, output_dir)
    if analysis_type_config.tert:
        OutputWriter.write_tert_cumulative_coverage_file(
            sample_with_coverage_info_list, panel.tert_site, output_dir)
        for min_coverage in min_coverages:
            OutputWriter.write_tert_count_min_coverage_file(
                sample_with_coverage_info_list, panel.tert_site, min_coverage, output_dir)


def create_depth_file(samtools: Path, bam: Path, depth_file: Path) -> None:
    cli_args = [str(samtools), "depth", "-s", str(bam)]
    with open(depth_file, "w") as depth_f:
        subprocess.run(cli_args, stdout=depth_f)
