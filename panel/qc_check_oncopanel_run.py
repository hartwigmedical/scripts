#!/usr/bin/env python3
import argparse
import gzip
import json
import logging
import subprocess
from dataclasses import dataclass
from decimal import Decimal, ROUND_DOWN
from enum import Enum, auto

from pathlib import Path
from typing import List, Optional, Tuple, Dict

RB1_GENE_NAME = "RB1"

GCP_URL_START = "gs:"

TSV_SEPARATOR = "\t"
DOMAIN_SEPARATOR = ";"

AMP_MANUAL_INTERPRETATION_THRESHOLD = 7
DEL_MIN_PURITY_THRESHOLD = Decimal("0.30")
AMP_MIN_PURITY_THRESHOLD = Decimal("0.20")
TMB_MIN_PURITY_THRESHOLD = Decimal("0.10")
MSI_MIN_PURITY_THRESHOLD = Decimal("0.10")
SOMATIC_VARIANT_MIN_PURITY_THRESHOLD = Decimal("0.10")

DELETION_COPY_NUMBER_THRESHOLD = Decimal("0.5")
RECURRENT_PARTIAL_DEL_GENES = {"BRCA2", "PTEN", "RASA1", "RB1", "EYS"}
GENES_EXCLUDED_FROM_DESIGN = {"SPATA31A7", "LINC01001", "U2AF1"}

BASES_PER_GBASE = 1000000000
TARGET_YIELD_IN_BASES = 50 * BASES_PER_GBASE

@dataclass(frozen=True, eq=True)
class Metadata:
    tumor_name: str
    tumor_barcode: str
    pathology_id: Optional[str]
    pipeline_version: str
    gcp_run_url: str
    yield_in_bases: int


class DriverType(Enum):
    AMP = auto()
    DEL = auto()
    MUTATION = auto()
    PARTIAL_AMP = auto()

    def __str__(self) -> str:
        return self.name


class GeneType(Enum):
    ONCO = auto()
    TSG = auto()

    def __str__(self) -> str:
        return self.name


class LikelihoodMethod(Enum):
    AMP = auto()
    BIALLELIC = auto()
    DEL = auto()
    DNDS = auto()
    HOTSPOT = auto()
    INFRAME = auto()

    def __str__(self) -> str:
        return self.name


@dataclass(frozen=True, eq=True)
class Driver:
    chromosome: str
    chromosome_band: str
    gene_name: str
    transcript: str
    is_canonical: bool
    driver_type: DriverType
    gene_type: GeneType
    likelihood_method: LikelihoodMethod
    driver_likelihood: Decimal
    missense: int
    nonsense: int
    splice: int
    inframe: int
    frameshift: int
    biallelic: bool
    min_copy_number: Decimal
    max_copy_number: Decimal

    def __str__(self) -> str:
        data = [
            self.chromosome,
            self.chromosome_band,
            self.gene_name,
            self.transcript,
            str(self.is_canonical),
            str(self.driver_type),
            str(self.gene_type),
            str(self.likelihood_method),
            str(self.driver_likelihood),
            str(self.missense),
            str(self.nonsense),
            str(self.splice),
            str(self.inframe),
            str(self.frameshift),
            str(self.biallelic),
            str(self.min_copy_number),
            str(self.max_copy_number),
        ]
        return "\t".join(data)


@dataclass(frozen=True, eq=True)
class PurpleQc:
    qc_status: str
    copy_number_segments: int
    unsupported_copy_number_segments: int
    purity: Decimal
    ploidy: Decimal
    amber_gender: str
    cobalt_gender: str
    deleted_genes: int
    amber_mean_depth: int
    ms_indels_per_mb: Decimal
    msi_status: str
    tmb_per_mb: Decimal
    tmb_status: str


@dataclass(frozen=True, eq=True)
class GeneCopyNumber:
    chromosome: str
    start: int
    end: int
    gene_name: str
    min_copy_number: Decimal
    max_copy_number: Decimal
    transcript: str
    is_canonical: bool
    depth_window_count: int
    min_region_start: int
    min_region_end: int


@dataclass(frozen=True, eq=True)
class SomaticCopyNumber:
    chromosome: str
    start: int
    end: int
    copy_number: Decimal
    depth_window_count: int


@dataclass(frozen=True, eq=True)
class ExonMedianCoverage:
    gene_name: str
    chromosome: str
    start: int
    end: int
    exon_rank: int
    median_depth: int


@dataclass(frozen=True, eq=True)
class WgsMetrics:
    mean_coverage: Decimal
    median_coverage: Decimal
    percent_excluded_adapter: Decimal
    percent_excluded_mapq: Decimal
    percent_excluded_duplicated: Decimal
    percent_excluded_unpaired: Decimal
    percent_excluded_baseq: Decimal
    percent_excluded_overlap: Decimal
    percent_excluded_capped: Decimal
    percent_excluded_total: Decimal
    percent_100x: Decimal


@dataclass(frozen=True, eq=True)
class Variant:
    chromosome: str
    position: int
    ref: str
    alt: str

    def __str__(self) -> str:
        return f"{self.chromosome}:{self.position}{self.ref}>{self.alt}"


@dataclass(frozen=True, eq=True)
class AnnotatedVariant:
    variant: Variant
    total_depth: int
    alt_read_count: int


@dataclass(frozen=True, eq=True)
class Fusion:
    name: str
    reported: bool
    reported_type: str
    phased: str
    driver_likelihood: str
    chain_links: int
    chain_terminated: bool
    domains_kept: Tuple[str, ...]
    domains_lost: Tuple[str, ...]
    fused_exon_up: str
    fused_exon_down: str
    gene_start: str
    gene_end: str
    junction_copy_number: Decimal


@dataclass(frozen=True, eq=True)
class RunData:
    drivers: Tuple[Driver, ...]
    purple_qc: PurpleQc
    gene_copy_numbers: Tuple[GeneCopyNumber, ...]
    somatic_copy_numbers: Tuple[SomaticCopyNumber, ...]
    exon_median_coverages: Tuple[ExonMedianCoverage, ...]
    sage_unfiltered_variants: Tuple[AnnotatedVariant, ...]
    # wgs_metrics: WgsMetrics
    reportable_fusions: Tuple[Fusion, ...]
    driver_catalog_header: str


@dataclass(frozen=True, eq=True)
class DriverGenePanelEntry:
    gene: str


@dataclass(frozen=True, eq=True)
class ExonExclusion:
    gene_name: str
    gene_id: str
    canonical_transcript_id: str
    excluded_exons: Tuple[int, ...]


ID_TO_RESISTANCE_VARIANT = {
    "EGFR T790M": Variant("chr7", 55181378, "C", "T"),
    "EGFR C797S": Variant("chr7", 55181398, "T", "A"),
}


def main(gcp_run_url: str, working_directory: Path, driver_gene_panel: Path, excluded_exons: Path, output_file: Optional[Path]) -> None:
    logging.info("Start QC checking OncoPanel run")

    logging.info("Perform sanity checks")
    if gcp_run_url[:len(GCP_URL_START)] != GCP_URL_START:
        raise ValueError(f"GCP run URL needs to start with '{GCP_URL_START}'")
    while gcp_run_url and gcp_run_url[-1] == "/":
        gcp_run_url = gcp_run_url[:-1]

    logging.info("Get metadata")
    metadata = get_metadata(gcp_run_url)

    logging.info(f"Handling run for {metadata.tumor_name} ({metadata.pathology_id}): {gcp_run_url}")

    if not working_directory.exists():
        logging.info(f"Create working directory: {working_directory}")
        working_directory.mkdir(parents=True)
    else:
        logging.info(f"Working directory exists: {working_directory}")

    local_directory = working_directory / metadata.tumor_name
    if local_directory.is_dir():
        raise ValueError(f"Local directory {local_directory} already exists!")

    copy_run_files_to_local(gcp_run_url, local_directory, metadata.tumor_name)

    run_data = load_run_data(local_directory, metadata.tumor_name)

    logging.info("Load driver gene panel")
    driver_gene_panel_entries = load_driver_gene_panel(driver_gene_panel)

    logging.info("Load excluded_exons")
    excluded_exons = load_excluded_exons(excluded_exons)

    logging.info(f"Do QC checks for {metadata.tumor_name}")
    if output_file is None:
        output_file = local_directory / f"{metadata.tumor_name}.qc_check.txt"

    do_qc_checks(run_data, driver_gene_panel_entries, excluded_exons, metadata, output_file)

    logging.info("Finished QC checking OncoPanel run")


def do_qc_checks(
        run_data: RunData,
        driver_gene_panel_entries: List[DriverGenePanelEntry],
        excluded_exons: List[ExonExclusion],
        metadata: Metadata,
        output_file: Path
) -> None:
    output_text = get_qc_check_output_text(run_data, driver_gene_panel_entries, excluded_exons, metadata)

    with open(output_file, "w") as out_f:
        out_f.write(output_text)


def get_qc_check_output_text(
        run_data: RunData, driver_gene_panel_entries: List[DriverGenePanelEntry], excluded_exons: List[ExonExclusion], metadata: Metadata,
) -> str:
    sections = [
        get_sample_info_text(metadata),
        get_sample_overview_stats_text(run_data.purple_qc, run_data.gene_copy_numbers, run_data.exon_median_coverages, driver_gene_panel_entries, excluded_exons),
        # get_percent_excluded_text(run_data.wgs_metrics),
        get_qc_warning_text(run_data, metadata),
        get_resistance_text(run_data.sage_unfiltered_variants, metadata),
        get_driver_catalog_text(run_data.drivers, run_data.driver_catalog_header),
        get_fusion_text(run_data.reportable_fusions),
    ]
    return "\n\n".join([section for section in sections if section]) + "\n"


def get_sample_info_text(metadata: Metadata) -> str:
    lines = [
        f"Pathology ID: {metadata.pathology_id}",
        f"Internal sample name: {metadata.tumor_name}",
        f"Pipeline version: {metadata.pipeline_version}",
        f"Yield: {convert_bases_to_rounded_gbases(metadata.yield_in_bases)}Gb"
    ]
    return "\n".join(lines)


def get_sample_overview_stats_text(
        purple_qc: PurpleQc,
        gene_copy_numbers: Tuple[GeneCopyNumber, ...],
        exon_median_coverages: Tuple[ExonMedianCoverage, ...],
        driver_gene_panel_entries: List[DriverGenePanelEntry],
        excluded_exons: List[ExonExclusion],
) -> str:
    driver_panel_genes = {entry.gene for entry in driver_gene_panel_entries}
    min_copy_numbers_in_driver_genes = [
        gene_cn.min_copy_number for gene_cn in gene_copy_numbers
        if gene_cn.gene_name in driver_panel_genes and gene_cn.is_canonical
    ]
    negative_min_copy_numbers_count = len(
        [
            min_copy_number for min_copy_number in min_copy_numbers_in_driver_genes if min_copy_number < 0
        ]
    )

    gene_to_excluded_exons = {exon_exclusion.gene_name: exon_exclusion.excluded_exons for exon_exclusion in excluded_exons}
    relevant_exons = [exon for exon in exon_median_coverages if is_relevant_for_coverage(exon, gene_to_excluded_exons)]
    if not relevant_exons:
        raise ValueError("No exons found relevant for coverage. Check excluded genes and exons!")
    exons_with_median_coverage_at_least_100 = [exon for exon in relevant_exons if exon.median_depth >= 100]
    percent_exons_median_coverage_at_least_100 = 100 * len(exons_with_median_coverage_at_least_100) / len(relevant_exons)

    lines = [
        f"Purple QC status: {purple_qc.qc_status}",
        f"Purity: {purple_qc.purity}",
        f"Ploidy: {purple_qc.ploidy}",
        f"MSI status: {purple_qc.msi_status}",
        f"Ms indels per MB: {purple_qc.ms_indels_per_mb}",
        f"TMB status: {purple_qc.tmb_status}",
        f"TMB per MB: {purple_qc.tmb_per_mb}",
        f"Amber gender: {purple_qc.amber_gender}",
        f"Cobalt gender: {purple_qc.cobalt_gender}",
        f"AMBER mean depth: {purple_qc.amber_mean_depth}",
        f"Exons where median coverage is checked: {len(relevant_exons)}",
        f"Exons with median coverage >= 100: {len(exons_with_median_coverage_at_least_100)}",
        f"Percent exons with median coverage >= 100: {percent_exons_median_coverage_at_least_100:.2f}%",
        # f"WgsMetrics mean depth: {wgs_metrics.mean_coverage}",
        # f"WgsMetrics median depth: {wgs_metrics.median_coverage}",
        # f"WgsMetrics percent 100X: {wgs_metrics.percent_100x}",
        f"Driver panel genes with negative minCN: {negative_min_copy_numbers_count} of {len(driver_gene_panel_entries)}",
        f"Lowest minCN in driver panel genes: {min(min_copy_numbers_in_driver_genes)}",
    ]
    return "\n".join(lines)


def is_relevant_for_coverage(exon: ExonMedianCoverage, gene_to_excluded_exons: Dict[str, Tuple[int, ...]]) -> bool:
    gene_explicitly_excluded_from_design = exon.gene_name in GENES_EXCLUDED_FROM_DESIGN
    gene_considered_for_excluded_exons = exon.gene_name in gene_to_excluded_exons.keys()
    exon_excluded = gene_considered_for_excluded_exons and exon.exon_rank in gene_to_excluded_exons[exon.gene_name]
    return not gene_explicitly_excluded_from_design and gene_considered_for_excluded_exons and not exon_excluded


# def get_percent_excluded_text(wgs_metrics: WgsMetrics) -> str:
#     lines = [
#         f"Percent excluded adapter: {wgs_metrics.percent_excluded_adapter}",
#         f"Percent excluded MAPQ: {wgs_metrics.percent_excluded_mapq}",
#         f"Percent excluded duplicated: {wgs_metrics.percent_excluded_duplicated}",
#         f"Percent excluded unpaired: {wgs_metrics.percent_excluded_unpaired}",
#         f"Percent excluded baseq: {wgs_metrics.percent_excluded_baseq}",
#         f"Percent excluded overlap: {wgs_metrics.percent_excluded_overlap}",
#         f"Percent excluded capped: {wgs_metrics.percent_excluded_capped}",
#         f"Percent excluded total: {wgs_metrics.percent_excluded_total}",
#     ]
#     return "\n".join(lines)


def get_qc_warning_text(run_data: RunData, metadata: Metadata) -> str:
    lines = []
    if metadata.yield_in_bases < TARGET_YIELD_IN_BASES:
        lines.append(f"WARN: Yield {convert_bases_to_rounded_gbases(metadata.yield_in_bases)}Gb less than target yield {convert_bases_to_rounded_gbases(TARGET_YIELD_IN_BASES)}Gb")

    if run_data.purple_qc.qc_status != "PASS":
        lines.append(f"WARN: Non-PASS Purple QC status: {run_data.purple_qc.qc_status}")

    if run_data.purple_qc.purity < TMB_MIN_PURITY_THRESHOLD:
        lines.append(
            f"WARN: TMB unreliable due to purity {run_data.purple_qc.purity} < {TMB_MIN_PURITY_THRESHOLD}"
        )

    if run_data.purple_qc.purity < MSI_MIN_PURITY_THRESHOLD:
        lines.append(
            f"WARN: MSI unreliable due to purity {run_data.purple_qc.purity} < {MSI_MIN_PURITY_THRESHOLD}"
        )

    if run_data.purple_qc.purity < SOMATIC_VARIANT_MIN_PURITY_THRESHOLD:
        mutation_drivers = [driver for driver in run_data.drivers if driver.driver_type == DriverType.MUTATION]
        if mutation_drivers:
            gene_name_string = ", ".join([driver.gene_name for driver in mutation_drivers])
            lines.append(
                f"WARN: MUTATION drivers unreliable due to purity {run_data.purple_qc.purity} < {SOMATIC_VARIANT_MIN_PURITY_THRESHOLD}: {gene_name_string}"
            )

    if run_data.purple_qc.purity < DEL_MIN_PURITY_THRESHOLD:
        del_drivers = [driver for driver in run_data.drivers if driver.driver_type == DriverType.DEL]
        if del_drivers:
            gene_name_string = ", ".join([driver.gene_name for driver in del_drivers])
            lines.append(
                f"WARN: DEL drivers unreliable due to purity {run_data.purple_qc.purity} < {DEL_MIN_PURITY_THRESHOLD}: {gene_name_string}"
            )

    if run_data.purple_qc.purity < AMP_MIN_PURITY_THRESHOLD:
        amp_drivers = [driver for driver in run_data.drivers if driver.driver_type == DriverType.AMP]
        if amp_drivers:
            gene_name_string = ", ".join([driver.gene_name for driver in amp_drivers])
            lines.append(
                f"WARN: AMP drivers unreliable due to purity {run_data.purple_qc.purity} < {AMP_MIN_PURITY_THRESHOLD}: {gene_name_string}"
            )

    partial_amp_drivers = [driver for driver in run_data.drivers if driver.driver_type == DriverType.PARTIAL_AMP]
    if partial_amp_drivers:
        gene_name_string = ", ".join([driver.gene_name for driver in partial_amp_drivers])
        lines.append(f"WARN: Partial AMP drivers called, but haven't been validated: {gene_name_string}")

    amp_drivers_needing_manual_curation = [
        driver for driver in run_data.drivers
        if (driver.driver_type == DriverType.AMP) and
           driver.min_copy_number <= AMP_MANUAL_INTERPRETATION_THRESHOLD
    ]
    if amp_drivers_needing_manual_curation:
        gene_name_string = ", ".join([driver.gene_name for driver in amp_drivers_needing_manual_curation])
        lines.append(f"WARN: AMP drivers need manual curation: {gene_name_string}")

    rb1_copy_numbers = [copy_number for copy_number in run_data.gene_copy_numbers if copy_number.gene_name == RB1_GENE_NAME]
    if len(rb1_copy_numbers) != 1:
        raise ValueError(f"Did not find precisely one RB1 gene copy number for sample")
    rb1_copy_number = rb1_copy_numbers[0]
    rb1_del_driver_calls = [
        driver for driver in run_data.drivers if driver.gene_name == RB1_GENE_NAME and driver.driver_type == DriverType.DEL
    ]
    if rb1_copy_number.min_copy_number < DELETION_COPY_NUMBER_THRESHOLD and not rb1_del_driver_calls:
        min_region_length_string = get_pretty_base_count(rb1_copy_number.min_region_end - rb1_copy_number.min_region_start + 1)
        lines.append(
            f"WARN: Potential RB1 deletion missed: "
            f"minCN={rb1_copy_number.min_copy_number}, "
            f"maxCN={rb1_copy_number.max_copy_number}, "
            f"depthWindowCount={rb1_copy_number.depth_window_count}, "
            f"minRegionLength={min_region_length_string}"
        )

    dels_in_suspicious_partial_del_gene = [
        driver for driver in run_data.drivers
        if driver.driver_type == DriverType.DEL
           and driver.gene_name in RECURRENT_PARTIAL_DEL_GENES
           and driver.max_copy_number < DELETION_COPY_NUMBER_THRESHOLD
    ]
    if dels_in_suspicious_partial_del_gene:
        gene_name_string = ", ".join([driver.gene_name for driver in dels_in_suspicious_partial_del_gene])
        lines.append(f"WARN: Full deletion called in gene suspicious for partial DELs: {gene_name_string}")

    suspicious_partial_dels = [
        driver for driver in run_data.drivers
        if driver.driver_type == DriverType.DEL
           and driver.gene_name in RECURRENT_PARTIAL_DEL_GENES
           and driver.max_copy_number >= DELETION_COPY_NUMBER_THRESHOLD
    ]
    if suspicious_partial_dels:
        per_driver_lines = []
        for partial_del in suspicious_partial_dels:
            per_driver_lines.append(f"{partial_del.gene_name}:")
            per_driver_lines.append(str(partial_del))

            per_driver_lines.append("\t".join(["gene", "exon_rank", "deletion_status", "median_depth"]))
            relevant_exon_coverages = [
                exon_coverage for exon_coverage in run_data.exon_median_coverages if
                exon_coverage.gene_name == partial_del.gene_name
            ]
            for exon_coverage in sorted(relevant_exon_coverages, key=lambda x: x.exon_rank):
                deletion_status = get_deletion_status(exon_coverage, run_data.somatic_copy_numbers)
                line = "\t".join(
                    [
                        exon_coverage.gene_name,
                        str(exon_coverage.exon_rank),
                        deletion_status,
                        str(exon_coverage.median_depth)
                    ]
                )
                per_driver_lines.append(line)
        driver_string = "\n".join(per_driver_lines)
        lines.append(f"WARN: Partial del called in suspicious gene:\n{driver_string}")

    return "\n".join(lines)


def get_resistance_text(sage_unfiltered_variants: Tuple[AnnotatedVariant, ...], metadata: Metadata) -> str:
    lines = [
        "Resistance checks:",
    ]
    for resistance_variant_id in sorted(ID_TO_RESISTANCE_VARIANT.keys()):
        variant = ID_TO_RESISTANCE_VARIANT[resistance_variant_id]
        matching_variants = [
            annotated_variant for annotated_variant in sage_unfiltered_variants if annotated_variant.variant == variant
        ]
        if len(matching_variants) == 1:
            matching_variant = matching_variants[0]
            lines.append(
                f"{resistance_variant_id} ({variant}): {matching_variant.alt_read_count} of {matching_variant.total_depth} reads"
            )
        elif not matching_variants:
            lines.append(
                f"{resistance_variant_id} ({variant}): ? of ? reads"
            )
        else:
            raise ValueError(f"Variant {variant} is present in SAGE unfiltered VCF more than once")

    lines.append("IGV URLs:")
    lines.append(get_signed_bam_urls(metadata))

    return "\n".join(lines)


def get_signed_bam_urls(metadata: Metadata) -> str:
    gcp_bam_url = f"{metadata.gcp_run_url}/{metadata.tumor_name}/aligner/{metadata.tumor_name}.bam"
    return run_bash_command(["sign_bam_url", gcp_bam_url, "24h"])


def get_driver_catalog_text(drivers: Tuple[Driver, ...], driver_catalog_header: str) -> str:
    lines = ["Driver catalog:", driver_catalog_header]
    for driver in drivers:
        lines.append(str(driver))
    return "\n".join(lines)


def get_fusion_text(reportable_fusions: Tuple[Fusion, ...]) -> str:
    lines = ["Reportable fusions:"]
    header = "\t".join(["Name", "DriverLikelihood", "JunctionCopyNumber", "ReportedType", "Phased", "FusedExonUp", "FusedExonDown"])
    lines.append(header)
    for fusion in reportable_fusions:
        line = "\t".join([
            fusion.name,
            fusion.driver_likelihood,
            str(fusion.junction_copy_number),
            fusion.reported_type,
            fusion.phased,
            str(fusion.fused_exon_up),
            str(fusion.fused_exon_down),
        ])
        lines.append(line)
    return "\n".join(lines)


def get_deletion_status(exon_coverage: ExonMedianCoverage, somatic_copy_numbers: Tuple[SomaticCopyNumber, ...]) -> str:
    overlapping_somatic_copy_number_regions = [
        somatic_copy_number for somatic_copy_number in somatic_copy_numbers
        if exon_coverage.chromosome == somatic_copy_number.chromosome
           and exon_coverage.start <= somatic_copy_number.end
           and somatic_copy_number.start <= exon_coverage.end
    ]
    some_deleted = any(region for region in overlapping_somatic_copy_number_regions if
                       region.copy_number < DELETION_COPY_NUMBER_THRESHOLD)
    some_not_deleted = any(region for region in overlapping_somatic_copy_number_regions if
                           region.copy_number >= DELETION_COPY_NUMBER_THRESHOLD)
    if some_deleted and some_not_deleted:
        deletion_status = "PARTIAL_DELETION"
    elif some_deleted and not some_not_deleted:
        deletion_status = "DELETION"
    elif not some_deleted and some_not_deleted:
        deletion_status = "PRESENT"
    else:
        raise ValueError(
            f"Could not determine deletion status of exon {exon_coverage.gene_name}:{exon_coverage.exon_rank}")
    return deletion_status


def load_run_data(local_directory: Path, tumor_name: str) -> RunData:
    logging.info("Load drivers")
    drivers = load_drivers(local_directory, tumor_name)
    logging.info("Load driver catalog header")
    driver_catalog_header = load_driver_catalog_header(local_directory, tumor_name)
    logging.info("Load sample qc info")
    purple_qc = load_purple_qc(local_directory, tumor_name)
    logging.info("Load gene copy numbers")
    gene_copy_numbers = load_gene_copy_numbers(local_directory, tumor_name)
    logging.info("Load somatic copy numbers")
    somatic_copy_numbers = load_somatic_copy_numbers(local_directory, tumor_name)
    logging.info("Load exon median coverages")
    exon_median_coverages = load_exon_median_coverages(local_directory, tumor_name)
    logging.info("Load SAGE unfiltered variants")
    sage_unfiltered_variants = load_sage_unfiltered_variants(local_directory, tumor_name)
    # logging.info("Load wgs metrics")
    # wgs_metrics = load_wgs_metrics(local_directory, tumor_name)
    logging.info("Load reported fusions")
    reported_fusions = load_reported_fusions(local_directory, tumor_name)
    run_data = RunData(
        tuple(drivers),
        purple_qc,
        tuple(gene_copy_numbers),
        tuple(somatic_copy_numbers),
        tuple(exon_median_coverages),
        tuple(sage_unfiltered_variants),
        # wgs_metrics,
        tuple(reported_fusions),
        driver_catalog_header,
    )
    return run_data


def load_driver_gene_panel(driver_gene_panel_path: Path) -> List[DriverGenePanelEntry]:
    driver_gene_panel_entries = []
    with open(driver_gene_panel_path, "r") as in_f:
        header = next(in_f)
        for line in in_f:
            entry = get_driver_gene_panel_entry_from_line(line, header)
            driver_gene_panel_entries.append(entry)
    return driver_gene_panel_entries


def get_driver_gene_panel_entry_from_line(line: str, header: str) -> DriverGenePanelEntry:
    split_header = header.replace("\n", "").split(TSV_SEPARATOR)

    split_line = line.replace("\n", "").split(TSV_SEPARATOR)

    gene = split_line[split_header.index("gene")]

    return DriverGenePanelEntry(gene)


def load_excluded_exons(excluded_exons_path: Path) -> List[ExonExclusion]:
    excluded_exons = []
    with open(excluded_exons_path, "r") as in_f:
        header = next(in_f)
        # skip subheader
        next(in_f)
        for line in in_f:
            entry = get_excluded_exons_entry_from_line(line, header)
            excluded_exons.append(entry)
    return excluded_exons


def get_excluded_exons_entry_from_line(line: str, header: str) -> ExonExclusion:
    split_header = header.replace("\n", "").split(TSV_SEPARATOR)
    split_line = line.replace("\n", "").split(TSV_SEPARATOR)
    gene_name = split_line[split_header.index("Gene (symbol)")]
    gene_id = split_line[split_header.index("Ensembl gene ID")]
    canonical_transcript_id = split_line[split_header.index("Ensembl Transcript ID (canonical)")]
    excluded_exons_string = split_line[split_header.index("Trancript ID exons not included")]
    excluded_exons = tuple(int(rank) for rank in excluded_exons_string.split(",") if rank)
    return ExonExclusion(gene_name, gene_id, canonical_transcript_id, excluded_exons)


def load_reported_fusions(local_directory: Path, tumor_name: str) -> List[Fusion]:
    fusion_file_path = local_directory / f"{tumor_name}.linx.fusion.tsv"
    return load_reported_fusions_from_file(fusion_file_path)


def load_reported_fusions_from_file(fusion_file_path: Path) -> List[Fusion]:
    reported_fusions = []
    with open(fusion_file_path, "r") as in_f:
        header = next(in_f)
        for line in in_f:
            fusion = get_fusion_from_line(line, header)
            if fusion.reported:
                reported_fusions.append(fusion)
    return reported_fusions


def get_fusion_from_line(line: str, header: str) -> Fusion:
    split_header = header.replace("\n", "").split(TSV_SEPARATOR)

    split_line = line.replace("\n", "").split(TSV_SEPARATOR)

    name = split_line[split_header.index("name")]
    reported = get_bool_from_string(split_line[split_header.index("reported")])
    reported_type = split_line[split_header.index("reportedType")]
    phased = split_line[split_header.index("phased")]
    driver_likelihood = split_line[split_header.index("likelihood")]
    chain_links = int(split_line[split_header.index("chainLinks")])
    chain_terminated = get_bool_from_string(split_line[split_header.index("chainTerminated")])
    domains_kept = tuple(split_line[split_header.index("domainsKept")].split(DOMAIN_SEPARATOR))
    domains_lost = tuple(split_line[split_header.index("domainsLost")].split(DOMAIN_SEPARATOR))
    fused_exon_up = split_line[split_header.index("fusedExonUp")]
    fused_exon_down = split_line[split_header.index("fusedExonDown")]
    gene_start = split_line[split_header.index("geneStart")]
    gene_end = split_line[split_header.index("geneEnd")]
    junction_copy_number = Decimal(split_line[split_header.index("junctionCopyNumber")])

    fusion = Fusion(
        name,
        reported,
        reported_type,
        phased,
        driver_likelihood,
        chain_links,
        chain_terminated,
        domains_kept,
        domains_lost,
        fused_exon_up,
        fused_exon_down,
        gene_start,
        gene_end,
        junction_copy_number,
    )
    return fusion


# def load_wgs_metrics(local_directory: Path, tumor_name: str) -> WgsMetrics:
#     wgs_metrics_path = local_directory / f"{tumor_name}.wgsmetrics"
#     return load_wgs_metrics_from_file(wgs_metrics_path)
#
#
# def load_wgs_metrics_from_file(wgs_metrics_path: Path) -> WgsMetrics:
#     with open(wgs_metrics_path, "r") as in_f:
#         # skip first line
#         next(in_f)
#         header = next(in_f)
#         data = next(in_f)
#     return get_wgs_metrics_from_line(data, header)
#
#
# def get_wgs_metrics_from_line(line: str, header: str) -> WgsMetrics:
#     split_header = header.replace("\n", "").split(TSV_SEPARATOR)
#
#     split_line = line.replace("\n", "").split(TSV_SEPARATOR)
#
#     mean_coverage = Decimal(split_line[split_header.index("MEAN_COVERAGE")])
#     median_coverage = Decimal(split_line[split_header.index("MEDIAN_COVERAGE")])
#     percent_excluded_adapter = Decimal(split_line[split_header.index("PCT_EXC_ADAPTER")])
#     percent_excluded_mapq = Decimal(split_line[split_header.index("PCT_EXC_MAPQ")])
#     percent_excluded_duplicated = Decimal(split_line[split_header.index("PCT_EXC_DUPE")])
#     percent_excluded_unpaired = Decimal(split_line[split_header.index("PCT_EXC_UNPAIRED")])
#     percent_excluded_baseq = Decimal(split_line[split_header.index("PCT_EXC_BASEQ")])
#     percent_excluded_overlap = Decimal(split_line[split_header.index("PCT_EXC_OVERLAP")])
#     percent_excluded_capped = Decimal(split_line[split_header.index("PCT_EXC_CAPPED")])
#     percent_excluded_total = Decimal(split_line[split_header.index("PCT_EXC_TOTAL")])
#     percent_100x = Decimal(split_line[split_header.index("PCT_100X")])
#
#     wgs_metrics = WgsMetrics(
#         mean_coverage,
#         median_coverage,
#         percent_excluded_adapter,
#         percent_excluded_mapq,
#         percent_excluded_duplicated,
#         percent_excluded_unpaired,
#         percent_excluded_baseq,
#         percent_excluded_overlap,
#         percent_excluded_capped,
#         percent_excluded_total,
#         percent_100x,
#     )
#     return wgs_metrics


def load_exon_median_coverages(local_directory: Path, tumor_name: str) -> List[ExonMedianCoverage]:
    sage_exon_medians_path = local_directory / f"{tumor_name}.sage.exon.medians.tsv"
    return load_exon_median_coverages_from_file(sage_exon_medians_path)


def load_exon_median_coverages_from_file(sage_exon_medians_path: Path) -> List[ExonMedianCoverage]:
    sage_exon_medians = []
    with open(sage_exon_medians_path, "r") as in_f:
        header = next(in_f)
        for line in in_f:
            sage_exon_median = get_sage_exon_median_from_line(line, header)
            sage_exon_medians.append(sage_exon_median)
    return sage_exon_medians


def get_sage_exon_median_from_line(line: str, header: str) -> ExonMedianCoverage:
    split_header = header.replace("\n", "").split(TSV_SEPARATOR)

    split_line = line.replace("\n", "").split(TSV_SEPARATOR)

    gene_name = split_line[split_header.index("gene")]
    chromosome = split_line[split_header.index("chromosome")]
    start = int(split_line[split_header.index("posStart")])
    end = int(split_line[split_header.index("posEnd")])
    exon_rank = int(split_line[split_header.index("exonRank")])
    median_depth = int(split_line[split_header.index("medianDepth")])

    exon_median_coverage = ExonMedianCoverage(
        gene_name,
        chromosome,
        start,
        end,
        exon_rank,
        median_depth,
    )
    return exon_median_coverage


def load_somatic_copy_numbers(local_directory: Path, tumor_name: str) -> List[SomaticCopyNumber]:
    purple_cnv_somatic_path = local_directory / f"{tumor_name}.purple.cnv.somatic.tsv"
    return load_somatic_copy_numbers_from_file(purple_cnv_somatic_path)


def load_somatic_copy_numbers_from_file(purple_cnv_somatic_path: Path) -> List[SomaticCopyNumber]:
    somatic_copy_numbers = []
    with open(purple_cnv_somatic_path, "r") as in_f:
        header = next(in_f)
        for line in in_f:
            somatic_copy_number = get_somatic_copy_number_from_line(line, header)
            somatic_copy_numbers.append(somatic_copy_number)
    return somatic_copy_numbers


def get_somatic_copy_number_from_line(line: str, header: str) -> SomaticCopyNumber:
    split_header = header.replace("\n", "").split(TSV_SEPARATOR)

    split_line = line.replace("\n", "").split(TSV_SEPARATOR)

    chromosome = split_line[split_header.index("chromosome")]
    start = int(split_line[split_header.index("start")])
    end = int(split_line[split_header.index("end")])
    copy_number = Decimal(split_line[split_header.index("copyNumber")])
    depth_window_count = int(split_line[split_header.index("depthWindowCount")])

    somatic_copy_number = SomaticCopyNumber(
        chromosome,
        start,
        end,
        copy_number,
        depth_window_count,
    )
    return somatic_copy_number


def load_gene_copy_numbers(local_directory: Path, tumor_name: str) -> List[GeneCopyNumber]:
    purple_cnv_gene_path = local_directory / f"{tumor_name}.purple.cnv.gene.tsv"
    return load_gene_copy_numbers_from_file(purple_cnv_gene_path)


def load_gene_copy_numbers_from_file(purple_cnv_gene_path: Path) -> List[GeneCopyNumber]:
    gene_copy_numbers = []
    with open(purple_cnv_gene_path, "r") as in_f:
        header = next(in_f)
        for line in in_f:
            gene_copy_number = get_gene_copy_number_from_line(line, header)
            gene_copy_numbers.append(gene_copy_number)
    return gene_copy_numbers


def get_gene_copy_number_from_line(line: str, header: str) -> GeneCopyNumber:
    split_header = header.replace("\n", "").split(TSV_SEPARATOR)

    split_line = line.replace("\n", "").split(TSV_SEPARATOR)

    chromosome = split_line[split_header.index("chromosome")]
    gene_name = split_line[split_header.index("gene")]
    transcript = split_line[split_header.index("transcriptId")]
    is_canonical = get_bool_from_string(split_line[split_header.index("isCanonical")])
    start = int(split_line[split_header.index("start")])
    end = int(split_line[split_header.index("end")])
    min_copy_number = Decimal(split_line[split_header.index("minCopyNumber")])
    max_copy_number = Decimal(split_line[split_header.index("maxCopyNumber")])
    depth_window_count = int(split_line[split_header.index("depthWindowCount")])
    min_region_start = int(split_line[split_header.index("minRegionStart")])
    min_region_end = int(split_line[split_header.index("minRegionEnd")])

    gene_copy_number = GeneCopyNumber(
        chromosome,
        start,
        end,
        gene_name,
        min_copy_number,
        max_copy_number,
        transcript,
        is_canonical,
        depth_window_count,
        min_region_start,
        min_region_end,
    )
    return gene_copy_number


def load_purple_qc(local_directory: Path, tumor_name: str) -> PurpleQc:
    purple_qc_path = local_directory / f"{tumor_name}.purple.qc"
    purple_purity_path = local_directory / f"{tumor_name}.purple.purity.tsv"
    return load_purple_qc_from_file(purple_qc_path, purple_purity_path)


def load_purple_qc_from_file(purple_qc_path: Path, purple_purity_path: Path) -> PurpleQc:
    field_to_string_value = {}
    with open(purple_qc_path, "r") as in_f:
        for line in in_f:
            split_line = line.replace("\n", "").split(TSV_SEPARATOR)
            field_to_string_value[split_line[0]] = split_line[1]

        qc_status = field_to_string_value["QCStatus"]
        copy_number_segments = int(field_to_string_value["CopyNumberSegments"])
        unsupported_copy_number_segments = int(field_to_string_value["UnsupportedCopyNumberSegments"])
        purity = Decimal(field_to_string_value["Purity"])
        amber_gender = field_to_string_value["AmberGender"]
        cobalt_gender = field_to_string_value["CobaltGender"]
        deleted_genes = int(field_to_string_value["DeletedGenes"])
        amber_mean_depth = int(field_to_string_value["AmberMeanDepth"])

    with open(purple_purity_path, "r") as in_f:
        header = next(in_f)
        line = next(in_f)
        split_header = header.replace("\n", "").split(TSV_SEPARATOR)
        split_line = line.replace("\n", "").split(TSV_SEPARATOR)

        ploidy = Decimal(split_line[split_header.index("ploidy")])
        ms_indels_per_mb = Decimal(split_line[split_header.index("msIndelsPerMb")])
        msi_status = split_line[split_header.index("msStatus")]
        tmb_per_mb = Decimal(split_line[split_header.index("tmbPerMb")])
        tmb_status = split_line[split_header.index("tmbStatus")]

    purple_qc = PurpleQc(
        qc_status,
        copy_number_segments,
        unsupported_copy_number_segments,
        purity,
        ploidy,
        amber_gender,
        cobalt_gender,
        deleted_genes,
        amber_mean_depth,
        ms_indels_per_mb,
        msi_status,
        tmb_per_mb,
        tmb_status
    )
    return purple_qc


def load_drivers(local_directory: Path, tumor_name: str) -> List[Driver]:
    return load_drivers_from_file(get_local_purple_driver_catalog(local_directory, tumor_name))


def load_drivers_from_file(driver_catalog_path: Path) -> List[Driver]:
    drivers = []
    with open(driver_catalog_path, "r") as in_f:
        header = next(in_f)
        for line in in_f:
            driver = get_driver_from_line(line, header)
            drivers.append(driver)
    return drivers


def load_driver_catalog_header(local_directory: Path, tumor_name: str) -> str:
    return load_driver_catalog_header_from_file(get_local_purple_driver_catalog(local_directory, tumor_name))


def load_driver_catalog_header_from_file(driver_catalog_path: Path) -> str:
    with open(driver_catalog_path, "r") as in_f:
        return next(in_f).replace("\n", "")


def get_local_purple_driver_catalog(local_directory, tumor_name):
    local_purple_driver_catalog = local_directory / f"{tumor_name}.purple.driver.catalog.somatic.tsv"
    if not local_purple_driver_catalog.exists():
        local_purple_driver_catalog = local_directory / f"{tumor_name}.driver.catalog.somatic.tsv"
    return local_purple_driver_catalog


def get_driver_from_line(line: str, header: str) -> Driver:
    split_header = header.replace("\n", "").split(TSV_SEPARATOR)
    split_line = line.replace("\n", "").split(TSV_SEPARATOR)

    chromosome = split_line[split_header.index("chromosome")]
    chromosome_band = split_line[split_header.index("chromosomeBand")]
    gene_name = split_line[split_header.index("gene")]
    transcript = split_line[split_header.index("transcript")]
    is_canonical = get_bool_from_string(split_line[split_header.index("isCanonical")])
    driver_type = DriverType[split_line[split_header.index("driver")]]
    gene_type = GeneType[split_line[split_header.index("category")]]
    likelihood_method = LikelihoodMethod[split_line[split_header.index("likelihoodMethod")]]
    driver_likelihood = Decimal(split_line[split_header.index("driverLikelihood")])
    missense = int(split_line[split_header.index("missense")])
    nonsense = int(split_line[split_header.index("nonsense")])
    splice = int(split_line[split_header.index("splice")])
    inframe = int(split_line[split_header.index("inframe")])
    frameshift = int(split_line[split_header.index("frameshift")])
    biallelic = get_bool_from_string(split_line[split_header.index("biallelic")])
    min_copy_number = Decimal(split_line[split_header.index("minCopyNumber")])
    max_copy_number = Decimal(split_line[split_header.index("maxCopyNumber")])

    driver = Driver(
        chromosome,
        chromosome_band,
        gene_name,
        transcript,
        is_canonical,
        driver_type,
        gene_type,
        likelihood_method,
        driver_likelihood,
        missense,
        nonsense,
        splice,
        inframe,
        frameshift,
        biallelic,
        min_copy_number,
        max_copy_number
    )
    return driver


def load_sage_unfiltered_variants(local_directory: Path, tumor_name: str) -> List[AnnotatedVariant]:
    sage_unfiltered_vcf_path = local_directory / f"{tumor_name}.sage.somatic.vcf.gz"
    return load_variants_from_file(sage_unfiltered_vcf_path, tumor_name)


def load_variants_from_file(vcf_path: Path, tumor_name: str) -> List[AnnotatedVariant]:
    variants = []
    sample_index = -1
    allelic_depth_index = -1
    depth_index = -1
    with gzip.open(vcf_path, "rt") as in_f:
        for line in in_f:
            if line[0:2] != "##":
                split_line = line.rstrip("\n").split("\t")
                if line[0] == "#":
                    sample_index = split_line.index(tumor_name)
                else:
                    chromosome = split_line[0]
                    position = int(split_line[1])
                    ref = split_line[3]
                    alt = split_line[4]
                    variant = Variant(chromosome, position, ref, alt)

                    if allelic_depth_index == -1:
                        info_format = split_line[8]
                        allelic_depth_index = info_format.split(":").index("AD")
                    if depth_index == -1:
                        info_format = split_line[8]
                        depth_index = info_format.split(":").index("DP")

                    sample_info = split_line[sample_index]
                    total_depth = int(sample_info.split(":")[depth_index])
                    allelic_depth = sample_info.split(":")[allelic_depth_index]
                    alt_read_count = int(allelic_depth.split(",")[1])

                    variants.append(AnnotatedVariant(variant, total_depth, alt_read_count))
    return variants

def get_bool_from_string(string: str) -> bool:
    if string == "true":
        return True
    elif string == "false":
        return False
    else:
        raise ValueError(f"Unrecognized boolean string: {string}")


def copy_run_files_to_local(gcp_run_url: str, local_directory: Path, tumor_name: str) -> None:
    logging.info("Copy Purple files to local")
    copy_purple_files_to_local(gcp_run_url, local_directory)
    logging.info("Copy SAGE files to local")
    copy_sage_somatic_files_to_local(gcp_run_url, local_directory)
    # logging.info("Copy WgsMetrics file to local")
    # copy_wgs_metrics_file_to_local(gcp_run_url, local_directory, tumor_name)
    logging.info("Copy Linx fusion file to local")
    copy_fusion_file_to_local(gcp_run_url, local_directory, tumor_name)
    logging.info("Copy Orange PDF file to local")
    copy_orange_pdf_file_to_local(gcp_run_url, local_directory, tumor_name)


def copy_purple_files_to_local(gcp_run_url: str, local_directory: Path) -> None:
    purple_url = f"{gcp_run_url}/purple"
    sync_directory_to_local(purple_url, local_directory)


def copy_sage_somatic_files_to_local(gcp_run_url: str, local_directory: Path) -> None:
    purple_url = f"{gcp_run_url}/sage_somatic"
    sync_directory_to_local(purple_url, local_directory)


def copy_wgs_metrics_file_to_local(gcp_run_url: str, local_directory: Path, tumor_name: str) -> None:
    wgs_metrics_url = f"{gcp_run_url}/{tumor_name}/bam_metrics/{tumor_name}.wgsmetrics"
    copy_file_to_local_directory(wgs_metrics_url, local_directory)


def copy_fusion_file_to_local(gcp_run_url: str, local_directory: Path, tumor_name: str) -> None:
    fusion_url = f"{gcp_run_url}/linx/{tumor_name}.linx.fusion.tsv"
    copy_file_to_local_directory(fusion_url, local_directory)


def copy_orange_pdf_file_to_local(gcp_run_url: str, local_directory: Path, tumor_name: str) -> None:
    orange_pdf_url = f"{gcp_run_url}/orange/{tumor_name}.orange.pdf"
    copy_file_to_local_directory(orange_pdf_url, local_directory)


def sync_directory_to_local(source_url: str, local_directory: Path) -> None:
    local_directory.mkdir(parents=True, exist_ok=True)
    subprocess.run(["gsutil", "-m", "rsync", "-r", source_url, local_directory], capture_output=True)


def copy_file_to_local_directory(source_url: str, local_directory: Path) -> None:
    local_directory.mkdir(parents=True, exist_ok=True)
    subprocess.run(["gsutil", "-m", "cp", source_url, local_directory], capture_output=True)


def get_metadata(gcp_run_url: str) -> Metadata:
    metadata_json = json.loads(run_bash_command(["gsutil", "cat", f"{gcp_run_url}/metadata.json"]))
    tumor_name = metadata_json["tumor"]["sampleName"]
    tumor_barcode = metadata_json["tumor"]["barcode"]

    patient_reporter_data_json = json.loads(run_bash_command(["lama_get_patient_reporter_data", tumor_barcode]))
    pathology_id = patient_reporter_data_json["pathologyNumber"] if "pathologyNumber" in patient_reporter_data_json.keys() else None
    set_name = metadata_json["set"]

    pipeline_version = get_pipeline_version_from_api(set_name)
    yield_in_bases = get_yield_in_bases_from_api(tumor_barcode)

    return Metadata(tumor_name, tumor_barcode, pathology_id, pipeline_version, gcp_run_url, yield_in_bases)


def get_pipeline_version_from_api(set_name: str) -> str:
    api_output = run_bash_command(["hmf_api_get", f"runs?set_name={set_name}"])

    if not api_output:
        logging.warning(f"Could not reach API for pipeline version of set {set_name}")
        return "API not found"

    api_json = json.loads(api_output)
    finished_runs = [run_info for run_info in api_json if run_info["status"] == "Finished"]
    if not finished_runs:
        logging.warning(f"Could not find finished pipeline run for set {set_name}")
        return "Run not found in API"

    if len(finished_runs) > 1:
        logging.warning(
            f"Multiple runs found for set {set_name}.\n"
            f"Chose latest finished run in API for pipeline version for '{set_name}'"
        )
    return finished_runs[-1]["version"]


def get_yield_in_bases_from_api(tumor_barcode: str) -> int:
    api_output = run_bash_command(["hmf_api_get", f"samples?barcode={tumor_barcode}"])
    if not api_output:
        logging.warning(f"Could not reach API for sample yield")
        return -BASES_PER_GBASE

    api_json = json.loads(api_output)
    if not api_json:
        logging.warning(f"Could not find sample in API with barcode {tumor_barcode}")
        return -BASES_PER_GBASE

    if len(api_json) > 1:
        logging.warning(
            f"Found multiple samples in API with barcode {tumor_barcode}"
        )
        return -BASES_PER_GBASE

    return api_json[0]["yld"]


def run_bash_command(command: List[str]) -> str:
    return subprocess.run(command, capture_output=True, text=True).stdout


def get_pretty_base_count(count: int) -> str:
    if count >= 1000000:
        return f"{count/1000000:.1f}MB"
    elif count >= 1000:
        return f"{count/1000:.1f}KB"
    else:
        return f"{count}B"


def convert_bases_to_rounded_gbases(base_count: int) -> Decimal:
    return (Decimal(base_count) / Decimal(BASES_PER_GBASE)).quantize(Decimal("0.1"), rounding=ROUND_DOWN)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="qc_check_oncopanel_run.py", description="Perform QC checks for OncoPanel run")
    parser.add_argument("--gcp-run-url", "-u", type=str, required=True, help="GCP path to run directory")
    parser.add_argument("--working-directory", "-w", type=Path, required=True, help="Local working directory")
    parser.add_argument("--driver-gene-panel", "-d", type=Path, required=True, help="Driver gene panel")
    parser.add_argument("--excluded-exons", "-e", type=Path, required=True, help="Exons excluded from coverage analysis")
    parser.add_argument("--output-file", "-o", type=Path, help="Output file with summary of qc check results")
    return parser.parse_args()


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s - [%(levelname)-8s] - %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
    )
    args = parse_args()
    main(args.gcp_run_url, args.working_directory, args.driver_gene_panel, args.excluded_exons, args.output_file)
