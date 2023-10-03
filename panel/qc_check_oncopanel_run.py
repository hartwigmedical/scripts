#!/usr/bin/env python3
import argparse
import json
import logging
import subprocess
from dataclasses import dataclass
from decimal import Decimal
from enum import Enum, auto

from pathlib import Path
from typing import Any, Dict, List, Optional

GCP_URL_START = "gs:"

TSV_SEPARATOR = "\t"

AMP_MANUAL_INTERPRETATION_THRESHOLD = 7
DEL_MIN_PURITY_THRESHOLD = Decimal("0.30")
DELETION_COPY_NUMBER_THRESHOLD = Decimal("0.5")
RECURRENT_PARTIAL_DEL_GENES = {"BRCA2", "PTEN", "RASA1", "RB1"}

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
    amber_gender: str
    cobalt_gender: str
    deleted_genes: int
    amber_mean_depth: int


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
class DriverGenePanelEntry:
    gene: str


def main(gcp_run_url: str, working_directory: Path, driver_gene_panel: Path, output_file: Optional[Path]) -> None:
    logging.info("Start QC checking OncoPanel run")

    logging.info("Perform sanity checks")
    if gcp_run_url[:len(GCP_URL_START)] != GCP_URL_START:
        raise ValueError(f"GCP run URL needs to start with '{GCP_URL_START}'")
    while gcp_run_url and gcp_run_url[-1] == "/":
        gcp_run_url = gcp_run_url[:-1]

    logging.info("Get tumor sample name")
    metadata_json = get_metadata_json(gcp_run_url)
    tumor_name = metadata_json["tumor"]["sampleName"]
    logging.info(f"Handling run for {tumor_name}: {gcp_run_url}")

    if not working_directory.exists():
        logging.info(f"Create working directory: {working_directory}")
        working_directory.mkdir(parents=True)
    else:
        logging.info(f"Working directory exists: {working_directory}")

    logging.info("Copy files to local")
    local_directory = working_directory / tumor_name
    copy_purple_files_to_local(gcp_run_url, local_directory)
    copy_sage_exons_median_tsv_to_local(gcp_run_url, local_directory, tumor_name)
    copy_wgs_metrics_file(gcp_run_url, local_directory, tumor_name)

    logging.info("Load drivers")
    drivers = load_drivers(local_directory, tumor_name)

    logging.info("Load sample qc info")
    purple_qc = load_purple_qc(local_directory, tumor_name)

    logging.info("Load gene copy numbers")
    gene_copy_numbers = load_gene_copy_numbers(local_directory, tumor_name)

    logging.info("Load somatic copy numbers")
    somatic_copy_numbers = load_somatic_copy_numbers(local_directory, tumor_name)

    logging.info("Load exon median coverages")
    exon_median_coverages = load_exon_median_coverages(local_directory, tumor_name)

    logging.info("Load wgs metrics")
    wgs_metrics = load_wgs_metrics(local_directory, tumor_name)

    logging.info("Load driver gene panel")
    driver_gene_panel_entries = load_driver_gene_panel(driver_gene_panel)

    logging.info(f"Do QC checks for {tumor_name}")
    if output_file is None:
        output_file = local_directory / f"{tumor_name}.qc_check.txt"

    do_qc_checks(
        drivers,
        purple_qc,
        gene_copy_numbers,
        somatic_copy_numbers,
        exon_median_coverages,
        wgs_metrics,
        driver_gene_panel_entries,
        tumor_name,
        output_file,
    )

    logging.info("Finished QC checking OncoPanel run")


def do_qc_checks(
        drivers: List[Driver],
        purple_qc: PurpleQc,
        gene_copy_numbers: List[GeneCopyNumber],
        somatic_copy_numbers: List[SomaticCopyNumber],
        exon_median_coverages: List[ExonMedianCoverage],
        wgs_metrics: WgsMetrics,
        driver_gene_panel_entries: List[DriverGenePanelEntry],
        tumor_name: str,
        output_file: Path,
) -> None:
    output_text = get_qc_check_output_text(
        drivers,
        purple_qc,
        gene_copy_numbers,
        somatic_copy_numbers,
        exon_median_coverages,
        wgs_metrics,
        driver_gene_panel_entries,
        tumor_name,
    )

    with open(output_file, "w") as out_f:
        out_f.write(output_text)

def get_qc_check_output_text(
        drivers: List[Driver],
        purple_qc: PurpleQc,
        gene_copy_numbers: List[GeneCopyNumber],
        somatic_copy_numbers: List[SomaticCopyNumber],
        exon_median_coverages: List[ExonMedianCoverage],
        wgs_metrics: WgsMetrics,
        driver_gene_panel_entries: List[DriverGenePanelEntry],
        tumor_name: str,
) -> str:
    lines = []
    lines.append(f"Sample name: {tumor_name}")
    lines.append(f"Amber gender: {purple_qc.amber_gender}")
    lines.append(f"Cobalt gender: {purple_qc.cobalt_gender}")
    lines.append(f"AMBER mean depth: {purple_qc.amber_mean_depth}")
    lines.append(f"WgsMetrics mean depth: {wgs_metrics.mean_coverage}")
    lines.append(f"WgsMetrics median depth: {wgs_metrics.median_coverage}")
    lines.append(f"WgsMetrics percent 100X: {wgs_metrics.percent_100x}")
    lines.append("")
    lines.append(f"Percent excluded adapter: {wgs_metrics.percent_excluded_adapter}")
    lines.append(f"Percent excluded MAPQ: {wgs_metrics.percent_excluded_mapq}")
    lines.append(f"Percent excluded duplicated: {wgs_metrics.percent_excluded_duplicated}")
    lines.append(f"Percent excluded unpaired: {wgs_metrics.percent_excluded_unpaired}")
    lines.append(f"Percent excluded baseq: {wgs_metrics.percent_excluded_baseq}")
    lines.append(f"Percent excluded overlap: {wgs_metrics.percent_excluded_overlap}")
    lines.append(f"Percent excluded capped: {wgs_metrics.percent_excluded_capped}")
    lines.append(f"Percent excluded total: {wgs_metrics.percent_excluded_total}")
    lines.append(f"")

    driver_panel_genes = {entry.gene for entry in driver_gene_panel_entries}
    min_copy_numbers_in_driver_genes = [
        gene_cn.min_copy_number for gene_cn in gene_copy_numbers
        if gene_cn.gene_name in driver_panel_genes and gene_cn.is_canonical
    ]
    negative_min_copy_numbers_in_driver_genes_count = len([
        min_copy_number for min_copy_number in min_copy_numbers_in_driver_genes if min_copy_number < 0
    ])
    lines.append(
        f"Driver panel genes with negative minCN: "
        f"{negative_min_copy_numbers_in_driver_genes_count} of {len(driver_gene_panel_entries)}"
    )
    lines.append(f"Lowest minCN in driver panel genes: {min(min_copy_numbers_in_driver_genes)}")
    if purple_qc.qc_status != "PASS":
        lines.append(f"WARN: Non-PASS Purple QC status: {purple_qc.qc_status}")
    amp_drivers_needing_manual_curation = [
        driver for driver in drivers
        if (driver.driver_type == DriverType.AMP or driver.driver_type == DriverType.PARTIAL_AMP) and
           driver.min_copy_number <= AMP_MANUAL_INTERPRETATION_THRESHOLD
    ]
    if amp_drivers_needing_manual_curation:
        driver_string = "\n".join([str(driver) for driver in amp_drivers_needing_manual_curation])
        lines.append(f"WARN: AMP driver needs manual curation:\n{driver_string}")
    if purple_qc.purity < DEL_MIN_PURITY_THRESHOLD:
        del_drivers = [driver for driver in drivers if driver.driver_type == DriverType.DEL]
        if del_drivers:
            driver_string = "\n".join([str(driver) for driver in del_drivers])
            lines.append(
                f"WARN: DEL driver unreliable due to purity {purple_qc.purity} < {DEL_MIN_PURITY_THRESHOLD}:\n{driver_string}"
            )
    suspicious_partial_dels = [
        driver for driver in drivers
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
                exon_coverage for exon_coverage in exon_median_coverages if
                exon_coverage.gene_name == partial_del.gene_name
            ]
            for exon_coverage in sorted(relevant_exon_coverages, key=lambda x: x.exon_rank):
                deletion_status = get_deletion_status(exon_coverage, somatic_copy_numbers)
                per_driver_lines.append("\t".join(
                    [exon_coverage.gene_name, str(exon_coverage.exon_rank), deletion_status,
                     str(exon_coverage.median_depth)]))
        driver_string = "\n".join(per_driver_lines)
        lines.append(f"WARN: Partial del called in suspicious gene:\n{driver_string}")
    dels_in_suspicious_partial_del_gene = [
        driver for driver in drivers
        if driver.driver_type == DriverType.DEL
           and driver.gene_name in RECURRENT_PARTIAL_DEL_GENES
           and driver.max_copy_number < DELETION_COPY_NUMBER_THRESHOLD
    ]
    if dels_in_suspicious_partial_del_gene:
        driver_string = "\n".join([str(driver) for driver in dels_in_suspicious_partial_del_gene])
        lines.append(f"WARN: Full deletion called in gene suspicious for partial DELs:\n{driver_string}")

    return "\n".join(lines)


def get_deletion_status(exon_coverage: ExonMedianCoverage, somatic_copy_numbers: List[SomaticCopyNumber]) -> str:
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


def load_wgs_metrics(local_directory: Path, tumor_name: str) -> WgsMetrics:
    wgs_metrics_path = local_directory / f"{tumor_name}.wgsmetrics"
    return load_wgs_metrics_from_file(wgs_metrics_path)


def load_wgs_metrics_from_file(wgs_metrics_path) -> WgsMetrics:
    with open(wgs_metrics_path, "r") as in_f:
        # skip first line
        next(in_f)
        header = next(in_f)
        data = next(in_f)
    return get_wgs_metrics_from_line(data, header)


def get_wgs_metrics_from_line(line: str, header: str) -> WgsMetrics:
    split_header = header.replace("\n", "").split(TSV_SEPARATOR)

    split_line = line.replace("\n", "").split(TSV_SEPARATOR)

    mean_coverage = Decimal(split_line[split_header.index("MEAN_COVERAGE")])
    median_coverage = Decimal(split_line[split_header.index("MEDIAN_COVERAGE")])
    percent_excluded_adapter = Decimal(split_line[split_header.index("PCT_EXC_ADAPTER")])
    percent_excluded_mapq = Decimal(split_line[split_header.index("PCT_EXC_MAPQ")])
    percent_excluded_duplicated = Decimal(split_line[split_header.index("PCT_EXC_DUPE")])
    percent_excluded_unpaired = Decimal(split_line[split_header.index("PCT_EXC_UNPAIRED")])
    percent_excluded_baseq = Decimal(split_line[split_header.index("PCT_EXC_BASEQ")])
    percent_excluded_overlap = Decimal(split_line[split_header.index("PCT_EXC_OVERLAP")])
    percent_excluded_capped = Decimal(split_line[split_header.index("PCT_EXC_CAPPED")])
    percent_excluded_total = Decimal(split_line[split_header.index("PCT_EXC_TOTAL")])
    percent_100x = Decimal(split_line[split_header.index("PCT_100X")])

    wgs_metrics = WgsMetrics(
        mean_coverage,
        median_coverage,
        percent_excluded_adapter,
        percent_excluded_mapq,
        percent_excluded_duplicated,
        percent_excluded_unpaired,
        percent_excluded_baseq,
        percent_excluded_overlap,
        percent_excluded_capped,
        percent_excluded_total,
        percent_100x,
    )
    return wgs_metrics



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
    )
    return gene_copy_number


def load_purple_qc(local_directory: Path, tumor_name: str) -> PurpleQc:
    purple_qc_path = local_directory / f"{tumor_name}.purple.qc"
    return load_purple_qc_from_file(purple_qc_path)


def load_purple_qc_from_file(purple_qc_path: Path) -> PurpleQc:
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

        purple_qc = PurpleQc(
            qc_status,
            copy_number_segments,
            unsupported_copy_number_segments,
            purity,
            amber_gender,
            cobalt_gender,
            deleted_genes,
            amber_mean_depth,
        )
        return purple_qc


def load_drivers(local_directory: Path, tumor_name: str) -> List[Driver]:
    local_purple_driver_catalog = local_directory / f"{tumor_name}.purple.driver.catalog.somatic.tsv"
    if not local_purple_driver_catalog.exists():
        local_purple_driver_catalog = local_directory / f"{tumor_name}.driver.catalog.somatic.tsv"
    return load_drivers_from_file(local_purple_driver_catalog)


def load_drivers_from_file(driver_catalog_path: Path) -> List[Driver]:
    drivers = []
    with open(driver_catalog_path, "r") as in_f:
        header = next(in_f)
        for line in in_f:
            driver = get_driver_from_line(line, header)
            drivers.append(driver)
    return drivers


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


def get_bool_from_string(string: str) -> bool:
    if string == "true":
        return True
    elif string == "false":
        return False
    else:
        raise ValueError(f"Unrecognized boolean string: {string}")


def copy_purple_files_to_local(gcp_run_url: str, local_directory: Path) -> None:
    purple_url = f"{gcp_run_url}/purple"
    sync_directory_to_local(purple_url, local_directory)


def copy_sage_exons_median_tsv_to_local(gcp_run_url: str, local_directory: Path, tumor_name: str) -> None:
    sage_exons_median_tsv_url = f"{gcp_run_url}/sage_somatic/{tumor_name}.sage.exon.medians.tsv"
    copy_file_to_local_directory(sage_exons_median_tsv_url, local_directory)


def copy_wgs_metrics_file(gcp_run_url: str, local_directory: Path, tumor_name: str) -> None:
    wgs_metrics_url = f"{gcp_run_url}/{tumor_name}/bam_metrics/{tumor_name}.wgsmetrics"
    copy_file_to_local_directory(wgs_metrics_url, local_directory)


def sync_directory_to_local(source_url: str, local_directory: Path) -> None:
    local_directory.mkdir(parents=True, exist_ok=True)
    subprocess.run(["gsutil", "-m", "rsync", "-r", source_url, local_directory], capture_output=True)


def copy_file_to_local_directory(source_url: str, local_directory: Path) -> None:
    local_directory.mkdir(parents=True, exist_ok=True)
    subprocess.run(["gsutil", "-m", "cp", source_url, local_directory], capture_output=True)


def get_metadata_json(gcp_run_url: str) -> Dict[str, Any]:
    metadata_json = json.loads(gsutil_cat(f"{gcp_run_url}/metadata.json"))
    return metadata_json


def gsutil_cat(gcp_url: str) -> str:
    return subprocess.run(["gsutil", "cat", gcp_url], capture_output=True, text=True).stdout


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="qc_check_oncopanel_run.py", description="Perform QC checks for OncoPanel run")
    parser.add_argument("--gcp-run-url", "-u", type=str, required=True, help="GCP path to run directory")
    parser.add_argument("--working-directory", "-w", type=Path, required=True, help="Local working directory")
    parser.add_argument("--driver-gene-panel", "-d", type=Path, required=True, help="Driver gene panel")
    parser.add_argument("--output-file", "-o", type=Path, help="Output file with summary of qc check results")
    return parser.parse_args()


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s - [%(levelname)-8s] - %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
    )
    args = parse_args()
    main(args.gcp_run_url, args.working_directory, args.driver_gene_panel, args.output_file)
