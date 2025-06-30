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
from typing import List, Tuple, Dict, Set

CRLF2_GENE_NAME = "CRLF2"

WARN_HIGH_COPY_NUMBER_NOISE = "WARN_HIGH_COPY_NUMBER_NOISE"
WARN_DELETED_GENES = "WARN_DELETED_GENES"
FAIL_NO_TUMOR = "FAIL_NO_TUMOR"
FAIL_CONTAMINATION = "FAIL_CONTAMINATION"

CHR_X = "chrX"
CHR_X_NON_PAR_START = 2781479
CHR_X_NON_PAR_END = 156030895

TSV_SEPARATOR = "\t"

AMP_MANUAL_INTERPRETATION_THRESHOLD = 7
VCHORD_MIN_PURITY_THRESHOLD = Decimal("0.30")
DEL_MIN_PURITY_THRESHOLD = Decimal("0.30")
AMP_MIN_PURITY_THRESHOLD = Decimal("0.20")
TMB_MIN_PURITY_THRESHOLD = Decimal("0.10")
MSI_MIN_PURITY_THRESHOLD = Decimal("0.10")
SOMATIC_VARIANT_MIN_PURITY_THRESHOLD = Decimal("0.10")

DELETION_COPY_NUMBER_THRESHOLD = Decimal("0.5")
GENES_EXCLUDED_FROM_DESIGN = {"SPATA31A7", "LINC01001", "U2AF1"}

BASES_PER_GBASE = 1000000000
TARGET_YIELD_IN_BASES = 50 * BASES_PER_GBASE
TARGET_PERCENT_EXONS_WITH_MEDIAN_COVERAGE_AT_LEAST_100 = 95.

AMPLIFICATIONS_DESCRIPTION = "amplificaties"


@dataclass(frozen=True, eq=True)
class Metadata:
    tumor_name: str
    tumor_barcode: str
    yield_in_bases: int


class DriverType(Enum):
    AMP = auto()
    DEL = auto()
    MUTATION = auto()
    PARTIAL_AMP = auto()

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
    contamination: Decimal


@dataclass(frozen=True, eq=True)
class ExonMedianCoverage:
    gene_name: str
    chromosome: str
    start: int
    end: int
    exon_rank: int
    median_depth: int


@dataclass(frozen=True, eq=True)
class AmberBafPoint:
    chromosome: str
    position: int


@dataclass(frozen=True, eq=True)
class RunData:
    drivers: Tuple[Driver, ...]
    purple_qc: PurpleQc
    exon_median_coverages: Tuple[ExonMedianCoverage, ...]
    amber_baf_points: Tuple[AmberBafPoint, ...]


@dataclass(frozen=True, eq=True)
class ExonExclusion:
    gene_name: str
    gene_id: str
    canonical_transcript_id: str
    excluded_exons: Tuple[int, ...]


@dataclass(frozen=True, eq=True)
class ExonCoverageStatistics:
    relevant_exon_count: int
    exons_with_median_coverage_at_least_100_count: int

    @classmethod
    def from_exon_coverages(
            cls, exon_median_coverages: Tuple[ExonMedianCoverage, ...], excluded_exons: List[ExonExclusion]
    ) -> "ExonCoverageStatistics":
        gene_to_excluded_exons = {exon_exclusion.gene_name: exon_exclusion.excluded_exons for exon_exclusion in excluded_exons}
        relevant_exons = [exon for exon in exon_median_coverages if is_relevant_for_coverage(exon, gene_to_excluded_exons)]
        if not relevant_exons:
            raise ValueError("No exons found relevant for coverage. Check excluded genes and exons!")
        exons_with_median_coverage_at_least_100 = [exon for exon in relevant_exons if exon.median_depth >= 100]
        return ExonCoverageStatistics(len(relevant_exons), len(exons_with_median_coverage_at_least_100))

    def percent_exons_median_coverage_at_least_100(self) -> float:
        return 100 * self.exons_with_median_coverage_at_least_100_count / self.relevant_exon_count


@dataclass(frozen=True, eq=True)
class SampleSummary:
    tumor_name: str
    remarks: List[str]

    def summary_text(self) -> str:
        return " ".join(self.remarks)


def main(run_names: str, input_directory: Path, excluded_exons_file_path: Path, output_file_path: Path) -> None:
    logging.info("Start QC checking OncoPanel runs")

    pipeline_directories = determine_relevant_pipeline_directories(input_directory, run_names)

    logging.info("Load excluded_exons")
    excluded_exons = load_excluded_exons(excluded_exons_file_path)

    sample_summaries = determine_sample_summaries(pipeline_directories, excluded_exons)

    write_summary_file(sample_summaries, output_file_path)

    logging.info("Finished QC checking OncoPanel run")


def determine_sample_summaries(
        pipeline_directories: Set[Path], excluded_exons: List[ExonExclusion]
) -> List[SampleSummary]:
    return [determine_sample_summary(pipeline_directory, excluded_exons) for pipeline_directory in pipeline_directories]


def determine_sample_summary(
        pipeline_directory: Path, excluded_exons: List[ExonExclusion]
) -> SampleSummary:
    logging.info(f"Get metadata from {pipeline_directory}")
    metadata = get_metadata(pipeline_directory)

    logging.info(f"Handling run for {metadata.tumor_name}: {pipeline_directory}")

    run_data = load_run_data(pipeline_directory, metadata.tumor_name)

    exon_coverage_statistics = ExonCoverageStatistics.from_exon_coverages(run_data.exon_median_coverages, excluded_exons)

    remarks = determine_remarks(run_data, exon_coverage_statistics, metadata)
    return SampleSummary(metadata.tumor_name, remarks)


def determine_remarks(
        run_data: RunData,
        exon_coverage_statistics: ExonCoverageStatistics,
        metadata: Metadata,
) -> List[str]:
    remarks = []

    percent_exons_median_coverage_at_least_100 = exon_coverage_statistics.percent_exons_median_coverage_at_least_100()
    if metadata.yield_in_bases < TARGET_YIELD_IN_BASES:
        if percent_exons_median_coverage_at_least_100 > TARGET_PERCENT_EXONS_WITH_MEDIAN_COVERAGE_AT_LEAST_100:
            remarks.append(
                f"Dit sample heeft {convert_bases_to_rounded_gbases(metadata.yield_in_bases)}Gb yield, "
                f"wat minder is dan ons eigenlijke minimum van {convert_bases_to_rounded_gbases(TARGET_YIELD_IN_BASES)}Gb yield. "
                f"De coverage van dit sample is wel goed: "
                f"{percent_exons_median_coverage_at_least_100:.2f}% van de exonen heeft median coverage minstens 100x."
            )
        else:
            remarks.append(
                f"Dit sample heeft {convert_bases_to_rounded_gbases(metadata.yield_in_bases)}Gb yield, "
                f"wat minder is dan ons eigenlijke minimum van {convert_bases_to_rounded_gbases(TARGET_YIELD_IN_BASES)}Gb yield."
            )

    if percent_exons_median_coverage_at_least_100 < TARGET_PERCENT_EXONS_WITH_MEDIAN_COVERAGE_AT_LEAST_100:
        remarks.append(
            f"De coverage van dit sample is te laag: "
            f"slechts {percent_exons_median_coverage_at_least_100:.2f}% van de exonen heeft median coverage minstens 100x, "
            f"terwijl dit minstens {TARGET_PERCENT_EXONS_WITH_MEDIAN_COVERAGE_AT_LEAST_100}% zou moeten zijn. "
            f"Hierdoor kunnen varianten gemist worden."
        )

    if run_data.purple_qc.qc_status != "PASS":
        remarks.append(
            f"De Purple QC status van dit sample is {run_data.purple_qc.qc_status}."
        )
        if FAIL_CONTAMINATION in run_data.purple_qc.qc_status:
            remarks.append(
                f"'{FAIL_CONTAMINATION}' komt doordat er volgens het algoritme {run_data.purple_qc.contamination}% contaminatie aanwezig is.")
        if FAIL_NO_TUMOR in run_data.purple_qc.qc_status:
            remarks.append(
                f"'{FAIL_NO_TUMOR}' betekent dat er volgens het algoritme geen of slechts zeer weinig tumor aanwezig is in het sample. "
                f"Dit zou kunnen komen doordat dit een normal sample of een zeer laag mTCP sample is, "
                f"of doordat dit een rustig type tumor is waarvoor onze mTCP schatter niet goed werkt.")
        if WARN_DELETED_GENES in run_data.purple_qc.qc_status or WARN_HIGH_COPY_NUMBER_NOISE in run_data.purple_qc.qc_status:
            if WARN_DELETED_GENES in run_data.purple_qc.qc_status and WARN_HIGH_COPY_NUMBER_NOISE in run_data.purple_qc.qc_status:
                sentence_start = f"'{WARN_DELETED_GENES}' en '{WARN_HIGH_COPY_NUMBER_NOISE}' betekenen"
            elif WARN_DELETED_GENES in run_data.purple_qc.qc_status:
                sentence_start = f"'{WARN_DELETED_GENES}' betekent"
            else:
                sentence_start = f"'{WARN_HIGH_COPY_NUMBER_NOISE}' betekent"
            remarks.append(
                f"{sentence_start} dat er indicaties zijn dat de gekozen mTCP en/of ploidy voor dit sample niet correct zijn."
            )
            if WARN_HIGH_COPY_NUMBER_NOISE in run_data.purple_qc.qc_status:
                remarks.append("Om deze reden worden amplificaties en deleties alleen gecalld met SV support.")
            else:
                remarks.append("Om deze reden worden deleties alleen gecalld met SV support.")

    baf_points_in_x_outside_par = [
        baf_point for baf_point in run_data.amber_baf_points
        if baf_point.chromosome == CHR_X and CHR_X_NON_PAR_START <= baf_point.position <= CHR_X_NON_PAR_END
    ]
    percent_baf_points_in_x_outside_par = (len(baf_points_in_x_outside_par) * 100) / max(len(run_data.amber_baf_points), 1)
    if run_data.purple_qc.amber_gender == "MALE" and percent_baf_points_in_x_outside_par > 0.5:
        remarks.append(
            f"Door het relatief hoge aantal BAF punten op chromosoom X (aantal: {len(baf_points_in_x_outside_par)}, percentage: {percent_baf_points_in_x_outside_par:.2f}%) "
            f"is de gender call van dit sample ({run_data.purple_qc.amber_gender}) onbetrouwbaar."
        )

    if run_data.purple_qc.amber_gender == "FEMALE" and percent_baf_points_in_x_outside_par < 1.5:
        remarks.append(
            f"Door het relatief lage aantal BAF punten op chromosoom X (aantal: {len(baf_points_in_x_outside_par)}, percentage: {percent_baf_points_in_x_outside_par:.2f}%) "
            f"is de gender call van dit sample ({run_data.purple_qc.amber_gender}) onbetrouwbaar."
        )

    purity_too_low = set()
    if run_data.purple_qc.purity < TMB_MIN_PURITY_THRESHOLD:
        purity_too_low.add("TMB")

    if run_data.purple_qc.purity < MSI_MIN_PURITY_THRESHOLD:
        purity_too_low.add("MSI")

    if run_data.purple_qc.purity < SOMATIC_VARIANT_MIN_PURITY_THRESHOLD:
        purity_too_low.add("kleine varianten")

    if run_data.purple_qc.purity < DEL_MIN_PURITY_THRESHOLD:
        purity_too_low.add("deleties")

    if run_data.purple_qc.purity < AMP_MIN_PURITY_THRESHOLD:
        purity_too_low.add(AMPLIFICATIONS_DESCRIPTION)

    if run_data.purple_qc.purity < VCHORD_MIN_PURITY_THRESHOLD:
        purity_too_low.add("HRD")

    if purity_too_low:
        remarks.append(
            f"De mTCP van dit sample ({run_data.purple_qc.purity*100:.2f}%) is te laag "
            f"voor het betrouwbaar callen van {pretty_listing(purity_too_low)}."
        )

    amp_drivers_needing_manual_curation = [
        driver for driver in run_data.drivers
        if driver.driver_type == DriverType.AMP and driver.min_copy_number <= AMP_MANUAL_INTERPRETATION_THRESHOLD and driver.gene_name != CRLF2_GENE_NAME
    ]
    if amp_drivers_needing_manual_curation:
        gene_name_string = ", ".join([driver.gene_name for driver in amp_drivers_needing_manual_curation])
        remarks.append(f"De volgende amplificaties moeten vanwege hun lage copy number met de hand beoordeeld worden: {gene_name_string}.")

    crlf2_amp_present = any(
        driver for driver in run_data.drivers
        if driver.driver_type == DriverType.AMP and driver.gene_name == CRLF2_GENE_NAME
    )
    if crlf2_amp_present and not AMPLIFICATIONS_DESCRIPTION in purity_too_low:
        remarks.append(f"De {CRLF2_GENE_NAME} amplificatie is onbetrouwbaar.")

    partial_amp_drivers = [driver for driver in run_data.drivers if driver.driver_type == DriverType.PARTIAL_AMP]
    if partial_amp_drivers:
        gene_name_string = ", ".join(sorted(driver.gene_name for driver in partial_amp_drivers))
        remarks.append(f"Partiële amplificaties zijn niet gevalideerd: {gene_name_string}.")

    if not remarks:
        remarks.append("Geen opmerkingen.")

    return remarks


def determine_relevant_pipeline_directories(input_directory: Path, run_names: str) -> Set[Path]:
    pipeline_directories = set()
    for given_run_name in run_names.split(";"):
        if "*" in given_run_name:
            for run_name in input_directory.glob(given_run_name):
                pipeline_directories.add(run_name)
        else:
            pipeline_directories.add(input_directory / given_run_name)
    return pipeline_directories


def write_summary_file(sample_summaries: List[SampleSummary], output_file: Path) -> None:
    with open(output_file, "w") as output_f:
        output_f.write(construct_summary_file_text(sample_summaries))


def construct_summary_file_text(sample_summaries: List[SampleSummary]) -> str:
    lines = ["Sample ID\tOpmerkingen"]
    for sample_summary in sorted(sample_summaries, key=lambda sample_summary: sample_summary.tumor_name):
        lines.append(f"{sample_summary.tumor_name}\t{sample_summary.summary_text()}")
    return "\n".join(lines) + "\n"


def is_relevant_for_coverage(exon: ExonMedianCoverage, gene_to_excluded_exons: Dict[str, Tuple[int, ...]]) -> bool:
    gene_explicitly_excluded_from_design = exon.gene_name in GENES_EXCLUDED_FROM_DESIGN
    gene_considered_for_excluded_exons = exon.gene_name in gene_to_excluded_exons.keys()
    exon_excluded = gene_considered_for_excluded_exons and exon.exon_rank in gene_to_excluded_exons[exon.gene_name]
    return not gene_explicitly_excluded_from_design and gene_considered_for_excluded_exons and not exon_excluded


def load_run_data(local_directory: Path, tumor_name: str) -> RunData:
    logging.info("Load drivers")
    drivers = load_drivers(local_directory, tumor_name)
    logging.info("Load sample qc info")
    purple_qc = load_purple_qc(local_directory, tumor_name)
    logging.info("Load exon median coverages")
    exon_median_coverages = load_exon_median_coverages(local_directory, tumor_name)
    logging.info("Load Amber BAF points")
    amber_baf_points = load_amber_baf_points(local_directory, tumor_name)
    run_data = RunData(
        tuple(drivers),
        purple_qc,
        tuple(exon_median_coverages),
        tuple(amber_baf_points),
    )
    return run_data


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


def load_amber_baf_points(local_directory: Path, tumor_name: str) -> List[AmberBafPoint]:
    amber_baf_file_path = local_directory / "amber" / f"{tumor_name}.amber.baf.tsv.gz"
    return load_amber_baf_points_from_file(amber_baf_file_path)


def load_amber_baf_points_from_file(amber_baf_file_path: Path) -> List[AmberBafPoint]:
    baf_points = []
    with gzip.open(amber_baf_file_path, "rt") as in_f:
        header = next(in_f)
        for line in in_f:
            baf_point = get_baf_point_from_line(line, header)
            baf_points.append(baf_point)
    return baf_points


def get_baf_point_from_line(line: str, header: str) -> AmberBafPoint:
    split_header = header.replace("\n", "").split(TSV_SEPARATOR)
    split_line = line.replace("\n", "").split(TSV_SEPARATOR)

    chromosome = split_line[split_header.index("chromosome")]
    position = int(split_line[split_header.index("position")])
    return AmberBafPoint(chromosome, position)


def load_exon_median_coverages(local_directory: Path, tumor_name: str) -> List[ExonMedianCoverage]:
    sage_exon_medians_path = local_directory / "sage_somatic" / f"{tumor_name}.sage.exon.medians.tsv"
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


def load_purple_qc(local_directory: Path, tumor_name: str) -> PurpleQc:
    purple_qc_path = local_directory / "purple" / f"{tumor_name}.purple.qc"
    purple_purity_path = local_directory / "purple" / f"{tumor_name}.purple.purity.tsv"
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
        contamination = Decimal(field_to_string_value["Contamination"])

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
        tmb_status,
        contamination
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


def get_local_purple_driver_catalog(local_directory: Path, tumor_name: str) -> Path:
    local_purple_driver_catalog = local_directory / "purple" / f"{tumor_name}.purple.driver.catalog.somatic.tsv"
    if not local_purple_driver_catalog.exists():
        local_purple_driver_catalog = local_directory / "purple" / f"{tumor_name}.driver.catalog.somatic.tsv"
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


def get_metadata(pipeline_directory: Path) -> Metadata:
    with open(pipeline_directory / "metadata.json") as in_f:
        metadata_json = json.load(in_f)
    tumor_name = metadata_json["tumor"]["sampleName"]
    tumor_barcode = metadata_json["tumor"]["barcode"]

    yield_in_bases = get_yield_in_bases_from_api(tumor_barcode)

    return Metadata(tumor_name, tumor_barcode, yield_in_bases)


def get_yield_in_bases_from_api(tumor_barcode: str) -> int:
    api_output = run_bash_command(
        ["curl", "--fail", "--silent", "--show-error", "-H", "Content-Type: application/json", "-X",
         "GET", f"http://api.prod-1/hmf/v1/samples?barcode={tumor_barcode}"]
    )
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


def pretty_listing(things: Set[str]) -> str:
    sorted_things = sorted(things)
    if len(sorted_things) >= 2:
        pretty_list_string = ", ".join(sorted_things[:-1])
        pretty_list_string += f" en {sorted_things[-1]}"
    elif len(sorted_things) == 1:
        pretty_list_string = f"{sorted_things[0]}"
    else:
        pretty_list_string = ""
    return pretty_list_string


def convert_bases_to_rounded_gbases(base_count: int) -> Decimal:
    return (Decimal(base_count) / Decimal(BASES_PER_GBASE)).quantize(Decimal("0.1"), rounding=ROUND_DOWN)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="oncoact_panel_remarks.py", description="Create automatic remarks for OncoAct Panel runs")
    parser.add_argument("--input-directory", "-i", type=Path, required=True, help="Input directory")
    parser.add_argument("--run-names", "-n", type=str, required=True, help="List of run names to include, separated by ';'. Use '*' as a wildcard.")
    parser.add_argument("--excluded-exons", "-e", type=Path, required=True, help="Exons excluded from coverage analysis")
    parser.add_argument("--output-file", "-o", type=Path, required=True, help="Output file with sample remarks")
    return parser.parse_args()


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s - [%(levelname)-8s] - %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
    )
    args = parse_args()
    main(args.run_names, args.input_directory, args.excluded_exons, args.output_file)
