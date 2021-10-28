import logging
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple, Dict, DefaultDict

from genome import BafSite, PgxSite, MsiSite, Position, Exon, Interval, FusionSite
from coverage_info import CoverageInfo


class OutputWriter(object):
    @classmethod
    def write_baf_coverage_file(
            cls,
            sample_with_coverage_info_list: List[Tuple[str, CoverageInfo]],
            baf_sites_list: Tuple[BafSite, ...],
            output_dir: Path,
    ) -> None:
        logging.info("Started BAF point coverage analysis")
        baf_site_output_file = f"{output_dir}/baf_coverage.tsv"
        if Path(baf_site_output_file).exists():
            error_msg = (
                f"At least one of the output files already exists: "
                f"{baf_site_output_file}."
            )
            raise FileExistsError(error_msg)

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        with open(baf_site_output_file, "w") as baf_o:
            base_header_elements = ["chrom", "position", "label", "probe_start", "probe_end"]
            header_elements = base_header_elements + [sample for sample, _ in sample_with_coverage_info_list]
            header = "\t".join(header_elements) + "\n"
            baf_o.write(header)
            for site in baf_sites_list:
                base_line_data = [
                    site.chromosome,
                    str(site.position),
                    site.label,
                    str(site.probe_start),
                    str(site.probe_end)
                ]
                samples_line_data = [
                    str(coverage_info.get_cumulative_coverage(site.get_site_interval()))
                    for sample, coverage_info in sample_with_coverage_info_list
                ]
                line_data = base_line_data + samples_line_data
                line = "\t".join(line_data) + "\n"
                baf_o.write(line)

        logging.info("Finished BAF point coverage analysis")

    @classmethod
    def write_exome_count_min_coverage_file(
            cls,
            sample_with_coverage_info_list: List[Tuple[str, CoverageInfo]],
            exons: Tuple[Exon, ...],
            min_coverage: int,
            output_dir: Path,
    ) -> None:
        logging.info(f"Started exome min coverage {min_coverage} count analysis")
        exon_min_coverage_count_output_file = f"{output_dir}/exon_min_coverage_count.{min_coverage}.tsv"
        gene_min_coverage_count_output_file = f"{output_dir}/gene_min_coverage_count.{min_coverage}.tsv"
        if Path(exon_min_coverage_count_output_file).exists() or Path(gene_min_coverage_count_output_file).exists():
            error_msg = (
                f"At least one of the output files already exists: "
                f"{exon_min_coverage_count_output_file} and/or {gene_min_coverage_count_output_file}."
            )
            raise FileExistsError(error_msg)

        sample_to_exon_to_min_coverage_count = {
            sample: {exon: coverage_info.get_count_with_min_coverage(exon.interval, min_coverage) for exon in exons}
            for sample, coverage_info in sample_with_coverage_info_list
        }

        sample_to_gene_to_min_coverage_count = cls.__sum_exon_data_to_gene_data(sample_to_exon_to_min_coverage_count)

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        samples = [sample for sample, _ in sample_with_coverage_info_list]
        cls.__write_exon_min_coverage_count_output_file(
            sample_to_exon_to_min_coverage_count, samples, exons, exon_min_coverage_count_output_file
        )
        cls.__write_gene_min_coverage_count_output_file(
            sample_to_gene_to_min_coverage_count, samples, exons, gene_min_coverage_count_output_file
        )

        logging.info(f"Finished exome min coverage {min_coverage} count analysis")

    @classmethod
    def write_exome_cumulative_coverage_analysis(
            cls,
            sample_with_coverage_info_list: List[Tuple[str, CoverageInfo]],
            exons: Tuple[Exon, ...],
            output_dir: Path,
    ) -> None:
        logging.info("Started exome cumulative coverage analysis")
        exon_coverage_output_file = f"{output_dir}/exon_cumulative_coverage.tsv"
        gene_coverage_output_file = f"{output_dir}/gene_cumulative_coverage.tsv"
        if Path(exon_coverage_output_file).exists() or Path(gene_coverage_output_file).exists():
            error_msg = (
                f"At least one of the output files already exists: "
                f"{exon_coverage_output_file} and/or {gene_coverage_output_file}."
            )
            raise FileExistsError(error_msg)

        sample_to_exon_to_cumulative_coverage = {
            sample: {exon: coverage_info.get_cumulative_coverage(exon.interval) for exon in exons}
            for sample, coverage_info in sample_with_coverage_info_list
        }

        sample_to_gene_to_cumulative_coverage = cls.__sum_exon_data_to_gene_data(sample_to_exon_to_cumulative_coverage)

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        samples = [sample for sample, _ in sample_with_coverage_info_list]
        cls.__write_exon_cumulative_coverage_output_file(
            sample_to_exon_to_cumulative_coverage, samples, exons, exon_coverage_output_file
        )
        cls.__write_gene_cumulative_coverage_output_file(
            sample_to_gene_to_cumulative_coverage, samples, exons, gene_coverage_output_file
        )

        logging.info("Finished exome cumulative coverage analysis")

    @classmethod
    def write_fusion_count_min_coverage_file(
            cls,
            sample_with_coverage_info_list: List[Tuple[str, CoverageInfo]],
            fusion_sites_list: Tuple[FusionSite, ...],
            min_coverage: int,
            output_dir: Path,
    ) -> None:
        logging.info(f"Started fusion min coverage {min_coverage} count analysis")
        fusion_site_min_coverage_count_output_file = f"{output_dir}/fusion_site_min_coverage_count.{min_coverage}.tsv"
        if Path(fusion_site_min_coverage_count_output_file).exists():
            error_msg = (
                f"At least one of the output files already exists: "
                f"{fusion_site_min_coverage_count_output_file}."
            )
            raise FileExistsError(error_msg)

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        with open(fusion_site_min_coverage_count_output_file, "w") as exon_o:
            base_header_elements = ["chrom", "start_position", "end_position", "gene", "intron_start", "intron_end"]
            header_elements = base_header_elements + [sample for sample, _ in sample_with_coverage_info_list]
            header = "\t".join(header_elements) + "\n"
            exon_o.write(header)

            for site in fusion_sites_list:
                base_line_data = [
                    site.interval.chromosome,
                    str(site.interval.start_position),
                    str(site.interval.end_position),
                    site.gene,
                    str(site.intron_start),
                    str(site.intron_end),
                ]
                samples_line_data = [
                    str(coverage_info.get_count_with_min_coverage(site.interval, min_coverage))
                    for sample, coverage_info in sample_with_coverage_info_list
                ]
                line_data = base_line_data + samples_line_data
                line = "\t".join(line_data) + "\n"
                exon_o.write(line)

        logging.info(f"Finished fusion min coverage {min_coverage} count analysis")

    @classmethod
    def write_fusion_cumulative_coverage_file(
            cls,
            sample_with_coverage_info_list: List[Tuple[str, CoverageInfo]],
            fusion_sites_list: Tuple[FusionSite, ...],
            output_dir: Path,
    ) -> None:
        logging.info("Started fusion cumulative coverage analysis")
        fusion_site_cumulative_coverage_output_file = f"{output_dir}/fusion_site_cumulative_coverage.tsv"
        if Path(fusion_site_cumulative_coverage_output_file).exists():
            error_msg = (
                f"At least one of the output files already exists: "
                f"{fusion_site_cumulative_coverage_output_file}."
            )
            raise FileExistsError(error_msg)

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        with open(fusion_site_cumulative_coverage_output_file, "w") as exon_o:
            base_header_elements = ["chrom", "start_position", "end_position", "gene", "intron_start", "intron_end"]
            header_elements = base_header_elements + [sample for sample, _ in sample_with_coverage_info_list]
            header = "\t".join(header_elements) + "\n"
            exon_o.write(header)

            for site in fusion_sites_list:
                base_line_data = [
                    site.interval.chromosome,
                    str(site.interval.start_position),
                    str(site.interval.end_position),
                    site.gene,
                    str(site.intron_start),
                    str(site.intron_end),
                ]
                samples_line_data = [
                    str(coverage_info.get_cumulative_coverage(site.interval))
                    for sample, coverage_info in sample_with_coverage_info_list
                ]
                line_data = base_line_data + samples_line_data
                line = "\t".join(line_data) + "\n"
                exon_o.write(line)

        logging.info("Finished fusion cumulative coverage analysis")

    @classmethod
    def write_hotspot_coverage_file(
            cls,
            sample_with_coverage_info_list: List[Tuple[str, CoverageInfo]],
            hotspot_list: Tuple[Position, ...],
            output_dir: Path,
    ) -> None:
        logging.info("Started hotspot coverage analysis")
        hotspot_output_file = f"{output_dir}/hotspot_coverage.tsv"
        if Path(hotspot_output_file).exists():
            error_msg = (
                f"At least one of the output files already exists: "
                f"{hotspot_output_file}."
            )
            raise FileExistsError(error_msg)

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        with open(hotspot_output_file, "w") as hotspot_o:
            base_header_elements = ["chrom", "position"]
            header_elements = base_header_elements + [sample for sample, _ in sample_with_coverage_info_list]
            header = "\t".join(header_elements) + "\n"
            hotspot_o.write(header)
            for hotspot in hotspot_list:
                base_line_data = [hotspot.chromosome, str(hotspot.position)]
                samples_line_data = [
                    str(coverage_info.get_cumulative_coverage(hotspot.get_interval()))
                    for sample, coverage_info in sample_with_coverage_info_list
                ]
                line_data = base_line_data + samples_line_data
                line = "\t".join(line_data) + "\n"
                hotspot_o.write(line)

        logging.info("Finished hotspot coverage analysis")

    @classmethod
    def write_msi_count_min_coverage_file(
            cls,
            sample_with_coverage_info_list: List[Tuple[str, CoverageInfo]],
            msi_sites_list: Tuple[MsiSite, ...],
            min_coverage: int,
            output_dir: Path,
    ) -> None:
        logging.info(f"Started msi min coverage {min_coverage} count analysis")
        msi_site_min_coverage_count_output_file = f"{output_dir}/msi_site_min_coverage_count.{min_coverage}.tsv"
        if Path(msi_site_min_coverage_count_output_file).exists():
            error_msg = (
                f"At least one of the output files already exists: "
                f"{msi_site_min_coverage_count_output_file}."
            )
            raise FileExistsError(error_msg)

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        with open(msi_site_min_coverage_count_output_file, "w") as exon_o:
            base_header_elements = [
                "chrom", "position", "repeat_count",
                "probe3_start", "probe3_end", "probe3_id",
                "probe5_start", "probe5_end", "probe5_id"
            ]
            header_elements = base_header_elements + [sample for sample, _ in sample_with_coverage_info_list]
            header = "\t".join(header_elements) + "\n"
            exon_o.write(header)

            for site in msi_sites_list:
                base_line_data = [
                    site.chromosome,
                    str(site.position),
                    str(site.repeat_count),
                    str(site.probe3_start),
                    str(site.probe3_end),
                    site.probe3_id,
                    str(site.probe5_start),
                    str(site.probe5_end),
                    site.probe5_id,
                ]
                samples_line_data = [
                    str(coverage_info.get_count_with_min_coverage(site.get_site_interval(), min_coverage))
                    for sample, coverage_info in sample_with_coverage_info_list
                ]
                line_data = base_line_data + samples_line_data
                line = "\t".join(line_data) + "\n"
                exon_o.write(line)

        logging.info(f"Finished msi min coverage {min_coverage} count analysis")

    @classmethod
    def write_msi_cumulative_coverage_file(
            cls,
            sample_with_coverage_info_list: List[Tuple[str, CoverageInfo]],
            msi_sites_list: Tuple[MsiSite, ...],
            output_dir: Path,
    ) -> None:
        logging.info("Started msi cumulative coverage analysis")
        msi_site_cumulative_coverage_output_file = f"{output_dir}/msi_site_cumulative_coverage.tsv"
        if Path(msi_site_cumulative_coverage_output_file).exists():
            error_msg = (
                f"At least one of the output files already exists: "
                f"{msi_site_cumulative_coverage_output_file}."
            )
            raise FileExistsError(error_msg)

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        with open(msi_site_cumulative_coverage_output_file, "w") as exon_o:
            base_header_elements = [
                "chrom", "position", "repeat_count",
                "probe3_start", "probe3_end", "probe3_id",
                "probe5_start", "probe5_end", "probe5_id"
            ]
            header_elements = base_header_elements + [sample for sample, _ in sample_with_coverage_info_list]
            header = "\t".join(header_elements) + "\n"
            exon_o.write(header)

            for site in msi_sites_list:
                base_line_data = [
                    site.chromosome,
                    str(site.position),
                    str(site.repeat_count),
                    str(site.probe3_start),
                    str(site.probe3_end),
                    site.probe3_id,
                    str(site.probe5_start),
                    str(site.probe5_end),
                    site.probe5_id,
                ]
                samples_line_data = [
                    str(coverage_info.get_cumulative_coverage(site.get_site_interval()))
                    for sample, coverage_info in sample_with_coverage_info_list
                ]
                line_data = base_line_data + samples_line_data
                line = "\t".join(line_data) + "\n"
                exon_o.write(line)

        logging.info("Finished msi cumulative coverage analysis")

    @classmethod
    def write_pgx_count_min_coverage_file(
            cls,
            sample_with_coverage_info_list: List[Tuple[str, CoverageInfo]],
            pgx_sites_list: Tuple[PgxSite, ...],
            min_coverage: int,
            output_dir: Path,
    ) -> None:
        logging.info(f"Started pgx min coverage {min_coverage} count analysis")
        pgx_site_min_coverage_count_output_file = f"{output_dir}/pgx_site_min_coverage_count.{min_coverage}.tsv"
        if Path(pgx_site_min_coverage_count_output_file).exists():
            error_msg = (
                f"At least one of the output files already exists: "
                f"{pgx_site_min_coverage_count_output_file}."
            )
            raise FileExistsError(error_msg)

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        with open(pgx_site_min_coverage_count_output_file, "w") as exon_o:
            base_header_elements = ["chrom", "start_position", "end_position", "gene", "label"]
            header_elements = base_header_elements + [sample for sample, _ in sample_with_coverage_info_list]
            header = "\t".join(header_elements) + "\n"
            exon_o.write(header)

            for site in pgx_sites_list:
                base_line_data = [
                    site.interval.chromosome,
                    str(site.interval.start_position),
                    str(site.interval.end_position),
                    site.gene,
                    site.label,
                ]
                samples_line_data = [
                    str(coverage_info.get_count_with_min_coverage(site.interval, min_coverage))
                    for sample, coverage_info in sample_with_coverage_info_list
                ]
                line_data = base_line_data + samples_line_data
                line = "\t".join(line_data) + "\n"
                exon_o.write(line)

        logging.info(f"Finished pgx min coverage {min_coverage} count analysis")

    @classmethod
    def write_pgx_cumulative_coverage_file(
            cls,
            sample_with_coverage_info_list: List[Tuple[str, CoverageInfo]],
            pgx_sites_list: Tuple[PgxSite, ...],
            output_dir: Path,
    ) -> None:
        logging.info("Started pgx cumulative coverage analysis")
        pgx_site_cumulative_coverage_output_file = f"{output_dir}/pgx_site_cumulative_coverage.tsv"
        if Path(pgx_site_cumulative_coverage_output_file).exists():
            error_msg = (
                f"At least one of the output files already exists: "
                f"{pgx_site_cumulative_coverage_output_file}."
            )
            raise FileExistsError(error_msg)

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        with open(pgx_site_cumulative_coverage_output_file, "w") as exon_o:
            base_header_elements = ["chrom", "start_position", "end_position", "gene", "label"]
            header_elements = base_header_elements + [sample for sample, _ in sample_with_coverage_info_list]
            header = "\t".join(header_elements) + "\n"
            exon_o.write(header)

            for site in pgx_sites_list:
                base_line_data = [
                    site.interval.chromosome,
                    str(site.interval.start_position),
                    str(site.interval.end_position),
                    site.gene,
                    site.label,
                ]
                samples_line_data = [
                    str(coverage_info.get_cumulative_coverage(site.interval))
                    for sample, coverage_info in sample_with_coverage_info_list
                ]
                line_data = base_line_data + samples_line_data
                line = "\t".join(line_data) + "\n"
                exon_o.write(line)

        logging.info("Finished pgx cumulative coverage analysis")

    @classmethod
    def write_tert_count_min_coverage_file(
            cls,
            sample_with_coverage_info_list: List[Tuple[str, CoverageInfo]],
            tert_site: Interval,
            min_coverage: int,
            output_dir: Path,
    ) -> None:
        logging.info(f"Started tert min coverage {min_coverage} count analysis")
        tert_site_min_coverage_count_output_file = f"{output_dir}/tert_site_min_coverage_count.{min_coverage}.tsv"
        if Path(tert_site_min_coverage_count_output_file).exists():
            error_msg = (
                f"At least one of the output files already exists: "
                f"{tert_site_min_coverage_count_output_file}."
            )
            raise FileExistsError(error_msg)

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        with open(tert_site_min_coverage_count_output_file, "w") as exon_o:
            base_header_elements = ["chrom", "start_position", "end_position"]
            header_elements = base_header_elements + [sample for sample, _ in sample_with_coverage_info_list]
            header = "\t".join(header_elements) + "\n"
            exon_o.write(header)

            base_line_data = [
                tert_site.chromosome,
                str(tert_site.start_position),
                str(tert_site.end_position),
            ]
            samples_line_data = [
                str(coverage_info.get_count_with_min_coverage(tert_site, min_coverage))
                for sample, coverage_info in sample_with_coverage_info_list
            ]
            line_data = base_line_data + samples_line_data
            line = "\t".join(line_data) + "\n"
            exon_o.write(line)

        logging.info(f"Finished tert min coverage {min_coverage} count analysis")

    @classmethod
    def write_tert_cumulative_coverage_file(
            cls,
            sample_with_coverage_info_list: List[Tuple[str, CoverageInfo]],
            tert_site: Interval,
            output_dir: Path,
    ) -> None:
        logging.info("Started tert cumulative coverage analysis")
        tert_site_cumulative_coverage_output_file = f"{output_dir}/tert_site_cumulative_coverage.tsv"
        if Path(tert_site_cumulative_coverage_output_file).exists():
            error_msg = (
                f"At least one of the output files already exists: "
                f"{tert_site_cumulative_coverage_output_file}."
            )
            raise FileExistsError(error_msg)

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        with open(tert_site_cumulative_coverage_output_file, "w") as exon_o:
            base_header_elements = ["chrom", "start_position", "end_position"]
            header_elements = base_header_elements + [sample for sample, _ in sample_with_coverage_info_list]
            header = "\t".join(header_elements) + "\n"
            exon_o.write(header)

            base_line_data = [
                tert_site.chromosome,
                str(tert_site.start_position),
                str(tert_site.end_position),
            ]
            samples_line_data = [
                str(coverage_info.get_cumulative_coverage(tert_site))
                for sample, coverage_info in sample_with_coverage_info_list
            ]
            line_data = base_line_data + samples_line_data
            line = "\t".join(line_data) + "\n"
            exon_o.write(line)

        logging.info("Finished tert cumulative coverage analysis")

    @classmethod
    def __write_exon_cumulative_coverage_output_file(
            cls,
            sample_to_exon_to_cumulative_coverage: Dict[str, Dict[Exon, int]],
            samples: List[str],
            exons: Tuple[Exon, ...],
            exon_coverage_output_file: str,
    ) -> None:
        sorted_exons = cls.__get_sorted_exons(exons)
        with open(exon_coverage_output_file, "w") as exon_o:
            base_header_elements = ["chrom", "gene", "gene_id", "exon_id", "start_position", "end_position"]
            header_elements = base_header_elements + samples
            header = "\t".join(header_elements) + "\n"
            exon_o.write(header)

            for exon in sorted_exons:
                base_line_data = [
                    exon.interval.chromosome,
                    exon.gene,
                    exon.gene_ensembl_id,
                    exon.exon_ensembl_id,
                    str(exon.interval.start_position),
                    str(exon.interval.end_position),
                ]
                samples_line_data = [str(sample_to_exon_to_cumulative_coverage[sample][exon]) for sample in samples]
                line_data = base_line_data + samples_line_data
                line = "\t".join(line_data) + "\n"
                exon_o.write(line)

    @classmethod
    def __write_exon_min_coverage_count_output_file(
            cls,
            sample_to_exon_to_min_coverage_count: Dict[str, Dict[Exon, int]],
            samples: List[str],
            exons: Tuple[Exon, ...],
            exon_coverage_output_file: str,
    ) -> None:
        sorted_exons = cls.__get_sorted_exons(exons)
        with open(exon_coverage_output_file, "w") as exon_o:
            base_header_elements = ["chrom", "gene", "gene_id", "exon_id", "start_position", "end_position"]
            header_elements = base_header_elements + samples
            header = "\t".join(header_elements) + "\n"
            exon_o.write(header)

            for exon in sorted_exons:
                base_line_data = [
                    exon.interval.chromosome,
                    exon.gene,
                    exon.gene_ensembl_id,
                    exon.exon_ensembl_id,
                    str(exon.interval.start_position),
                    str(exon.interval.end_position),
                ]
                samples_line_data = [str(sample_to_exon_to_min_coverage_count[sample][exon]) for sample in samples]
                line_data = base_line_data + samples_line_data
                line = "\t".join(line_data) + "\n"
                exon_o.write(line)

    @classmethod
    def __write_gene_cumulative_coverage_output_file(
            cls,
            sample_to_gene_to_cumulative_coverage: Dict[str, Dict[str, int]],
            samples: List[str],
            exons: Tuple[Exon, ...],
            gene_coverage_output_file: str,
    ) -> None:
        sorted_genes = cls.__get_sorted_genes(exons)
        gene_to_total_exons_length = cls.__get_gene_to_total_exons_length(exons)
        with open(gene_coverage_output_file, "w") as gene_o:
            base_header_elements = ["gene", "total_exons_length"]
            header_elements = base_header_elements + samples
            header = "\t".join(header_elements) + "\n"
            gene_o.write(header)
            for gene in sorted_genes:
                base_line_data = [gene, str(gene_to_total_exons_length[gene])]
                samples_line_data = [str(sample_to_gene_to_cumulative_coverage[sample][gene]) for sample in samples]
                line_data = base_line_data + samples_line_data
                line = "\t".join(line_data) + "\n"
                gene_o.write(line)

    @classmethod
    def __write_gene_min_coverage_count_output_file(
            cls,
            sample_to_gene_to_min_coverage_count: Dict[str, Dict[str, int]],
            samples: List[str],
            exons: Tuple[Exon, ...],
            gene_coverage_output_file: str,
    ) -> None:
        sorted_genes = cls.__get_sorted_genes(exons)
        gene_to_total_exons_length = cls.__get_gene_to_total_exons_length(exons)
        with open(gene_coverage_output_file, "w") as gene_o:
            base_header_elements = ["gene", "total_exons_length"]
            header_elements = base_header_elements + samples
            header = "\t".join(header_elements) + "\n"
            gene_o.write(header)
            for gene in sorted_genes:
                base_line_data = [gene, str(gene_to_total_exons_length[gene])]
                samples_line_data = [str(sample_to_gene_to_min_coverage_count[sample][gene]) for sample in samples]
                line_data = base_line_data + samples_line_data
                line = "\t".join(line_data) + "\n"
                gene_o.write(line)

    @classmethod
    def __sum_exon_data_to_gene_data(
            cls,
            sample_to_exon_to_integer: Dict[str, Dict[Exon, int]]
    ) -> Dict[str, Dict[str, int]]:
        sample_to_gene_to_integer: Dict[str, Dict[str, int]] = {}
        for sample, exon_to_integer in sample_to_exon_to_integer.items():
            gene_to_integer: DefaultDict[str, int] = defaultdict(int)
            for exon, integer in exon_to_integer.items():
                gene_to_integer[exon.gene] += integer

            sample_to_gene_to_integer[sample] = dict(gene_to_integer)
        return sample_to_gene_to_integer

    @classmethod
    def __get_sorted_genes(cls, exons: Tuple[Exon, ...]) -> List[str]:
        return sorted(list({exon.gene for exon in exons}))

    @classmethod
    def __get_gene_to_total_exons_length(cls, exons: Tuple[Exon, ...]) -> Dict[str, int]:
        gene_to_total_exons_length: DefaultDict[str, int] = defaultdict(int)
        for exon in exons:
            gene_to_total_exons_length[exon.gene] += exon.interval.end_position - exon.interval.start_position + 1
        return gene_to_total_exons_length

    @classmethod
    def __get_sorted_exons(cls, exons: Tuple[Exon, ...]) -> List[Exon]:
        return sorted(list(exons), key=lambda x: (x.interval.chromosome, x.interval.start_position))
