import argparse
import json
import os
import subprocess
import sys
from shutil import copyfile
from typing import List, Set, Tuple

from config.panel import Panel
from pgx_analysis import PgxAnalyser, PgxAnalysis
from pgx_reporter import HaplotypeReporter, GenotypeReporter
from vcf_reader import VcfReader


def main(vcf: str, sample_t_id: str, sample_r_id: str, version: str, panel_path: str, outputdir: str,
         recreate_bed: bool, vcftools: str, sourcedir: str) -> None:
    """ Run pharmacogenomics analysis on sample """
    print("\n[INFO] ## START PHARMACOGENOMICS ANALYSIS")

    # Check if output dir exists, create if it does not
    if not os.path.exists(outputdir):
        try:
            os.makedirs(outputdir)
        except FileExistsError:
            # Directory already exists
            pass

    # Get configuration
    panel = load_panel(panel_path)
    bed_file = get_bed_file(panel_path, recreate_bed, panel, sourcedir)

    if panel.is_empty():
        raise ValueError("No panel is given, so no analysis can be performed.")

    # Get data for patient
    filtered_vcf, filtered_vcf_log = get_filtered_vcf(vcf, bed_file, sample_r_id, sample_t_id, outputdir, vcftools)
    grch37_call_data = VcfReader.get_grch37_call_data(filtered_vcf, panel)

    # Compute output from input data
    pgx_analysis = PgxAnalyser.create_pgx_analysis(grch37_call_data, panel)

    # Output
    # TODO: maybe change file names of tsv files to .tsv
    print_calls_to_file(pgx_analysis, outputdir, sample_t_id)
    print_genotypes_to_file(pgx_analysis, panel, outputdir, sample_t_id, panel_path, version)
    # Also copy the bed-filtered VCF file for research purposes
    copy_filtered_vcf_file(filtered_vcf, outputdir, sample_t_id)

    # Clean up
    if os.path.exists(filtered_vcf):
        if os.path.exists(filtered_vcf):
            os.remove(filtered_vcf)
            print(f"[INFO] {filtered_vcf} removed.")
        if os.path.exists(filtered_vcf_log):
            os.remove(filtered_vcf_log)
            print(f"[INFO] {filtered_vcf_log} removed.")

    # TODO: add genes CYP2D6, CYP3A4, CYP3A5

    print("[INFO] ## PHARMACOGENOMICS ANALYSIS FINISHED\n")


def get_filtered_vcf(vcf: str, bed_file: str, sample_r_id: str,
                     sample_t_id: str, outputdir: str, vcftools: str) -> Tuple[str, str]:
    filtered_vcf_prefix = f"{outputdir}/{sample_t_id}_PGx"
    filtered_vcf = f"{filtered_vcf_prefix}.recode.vcf"
    filtered_vcf_log = f"{filtered_vcf_prefix}.log"
    # Check if output vcf does not already exist
    if os.path.exists(filtered_vcf) or os.path.exists(filtered_vcf_log):
        raise IOError(f"Temporary VCF file {filtered_vcf} already exists. Exiting.")

    subprocess.run([vcftools, '--gzvcf', vcf, '--bed', bed_file, '--out', filtered_vcf_prefix,
                    '--indv', sample_r_id, '--recode', '--recode-INFO-all'])
    print("[INFO] Subprocess completed.")

    if not os.path.exists(filtered_vcf) or not os.path.exists(filtered_vcf_log):
        error_msg = f"Filtered vcf or log of filtering does not exist"
        raise FileNotFoundError(error_msg)

    return filtered_vcf, filtered_vcf_log


def load_panel(panel_path: str) -> Panel:
    """ Load manually annotated JSON panel file """
    try:
        with open(panel_path, 'r+', encoding='utf-8') as json_file:
            data = json.load(json_file)
            return Panel.from_json(data)
    except IOError:
        raise FileNotFoundError(f"Panel file {panel_path} not found or cannot be opened.")


def get_bed_file(panel_path: str, recreate_bed: bool, panel: Panel, sourcedir: str) -> str:
    bed_file = replace_file_extension_of_path(panel_path, "bed")
    if recreate_bed:
        create_bed_file(panel.get_genes(), panel_path, sourcedir, bed_file)
    if not os.path.exists(bed_file):
        raise FileNotFoundError(
            "Could not locate bed-file. "
            "Could it be that it should be (re)created? "
            "Retry running with --recreate_bed."
        )
    return bed_file


def create_bed_file(genes_in_panel: Set[str], panel_path: str, sourcedir: str, bed_path: str) -> None:
    """ Generate bed file from gene panel and save as panel_path.bed """
    print("[INFO] Recreating bed-file...")
    header = (
        f'track name="{panel_path}" description="Bed file generated from {panel_path} with HMF_PGx main.py"\n'
    )
    bed_regions = []  # chrom, start, end, gene
    covered = []
    transcripts = open(f"{sourcedir}/all_genes.37.tsv", 'r')
    for line in transcripts:
        split_line = line.rstrip().split("\t")
        if split_line[4] in genes_in_panel and split_line[4] not in covered:
            bed_regions.append([split_line[0], split_line[1], split_line[2], split_line[4]])
            covered.append(split_line[4])
    if set(covered) != genes_in_panel:
        error_msg = (
            f"Missing genes from the gene panel in the transcript list. Please check:\n"
            f"Covered:\n"
            f"{covered}\n"
            f"Original gene panel:\n"
            f"{genes_in_panel}"
        )
        raise ValueError(error_msg)

    with open(bed_path, 'w') as bed:
        bed.write(header)
        for entry in bed_regions:
            bed.write("\t".join(entry) + "\n")

    print(f"[INFO] Created {bed_path}")


def print_calls_to_file(pgx_analysis: PgxAnalysis, outputdir: str, sample_t_id: str) -> None:
    calls_file = f"{outputdir}/{sample_t_id}_calls.txt"
    if os.path.exists(calls_file):
        raise IOError(f"Calls output file {calls_file} already exists. Exiting.")
    with open(calls_file, 'w') as f:
        f.write(GenotypeReporter.get_calls_tsv_text(pgx_analysis))
    if not os.path.exists(calls_file):
        raise FileNotFoundError(f"Failed to write calls output file {calls_file}")


def print_genotypes_to_file(pgx_analysis: PgxAnalysis, panel: Panel, outputdir: str, sample_t_id: str,
        panel_path: str, version: str) -> None:
    genotype_file = f"{outputdir}/{sample_t_id}_genotype.txt"
    if os.path.exists(genotype_file):
        raise IOError(f"Genotype output file {genotype_file} already exists. Exiting.")
    with open(genotype_file, 'w') as f:
        f.write(HaplotypeReporter.get_genotype_tsv_text(pgx_analysis, panel, panel_path, version))
    if not os.path.exists(genotype_file):
        raise FileNotFoundError(f"Failed to write calls output file {genotype_file}")


def copy_filtered_vcf_file(filtered_vcf: str, outputdir: str, sample_t_id: str) -> None:
    destination_path = f"{outputdir}/{sample_t_id}_PGx.vcf"
    if os.path.exists(destination_path):
        raise IOError(f"Copied filtered vcf file {destination_path} already exists. Exiting.")
    copyfile(filtered_vcf, destination_path)
    if not os.path.exists(destination_path):
        raise FileNotFoundError(f"Failed to copy filtered vcf file {destination_path}")


def parse_args(sys_args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=('Run pharmacogenomics panel on germline VCF file. The pharmacogenomic annotations are done on '
                     'GRCh38, so in the output both reference genome output is given where possible.')
    )
    parser.add_argument('vcf', type=str, help='VCF file to use for pharmacogenomics analysis')
    parser.add_argument('sample_t_id', type=str, help='The sample ID of the tumor')
    parser.add_argument('sample_r_id', type=str, help='The sample ID of the normal')
    parser.add_argument('version', type=str, help='The version of the tool')
    parser.add_argument('outputdir', type=str, help='Directory to store output of pharmacogenomic analysis')
    parser.add_argument('panel', type=str, help='Json file with the panel variants')
    parser.add_argument('vcftools', type=str, default='vcftools', help="Path to vcftools > 0.1.14 if not in $PATH")
    parser.add_argument(
        '--recreate_bed', default=False, action='store_true',
        help='Recreate bed-file from JSON files. If false, the panel file with extension .bed is searched for.'
    )
    parser.add_argument('--sourcedir', type=str, default='data', help="Optional path to location of source files")
    return parser.parse_args(sys_args)


def replace_file_extension_of_path(path: str, new_file_extension: str) -> str:
    split_path = path.split(".")
    new_path = ".".join(split_path[0:-1]) + "." + new_file_extension
    return new_path


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    main(args.vcf, args.sample_t_id, args.sample_r_id, args.version, args.panel,
         args.outputdir, args.recreate_bed, args.vcftools, args.sourcedir)
