#!/usr/bin/env python

import os
import gzip
import re
import gc
import argparse
import pandas as pd
from cyvcf2 import VCF

# =============================================================================
# Utility Functions
# =============================================================================

def safe_print(msg):
    """Utility for consistent logging."""
    print(msg)


# =============================================================================
# VCF Utility Functions
# =============================================================================

def parse_vcf_line(line):
    """
    Parse a VCF entry line into a structured dictionary.
    """
    fields = line.strip().split("\t")
    if len(fields) < 8:
        return "malformed_vcf_entry"

    parsed = {
        "CHROM": fields[0],
        "POS": fields[1],
        "ID": fields[2],
        "REF": fields[3],
        "ALT": fields[4],
        "QUAL": fields[5],
        "FILTER": fields[6],
        "INFO": fields[7],
        "FORMAT": fields[8] if len(fields) > 8 else "",
        "SAMPLES": fields[9:] if len(fields) > 9 else []
    }

    # Parse INFO field into a dictionary.
    parsed["INFO_PARSED"] = {}
    for info_field in parsed["INFO"].split(";"):
        if "=" in info_field:
            key, value = info_field.split("=", 1)
            parsed["INFO_PARSED"][key] = value

    return parsed


def format_vcf_entry(vcf_entry):
    """
    Format a VCF entry with newlines after each field for better readability.
    """
    if not vcf_entry or vcf_entry == "no_vcf_entry":
        return "no_vcf_entry"

    try:
        # If entry is already a string from VCF iteration.
        if isinstance(vcf_entry, str):
            fields = vcf_entry.strip().split("\t")
            return "\n".join(fields)
        # Otherwise, if it is a dict-like object.
        elif isinstance(vcf_entry, dict):
            fields = [
                f"CHROM={vcf_entry.get('CHROM', '')}",
                f"POS={vcf_entry.get('POS', '')}",
                f"ID={vcf_entry.get('ID', '')}",
                f"REF={vcf_entry.get('REF', '')}",
                f"ALT={','.join(vcf_entry.get('ALT', []))}",
                f"QUAL={vcf_entry.get('QUAL', '')}",
                f"FILTER={vcf_entry.get('FILTER', '')}",
                f"INFO={vcf_entry.get('INFO', '')}"
            ]
            return "\n".join(fields)
        else:
            return str(vcf_entry)
    except Exception as e:
        safe_print(f"Error formatting VCF entry: {e}")
        return "error_formatting_entry"


def search_vcf_by_chrom_pos(vcf_file_path, chrom, pos):
    """
    Search for a variant at a specific chromosome and position in a VCF file.
    Returns a formatted VCF entry string if found.
    """
    safe_print(f"Searching for variant at {chrom}:{pos} in VCF: {vcf_file_path}...")

    if not os.path.exists(vcf_file_path):
        safe_print(f"VCF file not found: {vcf_file_path}")
        return "no_vcf_file"

    try:
        vcf = VCF(vcf_file_path)
        region = f"{chrom}:{pos}-{pos}"
        for variant in vcf(region):
            if variant.CHROM == chrom and variant.POS == int(pos):
                safe_print(f"Found variant at {chrom}:{pos} in VCF.")
                return format_vcf_entry(str(variant))
    except Exception as e:
        safe_print(f"Error reading VCF file: {e}")
        return f"error_reading_vcf: {e}"
    finally:
        try:
            vcf.close()
        except Exception:
            pass
        gc.collect()

    safe_print(f"No variant found at {chrom}:{pos} in VCF.")
    return "no_vcf_entry"


def search_gripss_or_esvee_vcf(vcf_file_path, chrom, pos):
    """
    Search for a variant at a specific chromosome and position in a VCF file
    for GRIPSS or ESVEE.
    """
    safe_print(f"Searching for variant at {chrom}:{pos} in VCF: {vcf_file_path}...")

    if not os.path.exists(vcf_file_path):
        safe_print(f"VCF file not found: {vcf_file_path}")
        return "no_vcf_file"

    try:
        vcf = VCF(vcf_file_path)
        region = f"{chrom}:{pos}-{pos}"
        for variant in vcf(region):
            if variant.CHROM == chrom and variant.POS == int(pos):
                safe_print(f"Found variant at {chrom}:{pos} in VCF.")
                return format_vcf_entry(str(variant))
    except Exception as e:
        safe_print(f"Error reading VCF file: {e}")
        return f"error_reading_vcf: {e}"
    finally:
        try:
            vcf.close()
        except Exception:
            pass

    safe_print(f"No variant found at {chrom}:{pos} in VCF.")
    return "no_vcf_entry"


# =============================================================================
# COMPAR File Parsing Functions
# =============================================================================

def find_compar_files(compar_dir):
    """
    Find the .driver.tsv, .somatic_variant.tsv, and .disruption.tsv files in the COMPAR output directory.
    """
    safe_print(f"Searching for .driver.tsv, .somatic_variant.tsv, and .disruption.tsv files in {compar_dir}...")
    driver_file = None
    somatic_variant_file = None
    disruption_file = None

    for root, _, files in os.walk(compar_dir):
        for file in files:
            if file.endswith(".driver.tsv"):
                driver_file = os.path.join(root, file)
                safe_print(f"Found driver file: {driver_file}")
            elif file.endswith(".somatic_variant.tsv"):
                somatic_variant_file = os.path.join(root, file)
                safe_print(f"Found somatic_variant file: {somatic_variant_file}")
            elif file.endswith(".disruption.tsv"):
                disruption_file = os.path.join(root, file)
                safe_print(f"Found disruption file: {disruption_file}")

    if not driver_file or not somatic_variant_file or not disruption_file:
        raise FileNotFoundError("Could not find required .driver.tsv, .somatic_variant.tsv, or .disruption.tsv files in the COMPAR directory.")

    return driver_file, somatic_variant_file, disruption_file


def parse_compar_somatic(compar_somatic_file, sample_id, mismatch_type, gene_name):
    """
    Extract the Differences field from compar.somatic_variant.tsv for a given sample, mismatch type, and gene.
    """
    safe_print(f"Extracting Differences for SampleId={sample_id}, MismatchType={mismatch_type}, Gene={gene_name}...")
    if not os.path.exists(compar_somatic_file):
        safe_print("compar.somatic_variant.tsv not found.")
        return "no_compar_entry"

    compar_somatic_data = pd.read_csv(compar_somatic_file, sep="\t")
    if mismatch_type == "REF_ONLY":
        result = compar_somatic_data[
            (compar_somatic_data['SampleId'] == sample_id) &
            (compar_somatic_data['MismatchType'] == mismatch_type) &
            (compar_somatic_data['RefGene'] == gene_name)
        ]
    elif mismatch_type == "NEW_ONLY":
        result = compar_somatic_data[
            (compar_somatic_data['SampleId'] == sample_id) &
            (compar_somatic_data['MismatchType'] == mismatch_type) &
            (compar_somatic_data['NewGene'] == gene_name)
        ]
    else:  # VALUE
        result = compar_somatic_data[
            (compar_somatic_data['SampleId'] == sample_id) &
            (compar_somatic_data['MismatchType'] == mismatch_type) &
            ((compar_somatic_data['RefGene'] == gene_name) | (compar_somatic_data['NewGene'] == gene_name))
        ]

    if result.empty:
        safe_print("No matching entry found in somatic_variant.tsv.")
        return "no_compar_entry"

    differences = result['Differences'].values[0] if 'Differences' in result.columns else "no_differences_field"
    return differences


def parse_driver_file(driver_file_path):
    """
    Parse the compar.driver.tsv file to extract relevant rows for MUTATION, DISRUPTION, and CNV entries.
    """
    safe_print(f"Parsing driver file: {driver_file_path}")
    driver_data = pd.read_csv(driver_file_path, sep="\t")
    safe_print(f"Total rows in driver file: {len(driver_data)}")
    
    mutation_rows = driver_data[driver_data['Key'].str.startswith("MUTATION", na=False)]
    disruption_rows = driver_data[driver_data['Key'].str.startswith("DISRUPTION", na=False)]
    cnv_rows = driver_data[
        driver_data['Key'].str.startswith("DEL", na=False) |
        driver_data['Key'].str.startswith("AMP", na=False) |
        driver_data['Key'].str.startswith("PARTIAL_AMP", na=False) |
        driver_data['Key'].str.startswith("PARTIAL_DEL", na=False)
    ]
    
    safe_print(f"Found {len(mutation_rows)} MUTATION entries.")
    safe_print(f"Found {len(disruption_rows)} DISRUPTION entries.")
    safe_print(f"Found {len(cnv_rows)} CNV entries (DEL, AMP, PARTIAL_AMP, PARTIAL_DEL).")

    return mutation_rows, disruption_rows, cnv_rows


def parse_disruption_file(disruption_file_path, sample_id, mismatch_type, gene_name):
    """
    Parse the .disruption.tsv file to extract relevant rows for DISRUPTION entries.
    """
    safe_print(f"Parsing disruption file: {disruption_file_path} for SampleId={sample_id}, Gene={gene_name}...")
    
    disruption_data = pd.read_csv(disruption_file_path, sep="\t")
    
    def extract_gene_name(key):
        return re.match(r"^[^\s]+", key).group(0)

    disruption_data['ExtractedGeneName'] = disruption_data['Key'].apply(extract_gene_name)

    result = disruption_data[
        (disruption_data['SampleId'] == sample_id) &
        (disruption_data['MismatchType'] == mismatch_type) &
        (disruption_data['ExtractedGeneName'] == gene_name)
    ]

    disruption_data.drop(columns=['ExtractedGeneName'], inplace=True)

    if result.empty:
        safe_print("No matching entry found in disruption.tsv.")
        return "no_disruption_entry", None, None, None, None

    differences = result['Differences'].values[0] if 'Differences' in result.columns else "no_differences_field"
    
    ref_breakend_info = result['RefBreakendInfo'].values[0] if pd.notna(result['RefBreakendInfo'].values[0]) else None
    new_breakend_info = result['NewBreakendInfo'].values[0] if pd.notna(result['NewBreakendInfo'].values[0]) else None

    chrom, pos, chrom_new, pos_new = None, None, None, None

    if ref_breakend_info:
        parts = ref_breakend_info.split(" ")
        if len(parts) > 1:
            chrom, pos = parts[1].split(":")
    if new_breakend_info:
        parts = new_breakend_info.split(" ")
        if len(parts) > 1:
            chrom_new, pos_new = parts[1].split(":")
    elif chrom and pos:
        chrom_new, pos_new = chrom, pos

    return differences, chrom, int(pos) if pos else None, chrom_new, int(pos_new) if pos_new else None


# =============================================================================
# LINX & CNV Utility Functions
# =============================================================================

def search_linx_breakend(linx_dir, sample_id, gene_name):
    """
    Search for breakend information in the .breakend.tsv file located in the linx directory.
    """
    safe_print(f"Searching linx directory for breakend file for SampleId={sample_id}, Gene={gene_name}...")
    linx_file_path = os.path.join(linx_dir, sample_id, "linx", f"{sample_id}.linx.breakend.tsv")

    if not os.path.exists(linx_file_path):
        safe_print(f"No linx breakend file found for SampleId={sample_id}")
        return "no_linx_breakend_entry"

    linx_data = pd.read_csv(linx_file_path, sep="\t")
    result = linx_data[linx_data['gene'] == gene_name]

    if result.empty:
        safe_print(f"No breakend entry found in linx file for Gene={gene_name}.")
        return "no_linx_breakend_entry"

    return result.to_csv(sep="\t", index=False)


def reformat_linx_output(linx_data):
    """
    Reformat linx output to display values for each header on new lines.
    """
    lines = linx_data.strip().split("\n")
    header = lines[0].split("\t")
    values = [line.split("\t") for line in lines[1:]]
    reformatted_output = []
    for i, head in enumerate(header):
        values_list = [value[i] for value in values]
        reformatted_output.append(f"{head} " + " ".join(values_list))
    return "\n".join(reformatted_output)


# =============================================================================
# Row Processing Functions for MUTATION, DISRUPTION, and CNV
# =============================================================================

def process_row(row, ref_dir, new_dir):
    """
    Process a single row from the somatic_variant.tsv file and return the results.
    Uses coordinates as given without any liftover conversion.
    """
    sample_id = row['SampleId']
    mismatch_type = row['MismatchType']
    gene_name = row['RefGene'] if pd.notna(row['RefGene']) else row['NewGene']
    chrom, pos = row['Key'].split(' ')[0].split(':')

    # Use same coordinates for both REF and NEW
    chrom_new, pos_new = chrom, pos
    ref_purple_vcf_entry = "no_vcf_entry"
    new_purple_vcf_entry = "no_vcf_entry"
    ref_sage_vcf_entry = "no_vcf_entry"
    new_sage_vcf_entry = "no_vcf_entry"

    # Resolve sample directories in ref_dir and new_dir.
    directory_sample_id_ref = next((dir_name for dir_name in os.listdir(ref_dir) if sample_id in dir_name), None)
    if not directory_sample_id_ref:
        safe_print(f"Sample ID {sample_id} not found in {ref_dir}")
    directory_sample_id_new = next((dir_name for dir_name in os.listdir(new_dir) if sample_id in dir_name), None)
    if not directory_sample_id_new:
        safe_print(f"Sample ID {sample_id} not found in {new_dir}")

    try:
        # Process REF Purple VCF.
        ref_purple_vcf_file = os.path.join(ref_dir, directory_sample_id_ref, "purple", f"{sample_id}.purple.somatic.vcf.gz")
        ref_variant = search_vcf_by_chrom_pos(ref_purple_vcf_file, chrom, int(pos))
        if ref_variant:
            ref_purple_vcf_entry = format_vcf_entry(str(ref_variant))

        # Process NEW Purple VCF.
        new_purple_vcf_file = os.path.join(new_dir, directory_sample_id_new, "purple", f"{sample_id}.purple.somatic.vcf.gz")
        new_variant = search_vcf_by_chrom_pos(new_purple_vcf_file, chrom_new, int(pos_new))
        if new_variant:
            new_purple_vcf_entry = format_vcf_entry(str(new_variant))

        # Process SAGE VCFs.
        if chrom and pos:
            ref_sage_vcf_file = os.path.join(ref_dir, directory_sample_id_ref, "sage_somatic", f"{sample_id}.sage.somatic.vcf.gz")
            ref_sage_vcf_entry = search_vcf_by_chrom_pos(ref_sage_vcf_file, chrom, pos)
        if chrom_new and pos_new:
            new_sage_vcf_file = os.path.join(new_dir, directory_sample_id_new, "sage_somatic", f"{sample_id}.sage.somatic.vcf.gz")
            new_sage_vcf_entry = search_vcf_by_chrom_pos(new_sage_vcf_file, chrom_new, pos_new)
    except Exception as e:
        safe_print(f"Error processing row for SampleId={sample_id}, Gene={gene_name}: {e}")

    safe_print(f"Completed processing row for SampleId={sample_id}, Gene={gene_name}")

    return {
        "SampleId": sample_id,
        "MismatchType": mismatch_type,
        "Gene": gene_name,
        "Chrom": chrom,
        "Pos": pos,
        "ChromNew": chrom_new,
        "PosNew": pos_new,
        "RefPurpleVCF": ref_purple_vcf_entry,
        "NewPurpleVCF": new_purple_vcf_entry,
        "RefSageVCF": ref_sage_vcf_entry,
        "NewSageVCF": new_sage_vcf_entry,
    }


def process_disruption_row(row, disruption_file_path, linx_dir, ref_dir, new_dir, v6=False):
    """
    Process a single DISRUPTION row and return the results.
    Retrieves values from compar output files, linx breakend files, and VCFs.
    """
    sample_id = row['SampleId']
    mismatch_type = row['MismatchType']
    gene_name = row['Key'].split("_")[1]

    compar_differences, chrom, pos, chrom_new, pos_new = parse_disruption_file(
        disruption_file_path, sample_id, mismatch_type, gene_name
    )
    compar_differences = str(compar_differences) if pd.notna(compar_differences) else ""

    ref_linx_breakend_entry = search_linx_breakend(ref_dir, sample_id, gene_name)
    new_linx_breakend_entry = search_linx_breakend(new_dir, sample_id, gene_name)
    
    if ref_linx_breakend_entry != "no_linx_breakend_entry":
        ref_linx_breakend_entry = reformat_linx_output(ref_linx_breakend_entry)
    if new_linx_breakend_entry != "no_linx_breakend_entry":
        new_linx_breakend_entry = reformat_linx_output(new_linx_breakend_entry)

    ref_vcf_entry, new_vcf_entry = "no_vcf_entry", "no_vcf_entry"

    if v6:
        ref_gripss_vcf_file = os.path.join(ref_dir, sample_id, "esvee", f"{sample_id}.esvee.somatic.vcf.gz")
    else:
        ref_gripss_vcf_file = os.path.join(ref_dir, sample_id, "gripss_somatic", f"{sample_id}.gripss.somatic.vcf.gz")
    new_esvee_vcf_file = os.path.join(new_dir, sample_id, "esvee", f"{sample_id}.esvee.somatic.vcf.gz")
    
    if chrom and pos:
        ref_vcf_entry = search_gripss_or_esvee_vcf(ref_gripss_vcf_file, chrom, pos)
    if chrom_new and pos_new:
        new_vcf_entry = search_gripss_or_esvee_vcf(new_esvee_vcf_file, chrom_new, pos_new)

    safe_print(f"Completed processing DISRUPTION row for SampleId={sample_id}, Gene={gene_name}")

    return {
        "SampleId": sample_id,
        "MismatchType": mismatch_type,
        "Gene": gene_name,
        "ComparDifferences": compar_differences,
        "RefLinxBreakend": ref_linx_breakend_entry,
        "NewLinxBreakend": new_linx_breakend_entry,
        "RefVcf": ref_vcf_entry,
        "NewVcf": new_vcf_entry,
    }


def process_cnv_row(row, ref_dir, new_dir):
    """
    Process a single CNV row and return the results.
    Retrieves values from REF and NEW purple.cnv.gene.tsv files and formats them.
    """
    sample_id = row['SampleId']
    mismatch_type = row['MismatchType']

    key_parts = row['Key'].split("_")
    if key_parts[0] == "PARTIAL":
        cnv_type = f"{key_parts[0]}_{key_parts[1]}"
        gene_name = key_parts[2]
    else:
        cnv_type = key_parts[0]
        gene_name = key_parts[1]

    ref_cnv_entry = "no_cnv_entry"
    new_cnv_entry = "no_cnv_entry"

    ref_cnv_file = os.path.join(ref_dir, sample_id, "purple", f"{sample_id}.purple.cnv.gene.tsv")
    if os.path.exists(ref_cnv_file):
        ref_cnv_data = pd.read_csv(ref_cnv_file, sep="\t")
        ref_cnv_result = ref_cnv_data[ref_cnv_data['gene'] == gene_name]
        if not ref_cnv_result.empty:
            ref_cnv_entry = ref_cnv_result.to_csv(sep="\t", index=False)
            ref_cnv_entry = reformat_linx_output(ref_cnv_entry)

    new_cnv_file = os.path.join(new_dir, sample_id, "purple", f"{sample_id}.purple.cnv.gene.tsv")
    if os.path.exists(new_cnv_file):
        new_cnv_data = pd.read_csv(new_cnv_file, sep="\t")
        new_cnv_result = new_cnv_data[new_cnv_data['gene'] == gene_name]
        if not new_cnv_result.empty:
            new_cnv_entry = new_cnv_result.to_csv(sep="\t", index=False)
            new_cnv_entry = reformat_linx_output(new_cnv_entry)

    safe_print(f"Completed processing CNV row for SampleId={sample_id}, Gene={gene_name}")

    return {
        "SampleId": sample_id,
        "MismatchType": mismatch_type,
        "Gene": gene_name,
        "TypeCNV": cnv_type,
        "RefPurpleCNV": ref_cnv_entry,
        "NewPurpleCNV": new_cnv_entry,
    }


# =============================================================================
# Main Collection Functions
# =============================================================================

def collect_results(compar_dir, ref_dir, new_dir, output_file):
    """
    Main function to collect mutation (somatic) results and save them in a structured TSV file.
    """
    safe_print("Initializing collect_results...")

    _, somatic_variant_file, _ = find_compar_files(compar_dir)
    somatic_variant_data = pd.read_csv(somatic_variant_file, sep="\t")
    safe_print(f"Found {len(somatic_variant_data)} entries in somatic_variant file to process.")

    results = []
    for idx, row in somatic_variant_data.iterrows():
        safe_print(f"Processing row {idx + 1}/{len(somatic_variant_data)}: SampleId={row['SampleId']}, Gene={row['Key']}")
        try:
            result = process_row(row, ref_dir, new_dir)
            safe_print(f"Successfully processed row for SampleId={result['SampleId']}, Gene={result['Gene']}")
            results.append(result)
        except Exception as e:
            safe_print(f"Error processing row for SampleId={row['SampleId']} - {e}")

    if results:
        results_df = pd.DataFrame(results)
    else:
        safe_print("No mutation results to write. Creating an empty file.")
        results_df = pd.DataFrame(columns=["SampleId", "MismatchType", "Gene", "Chrom", "Pos",
                                           "ChromNew", "PosNew", "RefPurpleVCF", "NewPurpleVCF",
                                           "RefSageVCF", "NewSageVCF"])

    results_df.to_csv(output_file, sep="\t", index=False)
    safe_print("Results successfully written to output file.")


def collect_disruption_results(compar_dir, ref_dir, new_dir, output_file, v6=False):
    """
    Main function to collect results for DISRUPTION entries and save them in a structured TSV file.
    """
    safe_print("Initializing collect_disruption_results...")

    driver_file, _, disruption_file = find_compar_files(compar_dir)
    _, disruption_entries, _ = parse_driver_file(driver_file)
    safe_print(f"Found {len(disruption_entries)} disruption entries to process.")

    results = []
    for idx, row in disruption_entries.iterrows():
        safe_print(f"Processing DISRUPTION row {idx + 1}/{len(disruption_entries)}: SampleId={row['SampleId']}, Gene={row['Key']}")
        try:
            result = process_disruption_row(
                row=row,
                disruption_file_path=disruption_file,
                linx_dir=ref_dir,  # Assumes LINX directory under ref_dir.
                ref_dir=ref_dir,
                new_dir=new_dir,
                v6=v6
            )
            safe_print(f"Successfully processed DISRUPTION row for SampleId={result['SampleId']}, Gene={result['Gene']}")
            results.append(result)
        except Exception as e:
            safe_print(f"Error processing DISRUPTION row for SampleId={row['SampleId']} - {e}")

    if results:
        results_df = pd.DataFrame(results)
    else:
        safe_print("No disruption results to write. Creating an empty file.")
        results_df = pd.DataFrame(columns=["SampleId", "MismatchType", "Gene", "ComparDifferences",
                                           "RefLinxBreakend", "NewLinxBreakend", "RefVcf", "NewVcf"])

    results_df.to_csv(output_file, sep="\t", index=False)
    safe_print("DISRUPTION results successfully written to output file.")


def collect_cnv_results(compar_dir, ref_dir, new_dir, output_file):
    """
    Main function to collect results for CNV entries and save them in a structured TSV file.
    """
    safe_print("Initializing collect_cnv_results...")

    driver_file, _, _ = find_compar_files(compar_dir)
    _, _, cnv_entries = parse_driver_file(driver_file)
    safe_print(f"Found {len(cnv_entries)} CNV entries to process.")

    results = []
    for idx, row in cnv_entries.iterrows():
        safe_print(f"Processing CNV row {idx + 1}/{len(cnv_entries)}: SampleId={row['SampleId']}, Gene={row['Key']}")
        try:
            result = process_cnv_row(row, ref_dir, new_dir)
            safe_print(f"Successfully processed CNV row for SampleId={result['SampleId']}, Gene={result['Gene']}")
            results.append(result)
        except Exception as e:
            safe_print(f"Error processing CNV row for SampleId={row['SampleId']} - {e}")

    if results:
        results_df = pd.DataFrame(results)
    else:
        safe_print("No CNV results to write. Creating an empty file.")
        results_df = pd.DataFrame(columns=["SampleId", "MismatchType", "Gene", "TypeCNV", "RefPurpleCNV", "NewPurpleCNV"])

    results_df.to_csv(output_file, sep="\t", index=False)
    safe_print("CNV results successfully written to output file.")


# =============================================================================
# Main Execution Block
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process COMPAR output to collect variant information.")
    parser.add_argument("--compar_dir", required=True, help="Path to the COMPAR output directory. (output_reported/output_all)")
    parser.add_argument("--ref_dir", required=True, help="Path to the reference pipeline output directory.")
    parser.add_argument("--new_dir", required=True, help="Path to the new pipeline output directory.")
    parser.add_argument("--output_file", required=True, help="Path to save the output TSV file.")
    parser.add_argument("--v6", action="store_true", help="Indicate to use v6 format for finding esvee instead of gripss")
    
    args = parser.parse_args()

    # Collect results for MUTATION entries.
    collect_results(
        compar_dir=args.compar_dir,
        ref_dir=args.ref_dir,
        new_dir=args.new_dir,
        output_file=args.output_file.replace(".tsv", "_mutation.tsv")
    )

    # Collect results for DISRUPTION entries.
    collect_disruption_results(
        compar_dir=args.compar_dir,
        ref_dir=args.ref_dir,
        new_dir=args.new_dir,
        output_file=args.output_file.replace(".tsv", "_disruption.tsv"),
        v6=args.v6
    )

    # Collect results for CNV entries.
    collect_cnv_results(
        compar_dir=args.compar_dir,
        ref_dir=args.ref_dir,
        new_dir=args.new_dir,
        output_file=args.output_file.replace(".tsv", "_cnv.tsv")
    )