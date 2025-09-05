import argparse
import gzip
import sys
import pandas as pd
import pysam
import mysql.connector
import csv

# Genotype codes
GENOTYPE_MAP = {
    "DO_NOT_MATCH": 0,
    "HOM_REF": 1,
    "HET": 2,
    "HOM_ALT": 3,
}

def derive_genotype_from_record(sample_record):
    """
    Derives a genotype code (0, 1, 2, 3) from a VCF sample record
    based on the specified rules.
    """
    try:
        # Get Allelic Depth (AD) and Total Depth (DP)
        total_coverage = sample_record['DP']
        ref_count, alt_count = sample_record['AD']
    except (KeyError, TypeError, IndexError):
        # If format fields are missing, we cannot determine the genotype
        return GENOTYPE_MAP["DO_NOT_MATCH"]
    
    print("Total coverage: ", total_coverage, "Ref_count:", ref_count, "Alt_count:", alt_count)

    # Rule: Minimum depth must be >= 10
    if total_coverage < 10:
        return GENOTYPE_MAP["DO_NOT_MATCH"]

    # Rule: To call homozygous, one allele count must equal total coverage.
    if ref_count == total_coverage and alt_count == 0:
        return GENOTYPE_MAP["HOM_REF"]

    if alt_count == total_coverage and ref_count == 0:
        return GENOTYPE_MAP["HOM_ALT"]

    # Rule: To call heterozygous, allele fraction must be between 0.4 and 0.65
    try:
        allele_fraction = alt_count / total_coverage
        if 0.4 <= allele_fraction <= 0.65:
            return GENOTYPE_MAP["HET"]
    except ZeroDivisionError:
        # Should not happen if depth >= 10, but good to be safe
        return GENOTYPE_MAP["DO_NOT_MATCH"]

    # If none of the above rules are met, it's a no-match site
    return GENOTYPE_MAP["DO_NOT_MATCH"]

def get_100_snp_loci(snp_check_vcf_path):
    """
    Reads the SNPCheck VCF to get the exact 100 loci (and their order)
    that are used for patient matching. Returns a list of (chrom, pos) tuples.
    """
    loci = []
    # The snpcheck vcf may or may not be gzipped.
    _open = gzip.open if snp_check_vcf_path.endswith('.gz') else open
    with _open(snp_check_vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            chrom = parts[0]
            pos = int(parts[1])
            loci.append((chrom, pos))
    if len(loci) != 100:
        print(f"Warning: Expected 100 loci from SNPCheck VCF, but found {len(loci)}.", file=sys.stderr)
    return loci

def derive_genotypes_from_vcf(vcf_path, ordered_loci):
    """
    Parses the AMBER SNP VCF and derives genotypes for the 100 ordered loci.
    Returns a list of genotypes in the same order as ordered_loci.
    """
    vcf_file = pysam.VariantFile(vcf_path)

    # Create a dictionary of genotypes from the VCF for quick lookup
    vcf_genotypes = {}
    for record in vcf_file.fetch():
        pos_tuple = (record.chrom, record.pos)
        genotype = derive_genotype_from_record(record.samples[0])
        print(genotype)
        vcf_genotypes[pos_tuple] = genotype

    # Build the final, ordered list of 100 genotypes
    final_genotypes = []
    for locus in ordered_loci:
        # If a SNP from the snpcheck file is not in the amber output, it's a DO_NOT_MATCH
        final_genotypes.append(vcf_genotypes.get(locus, GENOTYPE_MAP["DO_NOT_MATCH"]))

    return final_genotypes

def find_matches_in_db(new_sample_genotypes, db_config):
    """
    Connects to the DB and returns a list of all samples with >= 80 genotype matches.
    Each item in the list is a tuple: (matched_sample_id, match_count).
    """
    all_matches = []
    try:
        cnx = mysql.connector.connect(
            user=db_config['user'],
            password=db_config['pass'],
            host=db_config['host'],
            port=db_config['port'], # Using the port here
            database=db_config['name']
        )
        db_df = pd.read_sql("SELECT * FROM amberSample", cnx)
        cnx.close()
    except mysql.connector.Error as err:
        print(f"Error connecting to database: {err}", file=sys.stderr)
        return []

    # Iterate through each sample in the database
    for index, db_sample in db_df.iterrows():
        db_genotypes = db_sample[f"site1":f"site100"].values.astype(int)

        match_count = sum(g1 == g2 for g1, g2 in zip(new_sample_genotypes, db_genotypes))

        # If the match count is 80 or more, add it to our list of results
        if match_count >= 80:
            all_matches.append((db_sample['sampleId'], match_count))

    # Sort results by match count, descending
    all_matches.sort(key=lambda x: x[1], reverse=True)
    return all_matches


def main():
    parser = argparse.ArgumentParser(description="Derive AMBER genotypes and match against a database.")
    parser.add_argument("--vcf", required=True, help="Input AMBER SNP VCF file (e.g., SAMPLE.amber.snp.vcf.gz)")
    parser.add_argument("--snp_check_vcf", required=True, help="VCF defining the 100 SNPs to check (e.g., Amber.snpcheck.37.vcf)")
    parser.add_argument("--sample_id", required=True, help="The original sample ID from the VCF file name (e.g., SAMPLE001_T)")
    parser.add_argument("--output_file", required=True, help="Path to save the output TSV file")
    parser.add_argument("--db_host", required=True, help="Database host")
    parser.add_argument("--db_user", required=True, help="Database user")
    parser.add_argument("--db_pass", required=True, help="Database password")
    parser.add_argument("--db_name", required=True, help="Database name")
    parser.add_argument("--db_port", required=True, help="Database port")
    args = parser.parse_args()

    # Step 1: Get the ordered list of 100 SNPs to check
    ordered_loci = get_100_snp_loci(args.snp_check_vcf)

    # Step 2: Derive the 100 genotypes from the input VCF
    new_genotypes = derive_genotypes_from_vcf(args.vcf, ordered_loci)

    # Step 3: Find the best match in the database
    db_config = {
        'host': args.db_host,
        'user': args.db_user,
        'pass': args.db_pass,
        'name': args.db_name,
        'port': args.db_port
    }

    found_matches = find_matches_in_db(new_genotypes, db_config)

    # Step 4: Write the results to the specified TSV file
    with open(args.output_file, 'w', newline='') as f_out:
        writer = csv.writer(f_out, delimiter='\t')

        # Write header
        writer.writerow(['original_sampleId', 'matched_sampleId', 'match_count'])

        if not found_matches:
            # If no matches were found, write a single line indicating this
            writer.writerow([args.sample_id, 'NO_MATCH', 0])
        else:
            # Write all found matches to the file
            for matched_id, match_count in found_matches:
                writer.writerow([args.sample_id, matched_id, match_count])

    print(f"Results successfully written to {args.output_file}")


if __name__ == "__main__":
    main()