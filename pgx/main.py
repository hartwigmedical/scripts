import argparse
import allel
import pandas as pd
import sys
import requests
import json
import time
import collections
import itertools
import subprocess
import os
from shutil import copyfile


def compare_collection(a, b):
    if collections.Counter(a) == collections.Counter(b):
        return True
    else:
        return False


def load_panel(path):
    """ Load in a TSV file that contains the variants in the panel """
    try:
        panel = pd.read_csv(path, sep='\t')
        panel = panel[~panel['chrom'].isnull()]
        panel = panel[~panel['hg19_start'].isnull()]
        panel[['chrom']] = panel[['chrom']].astype(int)
        panel[['hg19_start']] = panel[['hg19_start']].astype(int)
    except IOError:
        sys.exit("[ERROR] File " + path + " not found or cannot be opened.")

    genes = {}
    ids = {}

    for index, row in panel.iterrows():
        if ',' in row['gene']:
            sep_genes = row['gene'].split(',')
        else:
            sep_genes = [row['gene']]
        for gene in sep_genes:
            if gene not in genes:
                genes[gene] = []
            genes[gene].append(row['rsid'])
        if row['rsid'] not in ids:
            ids[row['rsid']] = str(row['chrom']) + ":" + str(row['hg19_start'])

    return genes, ids


def load_json_cache(gene_panel, sourcedir):
    genes = {}
    ids = {}
    skipped_variants = []

    for gene in gene_panel:
        genes[gene] = []
        with open(sourcedir + '/pharmgkb_cache/' + gene + ".json") as json_file:
            data = json.load(json_file)
            for p in data['data']:
                if 'buildVersion' not in p['location']:
                    continue
                if p['location']['buildVersion'] not in ['GRCh37.p13', 'GRCh37']:
                    skipped_variants.append([p['location']['buildVersion'], p['location']['displayName'],
                                             p['location']['chromosomeName'].replace('chr', '') + ':' +
                                             str(p['location']['gpPosition'])])
                genomic_location = p['location']['chromosomeName'].replace('chr', '') + ':' + \
                                   str(p['location']['gpPosition'])
                genes[gene].append(p['location']['displayName'])
                if p['location']['displayName'] not in ids:
                    ids[p['location']['displayName']] = genomic_location

    for var in skipped_variants:
        if var[1] not in ids:
            print("[INFO] Skipped variant:")
            print("[INFO] " + "\t".join(var))

    return genes, ids


def load_json(panel):
    """ Load manually annotated JSON panel file """
    genes = {}  # Nested dict: gene, haplotype, rsinfo
    ids = {}

    try:
        with open(panel, 'r+', encoding='utf-8') as json_file:
            data = json.load(json_file)
            for gene_info in data['genes']:
                if gene_info['genomeBuild'] != "GRCh37":
                    sys.exit("Exiting, we only support GRCh37, not " + str(data['orientation']))
                if gene_info['gene'] not in genes:
                    genes[gene_info['gene']] = gene_info
                for variant in gene_info['variants']:
                    if variant['rsid'] not in ids:
                        ids[variant['rsid']] = str(variant['chromosome']) + ":" + str(variant['position'])

    except IOError:
        sys.exit("[ERROR] File " + panel + " not found or cannot be opened.")

    return genes, ids


def parse_vcf(vcf, rs_ids, bed_file, outputdir, sampleId, vcftools):
    match_on_rsid = 0
    match_on_location = 0

    # Slice VCF on bed file
    runname = vcf.split("/")[-2]
    temp_vcf_prefix = outputdir + '/' + sampleId + '_PGx'
    temp_vcf = outputdir + '/' + sampleId + '_PGx.recode.vcf'

    # Check if output vcf does not already exist
    if os.path.exists(temp_vcf):
        raise IOError("Temporary VCF file " + temp_vcf + ".recode.vcf already exists. Exiting.")
    subprocess.run([vcftools, '--gzvcf', vcf, '--bed', bed_file, '--out',
                    temp_vcf_prefix, '--recode', '--recode-INFO-all'])
    print("[INFO] Subprocess completed.")

    # Read in VCF file
    try:
        variants = allel.read_vcf(temp_vcf, fields=['samples', 'calldata/GT', 'variants/ALT', 'variants/CHROM',
                                                    'variants/FILTER', 'variants/ID', 'variants/POS',
                                                    'variants/QUAL', 'variants/REF', 'variants/ANN'],
                                  transformers=allel.ANNTransformer())
    except IOError:
        sys.exit("[ERROR] File " + temp_vcf + " not found or cannot be opened.")

    ids_found_in_patient = pd.DataFrame(columns=['position_GRCh37', 'ref_GRCh37', 'alt_GRCh37', 'rsid',
                                                 'variant_annotation', 'gene', 'filter'])

    for i, rs_number in enumerate(variants['variants/ID']):
        chr = variants['variants/CHROM'][i]
        pos = variants['variants/POS'][i]

        if ";" in rs_number:
            rs_id_filt = []
            cur_rs = rs_number.split(";")
            for rs in cur_rs:
                if rs.startswith("rs"):
                    rs_id_filt.append(rs)
        else:
            rs_id_filt = [rs_number]

        if any(rs in rs_id_filt for rs in rs_ids) or str(chr) + ":" + str(pos) in rs_ids.values():
            if any(rs in rs_id_filt for rs in rs_ids):
                match_on_rsid += 1
            else:
                match_on_location += 1
            new_id = {}
            if variants['variants/FILTER_PASS'][i] == True:
                filter = "PASS"
            else:
                filter = "FILTERED"
            alt = variants['variants/ALT'][i]
            ref = variants['variants/REF'][i]
            genotype = ""
            for geno in variants['calldata/GT'][i][0]:
                if geno == 0:
                    genotype = genotype + ref
                elif geno == 1:
                    genotype = genotype + alt[0]
                elif geno == 2:
                    genotype = genotype + alt[1]
                else:
                    print(geno)
                    raise ValueError("Genotype looks weird")

            new_id['position_GRCh37'] = str(chr) + ":" + str(pos)
            new_id['rsid'] = ";".join(rs_id_filt)
            new_id['ref_GRCh37'] = genotype[0]
            new_id['alt_GRCh37'] = genotype[1]
            new_id['variant_annotation'] = variants['variants/ANN_HGVS_c'][i]
            new_id['filter'] = filter
            new_id['gene'] = variants['variants/ANN_Gene_Name'][i]
            ids_found_in_patient = ids_found_in_patient.append(new_id, ignore_index=True)

    print("[INFO] Matches on RS id: " + str(match_on_rsid))
    print("[INFO] Matches on location: " + str(match_on_location))
    return ids_found_in_patient, temp_vcf


def process_exceptions(ids_found, panel):
    try:
        panel_exceptions_file = panel.replace(".json", "") + "_exceptions.json"
        with open(panel_exceptions_file, 'r+', encoding='utf-8') as json_file:
            exceptions = json.load(json_file)
            exceptions = exceptions['exceptions']
            for variant in exceptions:
                if variant['rsid'] in ids_found['rsid'].tolist():
                    # get line and index from ids_found
                    found_var = ids_found[ids_found['rsid'].str.contains(variant['rsid'])]
                    # Delete found variant from results if the ref base and alt base are the correctedRefBase
                    if found_var['ref_GRCh37'].values == variant['referenceAlleleGRCh38'] and \
                            found_var['alt_GRCh37'].values == variant['referenceAlleleGRCh38']:
                        ids_found = ids_found.drop(found_var.index[0])
                    elif found_var['ref_GRCh37'].values == variant['referenceAlleleGRCh38'] and \
                            found_var['alt_GRCh37'].values == variant['altAlleleGRCh38']:
                        print("[WARNING] What should we do? ref = corRef, alt = corAlt")
                    elif found_var['ref_GRCh37'].values == variant['altAlleleGRCh38'] and \
                            found_var['alt_GRCh37'].values == variant['referenceAlleleGRCh38']:
                        # Change variant_annotation and ref and alt base
                        ids_found.at[found_var.index[0], 'variant_annotation'] = variant['annotationGRCh38']
                        ids_found.at[found_var.index[0], 'ref_GRCh38'] = variant['referenceAlleleGRCh38']
                        ids_found.at[found_var.index[0], 'alt_GRCh38'] = variant['altAlleleGRCh38']
                        ids_found.at[found_var.index[0], 'position_GRCh38'] = str(variant['chromosome']) + ":" + \
                                                                              str(variant['positionGRCh38'])
                    elif found_var['ref_GRCh37'].values == variant['altAlleleGRCh38'] and \
                            found_var['alt_GRCh37'].values == variant['altAlleleGRCh38']:
                        print("[WARNING] What should we do? ref = corAlt, alt = corAlt")
                    else:
                        print("[ERROR] Complete mismatch:")
                        print(found_var)
                        raise ValueError("[ERROR] Exceptions cannot be processed. Please check. Exiting.")
                else:
                    print("[INFO] Exception variant not found in this patient. This means ref/ref call, but should be "
                          "flipped. Add to table.")
                    new_id = {'position_GRCh37': str(variant['chromosome']) + ":" + str(variant['position']),
                              'rsid': variant['rsid'],
                              'ref_GRCh37': variant['altAlleleGRCh38'],
                              'alt_GRCh37': variant['altAlleleGRCh38'],
                              'variant_annotation': variant['annotationGRCh38'],
                              'filter': "INFERRED_REF_CALL",
                              'gene': variant['gene'],
                              'position_GRCh38': str(variant['chromosome']) + ":" + str(variant['positionGRCh38']),
                              'ref_GRCh38': variant['altAlleleGRCh38'],
                              'alt_GRCh38': variant['altAlleleGRCh38']}
                    ids_found = ids_found.append(new_id, ignore_index=True)
    except IOError:
        sys.exit("[ERROR] Exceptions file not found or cannot be opened.")

    return ids_found


def convert_results_into_haplotypes(haplotypes_info, ids_found_in_patient, rs_ids, panel):
    results = {}
    severity = {}
    drug_info = {}

    # Process the exceptions
    ids_found_in_patient = process_exceptions(ids_found_in_patient, panel)
    # Fill in the GRCh38 gaps, they should be the same as the GRCh37 equivalent
    if not ids_found_in_patient.empty:
        for index, row in ids_found_in_patient.iterrows():
            if 'ref_GRCh38' in row:
                if pd.isna(row['ref_GRCh38']):
                    ids_found_in_patient.at[index, 'ref_GRCh38'] = row['ref_GRCh37']
                    ids_found_in_patient.at[index, 'alt_GRCh38'] = row['alt_GRCh37']
            else:
                ids_found_in_patient.at[index, 'ref_GRCh38'] = row['ref_GRCh37']
                ids_found_in_patient.at[index, 'alt_GRCh38'] = row['alt_GRCh37']
                ids_found_in_patient.at[index, 'position_GRCh38'] = ""

    # Generate a list of ids not found in patient
    ids_not_found_in_patient = pd.DataFrame(columns=['position_GRCh37', 'ref_GRCh37', 'alt_GRCh37', 'position_GRCh38',
                                                     'ref_GRCh38', 'alt_GRCh38', 'rsid', 'variant_annotation', 'gene',
                                                     'filter'])
    for rs in ids_found_in_patient.rsid.tolist():
        rs_ids.pop(rs)
    rs_ids = {key: val for key, val in rs_ids.items() if val not in ids_found_in_patient.position_GRCh37.tolist()}
    for item in rs_ids:
        new_id = {}
        for gene in haplotypes_info:
            for variant in haplotypes_info[gene]['variants']:
                if variant['rsid'] == item:
                    new_id['position_GRCh37'] = rs_ids[item]
                    new_id['rsid'] = item
                    new_id['ref_GRCh37'] = variant['referenceAllele']
                    new_id['alt_GRCh37'] = variant['referenceAllele']  # Assuming REF/REF
                    new_id['variant_annotation'] = "REF_CALL"
                    new_id['filter'] = "NO_CALL"
                    new_id['gene'] = gene
                    new_id['ref_GRCh38'] = variant['referenceAlleleGRCh38']  # Again assuming REF/REF
                    new_id['alt_GRCh38'] = variant['referenceAlleleGRCh38']
                    new_id['position_GRCh38'] = str(variant['chromosome']) + ":" + str(variant['positionGRCh38'])
                    ids_not_found_in_patient = ids_not_found_in_patient.append(new_id, ignore_index=True)

    # Now we want to process all the variants in terms of the alleles
    all_ids_in_panel = pd.concat([ids_found_in_patient, ids_not_found_in_patient], sort=True)
    all_ids_in_panel = all_ids_in_panel.sort_values(by='position_GRCh37').reset_index(drop=True)
    all_ids_in_panel = all_ids_in_panel[['gene', 'position_GRCh37', 'ref_GRCh37', 'alt_GRCh37', 'position_GRCh38',
                                         'ref_GRCh38', 'alt_GRCh38', 'rsid', 'variant_annotation', 'filter']]

    for gene in haplotypes_info:
        print("[INFO] PROCESSING GENE " + gene)
        gene_info = haplotypes_info[gene]
        ids_found_in_gene = all_ids_in_panel[all_ids_in_panel['gene'].str.contains(gene)]
        perfect_match = False
        severity[gene_info['referenceAllele']] = "Normal Function"
        drug_info[gene] = [";".join([x['name'] for x in gene_info['drugs']]), ";".join(x['url_prescription_info'] for x in gene_info['drugs'])]

        # If all variants are assumed_ref, return reference allele
        if len(all_ids_in_panel.loc[all_ids_in_panel['filter'] == "NO_CALL"]) == len(all_ids_in_panel):
            print("[INFO] Found reference allele")
            results[gene] = [gene_info['referenceAllele'] + "_HOM"]
        else:
            results[gene] = []
            haplotypes_matching = []
            for allele in gene_info['alleles']:
                severity[allele['alleleName']] = allele['function']
                variants_sample = list(zip(all_ids_in_panel.rsid.tolist(), all_ids_in_panel.alt_GRCh38.tolist()))
                if set(variants_sample) == set([(x['rsid'], x['altAlleleGRCh38']) for x in allele['alleleVariants']]):
                    perfect_match = True
                    print("[INFO] Found 1:1 match with allele " + allele['alleleName'])
                    # Now we want to see if we have hetrozygous or homozygous calls
                    allele_status = []
                    for index, row in all_ids_in_panel.iterrows():
                        if row['ref_GRCh38'] == row['alt_GRCh38']:
                            allele_status.append("HOM")
                        else:
                            allele_status.append("HET")
                        # print(allele_status)
                        if all(x == allele_status[0] for x in allele_status):
                            allele_status = allele_status[0]
                        else:
                            allele_status = "HOMHET"
                        # Add to results
                        results[gene].append(allele['alleleName'] + "_" + str(allele_status))
                        if allele_status == "HET":
                            # Assume if perfect match with HET, we are also looking at reference allele
                            results[gene].append(gene_info['referenceAllele'] + "_HET")
                        break
                else:
                    if set([(x['rsid'], x['altAlleleGRCh38']) for x in allele['alleleVariants']]) <= \
                            set(variants_sample):
                        print("[INFO] A subset of rsids matches " + str(allele['alleleName']) + " in part")
                        haplotypes_matching.append(allele)

            if not perfect_match:
                if not haplotypes_matching:
                    sys.exit("[ERROR] No allele match found for " + str(gene))
                else:
                    print("[INFO] Test all possible combinations of haplotypes to see if a perfect match can be found")
                    optimal_set = []
                    for k in range(len(haplotypes_matching) + 1, 0, -1):
                        for subset in itertools.combinations(haplotypes_matching, k):
                            # See if this combination results in a perfect match, otherwise store the score
                            allele_variants = [x['alleleVariants'] for x in subset]
                            rs_ids_subset = []
                            for x in allele_variants:
                                for var in x:
                                    rs_ids_subset.append((var['rsid'], var['altAlleleGRCh38']))
                            if compare_collection(variants_sample, rs_ids_subset):
                                print("[INFO] Perfect haplotype combination found!")
                                perfect_match = True
                                for allele in subset:
                                    allele_status = []
                                    rs_ids_in_allele = [x['rsid'] for x in allele['alleleVariants']]
                                    found_vars = ids_found_in_gene[ids_found_in_gene['rsid'].isin(rs_ids_in_allele)]
                                    for index, row in found_vars.iterrows():
                                        if row['ref_GRCh38'] == row['alt_GRCh38']:
                                            allele_status.append("HOM")
                                        else:
                                            allele_status.append("HET")
                                    if all(x == allele_status[0] for x in allele_status):
                                        allele_status = allele_status[0]
                                    else:
                                        allele_status = "HOMHET"
                                    results[gene].append(allele['alleleName'] + "_" + str(allele_status))
                            else:
                                rs_matched = list(set(variants_sample) & set(rs_ids_subset))
                                rs_not_in_haplotype = list(set(variants_sample) - set(rs_ids_subset))
                                rs_not_found = list(set(rs_ids_subset) - set(variants_sample))
                                optimal_set.append([len(rs_matched), subset, rs_matched, rs_not_found,
                                                    rs_not_in_haplotype])
                    if not perfect_match:
                        # Get best scoring set and give as result
                        print("[INFO] No perfect match found. Start testing all possible combinations and choose best "
                              "one")
                        optimal_set.sort(key=lambda x: x[0], reverse=True)
                        for options in optimal_set:
                            print("[INFO] Option to be tested: " + str(options))
                        if optimal_set[0][0] >= 1:
                            subset = optimal_set[0][1]
                            for allele in subset:
                                allele_status = []
                                rs_ids_in_allele = [x['rsid'] for x in allele['alleleVariants']]
                                found_vars = ids_found_in_gene[ids_found_in_gene['rsid'].isin(rs_ids_in_allele)]
                                for index, row in found_vars.iterrows():
                                    if row['ref_GRCh38'] == row['alt_GRCh38']:
                                        allele_status.append("HOM")
                                    else:
                                        allele_status.append("HET")
                                if all(x == allele_status[0] for x in allele_status):
                                    allele_status = allele_status[0]
                                else:
                                    allele_status = "HOMHET"
                                results[gene].append(allele['alleleName'] + "_" + str(allele_status))
                        else:
                            sys.exit("[ERROR] No haplotype match was found. Exiting.")
        # If we only find one allele and it is HET, assume we're also dealing with reference allele
        if len(results[gene]) == 1:
            if results[gene][0].split("_")[-1] == "HET":
                results[gene].append(gene_info['referenceAllele'] + "_HET")

    return ids_found_in_patient, results, severity, all_ids_in_panel, drug_info


def get_bed_file(gene_panel, panel_path, sourcedir):
    """ Generate bed file from gene panel and save as panel_path.bed """
    print("[INFO] Recreating bed-file...")
    header = 'track name="' + panel_path + '" description="Bed file generated from ' + panel_path + \
             ' with HMF_PGx main.py"\n'
    bed_regions = []    # chrom, start, end, gene
    covered = []
    transcripts = open(sourcedir + "/all_genes.37.tsv", 'r')
    for line in transcripts:
        line = line.rstrip().split("\t")
        if line[4] in gene_panel:
            if line[4] not in covered:
                bed_regions.append([line[0], line[1], line[2], line[4]])
                covered.append(line[4])
    if set(covered) != set(gene_panel):
        raise ValueError("[ERROR] Missing genes from the gene panel in the transcript list. Please check:\nCovered:\n"
                         + str(covered) + "\nOriginal gene panel:\n" + str(gene_panel))

    if panel_path:
        bed_path = panel_path.split(".")
        bed_path = ".".join(bed_path[0:-1]) + ".bed"
    else:
        bed_path = sourcedir + "/pharmgkb_cache/cache_bed.bed"
    print("[INFO] Created " + bed_path)

    bed = open(bed_path, 'w')
    bed.write(header)
    for entry in bed_regions:
        bed.write("\t".join(entry) + "\n")
    bed.close()

    return bed_path


def main(vcf, sampleID, panel, requery, outputdir, recreate_bed, vcftools, sourcedir):
    """ Run pharmacogenomics analysis on sample """
    print("\n[INFO] ## START PHARMACOGENOMICS ANALYSIS")

    haplotypes_info = None
    genes = None
    
    # Check if output dir exists, create if it does not
    if not os.path.exists(outputdir):
        try:
            os.makedirs(outputdir)
        except FileExistsError:
            # Directory already exists
            pass

    # If panel file given as input, load panel
    if panel:
        if panel.endswith('.json'):
            # Load JSON panel file
            haplotypes_info, rs_ids = load_json(panel)
        else:
            # Give in a TSV file with the variants of the panel
            genes, rs_ids = load_panel(panel)
    else:
        # No panel is found, so load the JSON files from the pharmGKB
        gene_panel = ['CFTR', 'CYP2C19', 'CYP2C9', 'CYP2D6', 'CYP3A5', 'CYP4F2', 'DPYD', 'IFNL3', 'SLCO1B1', 'TPMT',
                      'UGT1A1', 'VKORC1']

        if requery:
            print("[INFO] Loading data from pharmGKB")
            for gene in gene_panel:
                time.sleep(1)
                url = "https://api.pharmgkb.org/v1/data/clinicalAnnotation?location.genes.symbol=" + gene
                print(url)
                r = requests.get(url)
                if r.status_code == 200:
                    with open(sourcedir + '/pharmgkb_cache/' + gene + ".json", 'w') as outfile:
                        json.dump(r.json(), outfile)
                elif r.status_code == 404:
                    print(gene + " not found in pharmGKB API")
                elif r.status_code == 429:
                    print("Too many requests per second, please delay querying.")
                else:
                    print("Other error: " + r.status_code)
        print("[INFO] Reading in stored pharmGKB data")
        genes, rs_ids = load_json_cache(gene_panel, sourcedir)

    if rs_ids != {}:
        if recreate_bed:
            if genes:
                gene_panel = list(genes)
            elif haplotypes_info:
                gene_panel = list(haplotypes_info)
            else:
                raise ValueError("[ERROR] No genes found for bed-file filtering. Exiting.")
            bed_file = get_bed_file(gene_panel, panel, sourcedir)
        else:
            if panel:
                bed_file = panel.split(".")
                bed_file = ".".join(bed_file[0:-1]) + ".bed"
            else:
                bed_file = sourcedir + "/pharmgkb_cache/cache_bed.bed"
            if not os.path.exists(bed_file):
                sys.exit("[ERROR] Could not locate bed-file. Could it be that it should be (re)created? Retry running with --recreate_bed.")
        ids_found_in_patient, temp_vcf = parse_vcf(vcf, rs_ids, bed_file, outputdir, sampleID, vcftools)
    else:
        sys.exit("[ERROR] No panel variants are given, no analysis is performed.")

    if haplotypes_info:
        ids_found_in_patient, results, severity, all_ids_in_panel, drug_info = \
            convert_results_into_haplotypes(haplotypes_info, ids_found_in_patient, rs_ids, panel)
        if outputdir:
            out = outputdir + "/" + sampleID
            all_ids_in_panel.to_csv(out + "_calls.txt", sep='\t', index=False)
            f = open(out + "_genotype.txt", 'w')
            f.write("gene\thaplotype\tfunction\tlinked_drugs\turl_prescription_info\tpanel_version\trepo_version\n")
            # Below solution for > Python 3.7
            # git_describe = subprocess.run(["git", "--git-dir=" +
            # os.path.dirname(os.path.realpath(__file__)) + "/.git", "describe", "--always"], capture_output=True)
            git_describe = subprocess.run(["git", "--git-dir=" + os.path.dirname(os.path.realpath(__file__)) + "/.git",
                                           "describe", "--always"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print("git version")
            print(git_describe)
            git_describe = git_describe.stdout.decode("utf-8").strip()
            print(git_describe)
            for gene in results:
                for haplotype in results[gene]:
                    f.write(gene + "\t" + haplotype + "\t" + severity[haplotype.split("_")[0]] + "\t" + drug_info[gene][0] + "\t" + drug_info[gene][1] + "\t" + panel + "\t" +
                            git_describe + "\n")
            f.close()
            # Also copy the bed-filtered VCF file for research purposes
            copyfile(temp_vcf, out + "_PGx.vcf")

        else:
            # Should not be executed, as outputdir is now required argument
            for gene in results:
                for haplotype in results[gene]:
                    print(gene + "\t" + haplotype + "\t" + severity[haplotype.split("_")[0]] + "\n")

    else:
        print(ids_found_in_patient)

    # Clean up temp_vcf
    if os.path.exists(temp_vcf):
        if os.path.exists(temp_vcf):
            os.remove(temp_vcf)
            print("[INFO] " + temp_vcf + " removed.")
        if os.path.exists(temp_vcf.replace(".recode.vcf", ".log")):
            os.remove(temp_vcf.replace(".recode.vcf", ".log"))
            print("[INFO] " + temp_vcf.replace(".recode.vcf", ".log") + " removed.")

    # TODO: add genes CYP2D6, CYP3A4, CYP3A5

    print("[INFO] ## PHARMACOGENOMICS ANALYSIS FINISHED\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run pharmacogenomics panel on germline VCF file. The pharmacogenomic '
                                                 'annotations are done on GRCh38, so in the output both reference '
                                                 'genome output is given where possible.')
    parser.add_argument('vcf', type=str, help='VCF file to use for pharmacogenomics analysis')
    parser.add_argument('sampleID', type=str, help='The sample ID of the run')
    parser.add_argument('outputdir', type=str, help='Directory to store output of pharmacogenomic analysis')
    parser.add_argument('--panel', type=str, help='TSV file with the panel variants')
    parser.add_argument('--requery', default=False, action='store_true', help='Requery genes in pharmGKB')
    parser.add_argument('--recreate_bed', default=False, action='store_true', help='Recreate bed-file from JSON files. '
                                                                                   'If false, the panel file with '
                                                                                   'extension .bed is searched for.')
    parser.add_argument('--vcftools', type=str, default='vcftools', help="Path to vcftools > 0.1.14 if not in $PATH")
    parser.add_argument('--sourcedir', type=str, default='data', help="Optional path to location of source files")

    args = parser.parse_args()

    main(args.vcf, args.sampleID, args.panel, args.requery, args.outputdir, args.recreate_bed, args.vcftools, args.sourcedir)
