#!/usr/bin/env python

"""
melt_somatic_vcf.py - Version 2
Melts / merges somatic vcf files coming from the IAP.
Supported somatic callers: Freebayes, Mutect, Strelka and Varscan.

Input:
    - vcf_file: merged somatic vcf file created by the IAP.
    - tumor_sample: tumor sample name, used to find somatic tumor columns and as the sample name in the output vcf.

Output:
    - Melted vcf to stdout:
        - Prints original header + extra header lines describing new CSA and CSP info fields.
        - New info fields:
            - CSA: Number of callers suporting each allele -> [ref, alt:]
            - CSP: Number of callers suporting the position
        - Melted GT format:
            - GT: Contains all suported allele calls, because of this the genotype can be 'non-diploid', for example 0/1/2
            - AD: Contains ad values for all suported alleles.
                - AD values are averages from callers with support for the allele.
                - Ref ad is set to 0 if there is no support from the callers.
            - DP:
                - Average DP from callers with support for position.


Known limitations:
- GATK combineVariants does not adjust freebayes AO field when changing the order of alt alleles. This will result in incorrect AD fields for some variants with multiple alt alleles.
This can not be fixed in this melt script.

"""

import sys
import argparse
from itertools import izip_longest
import re

def melt_somatic_vcf(vcf_file, remove_filtered, tumor_sample):
    try:
        f = open(vcf_file, 'r')
    except IOError:
        sys.exit("Error: Can't open vcf file: {0}".format(vcf_file))
    else:
        with f:
            vcf_header = []
            for line in f:
                line = line.strip('\n')
                if line.startswith('##'):
                    '''Print original vcf meta-information lines '''
                    print line

                elif line.startswith("#CHROM"):
                    '''Parse samples and print new header'''
                    header = line.split('\t')
                    samples = header[9:]

                    # Find tumor samples index
                    tumor_samples_index = {}
                    for index, sample in enumerate(samples):
                        if '{0}.freebayes'.format(tumor_sample) == sample:
                            tumor_samples_index['freebayes'] = index
                        elif '{0}.mutect'.format(tumor_sample) == sample:
                            tumor_samples_index['mutect'] = index
                        elif 'TUMOR.strelka' == sample:
                            tumor_samples_index['strelka'] = index
                        elif 'TUMOR.varscan' == sample:
                            tumor_samples_index['varscan'] = index

                    # Check somatic samples
                    if len(samples)/2 != len(tumor_samples_index.keys()):
                        sys.exit("Error: Found {0} sample somatic variant caller combinations, expected {1}. Please check tumor_sample name.".format(len(tumor_samples_index.keys()),len(samples)/2))

                    ## print meta-information lines with melter info to vcf
                    print "##source=IAP/scripts/melt_somatic_vcf.py"
                    print "##INFO=<ID=CSA,Number=R,Type=Integer,Description=\"Number of somatic variant callers with support for allele.\">"
                    print "##INFO=<ID=CSP,Number=R,Type=Integer,Description=\"Number of somatic variant callers with support for position.\">"
                    print "{header}\t{sample}".format(
                        header = '\t'.join(header[:9]),
                        sample = tumor_sample
                    )

                else:
                    '''Parse variant line and print'''
                    variant = line.split('\t')

                    #Skip variants with a filter flag other than PASS.
                    if remove_filtered:
                        if variant[6].upper() != 'PASS' and variant[6] != '.':
                            continue

                    #Setup variant variables
                    ref = variant[3]
                    alts = variant[4].split(',')
                    variant_calls = variant[9:]
                    variant_dp_counts = []
                    caller_count = 0 # used for PCS info field

                    if len(ref) == 1 and len(max(alts, key=len)) == 1:
                        type = 'snp'
                    else:
                        type = 'indel'

                    # Setup gt format indexes
                    variant_gt_format = variant[8].split(':')
                    dp_index = variant_gt_format.index('DP')
                    gt_index = variant_gt_format.index('GT')

                    ### Check ref support among callers
                    ref_support = 0 #used for ACS info field
                    ref_ad = []

                    for somatic_caller, tumor_sample_index in tumor_samples_index.iteritems():
                        # Skip no calls
                        if variant_calls[tumor_sample_index] == './.':
                            continue

                        variant_call = variant_calls[tumor_sample_index].split(':')

                        #Add one to caller count and save variant DP
                        #Variant DP, field is the same per caller
                        caller_count += 1
                        variant_dp_counts.append(float(variant_call[dp_index]))

                        # Check support for ref allele call
                        if '0' in variant_call[gt_index]:
                            ref_support += 1

                            ## Collect AD
                            if somatic_caller == 'freebayes':
                                # freebayes_example.vcf:##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">
                                ro_index = variant_gt_format.index('RO')
                                ref_ad.append(float(variant_call[ro_index]))
                            elif somatic_caller == 'mutect':
                                # mutect_example.vcf:##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
                                if 'AD' in variant_gt_format:
                                    ad_index = variant_gt_format.index('AD')
                                    ref_ad.append(float(variant_call[ad_index].split(',')[0]))
                                else: # AD NOT ALWAYS present.   Should use RD = (1-FA)*DP, AD = FA * DP  in this case
                                    freq = float(variant_call[variant_gt_format.index('FA')])
                                    variant_dp = float(variant_call[dp_index])
                                    ref_ad.append((1-freq) * variant_dp)
                            elif somatic_caller == 'varscan':
                                # varscan_example.vcf:##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
                                rd_index = variant_gt_format.index('RD')
                                ref_ad.append(float(variant_call[rd_index]))
                            elif somatic_caller == 'strelka':
                                if type == 'snp':
                                    # strelka_example.vcf:##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of 'A' alleles used in tiers 1,2">
                                    # strelka_example.vcf:##FORMAT=<ID=CU,Number=2,Type=Integer,Description="Number of 'C' alleles used in tiers 1,2">
                                    # strelka_example.vcf:##FORMAT=<ID=GU,Number=2,Type=Integer,Description="Number of 'G' alleles used in tiers 1,2">
                                    # strelka_example.vcf:##FORMAT=<ID=TU,Number=2,Type=Integer,Description="Number of 'T' alleles used in tiers 1,2">
                                    ref_ad.append(float(variant_call[variant_gt_format.index(ref+'U')].split(',')[0]))
                                elif type == 'indel':
                                    #strelka_example.vcf:##FORMAT=<ID=TAR,Number=2,Type=Integer,Description="Reads strongly supporting alternate allele for tiers 1,2"> == REF
                                    ref_ad.append(float(variant_call[variant_gt_format.index('TAR')].split(',')[0]))

                    ## Check support for each alternative allele among callers
                    alts_support = [] #used for ACS info field
                    alts_ad = []
                    for alt_index, alt in enumerate(alts):
                        alt_support = 0
                        alt_ad = []
                        alt_allele_num = alt_index + 1

                        for somatic_caller, tumor_sample_index in tumor_samples_index.iteritems():
                            # Skip no calls
                            if variant_calls[tumor_sample_index] == './.':
                                continue

                            variant_call = variant_calls[tumor_sample_index].split(':')

                            # Check support for alt allele call
                            variant_call_gt = variant_call[gt_index].split('/')
                            if str(alt_allele_num) in variant_call_gt:
                                alt_support += 1

                                ## Collect AD
                                if somatic_caller == 'freebayes':
                                    # freebayes_example.vcf:##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
                                    # GATK combineVariants does not change AO field if the order of alt alleles is changed. Therefor in some cases an incorrect AO value is selected.
                                    variant_ao_index = variant_gt_format.index('AO')
                                    variant_ao =  variant_call[variant_ao_index].split(',')
                                    gt_alt_index = min(alt_index,(len(variant_ao)-1))
                                    alt_ad.append(float(variant_ao[gt_alt_index]))
                                elif somatic_caller == 'mutect':
                                    # mutect_example.vcf:##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
                                    if 'AD' in variant_gt_format:
                                        ad_index = variant_gt_format.index('AD')
                                        alt_ad.append(float(variant_call[ad_index].split(',')[alt_allele_num]))
                                    else: # AD NOT ALWAYS present.   Should use RD = (1-FA)*DP, AD = FA * DP  in this case
                                        freq = float(variant_call[variant_gt_format.index('FA')])
                                        variant_dp = float(variant_call[dp_index])
                                        alt_ad.append(freq * variant_dp)
                                elif somatic_caller == 'varscan':
                                    # varscan_example.vcf:##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
                                    if 'AD' in variant_gt_format:
                                        ad_index = variant_gt_format.index('AD')
                                        alt_ad.append(float(variant_call[ad_index]))
                                    else:
                                        # AD NOT ALWAYS present for varscan.  use FREQ * DP in this case
                                        freq = float(variant_call[variant_gt_format.index('FREQ')].rstrip('%'))/100
                                        variant_dp = float(variant_call[dp_index])
                                        alt_ad.append(freq * variant_dp)
                                elif somatic_caller == 'strelka':
                                    if type == 'snp':
                                        # strelka_example.vcf:##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of 'A' alleles used in tiers 1,2">
                                        # strelka_example.vcf:##FORMAT=<ID=CU,Number=2,Type=Integer,Description="Number of 'C' alleles used in tiers 1,2">
                                        # strelka_example.vcf:##FORMAT=<ID=GU,Number=2,Type=Integer,Description="Number of 'G' alleles used in tiers 1,2">
                                        # strelka_example.vcf:##FORMAT=<ID=TU,Number=2,Type=Integer,Description="Number of 'T' alleles used in tiers 1,2">
                                        alt_ad.append(float(variant_call[variant_gt_format.index(alt+'U')].split(',')[0]))
                                    else:
                                        #strelka_example.vcf:##FORMAT=<ID=TIR,Number=2,Type=Integer,Description="Reads strongly supporting indel allele for tiers 1,2"> == ALT
                                        alt_ad.append(float(variant_call[variant_gt_format.index('TIR')].split(',')[0]))

                        # Store results
                        alts_support.append(alt_support)
                        alts_ad.append(alt_ad)

                    # Create gt and ad lists
                    gt = []
                    ad = []
                    if ref_support:
                        gt.append('0')
                        ad.append(int(round(sum(ref_ad)/ref_support)))
                    else: # If no support for ref allele, set ref ad to 0.
                        ad.append(0)

                    for alt_index, alt_support in enumerate(alts_support):
                        if alt_support:
                            gt.append(str(alt_index + 1))
                            ad.append(int(round(sum(alts_ad[alt_index])/alt_support)))

                    # Correct hom calls
                    if len(gt) == 1:
                        gt.append(gt[0])

                    #Calculate mean DP among callers with support for position
                    variant_dp = int(round(sum(variant_dp_counts)/len(variant_dp_counts)))

                    # Print variant
                    print "{var_data};CSA={csa};CSP={csp}\t{gt_format}\t{gt}:{ad}:{dp}".format(
                        var_data = "\t".join(variant[:8]),
                        csa = ','.join(map(str,[ref_support]+alts_support)),
                        csp = caller_count,
                        gt_format = "GT:AD:DP",
                        gt = '/'.join(gt),
                        ad = ','.join(map(str,ad)),
                        dp = str(variant_dp),
                    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=100, width=200))

    parser.add_argument('--remove_filtered', action='store_true', help='Skip variants with a filter flag other than PASS.')

    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument('-v', '--vcf_file', help='path/to/file.vcf', required=True)
    required_named.add_argument('-t', '--tumor_sample', help='Tumor sample name', required=True)

    args = parser.parse_args()
    melt_somatic_vcf(args.vcf_file, args.remove_filtered, args.tumor_sample)
