#!/usr/bin/env python

"""
melt_somatic_vcfNEW.py
Melts / merges somatic vcf files coming from the IAP.
Supported somatic callers: Freebayes, Mutect, Strelka and Varscan.

KEY CHANGES:
* Select the allele with the most allelic depth support rather than always choosing the first
* Only average AD, RD & DP for callers which AGREE on variant allele
* Handle unhandled corner cases / bugs in calculation of AD/RD
    * Mutect - when AD does not exist use FA*DP
    * Varscan - when AD does not exist use FREQ*DP
    * Freebayes - chose the AO for the selected allele instead of the 1st
    * Strelka - choose the AD for the selected allele only instead of all alternates
* Add 'diffAllele' field to "INFO" to indicate if any callers disagree with the allele

"""

import argparse
from itertools import izip_longest
import re

def melt_somatic_vcf(vcf_file, remove_filtered, tumor_sample):
    try:
        f = open(vcf_file, 'r')
    except IOError:
        print "Can't open vcf file: {0}".format(vcf_file)
    else:
        with f:
            for line in f:
                line = line.strip('\n')

                if line.startswith('##'):
                    '''Print original vcf meta-information lines '''
                    print line
                    continue

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

                    # sample name == file_name
                    sample_name = vcf_file.split('.')[0]

                    ## Add meta-information lines with melter info to vcf
                    print "##source=IAP/scripts/melt_somatic_vcf.py"
                    print "##INFO=<ID=CC,Number=1,Type=Integer,Description=\"Number of somatic variant callers with call.\">"
                    print "##INFO=<ID=diffAllele,Number=1,Type=String,Description=\"Callers with different Alleles\">"

                    ## print header
                    print "{header}\t{sample}".format(
                        header = '\t'.join(header[:9]),
                        sample = sample_name
                    )

                else:

                    '''Parse variant line and print'''
                    variant = line.split('\t')
                    variant_gt_format = variant[8].split(':')
                    variant_calls = variant[9:]
                    caller_count = 0

                    #Skip variants with a filter flag other than PASS.
                    if remove_filtered:
                        if variant[6].upper() != 'PASS' and variant[6] != '.':
                            continue

                    calledVariants = []  #list of [somaticCaller,altNum,DP,RD,AD]
                    dp_index = variant_gt_format.index('DP')
                    gt_index = variant_gt_format.index('GT')
                    ADPerAllele = {}

                    for somatic_caller, tumor_sample_index in tumor_samples_index.iteritems():
                        # Skip no calls
                        if variant_calls[tumor_sample_index] == './.':
                            continue

                        variant_call = variant_calls[tumor_sample_index].split(':')

                        caller_count += 1

                        altNum = int(variant_call[gt_index].split('/')[1])   # Always choose allele B for genotype A/B
                        variant_dp = float(variant_call[dp_index])

                        #Populate Variant info for each caller
                        if somatic_caller == 'freebayes':
                            ro_index = variant_gt_format.index('RO')
                            ao_index = variant_gt_format.index('AO')
                            # freebayes_example.vcf:##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
                            # freebayes_example.vcf:##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">
                            variant_rd = float(variant_call[ro_index])
                            ao_split = variant_call[ao_index].split(',')
                            variant_ad = float(ao_split[min(altNum,len(ao_split))-1])

                        elif somatic_caller == 'mutect':
                            if 'AD' in variant_gt_format:
                                ad_index = variant_gt_format.index('AD')
                                # mutect_example.vcf:##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
                                variant_rd = float(variant_call[ad_index].split(',')[0])
                                variant_ad = float(variant_call[ad_index].split(',')[altNum])
                            else: # PAD NOT ALWAYS present.   Should use RD = (1-FA)*DP, AD = FA * DP  in this case
                                freq = float(variant_call[variant_gt_format.index('FA')])
                                variant_rd = (1-freq) * variant_dp
                                variant_ad = freq * variant_dp

                        elif somatic_caller == 'varscan':
                            rd_index = variant_gt_format.index('RD')
                            variant_rd = float(variant_call[rd_index])
                            if 'AD' in variant_gt_format:
                                ad_index = variant_gt_format.index('AD')
                                # varscan_example.vcf:##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
                                # varscan_example.vcf:##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
                                variant_ad = float(variant_call[ad_index].split(',')[altNum - 1])
                            else:  # AD NOT ALWAYS present for varscan.  use FREQ * DP in this case
                                freq = float(variant_call[variant_gt_format.index('FREQ')].split('%')[0])/100
                                variant_ad = freq*variant_dp

                        elif somatic_caller == 'strelka':
                            ref = variant[3]
                            altAllele = variant[4].split(',')[altNum-1]
                            # Parse snps
                            if len(ref) == 1 and len(altAllele) == 1:
                                variant_rd = float(variant_call[variant_gt_format.index(ref+'U')].split(',')[0])
                                variant_ad = float(variant_call[variant_gt_format.index(altAllele + 'U')].split(',')[0])
                                # strelka_example.vcf:##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of 'A' alleles used in tiers 1,2">
                                # strelka_example.vcf:##FORMAT=<ID=CU,Number=2,Type=Integer,Description="Number of 'C' alleles used in tiers 1,2">
                                # strelka_example.vcf:##FORMAT=<ID=GU,Number=2,Type=Integer,Description="Number of 'G' alleles used in tiers 1,2">
                                # strelka_example.vcf:##FORMAT=<ID=TU,Number=2,Type=Integer,Description="Number of 'T' alleles used in tiers 1,2">
                            # Parse Indels
                            else:
                                variant_rd = float(variant_call[variant_gt_format.index('TAR')].split(',')[0])
                                variant_ad = float(variant_call[variant_gt_format.index('TIR')].split(',')[0])
                                #strelka_example.vcf:##FORMAT=<ID=TAR,Number=2,Type=Integer,Description="Reads strongly supporting alternate allele for tiers 1,2"> REFERENCE!!!
                                #strelka_example.vcf:##FORMAT=<ID=TIR,Number=2,Type=Integer,Description="Reads strongly supporting indel allele for tiers 1,2">

                        # Add to array and allelic depth tracker
                        calledVariants.append([somatic_caller, altNum, variant_dp, variant_rd, variant_ad])
                        if ADPerAllele.has_key(altNum):
                            ADPerAllele[altNum] += variant_ad
                        else:
                            ADPerAllele[altNum] = variant_ad

                    # Filter for allele with the most AD and calculate averages
                    selectedAllele = max(ADPerAllele, key=ADPerAllele.get)
                    diffAlleleVariants = [x for x in calledVariants if x[1] != selectedAllele]
                    calledVariants = [x for x in calledVariants if x[1] == selectedAllele]
                    DP = int(round(sum(calledVariant[2] for calledVariant in calledVariants) / len(calledVariants)))
                    RD = int(round(sum(calledVariant[3] for calledVariant in calledVariants) / len(calledVariants)))
                    AD = int(round(sum(calledVariant[4] for calledVariant in calledVariants) / len(calledVariants)))

                    # Create "INFO" string if there is a caller with a different allele
                    if diffAlleleVariants:
                        diffAlleleString = ";diffAllele="
                        for diffAlleleVariant in diffAlleleVariants:
                            if diffAlleleString[-1] != '=':
                                diffAlleleString += "-"
                            diffAlleleString += diffAlleleVariant[0]
                    else:
                        diffAlleleString = ""

                    print "{var_data};CC={cc}{diffAllele}\t{gt_format}\t{gt}:{ad}:{dp}".format(
                        var_data = "\t".join(variant[:8]),
                        cc = caller_count,
                        diffAllele = diffAlleleString,
                        gt_format = "GT:AD:DP",
                        gt = "0/" + str(selectedAllele),
                        ad = str(RD) + "," + str(AD),
                        dp = str(DP),
                    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=100, width=200))

    parser.add_argument('--remove_filtered', action='store_true', help='Skip variants with a filter flag other than PASS.')

    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument('-v', '--vcf_file', help='path/to/file.vcf', required=True)
    required_named.add_argument('-t', '--tumor_sample', help='Tumor sample name', required=True)

    args = parser.parse_args()
    melt_somatic_vcf(args.vcf_file, args.remove_filtered, args.tumor_sample)
