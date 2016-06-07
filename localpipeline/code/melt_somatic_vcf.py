#!/usr/bin/env python

"""
melt_somatic_vcf.py
Melts / merges somatic vcf files coming from the IAP.
Supported somatic callers: Freebayes, Mutect, Strelka and Varscan.
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

                    dp_index = variant_gt_format.index('DP')
                    variant_dp = []
                    freebayes_ad, mutect_ad, varscan_ad, strelka_ad = ([],[],[],[])
                    variant_ad = []

                    for somatic_caller, tumor_sample_index in tumor_samples_index.iteritems():
                        # Skip no calls
                        if variant_calls[tumor_sample_index] == './.':
                            continue

                        variant_call = variant_calls[tumor_sample_index].split(':')
                        caller_count += 1
                        #Variant DP, field is the same per caller
                        variant_dp.append(float(variant_call[dp_index]))

                        #Variant AD, field differs per caller.
                        if somatic_caller == 'freebayes':
                            try:
                                ro_index = variant_gt_format.index('RO')
                                ao_index = variant_gt_format.index('AO')
                            except ValueError:
                                continue
                            freebayes_ad = map(float,[variant_call[ro_index]] + variant_call[ao_index].split(','))
                            #freebayes_example.vcf:##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
                            #freebayes_example.vcf:##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">

                        elif somatic_caller == 'mutect':
                            try:
                                ad_index = variant_gt_format.index('AD')
                            except ValueError:
                                continue
                            mutect_ad = map(float, variant_call[ad_index].split(','))
                            #mutect_example.vcf:##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">

                        elif somatic_caller == 'varscan':
                            try:
                                rd_index = variant_gt_format.index('RD')
                                ad_index = variant_gt_format.index('AD')
                            except ValueError:
                                continue
                            varscan_ad = map(float,[variant_call[rd_index]] + variant_call[ad_index].split(','))
                            #varscan_example.vcf:##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
                            #varscan_example.vcf:##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">

                        elif somatic_caller == 'strelka':
                            gt = [variant[3]] + variant[4].split(',') #[ref, alt, alt]

                            # Parse snps
                            if len(gt[0]) == 1 and len(gt[1]) == 1:
                                for nucl in gt:
                                    if nucl == "A":
                                        strelka_ad.append(float(variant_call[variant_gt_format.index('AU')].split(',')[0]))
                                    elif nucl == "C":
                                        strelka_ad.append(float(variant_call[variant_gt_format.index('CU')].split(',')[0]))
                                    elif nucl == "G":
                                        strelka_ad.append(float(variant_call[variant_gt_format.index('GU')].split(',')[0]))
                                    elif nucl == "T":
                                        strelka_ad.append(float(variant_call[variant_gt_format.index('TU')].split(',')[0]))
                                        #strelka_example.vcf:##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of 'A' alleles used in tiers 1,2">
                                        #strelka_example.vcf:##FORMAT=<ID=CU,Number=2,Type=Integer,Description="Number of 'C' alleles used in tiers 1,2">
                                        #strelka_example.vcf:##FORMAT=<ID=GU,Number=2,Type=Integer,Description="Number of 'G' alleles used in tiers 1,2">
                                        #strelka_example.vcf:##FORMAT=<ID=TU,Number=2,Type=Integer,Description="Number of 'T' alleles used in tiers 1,2">

                            # Parse Indels
                            else:
                                strelka_ad.append(float(variant_call[variant_gt_format.index('TAR')].split(',')[0]))
                                strelka_ad.append(float(variant_call[variant_gt_format.index('TIR')].split(',')[0]))
                                #strelka_example.vcf:##FORMAT=<ID=TAR,Number=2,Type=Integer,Description="Reads strongly supporting alternate allele for tiers 1,2"> REFERENCE!!!
                                #strelka_example.vcf:##FORMAT=<ID=TIR,Number=2,Type=Integer,Description="Reads strongly supporting indel allele for tiers 1,2">

                    '''Calculate mean dp and ad'''
                    variant_dp = int(round(sum(variant_dp)/len(variant_dp)))
                    ad_callers = izip_longest(freebayes_ad, mutect_ad, varscan_ad, strelka_ad)
                    for allele in ad_callers:
                        allele_ad = [ad for ad in allele if ad is not None]
                        variant_ad.append(str(int(round(sum(allele_ad)/len(allele_ad)))))

                    print "{var_data};CC={cc}\t{gt_format}\t{gt}:{ad}:{dp}".format(
                        var_data = "\t".join(variant[:8]),
                        cc = caller_count,
                        gt_format = "GT:AD:DP",
                        gt = "0/1",
                        ad = ','.join(variant_ad),
                        dp = variant_dp,
                    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=100, width=200))

    parser.add_argument('--remove_filtered', action='store_true', help='Skip variants with a filter flag other than PASS.')

    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument('-v', '--vcf_file', help='path/to/file.vcf', required=True)
    required_named.add_argument('-t', '--tumor_sample', help='Tumor sample name', required=True)

    args = parser.parse_args()
    melt_somatic_vcf(args.vcf_file, args.remove_filtered, args.tumor_sample)
