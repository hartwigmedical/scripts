#!/usr/bin/env python

"""
filterFreebayes.py
Applies BCBIO developed filters to the freebayes somatic output.
2 main filters applied:
1. LOD Tumor, LOD Normal > 3.5
2. Freq of Tumor / Freq Normal > 2.7
"""

import argparse
from itertools import izip_longest
import re

TUMOR_PARTS_INDEX = 10
NORMAL_PARTS_INDEX = 9

def customFilterFreebayes(vcf_file):

    try:
        f = open(vcf_file, 'r')
    except IOError:
        print "Can't open vcf file: {0}".format(vcf_file)
    else:
        with f:
            for line in f:
                line = line.strip('\n')
                # Code copied from BCBIO Freebayes pipeline.
                # Thresholds are like phred scores, so 3.5 = phred35
                tumor_thresh, normal_thresh = 3.5, 3.5
                if line.startswith("#CHROM"):
                    headers = [('##FILTER=<ID=REJECT,Description="Not somatic due to normal call frequency '
                                'or phred likelihoods: tumor: %s, normal %s.">')
                               % (int(tumor_thresh * 10), int(normal_thresh * 10))]
                    print "\n".join(headers) + "\n" + line
                elif line.startswith("#"):
                    print line
                else:
                    parts = line.split("\t")
                    if not (_check_lods(parts, tumor_thresh, normal_thresh) and _check_freqs(parts)):
                        if parts[6] in set([".", "PASS"]):
                            parts[6] = "REJECT"
                        else:
                            parts[6] += ";REJECT"
                    else:
                        line = "\t".join(parts)
                        print line

### Filtering

def _check_lods(parts, tumor_thresh, normal_thresh):
    """Ensure likelihoods for tumor and normal pass thresholds.

    Skipped if no FreeBayes GL annotations available.
    """
    try:
        gl_index = parts[8].split(":").index("GL")
    except ValueError:
        return True
    try:
        tumor_gls = [float(x) for x in parts[TUMOR_PARTS_INDEX].split(":")[gl_index].split(",") if x != "."]
        if tumor_gls:
            tumor_lod = max(tumor_gls[i] - tumor_gls[0] for i in range(1, len(tumor_gls)))
        else:
            tumor_lod = -1.0
    # No GL information, no tumor call (so fail it)
    except IndexError:
        tumor_lod = -1.0
    try:
        normal_gls = [float(x) for x in parts[NORMAL_PARTS_INDEX].split(":")[gl_index].split(",") if x != "."]
        if normal_gls:
            normal_lod = min(normal_gls[0] - normal_gls[i] for i in range(1, len(normal_gls)))
        else:
            normal_lod = normal_thresh
    # No GL inofmration, no normal call (so pass it)
    except IndexError:
        normal_lod = normal_thresh
    return normal_lod >= normal_thresh and tumor_lod >= tumor_thresh

def _check_freqs(parts):
    """Ensure frequency of tumor to normal passes a reasonable threshold.

    Avoids calling low frequency tumors also present at low frequency in normals,
    which indicates a contamination or persistent error.
    """
    thresh_ratio = 2.7
    try:  # FreeBayes
        ao_index = parts[8].split(":").index("AO")
        ro_index = parts[8].split(":").index("RO")
    except ValueError:
        ao_index, ro_index = None, None
    try:  # VarDict
        af_index = parts[8].split(":").index("AF")
    except ValueError:
        af_index = None
    if af_index is None and ao_index is None:
        raise NotImplementedError("Unexpected format annotations: %s" % parts[0])
    def _calc_freq(item):
        try:
            if ao_index is not None and ro_index is not None:
                ao = sum([int(x) for x in item.split(":")[ao_index].split(",")])
                ro = int(item.split(":")[ro_index])
                freq = ao / float(ao + ro)
            elif af_index is not None:
                freq = float(item.split(":")[af_index])
        except (IndexError, ValueError, ZeroDivisionError):
            freq = 0.0
        return freq
    tumor_freq, normal_freq = _calc_freq(parts[TUMOR_PARTS_INDEX]), _calc_freq(parts[NORMAL_PARTS_INDEX])
    return normal_freq <= 0.001 or normal_freq <= tumor_freq / thresh_ratio



if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=100, width=200))

    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument('-v', '--vcf_file', help='path/to/file.vcf', required=True)

    args = parser.parse_args()
    customFilterFreebayes(args.vcf_file)