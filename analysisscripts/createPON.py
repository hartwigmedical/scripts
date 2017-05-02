#!/usr/bin/env python

import heapq
import re

TAB = '\t'
FORMAT_SEP = ':'

REFERENCE_SAMPLE_SUFFIXES = ( 'R', 'BL' ) # TODO maybe make this command line argument

# TODO perhaps a bit heavy handed
GT_MATCH = re.compile('[1-9]+')
def GenotypeFilter(vcf, variant):
    assert variant.FORMAT.startswith('GT')
    gt = vcf.getReferenceSampleFromVariant(variant).split(FORMAT_SEP, 1)[0]
    return GT_MATCH.search(gt) is not None

# basic helper for reading and validating a VCF file
class VCFReader():
    def __init__(self, file):
        self.file = file
        self.headers = []
        self.samples = []
        self.processHeaders()

    def processHeaders(self):
        assert self.file.readline().startswith('##fileformat=VCFv4.')
        while True:
            line = self.file.readline().rstrip()
            assert len(line)
            if line[0:2] == "##":
                self.processMetaInformation(line)
            elif line[0] == "#":
                self.processVariantHeader(line)
                return
            else:
                raise Exception("VCF format derailment")

    def processMetaInformation(self, line):
        pass

    def processVariantHeader(self, line):
        self.headers = line[1:].split(TAB)
        assert self.headers[:8] == ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        if len(self.headers) > 8:
            assert self.headers[8] == 'FORMAT'
            self.samples = self.headers[9:]
        from collections import namedtuple
        self.tuple = namedtuple('Variant', self.headers)

    def getSamples(self):
        return self.samples

    def setReferenceSample(self, sample):
        self.reference_idx = self.headers.index(sample)

    def getReferenceSampleFromVariant(self, variant):
        return variant[self.reference_idx]

    def readVariant(self):
        line = self.file.readline()
        return self.tuple._make(line.split(TAB)) if line else None

    def readVariantMatchingFilter(self, filter):
        variant = self.readVariant()
        while variant and not filter(variant):
            variant = self.readVariant()
        return variant

class PONGenerator():

    def __init__(self, outputFile, minCountThreshold):
        self._heap = []
        self._outputFile = outputFile
        self._minCountThreshold = minCountThreshold

    def merge(self, vcf_readers):

        for vcf in vcf_readers:
            # find the reference sample
            reference = next(sample for sample in vcf.getSamples() for suffix in REFERENCE_SAMPLE_SUFFIXES if sample.endswith(suffix))
            vcf.setReferenceSample(reference)
            # prime the heap
            self.pushVariantToHeap(vcf.readVariantMatchingFilter(lambda x : GenotypeFilter(vcf, x)), vcf)

        previousCount = 0
        currentIdx, variant, vcf = self._heap[0]

        while self._heap:

            previousIdx, previousVariant = currentIdx, variant
            currentIdx, variant, vcf = heapq.heappop(self._heap)

            if previousIdx < currentIdx:
                self.writeToOutput(previousVariant, previousCount)
                previousCount = 1
            else:
                previousCount += 1

            self.pushVariantToHeap(vcf.readVariantMatchingFilter(lambda x : GenotypeFilter(vcf, x)), vcf)

        self.writeToOutput(variant, previousCount)

    def pushVariantToHeap(self, variant, vcf):
        def chromosomeToNumber(chromosome):
            if chromosome == 'X':
                return 23
            elif chromosome == 'Y':
                return 24
            elif chromosome == 'MT':
                return 25
            else:
                return int(chromosome)
        if variant:
            heapq.heappush(self._heap,
                (
                    (chromosomeToNumber(variant.CHROM), int(variant.POS)), # sorted on this field
                    variant,
                    vcf
                )
            )

    def writeToOutput(self, variant, count):
        if count < self._minCountThreshold:
            return
        self._outputFile.write(
            TAB.join((
                variant.CHROM,
                variant.POS,
                variant.REF,
                variant.ALT[0], # NB - only pushing The 1st ALT into PON file
                str(count)
            )) + '\n'
        )

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(
        description="Generates a Panel of Normals (PON) file",
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=100, width=200)
    )
    required = parser.add_argument_group('required arguments')
    required.add_argument('-m', '--minCountThreshold', help='minCount to add to PON output.  eg: 2', required=True, type=int)
    required.add_argument('-o', '--outputFile', help='output file name', required=True, type=argparse.FileType('w'))
    required.add_argument('-i', '--inputFiles', nargs='+', help='list of vcf files to merge, use bash ', required=True, type=argparse.FileType('r'))
    args = parser.parse_args()

    try:
        generator = PONGenerator(args.outputFile, args.minCountThreshold)
        generator.merge([ VCFReader(f) for f in args.inputFiles ])
    finally: # be a good citizen
        args.outputFile.close()
        for f in args.inputFiles: f.close()
