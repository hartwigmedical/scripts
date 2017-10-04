#!/usr/bin/env python

import heapq
import pprint
import warnings
from collections import namedtuple

TAB = '\t'
COLON = ':'
COMMA = ','
EQUALS = '='
FORWARDSLASH = '/'

REFERENCE_SAMPLE_SUFFIXES = ( 'R', 'BL' ) # TODO maybe make this command line argument

REQUIRED_VARIANT_FIELDS = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')

def ScriptName():
    from os import path
    return path.basename(__file__)

# basic helper for writing valid VCF files with no samples
class VCFWriter():

    VARIANT_MISSING_FIELD = '.'
    REQUIRED_INFO_DESCRIPTORS = ('ID', 'Number', 'Type', 'Description')
    VALID_INFO_TYPES = ('Integer', 'Float', 'Flag', 'Character', 'String')

    InformationField = namedtuple('InformationField', REQUIRED_INFO_DESCRIPTORS)
    Variant = namedtuple('Variant', REQUIRED_VARIANT_FIELDS)

    def __init__(self, file):
        self.file = file
        self.writeFileHeader()
        self.writeFileDateHeader()
        self.writeSourceHeader()

    def writeFileHeader(self):
        self.file.write('##fileformat=VCFv4.1\n')

    def writeFileDateHeader(self):
        from datetime import date
        self.file.write('##fileDate=%s\n' % date.today().strftime("%Y%m%d"))

    def writeSourceHeader(self):
        self.file.write('##source=%s\n' % ScriptName())

    def writeInfoHeader(self, **kwargs):
        kwargs['Description'] = '"%s"' % kwargs['Description']
        assert kwargs['Type'] in self.VALID_INFO_TYPES
        self.file.write( '##INFO=<%s>\n' % COMMA.join(EQUALS.join((k, str(v))) for k, v in zip(self.InformationField._fields, self.InformationField(**kwargs))) )

    def writeCustomHeader(self, name, value):
        self.file.write( '##%s=%s\n' % (name, value) )

    def writeVariantHeader(self):
        self.file.write('#%s\n' % TAB.join(REQUIRED_VARIANT_FIELDS))

    def writeVariant(self, **kwargs):
        for field in REQUIRED_VARIANT_FIELDS: kwargs.setdefault(field, self.VARIANT_MISSING_FIELD)
        self.file.write(TAB.join(self.Variant(**kwargs)) + '\n')

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
        assert self.headers[:8] == list(REQUIRED_VARIANT_FIELDS)
        if len(self.headers) > 8:
            assert self.headers[8] == 'FORMAT'
            self.samples = self.headers[9:]
        self.tuple = namedtuple('VCFReaderVariant', self.headers)

    def readVariant(self):
        line = self.file.readline()
        return self.tuple._make(line.rstrip().split(TAB)) if line else None

    def readVariantMatchingFilter(self, filter):
        variant = self.readVariant()
        while variant and not filter(variant):
            variant = self.readVariant()
        return variant

    def getSamples(self):
        return self.samples

    def setReferenceSample(self, sample):
        self.reference_idx = self.headers.index(sample)

    def getReferenceSampleFromVariant(self, variant):
        return variant[self.reference_idx]

class PONGenerator():

    def __init__(self, outputFile, minCountThreshold):
        self._heap = []
        self._outputFile = outputFile
        self._minCountThreshold = minCountThreshold
        self._outputFile.writeInfoHeader(ID="PON_COUNT", Number=1, Type="Integer", Description="how many samples had the variant")
        self._outputFile.writeInfoHeader(ID="PON_HET_COUNT", Number=1, Type="Integer", Description="how many samples had the variant in heterozygous genotype")
        self._outputFile.writeInfoHeader(ID="PON_HOM_COUNT", Number=1, Type="Integer", Description="how many samples had the variant in homozygous genotype")

    def merge(self, vcf_readers):
        def readAndPushVariant(vcf):
            done = False
            while not done:
                var = vcf.readVariant()
                if var is None: return
                done = self.pushVariantToHeap(var, vcf)

        samples = set()
        for vcf in vcf_readers:
            # find the reference sample
            sample = next(s for s in vcf.getSamples() for suffix in REFERENCE_SAMPLE_SUFFIXES if s.endswith(suffix))
            if sample in samples: # check it is unique
                continue
            samples.add(sample)
            vcf.setReferenceSample( sample )
            # prime the heap
            print "[DEBUG] Start pushing variants for vcf:" + str(vcf.file.name)
            readAndPushVariant(vcf)

        # finalise headers
        self._outputFile.writeCustomHeader("PonInputSamples", COMMA.join(samples))
        self._outputFile.writeVariantHeader()

        previousCount = 0
        previousHetCount = 0
        previousHomCount = 0
        
        location, variant, alt, genotype_of_alt, vcf = self._heap[0]

        while self._heap:

            previousLocation, previousVariant, previousAlt = location, variant, alt
            location, variant, alt, genotype_of_alt, vcf = heapq.heappop(self._heap)

            if location > previousLocation:
                self.writeToOutput(previousVariant, previousAlt, previousCount, previousHetCount, previousHomCount)
                previousCount = 1
                if genotype_of_alt == 'HET':
                    previousHetCount = 1
                    previousHomCount = 0
                elif genotype_of_alt == 'HOM':
                    previousHetCount = 0
                    previousHomCount = 1
                else:
                    previousHetCount = 0
                    previousHomCount = 0
            else:
                previousCount += 1
                if genotype_of_alt == 'HET':
                    previousHetCount += 1
                elif genotype_of_alt == 'HOM':
                    previousHomCount += 1
            
            if vcf:
                readAndPushVariant(vcf)

        self.writeToOutput(variant, alt, previousCount, previousHetCount, previousHomCount)

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

        # check that reference sample shows the alt in GT field
        #print "[DEBUG] Variant: " + str(variant.POS)

        alleles = [variant.REF]
        alleles.extend( variant.ALT.split(COMMA) )
        
        refSampleGenotypeString = vcf.getReferenceSampleFromVariant(variant).split(COLON, 1)[0]
        alts = [ alt for alt_idx, alt in enumerate(variant.ALT.split(COMMA), start=1) ]
            
        # 117157372
        for idx, alt in enumerate(alts):
            
            alt_int = idx+1
            alt_count = refSampleGenotypeString.count( str(alt_int) )
            genotype_of_alt = 'NA'
                
            if alt_count == 2:
                genotype_of_alt = 'HOM'
            elif alt_count == 1:
                genotype_of_alt = 'HET'
            elif alt_count == 0:
                genotype_of_alt = 'NotPresent'
                continue
            else:
                raise Exception('Unsupported alt count (should be one of 0/1/2): ' + str(alt_count) )
            
            heapq.heappush(self._heap,
                (
                    (chromosomeToNumber(variant.CHROM), int(variant.POS), hash(variant.REF), hash(alt)), # location tuple, sorted on this field
                    variant,
                    alt,
                    genotype_of_alt,
                    vcf if idx==len(alts)-1 else None
                )
            )
            
        return False

    def writeToOutput(self, variant, alt, totalcount, hetcount, homcount):
        
        chrom=variant.CHROM
        pos=variant.POS
        ref=variant.REF
        
        # skip if not enough samples have the alt
        if totalcount < self._minCountThreshold:
            return
            
        if totalcount != hetcount+homcount:
            warnings.warn( '[ERR] Somehow total is not het+hom: TOTAL=%s HET=%s HOM=%s at CHROM=%s POS=%s' % ( str(totalcount), str(hetcount), str(homcount), chrom, pos ) )
            
        self._outputFile.writeVariant(
            CHROM = variant.CHROM,
            POS = variant.POS,
            REF = variant.REF,
            ALT = alt,
            FILTER = 'PASS',
            INFO = 'PON_COUNT=%i;PON_HET_COUNT=%i;PON_HOM_COUNT=%i' % (totalcount, hetcount, homcount)
        )

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(
        description="Generates a Panel of Normals (PON) VCF file",
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=100, width=200)
    )
    required = parser.add_argument_group('required arguments')
    required.add_argument('-m', '--minCountThreshold', help='minCount to add to PON output. eg: 2', required=True, type=int)
    required.add_argument('-o', '--outputFile', help='output file name', required=True, type=argparse.FileType('w'))
    required.add_argument('-i', '--inputFiles', nargs='+', help='list of vcf files to merge', required=True, type=argparse.FileType('r'))
    args = parser.parse_args()

    try:
        generator = PONGenerator( VCFWriter(args.outputFile), args.minCountThreshold )
        generator.merge( VCFReader(f) for f in args.inputFiles )
    finally: # be a good citizen
        args.outputFile.close()
        for f in args.inputFiles: f.close()
