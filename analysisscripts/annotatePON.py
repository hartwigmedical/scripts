#!/usr/bin/env python

from collections import namedtuple

TAB = '\t'
COLON = ':'
COMMA = ','
EQUALS = '='

REQUIRED_VARIANT_FIELDS = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')

REQUIRED_INFO_DESCRIPTORS = ('ID', 'Number', 'Type', 'Description')
VALID_INFO_TYPES = ('Integer', 'Float', 'Flag', 'Character', 'String')
InformationField = namedtuple('InformationField', REQUIRED_INFO_DESCRIPTORS)
def generateInfoHeader(**kwargs):
    kwargs['Description'] = '"%s"' % kwargs['Description']
    assert kwargs['Type'] in VALID_INFO_TYPES
    return '##INFO=<%s>' % COMMA.join(EQUALS.join((k, str(v))) for k, v in zip(InformationField._fields, InformationField(**kwargs)))

# basic helper for writing valid VCF files with no samples
class VCFWriter():

    def __init__(self, file):
        self.file = file

    def writeMetaInformation(self, metaInformationLines):
        for l in metaInformationLines:
            self.file.write(l)
            self.file.write('\n')

    def writeTuple(self, t):
        self.file.write(TAB.join(t))
        self.file.write('\n')

# basic helper for reading and validating a VCF file
class VCFReader():
    def __init__(self, file):
        self.file = file
        self.headers = []
        self.samples = []
        self.metaInformationLines = []
        self.processHeaders()

    def processHeaders(self):
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
        self.metaInformationLines.append(line)

    def processVariantHeader(self, line):
        self.metaInformationLines.append(line)
        self.headers = line[1:].split(TAB)
        assert self.headers[:8] == list(REQUIRED_VARIANT_FIELDS)
        if len(self.headers) > 8:
            assert self.headers[8] == 'FORMAT'
            self.samples = self.headers[9:]
            # normalise sample names
            for i, h in enumerate(self.headers[9:]):
                self.headers[9+i] = "SAMPLE%i" % i
        self.tuple = namedtuple('VCFReaderVariant', self.headers)

    def readVariant(self):
        line = self.file.readline()
        return self.tuple._make(line.split(TAB)) if line else None

    def getSamples(self):
        return self.samples

    def setReferenceSample(self, sample):
        self.reference_idx = self.headers.index(sample)

    def getReferenceSampleFromVariant(self, variant):
        return variant[self.reference_idx]

def position(variant):
    def chromosomeToNumber(chromosome):
        if chromosome == 'X':
            return 23
        elif chromosome == 'Y':
            return 24
        elif chromosome == 'MT':
            return 25
        else:
            return int(chromosome)
    return ( chromosomeToNumber(variant.CHROM), int(variant.POS) ) if variant is not None else None

class PONAnnotator():

    def __init__(self, outputFile, inputFile, ponFile):
        self._outputFile = outputFile
        self._inputFile = inputFile
        self._ponFile = ponFile

    def updateInfoHeader(self):
        insertIndex = len(self._inputFile.metaInformationLines) - 1 # before #CHROM
        self._inputFile.metaInformationLines.insert(insertIndex,
            generateInfoHeader(
                ID="PON_COUNT",
                Number='A', # one per ALT
                Type="Integer",
                Description="how many samples in the PON had this variant"
            )
        )
        self._outputFile.writeMetaInformation(self._inputFile.metaInformationLines)

    def annotate(self):
        def advanceInput():
            return self._inputFile.readVariant()
        def advancePON():
            return self._ponFile.readVariant()
        
        self.updateInfoHeader()
        inputVariant = advanceInput()
        ponVariant = advancePON()

        while inputVariant is not None:
            inputPos, ponPos = position(inputVariant), position(ponVariant)
            if ponPos is None or inputPos < ponPos:
                self._outputFile.writeTuple(inputVariant)
                inputVariant = advanceInput()
            elif ponPos < inputPos:
                ponVariant = advancePON()
            else:
                # retrieve all possibilities at this position
                altList = []
                while ponVariant is not None and position(ponVariant) == inputPos:
                    altList.append(ponVariant)
                    ponVariant = advancePON()
                # process all inputs at this position
                while inputVariant is not None and position(inputVariant) == inputPos:
                    matches = []
                    for alt in inputVariant.ALT.split(COMMA):
                        matches.append(
                            next(
                                (a.INFO.split(EQUALS)[1] for a in altList if inputVariant.REF == a.REF and alt == a.ALT),
                                '0'
                            )
                        )
                    # write to output
                    if all(m=='0' for m in matches):
                        self._outputFile.writeTuple(inputVariant)
                    else:
                        self._outputFile.writeTuple(
                            inputVariant._replace(INFO=inputVariant.INFO+";PON_COUNT=" + COMMA.join(matches))
                        )
                    inputVariant = advanceInput()

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(
        description="Annotates with a Panel of Normals (PON) VCF file",
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=100, width=200)
    )
    required = parser.add_argument_group('required arguments')
    required.add_argument('-o', '--outputFile', help='output VCF', required=True, type=argparse.FileType('w'))
    required.add_argument('-i', '--inputFile', help='input VCF to annotate', required=True, type=argparse.FileType('r'))
    required.add_argument('-p', '--ponFile', help='pon VCF', required=True, type=argparse.FileType('r'))
    args = parser.parse_args()

    try:
        annotator = PONAnnotator( VCFWriter(args.outputFile), VCFReader(args.inputFile), VCFReader(args.ponFile) )
        annotator.annotate()
    finally: # be a good citizen
        args.outputFile.close()
        args.inputFile.close()
        args.ponFile.close()
