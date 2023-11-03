#!/usr/bin/python3

import sys
import gzip
import math
from dataclasses import dataclass
import logging
import argparse

vcf_header = '''##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
##INFO=<ID=SAF,Number=A,Type=Integer,Description="Number of alternate observations on the forward strand">
##INFO=<ID=SAR,Number=A,Type=Integer,Description="Number of alternate observations on the reverse strand">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SampleName'''

AmpsDels=dict()

def readDriverCatalog(purple_driver):
    first = True
    with open(purple_driver,'rt') as fh:
        for lines in fh:
            lines = lines.rstrip('\n')
            if first == True:
                header= lines.split(tsvSplit)
                first = False
                continue
            arr = lines.split(tsvSplit)
            if arr[5] == 'AMP' or arr[5] == 'DEL' or arr[5] == 'PARTIAL_AMP':
                AmpsDels[arr[2]] = arr[5]



binSize = 1000

tsvSplit = '\t'



def getPurityPloidy(purple_purity):
    first = True
    header = []
    res = dict()
    with open(purple_purity,'r') as fh:
        for lines in fh:
            lines = lines.rstrip('\n')
            if first == True:
                header = lines.split(tsvSplit)
                first = False
                continue
            arr = lines.split(tsvSplit)
            res['purity'] = float(arr[header.index('purity')])
            res['ploidy'] = float(arr[header.index('ploidy')])
    return res
def transformRatioFile(cobalt, sampleId):
## xx
    ft=open('transformed_data/'+sampleId+'.transformed.cnr','w')
    first = True
    header=[]
    with gzip.open(cobalt,'rt') as fh:
        for lines in fh:
            lines = lines.rstrip("\n")
            if first == True:
                header= lines.split(tsvSplit)
                header.append("end")
                first = False
                ft.write(lines+tsvSplit+'start' + tsvSplit + 'end' + tsvSplit + 'log2\n')
                continue
            else:
                data = lines.split(tsvSplit)
                lines = lines[3:]
                start = data[header.index('position')]
                tumorGCRatio = float(data[header.index('tumorGCRatio')])
                if tumorGCRatio < 0:
                    continue
                else:
                    logR = math.log2(tumorGCRatio)
                end = str(int (start) + binSize -1)
                ft.write(lines+tsvSplit+start+tsvSplit+end+tsvSplit+str(logR)+'\n')

def transformGeneFile(purple_gene, sampleId):
    first=True
    ft=open('transformed_data/'+sampleId+'.transformed.genemetrics.cns','w')
    with open(purple_gene,'rt') as fh:
        for lines in fh:
            lines = lines.rstrip('\n')
            if first == True:
                ft.write(lines+tsvSplit+'log2\n')
                first = False
            lines = lines[3:]
            arr = lines.split(tsvSplit)
            if arr[3] in AmpsDels:
                if AmpsDels[arr[3]] == 'DEL':
                    ft.write(lines+tsvSplit+'-1\n') 
                else:
                    ft.write(lines+tsvSplit+'1\n')
def genVCFfile(amber, sampleId):
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SampleName
    ft = open('transformed_data/'+sampleId+'.transformed.vcf','w')
    ft.write(vcf_header+"\n")
    first = True
    header = []
    with gzip.open(amber,'rt') as fh:
        for lines in fh:
            lines = lines.rstrip('\n')
            if first == True:
                header = lines.split(tsvSplit)
                first = False
                continue
            arr = lines.split(tsvSplit)
            Depth = arr[header.index("tumorDepth")]

            BAF =  arr[header.index("tumorBAF")]
            SAF = str(int(round(int(Depth)*float(BAF)*0.5)))
            SAR = str(int(float(Depth)*float(BAF)) - int(SAF))



            CHROM = arr[header.index('chromosome')]
            CHROM = CHROM[3:]
            POS = arr[header.index('position')]
            ID = '.'
            REF = 'N'
            ALT = 'N'
            QUAL = '10000'
            FILTER = '.'
            INFO = 'DP=' + arr[header.index("tumorDepth")]+';'+"SAF=" + SAF+";SAR="+SAR

            FORMAT = 'DP:AD 1283:880,403'

            SampleName = 'NA'
            ft.write(CHROM+tsvSplit+POS + tsvSplit+ID+tsvSplit+REF+tsvSplit+ALT+tsvSplit+QUAL+tsvSplit+FILTER+tsvSplit+INFO+tsvSplit+FORMAT+tsvSplit+SampleName+'\n')

def transformSegmentedFile(purple_segmented,purple_purity, sampleId):

# xx
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SampleName
    ft=open('transformed_data/' + sampleId + '.transformed.cns','w')
    first = True
    header=[]
    with open(purple_segmented,'r') as fh:
        for lines in fh:
            lines=lines.rstrip('\n')
            if first == True:
                header=lines.split(tsvSplit)
                header.append('gene')
                header.append('cell_frac')
                header.append('log2')
                header[header.index("minorAlleleCopyNumber")] = 'minor_cn'
                header[header.index("majorAlleleCopyNumber")] = 'major_cn'

                ft.write(tsvSplit.join(header)+"\n")
                first = False
                continue
            lines = lines[3:]
            arr = lines.split(tsvSplit)
            arr.append("NA") #gene
            arr.append("1") #cell_frac 
            res = getPurityPloidy(purple_purity)
            purity = res['purity']
            ploidy = res['ploidy']
            copyNumber = float(arr[header.index('copyNumber')])
            FClog = math.log2((purity * copyNumber +2*(1-purity))/(purity*ploidy+2*(1-purity)))

            arr.append(str(FClog))
            ft.write(tsvSplit.join(arr)+'\n')

@dataclass(frozen=True)
class Config:
    cobaltRatio: str
    amber: str
    purpledriverCatalog: str
    purpleGene: str
    purpleSomatic: str
    purplePurity: str
    sampleId: str

def parse_args(sys_args):
    parser = argparse.ArgumentParser(prog="TransformData", description="Transform hmftools output to input for reconCNV")

    parser.add_argument("--cobaltRatio", "-r", type=str, required=True,  help="ratio file - ends with cobalt.ratio.tsv.gz")
    parser.add_argument("--amber", "-a", type=str, required=True, help="amber file - end with amber.baf.tsv.gz")
    parser.add_argument("--purpledriverCatalog", "-d", type=str, required=True, help="driver catalog - ends with driver.catalog.somatic.tsv")
    parser.add_argument("--purpleGene", "-g", type=str, required=True, help="gene file - ends with cnv.gene.tsv")
    parser.add_argument("--purpleSomatic", "-s", type=str, required=True, help="somatic file - ends with cnv.somatic.tsv")
    parser.add_argument("--purplePurity", "-p", type=str, required=True, help="purity file - ends with purple.purity.tsv")
    parser.add_argument("--sampleId","-i", type = str, required =True, help="sample Id - used for output file names")
    args = parser.parse_args(sys_args)
    return Config(args.cobaltRatio, args.amber, args.purpledriverCatalog, args.purpleGene, args.purpleSomatic, args.purplePurity, args.sampleId)


def main(args):
    print(args)
    sampleId=args.sampleId

    transformRatioFile(args.cobaltRatio, sampleId)
    genVCFfile(args.amber, sampleId)
    readDriverCatalog(args.purpledriverCatalog)
    transformGeneFile(args.purpleGene, sampleId)
    transformSegmentedFile(args.purpleSomatic,args.purplePurity,sampleId)
if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s - [%(levelname)-8s] - %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
    )

    main(parse_args(sys.argv[1:]))



