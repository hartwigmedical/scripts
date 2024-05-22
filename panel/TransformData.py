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

def overlaps(x, y):
    return max(x[0],y[0]) < min(x[1],y[1])

@dataclass(frozen=True)
class geneLoc:
    gene : str
    chr : str
    start : int
    end : int

def weightedAverage(values,weights):
    weightedSum=0
    for i in range(0,len(values)):
        weightedSum+= values[i]*weights[i]
    return weightedSum/sum(weights)


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
            res['gender'] =str(arr[header.index('gender')])
    return res

def writePurplePurity(purple_purity, sampleId, output_dir='transformed_data'):
  res = getPurityPloidy(purple_purity)

  with open(output_dir + '/' + sampleId + '.purity', 'w') as ft:
    ft.write(str(res['purity']))
  with open(output_dir + '/' + sampleId + '.gender', 'w') as ft:
    ft.write(res['gender'])
  with open(output_dir + '/' + sampleId + '.ploidy', 'w') as ft:
    ft.write(str(res['ploidy']))


def transformRatioFile(cobalt, sampleId,geneLocation, output_dir='transformed_data'):
## xx

    minTumorGCRatio=10
    zeroGCoffset = 0.2

    ft=open(output_dir + '/' + sampleId + '.transformed.cnr', 'w')
    first = True
    header=[]

    linesList = []



    with gzip.open(cobalt,'rt') as fh:

        for lines in fh:
            if first == True:
               lines=lines.rstrip("\n")
               header= lines.split(tsvSplit)
               header.append("end")
               first = False
               ft.write(lines+tsvSplit+'start' + tsvSplit + 'end' + tsvSplit + 'log2' + tsvSplit+'gene\n')
               continue
            linesList.append(lines.rstrip("\n"))
    for lines in linesList:
        lines = lines[3:]
        data = lines.split(tsvSplit)
        tumorGCRatio = float(data[header.index('tumorGCRatio')])
        if (tumorGCRatio < minTumorGCRatio) and (tumorGCRatio > 0):
            minTumorGCRatio = tumorGCRatio
            lines = lines.rstrip("\n")
            if first == True:
                header= lines.split(tsvSplit)
                header.append("end")
                first = False


    for lines in linesList:
        lines = lines[3:]
        data = lines.split(tsvSplit)
        start = data[header.index('position')]
        tumorGCRatio = float(data[header.index('tumorGCRatio')])
        if tumorGCRatio < 0:
            continue
        else:
            if tumorGCRatio == 0:
                logR = math.log2(minTumorGCRatio) - zeroGCoffset
            else:
                logR = math.log2(tumorGCRatio)
            end = str(int (start) + binSize -1)
            chr = data[0]
            gene = [geneLoc.gene for geneLoc in geneLocation if geneLoc.chr==data[0] and overlaps([geneLoc.start,geneLoc.end],[int(start),int(end)])]
            if gene:
                ft.write(lines+tsvSplit+start+tsvSplit+end+tsvSplit+str(logR)+tsvSplit+';'.join(gene)+ '\n')
            else:
                ft.write(lines+tsvSplit+start+tsvSplit+end+tsvSplit+str(logR)+tsvSplit+'NA'+'\n')


def transformGeneFile(purple_gene, sampleId, panelGenes, output_dir='transformed_data'):
    first=True
    geneLocation = []

    ft=open(output_dir + '/' + sampleId + '.transformed.genemetrics.cns', 'w')
    with open(purple_gene,'rt') as fh:
        for lines in fh:
            lines = lines.rstrip('\n')
            if first == True:
                ft.write(lines+tsvSplit+'log2\n')
                first = False
                continue
            lines = lines[3:]
            arr = lines.split(tsvSplit)
            if arr[3] in panelGenes:
                geneLocation.append(geneLoc(arr[3],arr[0],int(arr[1]),int(arr[2])))

            if arr[3] in AmpsDels:
                if AmpsDels[arr[3]] == 'DEL':
                    ft.write(lines+tsvSplit+'-1\n')
                else:
                    ft.write(lines+tsvSplit+'1\n')
    return geneLocation


def genVCFfile(amber, sampleId, output_dir='transformed_data'):
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SampleName
    ft = open(output_dir + '/' + sampleId + '.transformed.vcf', 'w')
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
            SAR = str(round(float(Depth)*float(BAF)) - int(SAF))


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

def getNormCorrection(cobalt_segmented):
    nProbes = []
    sValues = []
    first = True
    with open(cobalt_segmented,'r') as fh:
        for lines in fh:
            lines = lines.rstrip('\n')
            if first == True:
                header = lines.split(tsvSplit)
                first = False
                continue
            arr = lines.split(tsvSplit)
            nProbes.append(float(arr[header.index('n.probes')]))
            sValues.append(float(arr[header.index('mean')]))
    return weightedAverage(sValues,nProbes)


def getNormPurple(purple_segmented, purple_purity):
    first = True
    header=[]

    FClogs = []
    Weights = []

    with open(purple_segmented,'r') as fh:
        for lines in fh:
            lines=lines.rstrip('\n')
            if first == True:
                header=lines.split(tsvSplit)
                first = False
                continue
            arr = lines.split(tsvSplit)
            res = getPurityPloidy(purple_purity)
            purity = res['purity']
            ploidy = res['ploidy']
            copyNumber = float(arr[header.index('copyNumber')])
            FClog = math.log2(max((purity * copyNumber +2*(1-purity))/(purity*ploidy+2*(1-purity)),0.05))
            FClogs.append(FClog)
            Weights.append(int(arr[header.index('depthWindowCount')]))
    return weightedAverage(FClogs,Weights)

def transformSegmentedFile(purple_segmented,purple_purity, sampleId, normCorrection, output_dir='transformed_data'):

# xx
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SampleName
    ft=open(output_dir + '/' + sampleId + '.transformed.cns', 'w')
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
            FClog = math.log2(max((purity * copyNumber +2*(1-purity))/(purity*ploidy+2*(1-purity)),0.05))+normCorrection

            arr.append(str(FClog))
            ft.write(tsvSplit.join(arr)+'\n')
    ft.close()

def getPanelGenes(panel_file):
    panelGenes = dict()
    with open (panel_file,'r') as fh:
        for lines in fh:
            lines=lines.rstrip("\n")
            panelGenes[lines]=1
    return(panelGenes)

@dataclass(frozen=True)
class Config:
    cobaltRatio: str
    cobaltSegmented: str
    amber: str
    purpledriverCatalog: str
    purpleGene: str
    purpleSomatic: str
    purplePurity: str
    sampleId: str
    panelGenes: str
    outputDir: str = "transformed_data"

def parse_args(sys_args):
    parser = argparse.ArgumentParser(prog="TransformData", description="Transform hmftools output to input for reconCNV")

    parser.add_argument("--cobaltRatio", "-r", type=str, required=True,  help="ratio file - ends with cobalt.ratio.tsv.gz")
    parser.add_argument("--cobaltSegmented","-c",type=str,required=True, help="cobalt segment file -ends with cobalt.ratio.pcf")
    parser.add_argument("--amber", "-a", type=str, required=True, help="amber file - end with amber.baf.tsv.gz")
    parser.add_argument("--purpledriverCatalog", "-d", type=str, required=True, help="driver catalog - ends with driver.catalog.somatic.tsv")
    parser.add_argument("--purpleGene", "-g", type=str, required=True, help="gene file - ends with cnv.gene.tsv")
    parser.add_argument("--purpleSomatic", "-s", type=str, required=True, help="somatic file - ends with cnv.somatic.tsv")
    parser.add_argument("--purplePurity", "-p", type=str, required=True, help="purity file - ends with purple.purity.tsv")
    parser.add_argument("--sampleId","-i", type = str, required =True, help="sample Id - used for output file names")
    parser.add_argument("--panelGenes","-t", type = str, required = True, help="panel genes - used to annotate genes")
    parser.add_argument("--outputDir","-o", type = str, required = False, help="output directory",default = 'transformed_data')
    args = parser.parse_args(sys_args)
    return Config(args.cobaltRatio, args.cobaltSegmented, args.amber, args.purpledriverCatalog, args.purpleGene, args.purpleSomatic, args.purplePurity, args.sampleId, args.panelGenes, args.outputDir)


def main(args):
    print(args)
    sampleId=args.sampleId
    readDriverCatalog(args.purpledriverCatalog)

    panelGenes = getPanelGenes(args.panelGenes)
    normCorrection=getNormCorrection(args.cobaltSegmented)
    geneLocation=transformGeneFile(args.purpleGene,sampleId, panelGenes, args.outputDir)
    transformRatioFile(args.cobaltRatio, sampleId, geneLocation, args.outputDir)
    genVCFfile(args.amber, sampleId, args.outputDir)
    normPurple=getNormPurple(args.purpleSomatic, args.purplePurity)
    transformSegmentedFile(args.purpleSomatic,args.purplePurity,sampleId,-normPurple+normCorrection, args.outputDir)
    writePurplePurity(args.purplePurity,sampleId,args.outputDir)

if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s - [%(levelname)-8s] - %(message)s", level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"
    )

    main(parse_args(sys.argv[1:]))
