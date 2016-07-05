#!/usr/local/bin/python
import pandas as pd
import numpy as np

###############################################
# VCF CONFIG
VCF_SAMPLE = "CPCT11111111"
VCF_PATH = "/Users/peterpriestley/hmf/analyses/70-30sample/160524/"
VCF_FILE_NAME = VCF_SAMPLE + "R_" + VCF_SAMPLE + "T_merged_somatics.vcf"
SAMPLE_NAMES = {VCF_SAMPLE + 'T.mutect': 'mutect',
                VCF_SAMPLE + 'T.freebayes': 'freebayes',
                'TUMOR.strelka': 'strelka',
                'TUMOR.varscan': 'varscan'}

###############################################

#DEFINE CHR LENGTH
chromosomeLength = {}
chromosomeLength['1'] = 249250622
chromosomeLength['10'] = 135534749
chromosomeLength['11'] = 135006518
chromosomeLength['12'] = 133851897
chromosomeLength['13'] = 115169880
chromosomeLength['14'] = 107349541
chromosomeLength['15'] = 102531394
chromosomeLength['16'] = 90354755
chromosomeLength['17'] = 81195212
chromosomeLength['18'] = 78077250
chromosomeLength['19'] = 59128985
chromosomeLength['2'] = 243199374
chromosomeLength['20'] = 63025522
chromosomeLength['21'] = 48129897
chromosomeLength['22'] = 51304568
chromosomeLength['3'] = 198022431
chromosomeLength['4'] = 191154277
chromosomeLength['5'] = 180915261
chromosomeLength['6'] = 171115068
chromosomeLength['7'] = 159138664
chromosomeLength['8'] = 146364023
chromosomeLength['9'] = 141213432
chromosomeLength['MT'] = 16571
chromosomeLength['X'] = 155270561
chromosomeLength['Y'] = 59373567

class variantType():
    sameAsRef = "Same as Ref"
    missingGenotype = "Missing Genotype"
    indel = "INDEL"
    SNP = "SNP"
    mixed = "MIXED"

class subVariantType():
    none = ""
    insert = "INSERT"
    delete = "DELETE"
    indel = "INDEL"

def intChrom(chrom):
    if chrom == 'X':
        return 23
    elif chrom == 'Y':
        return 24
    elif chrom == 'MT':
        return 25
    else:
        return int(chrom)

def calculateReadDepth(format,genotype):
    formatSplit = format.split(':')
    genotypeSplit = genotype.split(':')
    if 'DP' in formatSplit:
        return int(genotypeSplit[formatSplit.index('DP')])
    else:
        return 1

def calculateSomaticGenotype(infoSplit,caller,aVariantType):
    #Ideally both Normal (ref, hom, het) and somatic (ref,hom,het).
    infoHeaders = [x.split('=')[0] for x in infoSplit]
    if caller == 'strelka' and aVariantType == variantType.indel:
        return infoSplit[infoHeaders.index("SGT")].split('=')[1]
    elif caller == 'strelka' and aVariantType == variantType.SNP:
        return infoSplit[infoHeaders.index("NT")].split('=')[1]
    elif caller == 'varscan':  # VARSCAN - 'SS=1,2,3,4' germline,somatic,LOH, unknown
        return infoSplit[infoHeaders.index("SS")].split('=')[1]
    elif caller == 'freebayes' and 'VT' in infoHeaders :  # FREEBAYES - Better would be GT(Normal)  0/0 = ref;  X/X = het;  X/Y = hom;
        return infoSplit[infoHeaders.index("VT")].split('=')[1]
    elif caller =='mutect':  # Mutect is always het ?
        return "ref-het"
    else:
        return 'unknown'


def calculateQualityScore(infoSplit,caller,aVariantType):
    infoHeaders = [x.split('=')[0] for x in infoSplit]
    if caller == 'strelka' and aVariantType == variantType.indel:
        return infoSplit[infoHeaders.index("QSI_NT")].split('=')[1]
    elif caller == 'strelka' and aVariantType == variantType.SNP:
        return infoSplit[infoHeaders.index("QSS_NT")].split('=')[1]
    elif caller == 'varscan':
        if "VS_SSC" in infoHeaders:
            return infoSplit[infoHeaders.index("VS_SSC")].split('=')[1]
        else:
            return infoSplit[infoHeaders.index("SSC")].split('=')[1]
    elif caller == 'freebayes':
        if "FB_SSC" in infoHeaders:
            return infoSplit[infoHeaders.index("FB_SSC")].split('=')[1]
        else:
            return infoSplit[infoHeaders.index("SSC")].split('=')[1]
    else:  # Mutect has no quality score
        return -1


def calculateAllelicFreq(format,genotype,caller,tumorVariantType,ref,alleleTumor2):
    formatSplit = format.split(':')
    genotypeSplit = genotype.split(':')
    if caller == 'mutect':
        return float(genotypeSplit[formatSplit.index('FA')])
    elif caller == 'varscan':
        return float(genotypeSplit[formatSplit.index('FREQ')].split('%')[0])/100
    else:
        if caller == 'freebayes':
            ao_split = genotypeSplit[formatSplit.index('AO')].split(',')
            ad = ao_split[min(int(genotype[2]),len(ao_split))-1]  #NB - Assume B if GT = A/B
            rd = genotypeSplit[formatSplit.index('RO')].split(',')[0]
        elif caller == 'strelka' and tumorVariantType == variantType.SNP:
            ad = genotypeSplit[formatSplit.index(alleleTumor2 + 'U')].split(',')[0]  #NB - Assume B if GT = A/B
            rd = genotypeSplit[formatSplit.index(ref + 'U')].split(',')[0]
        elif caller == 'strelka' and tumorVariantType == variantType.indel:
            ad = genotypeSplit[formatSplit.index('TIR')].split(',')[0]
            rd = genotypeSplit[formatSplit.index('TAR')].split(',')[0]
        elif caller == 'melted':
            ad, rd = genotypeSplit[formatSplit.index('AD')].split(',')[:2]
        else:
            return -1

        if float(ad) == 0:
            return 0
        else:
            return float(ad) / (float(rd) + float(ad))   #is this correct, or should it be /DP?

class genotype:

    def __init__(self,chrom,pos,caller,ref,alt,info,format,vennSegment,inputGenotype):
        altsplit = (ref + ","+ alt).split(',')
        self.tumorVariantSubType = subVariantType.none

        if inputGenotype[:3] == "./.":
            self.tumorVariantType = variantType.missingGenotype
        elif inputGenotype[:3] == "0/0" or alt == ".":   #STRELKA unfiltered
            self.tumorVariantType = variantType.sameAsRef
        else:
            if inputGenotype[1] != '/':  #STRELKA unfiltered case
                alleleTumor1 = altsplit[1]
                alleleTumor2 = altsplit[1]
            else:
                alleleTumor1 = altsplit[int(inputGenotype[0])]
                alleleTumor2 = altsplit[int(inputGenotype[2])]
            if len(alleleTumor1) == len(alleleTumor2) and len(alleleTumor1) == len(ref):
                self.tumorVariantType = variantType.SNP
            else:
                self.tumorVariantType = variantType.indel
                if len(alleleTumor1) <= len(ref) and len(alleleTumor2) <= len(ref):
                    self.tumorVariantSubType = subVariantType.delete
                elif len(alleleTumor1) >= len(ref) and len(alleleTumor2) >= len(ref):
                    self.tumorVariantSubType = subVariantType.insert
                else:
                    self.tumorVariantSubType = subVariantType.indel

            self.allelicFreq = calculateAllelicFreq(format, inputGenotype, caller, self.tumorVariantType, ref,alleleTumor2)
            self.readDepth = calculateReadDepth(format,inputGenotype)
            infoSplit = info.split(';')
            self.qualityScore = float(calculateQualityScore(infoSplit,caller,self.tumorVariantType))
            self.somaticGenotype = calculateSomaticGenotype(infoSplit,caller,self.tumorVariantType)
            self.allele = alleleTumor2

class somaticVariant:

    variantInfo = []
    bedItem = []

    def __init__(self, chrom, pos, id, ref, alt, filter, info,format,inputGenotypes,useFilter,useBed,aBedReverse):

        #Find the 1st Bed region with maxPos > variantPos
        if aBedReverse:
            if not somaticVariant.bedItem:
                somaticVariant.bedItem = aBedReverse.pop()
            while intChrom(chrom) > intChrom(somaticVariant.bedItem[0]) or (intChrom(chrom) == intChrom(somaticVariant.bedItem[0]) and int(pos) > int(somaticVariant.bedItem[2])):
                somaticVariant.bedItem = aBedReverse.pop()
        else:
            somaticVariant.bedItem = []

        #Only use if inside the next BED region
        if (somaticVariant.bedItem and int(somaticVariant.bedItem[1])<int(pos) and somaticVariant.bedItem[0]==chrom) or not useBed:
            if filter == "PASS" or filter == "." or useFilter == False:

                tumorCallerCountSNP = 0
                tumorCallerCountIndel = 0
                tumorCallerCountSubTypeIndel = 0
                tumorCallerCountSubTypeDelete = 0
                tumorCallerCountSubTypeInsert = 0
                variantGenotypes = {}
                vennSegment = ""

                formatSplit = format.split(';')
                for i in range(len(formatSplit)):
                    formatItem = formatSplit[i].split('=')
                    if formatItem[0] == "set":
                        vennSegment = formatItem[1]

                for key in inputGenotypes.iterkeys():
                    variantGenotypes[key] = genotype(chrom, pos, key, ref, alt, info,format, vennSegment, inputGenotypes[key])

                for key, value in variantGenotypes.items():
                    if value.tumorVariantType == variantType.SNP:
                        tumorCallerCountSNP += 1
                    if value.tumorVariantType == variantType.indel:
                        tumorCallerCountIndel += 1
                        if value.tumorVariantSubType == subVariantType.delete:
                            tumorCallerCountSubTypeDelete += 1
                        if value.tumorVariantSubType == subVariantType.insert:
                            tumorCallerCountSubTypeInsert += 1
                        if value.tumorVariantSubType == subVariantType.indel:
                            tumorCallerCountSubTypeIndel += 1

                ############### Pandas Prep ####################
                #PREPARE FIELDS:
                if chrom[:3] == 'chr':
                    chrom = chrom[3:]
                posPercent = float(pos) / chromosomeLength[chrom]
                numCallers = tumorCallerCountSNP + tumorCallerCountIndel
                mySubVariantType = ""
                if tumorCallerCountIndel > 0:
                    if tumorCallerCountSNP == 0:
                        myVariantType = variantType.indel
                    else:
                        myVariantType = variantType.mixed
                    if tumorCallerCountSubTypeDelete > 0:  # this logic still needs work
                        if tumorCallerCountSubTypeInsert > 0:
                            mySubVariantType = subVariantType.indel
                        else:
                            mySubVariantType = subVariantType.delete
                    elif tumorCallerCountSubTypeInsert > 0:
                        mySubVariantType = subVariantType.insert
                    elif tumorCallerCountIndel > 0:
                        mySubVariantType = subVariantType.indel
                elif tumorCallerCountSNP > 0:
                    myVariantType = variantType.SNP
                else:
                    myVariantType = variantType.missingGenotype

                # Append to list
                somaticVariant.variantInfo.append(
                    [chrom, pos, chrom + ':' + pos, intChrom(chrom) + posPercent, ref, vennSegment, numCallers,
                     myVariantType, mySubVariantType,filter])
                for caller, variantGenotype in variantGenotypes.items():
                    if variantGenotype.tumorVariantType == variantType.indel or variantGenotype.tumorVariantType == variantType.SNP:
                        callerSpecificFields = [variantGenotype.allele, variantGenotype.allelicFreq, variantGenotype.readDepth,
                                                variantGenotype.qualityScore,variantGenotype.somaticGenotype]
                    else:
                        callerSpecificFields = ['', '', '','','']
                    somaticVariant.variantInfo[-1] = somaticVariant.variantInfo[-1] + callerSpecificFields
                #########################################


def loadVaraintsFromVCF(aPath, aVCFFile,sampleNames,aPatientName,useFilter,useBed=False,aBed=[]):
    print "reading vcf file:",aVCFFile

    if useBed == True:
        aBed.reverse()

    variants = []
    with open(aPath + aVCFFile, 'r') as f:
        i=0
        for line in f:
            line = line.strip('\n')
            myGenotypes = {}
            a = [x for x in line.split('\t')]
            if a[0] == '#CHROM':
                headers = a[9:]
                header_index = {}
                for sampleName,sampleLabel in sampleNames.iteritems():
                    for index, header in enumerate(headers):
                        if sampleName == header:
                            header_index[sampleLabel] = index
                            break
                    if not header_index.has_key(sampleLabel):
                        print 'Error - missing sample inputs'
                        return -1
            if a[0][:1] != '#':
                variant_calls = a[9:]
                for caller,index in header_index.iteritems():
                    myGenotypes[caller] = variant_calls[index]
                variants.append(somaticVariant(a[0], a[1], a[2], a[3], a[4], a[6], a[7], a[8],myGenotypes, useFilter, useBed,aBed))
                i += 1
                if i% 100000 == 0:
                    print "reading VCF File line:",i

    #Reset bed item
    somaticVariant.bedItem = []

    print "Number variants loaded:",len(somaticVariant.variantInfo)

    df = pd.DataFrame(somaticVariant.variantInfo)
    myColumnList = ['chrom', 'pos', 'chromPos','chromFrac', 'ref', 'vennSegment','numCallers','variantType',
                    'variantSubType','filter']
    for caller in header_index.iterkeys():
        myColumnList = myColumnList + [caller + 'allele',caller+ 'allelicFreq',caller+'readDepth',caller+'qualityScore',caller+'somaticGenotype']
    df.columns = (myColumnList)
    df['patientName'] = aPatientName
    # Need to empty genotype.variantInfo in case we need to load multiple files
    del somaticVariant.variantInfo[:]
    return df

def loadBEDFile(aPath, aBEDFile):
    print "reading BED file"
    myBed = []
    with open(aPath + aBEDFile, 'r') as f:
        for line in f:
            line = line.strip('\n')
            splitLine = line.split('\t')
            if splitLine[0] != 'chrom':
                myBed.append(splitLine)
    print "Bed File Loaded"
    return myBed

def printStatistics(df):
    #Calculate 2+_caller precision and sensitivity
    outputdata = []
    for columnName in list(df):
        if columnName.endswith('allele'):
            myCaller = columnName[:-6]
            variantTypes = df[(df[myCaller + 'allele'] != '')].variantType.unique()
            for variantType in variantTypes:
                truePositives = len(
                    df[(df[myCaller + 'allele'] != '') & (df['numCallers'] >= 2) & (df['variantType'] == variantType)])
                falseNegatives = len(
                    df[(df[myCaller + 'allele'] == '') & (df['numCallers'] >= 2) & (df['variantType'] == variantType)])
                positives = len(df[(df[myCaller + 'allele'] != '') & (df['variantType'] == variantType)])
                truthSet = truePositives + falseNegatives
                if positives > 0:
                    outputdata.append(
                        [variantType, myCaller, truthSet, truePositives, positives - truePositives, falseNegatives,
                         round(truePositives / float(positives), 4), round(truePositives / float(truthSet), 4)])

    outputDF = pd.DataFrame(outputdata)
    outputDF.columns = (['variantType', 'caller', 'truthSet', 'truePositives', 'falsePositives', 'falseNegatives',
                         'precision_2+_callers','sensitivity_2+_callers'])
    print '\n2+ caller precision and sensitivity\n'
    print outputDF.sort_values(['variantType', 'caller'])

    # calculate # of variants by variant type
    print '\n# of variants by sub type\n'
    print df[['pos', 'variantSubType']].groupby(['variantSubType']).agg('count')

    # calculate # of callers
    print '\n# of variants by number of callers and variant type\n'
    df_pivot = df[['numCallers', 'pos', 'variantType']].groupby(['variantType', 'numCallers' ]).agg('count')
    print df_pivot.groupby(level=0).transform(lambda x: x / x.sum())


if __name__ == "__main__":
    df = loadVaraintsFromVCF(VCF_PATH,VCF_FILE_NAME,SAMPLE_NAMES,True,VCF_SAMPLE,False)
    printStatistics(df)




