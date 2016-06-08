#!/usr/local/bin/python
import pandas as pd

ROUNDING_PRECISION = 4
POS_PERCENTILE_BUCKET = 0.1
RUN_PANDAS = True
Path = "/Users/peterpriestley/hmf/analyses/70-30sample/160401Schuberg/"
VCFFile = "combined.vcf"
vcfMelted = 'False'
SAMPLE_NAMES = {'CPCT11111111T.mutect': 'mutect', \
                'CPCT11111111T.freebayes': 'freebayes', \
                'TUMOR.strelka': 'strelka', \
                'TUMOR.varscan': 'varscan'}
BED_PATH = "/Users/peterpriestley/hmf/analyses/70-30sample/truthset/"
BED_FILE_NAME = "union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio.bed"


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

def intChrom(chrom):
    if chrom == 'X':
        return 23
    elif chrom == 'Y':
        return 24
    elif chrom == 'MT':
        return 25
    else:
        return int(chrom)

def calculateAllelicFreq(info,genotype,caller,tumorVariantType,alt):
    infoSplit = info.split(':')
    genotypeSplit = genotype.split(':')
    if caller == 'mutect':
        return float(genotypeSplit[infoSplit.index('FA')])
    elif caller == 'varscan':
        return float(genotypeSplit[infoSplit.index('FREQ')].split('%')[0])/100
    else:
        if caller == 'freebayes':
            ad = genotypeSplit[infoSplit.index('AO')].split(',')[0]  #NB - does not take into account 2nd allele if exists
            rd = genotypeSplit[infoSplit.index('RO')].split(',')[0]  # NB - does not take into account 2nd allele if exists
        elif caller == 'strelka' and tumorVariantType == variantType.SNP:
            ad, rd = genotypeSplit[infoSplit.index(alt.split(',')[int(genotype[0])] + 'U')].split(',')  # NB - does not take into account genotype
        elif caller == 'strelka' and tumorVariantType == variantType.indel:
            ad = genotypeSplit[infoSplit.index('TIR')].split(',')[0]  # NB - does not take into account 2nd allele if exists
            rd = genotypeSplit[infoSplit.index('TAR')].split(',')[0]  # NB - does not take into account 2nd allele if exists
            return float(ad) / (float(rd) + float(ad))
        elif caller == 'melted':
            #print genotypeSplit
            ad, rd = genotypeSplit[infoSplit.index('AD')].split(',')[:2]
        else:
            return -1

        if float(ad) == 0:
            return 0
        else:
            return float(ad) / (float(rd) + float(ad))



class variantType():
    sameAsRef = "Same as Ref"
    missingGenotype = "Missing Genotype"
    indel = "INDEL"
    SNP = "SNP"

class subVariantType():
    none = ""
    insert = "INSERT"
    delete = "DELETE"

class genotype:

    variantCountTumor = {}
    variantCountTumorPrivate = {}
    variantCountSubTypeTumor = {}
    variantCountSubTypeTumorPrivate = {}
    variantFreebayesSNP=0
    df = pd.DataFrame()
    variantInfo = []


    def __init__(self,chrom,pos,caller,ref,alt,info,vennSegment,inputGenotype):
        altsplit = (ref + ","+ alt).split(',')
        self.tumorVariantSubType = subVariantType.none

        if inputGenotype[:3] == "./.":
            self.tumorVariantType = variantType.missingGenotype
        elif inputGenotype[:3] == "0/0":
            self.tumorVariantType = variantType.sameAsRef
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

            ############### Pandas ##################
            if RUN_PANDAS == True:
                allelicFreq = calculateAllelicFreq(info,inputGenotype,caller,self.tumorVariantType,alt)
                if chrom[:3] == 'chr':
                    chrom = chrom[3:]
                posPercent = float(pos) / chromosomeLength[chrom]
                genotype.variantInfo.append((chrom, pos, intChrom(chrom)+posPercent,caller, ref, alleleTumor1, alleleTumor2, \
                                             vennSegment,self.tumorVariantType,self.tumorVariantSubType,allelicFreq))
            #########################################

        if genotype.variantCountSubTypeTumor.has_key(str(self.tumorVariantType) + str(self.tumorVariantSubType) + " " + caller):
            genotype.variantCountSubTypeTumor[str(self.tumorVariantType) + str(self.tumorVariantSubType) + " " + caller] += 1
        else:
            genotype.variantCountSubTypeTumor[str(self.tumorVariantType) + str(self.tumorVariantSubType) + " " + caller] = 1

        if genotype.variantCountTumor.has_key(str(self.tumorVariantType) + " " + caller):
            genotype.variantCountTumor[str(self.tumorVariantType) + " " + caller] += 1
        else:
            genotype.variantCountTumor[str(self.tumorVariantType) + " " + caller] = 1

    def markPrivate(self,caller):
        if genotype.variantCountSubTypeTumorPrivate.has_key(str(self.tumorVariantType) + str(self.tumorVariantSubType) + " " + caller):
            genotype.variantCountSubTypeTumorPrivate[str(self.tumorVariantType) + str(self.tumorVariantSubType) + " " + caller] += 1
        else:
            genotype.variantCountSubTypeTumorPrivate[str(self.tumorVariantType) + str(self.tumorVariantSubType) + " " + caller] = 1

        if genotype.variantCountTumorPrivate.has_key(str(self.tumorVariantType) + " " + caller):
            genotype.variantCountTumorPrivate[str(self.tumorVariantType) + " " + caller] += 1
        else:
            genotype.variantCountTumorPrivate[str(self.tumorVariantType) + " " + caller] = 1



class somaticVariant:
    variantCountSNPTotal = 0
    variantCountSNPNumberCallers = {}
    variantCountIndelNumberCallers = {}
    variantCountIndelTotal = 0
    bedItem = []

    def __init__(self, chrom, pos, id, ref, alt, filter, format,info,inputGenotypes,useBed,aBedReverse):

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
            if filter == "PASS" or filter == ".":
                tumorCallerCountSNP = 0
                tumorCallerCountIndel = 0
                variantGenotypes = {}
                vennSegment = ""

                formatSplit = format.split(';')
                for i in range(len(formatSplit)):
                    formatItem = formatSplit[i].split('=')
                    if formatItem[0] == "set":
                        vennSegment = formatItem[1]

                for key in inputGenotypes.iterkeys():
                    variantGenotypes[key] = genotype(chrom, pos, key, ref, alt, info, vennSegment, inputGenotypes[key])

                for key, value in variantGenotypes.items():
                    if value.tumorVariantType == variantType.SNP:
                        tumorCallerCountSNP += 1
                    if value.tumorVariantType == variantType.indel:
                        tumorCallerCountIndel += 1

                if tumorCallerCountSNP > 0:
                    if somaticVariant.variantCountSNPNumberCallers.has_key(tumorCallerCountSNP):
                        somaticVariant.variantCountSNPNumberCallers[tumorCallerCountSNP] += 1
                    else:
                        somaticVariant.variantCountSNPNumberCallers[tumorCallerCountSNP] = 1

                if tumorCallerCountIndel > 0:
                    if somaticVariant.variantCountIndelNumberCallers.has_key(tumorCallerCountIndel):
                        somaticVariant.variantCountIndelNumberCallers[tumorCallerCountIndel] += 1
                    else:
                        somaticVariant.variantCountIndelNumberCallers[tumorCallerCountIndel] = 1

                for caller, variantGenotype in variantGenotypes.items():
                    if variantGenotype.tumorVariantType == variantType.SNP and tumorCallerCountSNP == 1:
                        variantGenotype.markPrivate(caller)
                    if variantGenotype.tumorVariantType == variantType.indel and tumorCallerCountIndel == 1:
                        variantGenotype.markPrivate(caller)

                if tumorCallerCountSNP > 0:
                    somaticVariant.variantCountSNPTotal += 1

                if tumorCallerCountIndel > 0:
                    somaticVariant.variantCountIndelTotal += 1



def loadVaraintsFromVCF(aPath, aVCFFile,sampleNames,useBed=False,aBedReverse=[]):
    print "reading vcf file. . .\n"

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
                variants.append(somaticVariant(a[0], a[1], a[2], a[3], a[4], a[6], a[7], a[8],myGenotypes,useBed,aBedReverse))
                i += 1
                if i% 100000 == 0:
                    print "reading VCF File line:",i
                #if i > 1000000000:
                #    break
        #Need to reset bed item
        somaticVariant.bedItem = []
    print "data frame loaded\n"
    if RUN_PANDAS == True:
        df = pd.DataFrame(genotype.variantInfo)
        df.columns = (['chrom', 'pos', 'chromFrac', 'caller','ref', 'alleleTumor1', 'alleleTumor2','vennSegment','variantType', \
                   'variantSubType','allelicFreq'])
        return df
    else:
        return 0

def loadBEDFile(aPath, aBEDFile):
    print "reading BED file. . .\n"
    myBed = []
    with open(aPath + aBEDFile, 'r') as f:
        for line in f:
            line = line.strip('\n')
            splitLine = line.split('\t')
            if splitLine[0] != 'chrom':
                myBed.append(splitLine)
    return myBed

def printStatistics():
    print "\nTUMOR VARIANTS"
    for key, value in sorted(genotype.variantCountSubTypeTumor.items()):
        print "%s: %s" % (key, value)
    indelTruthSet = somaticVariant.variantCountIndelTotal - somaticVariant.variantCountIndelNumberCallers[1]
    snpTruthSet = somaticVariant.variantCountSNPTotal - somaticVariant.variantCountSNPNumberCallers[1]
    print "\nIndel Total: ", somaticVariant.variantCountIndelTotal
    print "SNP Total: ", somaticVariant.variantCountSNPTotal
    print "\nIndel 'Truth Set': ", indelTruthSet
    print "SNP 'Truth Set': ", snpTruthSet
    print "\nSENSITIVITY"
    for myVariantType, myTumorCount in sorted(genotype.variantCountTumor.items()):
        if myVariantType[:3] == 'SNP':
            print myVariantType, ":", round(
                float(myTumorCount - genotype.variantCountTumorPrivate[myVariantType]) / snpTruthSet,
                ROUNDING_PRECISION)
        elif myVariantType[:5] == 'INDEL':
            print myVariantType, ":", round(
                float(myTumorCount - genotype.variantCountTumorPrivate[myVariantType]) / indelTruthSet,
                ROUNDING_PRECISION)
    print "\nPRECISION"
    for myVariantType, myTumorCount in sorted(genotype.variantCountSubTypeTumorPrivate.items()):
        print myVariantType, ":", round(
            1 - float(myTumorCount) / float(genotype.variantCountSubTypeTumor[myVariantType]), ROUNDING_PRECISION)
    print "\nNumber of Callers SNP: ", somaticVariant.variantCountSNPNumberCallers
    print "Number of Callers Indel: ", somaticVariant.variantCountIndelNumberCallers

if __name__ == "__main__":
    RUN_PANDAS = False
    bed = loadBEDFile(BED_PATH,BED_FILE_NAME)
    bed.reverse()
    if loadVaraintsFromVCF(Path,VCFFile,SAMPLE_NAMES,True,bed) != -1:
        printStatistics()




