#!/usr/local/bin/python
import pandas as pd
import numpy as np
import difflib as dl
import chromosomeDefinition as cd

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

class variantType():
    sameAsRef = "Same as Ref"
    missingGenotype = "Missing Genotype"
    indel = "INDEL"
    SNP = "SNP"
    mixed = "MIXED"
    SV = "SV"

class subVariantType():
    none = ""
    insert = "INSERT"
    delete = "DELETE"
    indel = "INDEL"



def indelDiff(ref,variantAllele):
    #LOGIC - MATCH 1st char and then use strDiff in the reverse direction
    if ref[0]==variantAllele[0]: i=1
    else: i = 0
    reverseRef = ref[i:][::-1]
    reverseVariantAllele = variantAllele[i:][::-1]
    strdiff = dl.SequenceMatcher(None, reverseRef, reverseVariantAllele)
    myIndelString = ""
    for item in strdiff.get_opcodes():
        if item[0] == 'delete':
            myIndelString = myIndelString + "-" + reverseRef[item[1]:item[2]]
        elif item[0] == 'insert':
            myIndelString = myIndelString + "+" + reverseVariantAllele[item[3]:item[4]]
    myIndelString = myIndelString[::-1]
    return myIndelString

def calculateReadDepth(format,genotype):
    formatSplit = format.split(':')
    genotypeSplit = genotype.split(':')
    if 'DP' in formatSplit:
        try:
            return int(genotypeSplit[formatSplit.index('DP')])
        except IndexError:
            return -1
        except ValueError:
            return -1
    else:
        return 1

def calculateNumCallers(infoSplit,infoHeaders,internalCount):
    if "CSP" in infoHeaders:
        return int(infoSplit[infoHeaders.index("CSP")].split('=')[1])
    elif "CC" in infoHeaders:
        return int(infoSplit[infoHeaders.index("CC")].split('=')[1])
    else:
        return internalCount

def calculateSomaticGenotype(infoSplit,infoHeaders,caller,aVariantType):
    #Ideally both Normal (ref, hom, het) and somatic (ref,hom,het).
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

def calculateSVLengthAndStart(pos,infoSplit,infoHeaders):
    # NOT SURE WHAT TO DO WITH BND => diff chromosome
    # ALSO CHECK INS.
    svLen = [0,0]
    svStart = [pos,pos]
    for i in range(2):
        if "SVLEN" in infoHeaders:
            svLen[i] = int(infoSplit[infoHeaders.index("SVLEN")].split('=')[1])
        elif "END" in infoHeaders:
            svLen[i] = int(infoSplit[infoHeaders.index("END")].split('=')[1]) - pos
        else:
            return svLen,svStart
        if "CIPOS" in infoHeaders:
            svLen[i] = svLen[i] - int(infoSplit[infoHeaders.index("CIPOS")].split('=')[1].split(",")[i])
            svStart[i] = svStart[i] + int(infoSplit[infoHeaders.index("CIPOS")].split('=')[1].split(",")[1-i])
        if "CIEND" in infoHeaders:
            svLen[i] = svLen[i] + int(infoSplit[infoHeaders.index("CIEND")].split('=')[1].split(",")[1-i])
    return svLen,svStart

def calculateQualityScore(infoSplit,infoHeaders,caller,qual,aVariantType,format,genotype):
    if caller == 'strelka' and aVariantType == variantType.indel:
        try:
            return infoSplit[infoHeaders.index("QSI_NT")].split('=')[1]
        except ValueError:
            return -1
    elif caller == 'strelka' and aVariantType == variantType.SNP:
        return infoSplit[infoHeaders.index("QSS_NT")].split('=')[1]
    elif caller == 'varscan':
        if "VS_SSC" in infoHeaders:
            return infoSplit[infoHeaders.index("VS_SSC")].split('=')[1]
        else:
            return infoSplit[infoHeaders.index("SSC")].split('=')[1]
    elif caller == 'freebayes':
        return qual
    elif caller == 'Set1GIAB12878':
        try:
            return infoSplit[infoHeaders.index("MQ")].split('=')[1]
        except ValueError:
            return -1
    else:  # Mutect has no quality score
        return -1

def calculateConsensusVariantType(tumorCallerCountSNP,tumorCallerCountSV,tumorCallerCountIndel,tumorCallerCountSubTypeDelete,tumorCallerCountSubTypeInsert,tumorCallerCountSubTypeIndel):
    if tumorCallerCountIndel > 0 and tumorCallerCountSNP > 0:
        return variantType.mixed, ""
    elif tumorCallerCountIndel > 0:
        if tumorCallerCountSubTypeDelete > 0 and tumorCallerCountSubTypeIndel == 0 and tumorCallerCountSubTypeInsert == 0:
            return variantType.indel, subVariantType.delete
        elif tumorCallerCountSubTypeInsert > 0 and tumorCallerCountSubTypeIndel == 0 and tumorCallerCountSubTypeDelete == 0:
            return variantType.indel, subVariantType.insert
        else:
            return variantType.indel,subVariantType.indel
    elif tumorCallerCountSNP > 0:
        return variantType.SNP,""
    elif tumorCallerCountSV > 0:
        return variantType.SV, ""
    else:
        return variantType.missingGenotype,""

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
        else:
            try:
                rd, ad = genotypeSplit[formatSplit.index('AD')].split(',')[:2]
            except ValueError:
                return -1
        if float(ad) == 0:
            return 0
        else:
            return float(ad) / (float(rd) + float(ad))   #is this correct, or should it be /DP?

class genotype:

    def __init__(self,caller,pos,ref,alt,qual,infoSplit,infoHeaders,format,inputGenotype):
        altsplit = (ref + ","+ alt).split(',')

        #self.tumorVariantSubType = subVariantType.none
        self.indelDiff = ""
        self.svStart = ["",""]
        self.svLen = ["", ""]

        if inputGenotype[:3] == "./." or inputGenotype[1]<>"/":
            if alt[0]== "<":   #SV case
                self.tumorVariantType = variantType.SV
                self.allele = alt[1:-1]
                self.svLen,self.svStart = calculateSVLengthAndStart(int(pos),infoSplit,infoHeaders)
            else:
                self.tumorVariantType = variantType.missingGenotype
                self.allele = ""
        elif inputGenotype[:3] == "0/0" or alt == ".":   #STRELKA unfiltered
            self.tumorVariantType = variantType.sameAsRef
            self.allele = ref
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

            self.allele = alleleTumor2
            self.indelDiff = indelDiff(ref, self.allele)

        self.allelicFreq = calculateAllelicFreq(format, inputGenotype, caller, self.tumorVariantType, ref,self.allele)
        self.readDepth = calculateReadDepth(format,inputGenotype)
        self.qualityScore = float(calculateQualityScore(infoSplit,infoHeaders,caller,qual,self.tumorVariantType,format,inputGenotype))
        self.somaticGenotype = calculateSomaticGenotype(infoSplit,infoHeaders,caller,self.tumorVariantType)
        if self.somaticGenotype == 'unknown':
            self.somaticGenotype = inputGenotype[:3]



class somaticVariant:

    variantInfo = []
    bedItem = []

    def __init__(self, chrom, pos, id, ref, alt, qual, filter, info,format,inputGenotypes,useFilter,useBed,aBedReverse,loadRegionsOutsideBed):

        #Find the 1st Bed region with maxPos > variantPos
        if aBedReverse:
            if not somaticVariant.bedItem:
                somaticVariant.bedItem = aBedReverse.pop()
            while cd.intChrom(chrom) > cd.intChrom(somaticVariant.bedItem[0]) or (cd.intChrom(chrom) == cd.intChrom(somaticVariant.bedItem[0]) and int(pos) > int(somaticVariant.bedItem[2])) and aBedReverse:
                somaticVariant.bedItem = aBedReverse.pop()
        else:
            somaticVariant.bedItem = []

        #Label the BED region
        bedRegion = ""
        if (somaticVariant.bedItem and int(somaticVariant.bedItem[1])<int(pos) and int(somaticVariant.bedItem[2])>=int(pos) and somaticVariant.bedItem[0]==chrom):
            try:
                bedRegion = somaticVariant.bedItem[3]
            except IndexError:
                bedRegion = "Default"

        #Process if in Bed region or not using BED or if loading whole file
        if bedRegion <> "" or not useBed or loadRegionsOutsideBed:
            if filter == "PASS" or filter == "." or useFilter == False:

                tumorCallerCountSNP = 0
                tumorCallerCountSV = 0
                tumorCallerCountIndel = 0
                tumorCallerCountSubTypeIndel = 0
                tumorCallerCountSubTypeDelete = 0
                tumorCallerCountSubTypeInsert = 0
                variantGenotypes = {}

                #Split info fields
                infoSplit = info.split(';')
                infoHeaders = [x.split('=')[0] for x in infoSplit]

                #CALLER SPECIFIC FIELDS
                for key in inputGenotypes.iterkeys():
                    variantGenotypes[key] = genotype(key, pos, ref, alt, qual,infoSplit,infoHeaders,format,inputGenotypes[key])

                #CALLER COUNTS
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
                    if value.tumorVariantType == variantType.SV:
                        tumorCallerCountSV += 1

                # META DATA
                if "set" in infoHeaders:
                    vennSegment = infoSplit[infoHeaders.index("set")].split('=')[1]
                else:
                    vennSegment = "test"   # to do - calculate somatic, LOH, or germline
                numCallers = calculateNumCallers(infoSplit,infoHeaders,tumorCallerCountSNP + tumorCallerCountIndel)
                if "filter" in vennSegment:
                    numCallers=numCallers -1
                myVariantType,mySubVariantType = calculateConsensusVariantType(tumorCallerCountSNP,tumorCallerCountSV,tumorCallerCountIndel,\
                                                tumorCallerCountSubTypeDelete,tumorCallerCountSubTypeInsert,tumorCallerCountSubTypeIndel)
                inDBSNP = any(['rs' in x for x in id.split(';')])
                inCOSMIC = any(['COSM' in x for x in id.split(';')])

                #ANNOTATIONS
                annWorstEffect = ""
                annAllEffects = ""
                annWorstImpact = ""
                annGene = ""
                if 'ANN' in infoHeaders:
                    annSplit = infoSplit[infoHeaders.index("ANN")].split('=')[1].split(',')
                    annWorstImpact = annSplit[0].split('|')[2]
                    annWorstEffect = annSplit[0].split('|')[1]
                    annAllEffects = '|'.join([annAllEffects + x.split('|')[1] for x in annSplit])
                    annGene = annSplit[0].split('|')[3]

                #CONSENSUS RULE
                if myVariantType == variantType.SNP:
                    consensus = int(numCallers) >= 3 or (int(numCallers) == 2 and bedRegion <> "" and (not inDBSNP or inCOSMIC))
                elif myVariantType == variantType.indel:
                    consensus = (int(numCallers) >= 2)
                else:
                    consensus = False

                ############### Pandas Prep ####################

                if(cd.intChrom(chrom)) < 25:
                    # APPEND NORMAL FIELDS
                    somaticVariant.variantInfo.append(
                        [chrom, pos, chrom + ':' + pos, cd.intChrom(chrom) + float(pos) / cd.chromosomeLength[chrom], id, ref, vennSegment, numCallers,
                        myVariantType, mySubVariantType,filter,bedRegion,inDBSNP,inCOSMIC,annGene,annWorstImpact,annWorstEffect,annAllEffects,consensus])

                    #APPEND CALLER SPECIFIC FIELDS
                    for caller, variantGenotype in variantGenotypes.items():
                        if variantGenotype.tumorVariantType == variantType.indel or variantGenotype.tumorVariantType == variantType.SNP or True:
                            callerSpecificFields = [variantGenotype.allele, variantGenotype.allelicFreq, variantGenotype.readDepth,
                                                variantGenotype.qualityScore,variantGenotype.somaticGenotype,
                                                variantGenotype.indelDiff,variantGenotype.svLen[1],variantGenotype.svLen[0],variantGenotype.svStart[1],variantGenotype.svStart[0]]
                        else:
                            callerSpecificFields = ['', '', '','','','']
                        somaticVariant.variantInfo[-1] = somaticVariant.variantInfo[-1] + callerSpecificFields
                #########################################


def loadVariantsFromVCF(aPath, aVCFFile,sampleNames,aPatientName,useFilter,useBed=False,aBed=[],loadRegionsOutsideBed=False):

    variants = []
    i = 0

    if useBed == True:
        aBed.reverse()

    print "reading vcf file:", aVCFFile
    with open(aPath + aVCFFile, 'r') as f:
        for line in f:

            line = line.strip('\n')
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
                        print 'Error - missing sample input: ',sampleLabel
                        return -1

            if a[0][:1] != '#':
                variant_calls = a[9:]
                myGenotypes = {}
                for caller,index in header_index.iteritems():
                    myGenotypes[caller] = variant_calls[index]
                variants.append(somaticVariant(a[0].lstrip("chr"), a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8],myGenotypes, useFilter, useBed,aBed,loadRegionsOutsideBed))
                i += 1
                if i% 100000 == 0:
                    print "reading VCF File line:",i
                    #break

    #Reset bed item
    somaticVariant.bedItem = []

    print "Number variants loaded:",len(somaticVariant.variantInfo)

    ###### PANDAS ##############
    df = pd.DataFrame(somaticVariant.variantInfo)
    if len(df)>0:
        myColumnList = ['chrom', 'pos', 'chromPos','chromFrac','id', 'ref', 'vennSegment','numCallers','variantType','variantSubType','filter',
                        'bedRegion','inDBSNP','inCOSMIC','annGene','annWorstImpact','annWorstEffect','annAllEffects','consensus']
        for caller in header_index.iterkeys():
            myColumnList = myColumnList + [caller + 'allele',caller+ 'AF',caller+'DP',caller+'QS',caller+'SGT',caller+'indelDiff',caller+'SVLenMin',caller+'SVLenMax',caller+'SVStartMin',caller+'SVStartMax']
        df.columns = (myColumnList)
        df['patientName'] = aPatientName
    ###### END PANDAS ###########

    # Need to empty genotype.variantInfo in case we need to load multiple files
    del somaticVariant.variantInfo[:]
    return df

def loadBEDFile(aPath, aBEDFile):
    myBed = []
    with open(aPath + aBEDFile, 'r') as f:
        for line in f:
            line = line.strip('\n')
            line = line.strip('\r')
            splitLine = line.split('\t')
            if splitLine[0] != 'chrom':
                myBed.append(splitLine)
    return myBed



