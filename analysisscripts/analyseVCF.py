#!/usr/local/bin/python

ROUNDING_PRECISION = 4
Path = "/Users/peterpriestley/hmf/70-30sample/"
VCFFile = "combined.vcf"

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

    variantCountSubTypeNormal = {}
    variantCountTumor = {}
    variantCountTumorPrivate = {}
    variantCountSubTypeTumor = {}
    variantCountSubTypeTumorPrivate = {}

    def __init__(self,caller,ref,alt,inputGenotype):
        altsplit = (ref + ","+ alt).split(',')
        self.tumorVariantSubType = subVariantType.none
        self.normalVariantSubType = subVariantType.none

        #TUMOR SAMPLE
        if inputGenotype[1] == "./.":
            self.tumorVariantType = variantType.missingGenotype
        elif inputGenotype[1] == "0/0":
            self.tumorVariantType = variantType.sameAsRef
        else:
            alleleTumor1 = altsplit[int(inputGenotype[1][0])]
            alleleTumor2 = altsplit[int(inputGenotype[1][2])]
            if len(alleleTumor1) == len(alleleTumor2) and len(alleleTumor1) == len(ref):
                self.tumorVariantType = variantType.SNP
            else:
                self.tumorVariantType = variantType.indel
                if len(alleleTumor1) <= len(ref) and len(alleleTumor2) <= len(ref):
                    self.tumorVariantSubType = subVariantType.delete
                elif len(alleleTumor1) >= len(ref) and len(alleleTumor2) >= len(ref):
                    self.tumorVariantSubType = subVariantType.insert

        if genotype.variantCountSubTypeTumor.has_key(str(self.tumorVariantType) + str(self.tumorVariantSubType) + " " + caller):
            genotype.variantCountSubTypeTumor[str(self.tumorVariantType) + str(self.tumorVariantSubType) + " " + caller] += 1
        else:
            genotype.variantCountSubTypeTumor[str(self.tumorVariantType) + str(self.tumorVariantSubType) + " " + caller] = 1

        if genotype.variantCountTumor.has_key(str(self.tumorVariantType) + " " + caller):
            genotype.variantCountTumor[str(self.tumorVariantType) + " " + caller] += 1
        else:
            genotype.variantCountTumor[str(self.tumorVariantType) + " " + caller] = 1


        # NORMAL SAMPLE
        if inputGenotype[0] == "./.":
            self.normalVariantType = variantType.missingGenotype
        elif inputGenotype[0] == "0/0" or inputGenotype[0][1] != "/":  # Mutect Normal Case
            self.normalVariantType = variantType.sameAsRef
        else:
            alleleNormal1 = altsplit[int(inputGenotype[0][0])]
            alleleNormal2 = altsplit[int(inputGenotype[0][2])]
            if len(alleleNormal1) == len(alleleNormal2) and len(alleleNormal1) == len(ref):
                self.normalVariantType = variantType.SNP
            else:
                self.normalVariantType = variantType.indel
                if len(alleleNormal1) <= len(ref) and len(alleleNormal2) <= len(ref):
                    self.normalVariantSubType = subVariantType.delete
                elif len(alleleNormal1) >= len(ref) and len(alleleNormal2) >= len(ref):
                    self.normalVariantSubType = subVariantType.insert

        if genotype.variantCountSubTypeNormal.has_key(str(self.normalVariantType) + str(self.normalVariantSubType) + " " + caller):
            genotype.variantCountSubTypeNormal[str(self.normalVariantType) + str(self.normalVariantSubType) + " " + caller] += 1
        else:
            genotype.variantCountSubTypeNormal[str(self.normalVariantType) + str(self.normalVariantSubType) + " " + caller] = 1


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

    def __init__(self, chrom, pos, id, ref, alt, filter, inputGenotypes):

        if filter == "PASS" or filter == ".":
            tumorCallerCountSNP = 0
            tumorCallerCountIndel = 0
            #self.chrom = chrom
            #self.pos = pos
            #self.id = id
            variantGenotypes = {}

            for key in inputGenotypes.iterkeys():
                variantGenotypes[key] = genotype(key, ref, alt, inputGenotypes[key])

            #####DEBUG#########
            #if self.chrom == "1" and self.pos == "170993444":
            #    print filter, chrom, pos, id, ref, alt, inputGenotypes
            #####DEBUG END#####

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


def loadVaraintsFromVCF(aPath, aVCFFile,aVariants):
    print "reading vcf file. . .\n"
    with open(aPath + aVCFFile, 'r') as f:
        for line in f:
            myGenotypes = {}
            a = [x for x in line.split('\t')]
            if a[0][:1] != '#':
                myGenotypes['freebayes']=[a[9][:3],a[11][:3]]
                myGenotypes['mutect'] = [a[10][:3], a[12][:3]]
                myGenotypes['strelka'] = [a[13][:3], a[15][:3]]
                myGenotypes['varscan'] = [a[14][:3], a[16][:3]]
                aVariants.append(somaticVariant(a[0], a[1], a[2], a[3], a[4], a[6], myGenotypes))

    return 1


variants = []
loadVaraintsFromVCF(Path,VCFFile,variants)


print "NORMAL VARIANTS"
for key, value in sorted(genotype.variantCountSubTypeNormal.items()):
    print "%s: %s" % (key, value)
print "\nTUMOR VARIANTS"
for key, value in sorted(genotype.variantCountSubTypeTumor.items()):
    print "%s: %s" % (key, value)
indelTruthSet = somaticVariant.variantCountIndelTotal - somaticVariant.variantCountIndelNumberCallers[1]
snpTruthSet = somaticVariant.variantCountSNPTotal - somaticVariant.variantCountSNPNumberCallers[1]
print "\nIndel Total: ",somaticVariant.variantCountIndelTotal
print "SNP Total: ",somaticVariant.variantCountSNPTotal
print "\nIndel 'Truth Set': ",indelTruthSet
print "SNP 'Truth Set': ",snpTruthSet
print "\nSENSITIVITY"
for myVariantType,myTumorCount in sorted(genotype.variantCountTumor.items()):
            if myVariantType[:3] == 'SNP':
                print myVariantType,":",round(float(myTumorCount-genotype.variantCountTumorPrivate[myVariantType])/snpTruthSet,ROUNDING_PRECISION)
            elif myVariantType[:5]  == 'INDEL':
                print myVariantType,":",round(float(myTumorCount-genotype.variantCountTumorPrivate[myVariantType])/indelTruthSet,ROUNDING_PRECISION)
print "\nPRECISION"
for myVariantType,myTumorCount in sorted(genotype.variantCountSubTypeTumorPrivate.items()):
            print myVariantType,":",round(1-float(myTumorCount)/float(genotype.variantCountSubTypeTumor[myVariantType]),ROUNDING_PRECISION)
print "\nNumber of Callers SNP: ",somaticVariant.variantCountSNPNumberCallers
print "Number of Callers Indel: ", somaticVariant.variantCountIndelNumberCallers
