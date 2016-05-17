#!/usr/local/bin/python
import math
import numpy as np
from array import array

ROUNDING_PRECISION = 4
Path = "/Users/peterpriestley/hmf/70-30sample/"
VCFFile = "combined.vcf"

class genotype:

    variantCountNormal = {}
    variantCountTumor = {}
    variantCountTumorPrivate = {}

    def __init__(self,caller,ref,alt,genotypeTumor,genotypeNormal):
        self.private = False
        self.caller = caller
        self.genotypeTumor = genotypeTumor
        self.genotypeNormal = genotypeNormal
        self.altsplit = (ref + ","+ alt).split(',')

        if self.genotypeTumor == "./.":
            self.tumorVariantType = "Missing Genotype"
        elif self.genotypeTumor == "0/0":
            self.tumorVariantType = "Same As Ref"
        else:
            self.alleleTumor1 = self.altsplit[int(self.genotypeTumor[0])]
            self.alleleTumor2 = self.altsplit[int(self.genotypeTumor[2])]
            if len(self.alleleTumor1) == len(self.alleleTumor2) and len(self.alleleTumor1) == len(ref):
                self.tumorVariantType = "SNP"
            else:
                self.tumorVariantType = "Indel"

            if genotype.variantCountTumor.has_key(caller + self.tumorVariantType):
                genotype.variantCountTumor[caller + self.tumorVariantType] += 1
            else:
                genotype.variantCountTumor[caller + self.tumorVariantType] = 1

        if self.genotypeNormal == "./.":
            self.normalVariantType = "Missing Genotype"
        elif self.genotypeNormal == "0/0" or self.genotypeNormal[1] != "/":  # Mutect Normal Case
            self.normalVariantType = "Same As Ref"
        else:
            self.alleleNormal1 = self.altsplit[int(self.genotypeNormal[0])]
            self.alleleNormal2 = self.altsplit[int(self.genotypeNormal[2])]
            if len(self.alleleNormal1) == len(self.alleleNormal2) and len(self.alleleNormal1) == len(ref):
                self.normalVariantType = "SNP"
            else:
                self.normalVariantType = "Indel"
            if genotype.variantCountNormal.has_key(caller + self.normalVariantType):
                genotype.variantCountNormal[caller + self.normalVariantType] += 1
            else:
                genotype.variantCountNormal[caller + self.normalVariantType] = 1

    def markPrivate(self):
        self.private = True
        if genotype.variantCountTumorPrivate.has_key(self.caller + self.tumorVariantType):
            genotype.variantCountTumorPrivate[self.caller + self.tumorVariantType] += 1
        else:
            genotype.variantCountTumorPrivate[self.caller + self.tumorVariantType] = 1




# TODO
# Define an SNP vs an INDEL
# KEEP SEPARATE COUNTS OF BOTH!
#

class somaticVariant:
    variantCountSNPTotal = 0
    variantCountSNPNumberCallers = {}
    variantCountIndelNumberCallers = {}
    variantCountIndelTotal = 0


    def __init__(self, fullstring, chrom,pos, id, ref,alt,format,freebayesNormalGenotype,mutectNormalGenotype, \
                 freebayesTumorGenotype,mutectTumorGenotype,strelkaNormalGenotype,varscanNormalGenotype, \
                 strelkaTumorGenotype,varscanTumorGenotype):
        self.tumorCallerCount = 0
        self.tumorCallerCountIndel = 0
        self.fullstring = fullstring
        self.chrom = chrom
        self.pos = pos
        self.id = id
        self.format = format
        self.variantGenotypes = {}

        self.variantGenotypes['mutect'] = genotype("mutect",ref,alt,mutectTumorGenotype,mutectNormalGenotype)
        self.variantGenotypes['freebayes'] = genotype("freebayes", ref, alt, freebayesTumorGenotype, freebayesNormalGenotype)
        self.variantGenotypes['varscan'] = genotype("varscan", ref, alt, varscanTumorGenotype, varscanNormalGenotype)
        self.variantGenotypes['strelka'] = genotype("strelka", ref, alt, strelkaTumorGenotype, strelkaNormalGenotype)

        for key,value in self.variantGenotypes.items():
            if value.tumorVariantType == "SNP":
                self.tumorCallerCount += 1
            if value.tumorVariantType == "Indel":
                self.tumorCallerCountIndel += 1

        if self.tumorCallerCount > 0:
            if somaticVariant.variantCountSNPNumberCallers.has_key(self.tumorCallerCount):
                somaticVariant.variantCountSNPNumberCallers[self.tumorCallerCount] += 1
            else:
                somaticVariant.variantCountSNPNumberCallers[self.tumorCallerCount] = 1

        if self.tumorCallerCountIndel > 0:
            if somaticVariant.variantCountIndelNumberCallers.has_key(self.tumorCallerCountIndel):
                somaticVariant.variantCountIndelNumberCallers[self.tumorCallerCountIndel] += 1
            else:
                somaticVariant.variantCountIndelNumberCallers[self.tumorCallerCountIndel] = 1

        # self.mutectGenotype = genotype("mutect",ref,alt,mutectTumorGenotype,mutectNormalGenotype)
        # self.freebayesGenotype = genotype("freebayes", ref, alt, freebayesTumorGenotype, freebayesNormalGenotype)
        # self.varscanGenotype = genotype("varscan", ref, alt, varscanTumorGenotype, varscanNormalGenotype)
        # self.strelkaGenotype = genotype("strelka", ref, alt, strelkaTumorGenotype, strelkaNormalGenotype)
        #
        # if self.mutectGenotype.tumorVariantType = "SNP":
        #     self.tumorCallerCount += 1
        # if self.mutectGenotype.tumorVariantType = "Indel":
        #     self.tumorCallerCountIndel + = 1
        # if self.mutectGenotype.tumorVariantType = "SNP":
        #     self.tumorCallerCount += 1
        # if self.mutectGenotype.tumorVariantType = "Indel":
        #     self.tumorCallerCountIndel + = 1
        # if self.mutectGenotype.tumorVariantType = "SNP":
        #     self.tumorCallerCount += 1
        # if self.Genotype.tumorVariantType = "Indel":
        #     self.tumorCallerCountIndel + = 1
        # if self.strelkaGenotype.tumorVariantType = "SNP":
        #     self.tumorCallerCount += 1
        # if self.strelkaGenotype.tumorVariantType = "Indel":
        #     self.tumorCallerCountIndel + = 1

        for caller, variantGenotype in self.variantGenotypes.items():
            if variantGenotype.tumorVariantType == "SNP" and self.tumorCallerCount == 1:
                value.markPrivate()
            if variantGenotype.tumorVariantType == "Indel" and self.tumorCallerCountIndel == 1:
                value.markPrivate()

        # if self.tumorCallerCountIndel == 1:
        #     self.mutectGenotype.markPrivate("Indel")
        #     self.freebayesGenotype.markPrivate("Indel")
        #     self.strelkaGenotype.markPrivate("Indel")
        #     self.varscanGenotype.markPrivate("Indel")

        if self.tumorCallerCount > 0:
            somaticVariant.variantCountSNPTotal += 1

        if self.tumorCallerCountIndel > 0:
            somaticVariant.variantCountIndelTotal += 1


# 
# class somaticVariant:
#     variantCountTotal = 0
#     variantCountPrivate = 0
#     variantCountOneCaller = 0
#     variantCountTwoCallers = 0
#     variantCountThreeCallers = 0
#     variantCountFourCallers = 0
#     variantCountMutectNormal = 0
#     variantCountMutectTumor = 0
#     variantCountVarscanNormal = 0
#     variantCountVarscanTumor = 0
#     variantCountFreebayesNormal = 0
#     variantCountFreebayesTumor = 0
#     variantCountStrelkaNormal = 0
#     variantCountStrelkaTumor = 0
#     variantCountMutectTumorPrivate = 0
#     variantCountVarscanTumorPrivate = 0
#     variantCountFreebayesTumorPrivate = 0
#     variantCountStrelkaTumorPrivate = 0
# 
#     variantCountIndelTotal = 0
#     variantCountIndelPrivate = 0
#     variantCountIndelOneCaller = 0
#     variantCountIndelTwoCallers = 0
#     variantCountIndelThreeCallers = 0
#     variantCountIndelFourCallers = 0
#     variantCountIndelMutectNormal = 0
#     variantCountIndelMutectTumor = 0
#     variantCountIndelVarscanNormal = 0
#     variantCountIndelVarscanTumor = 0
#     variantCountIndelFreebayesNormal = 0
#     variantCountIndelFreebayesTumor = 0
#     variantCountIndelStrelkaNormal = 0
#     variantCountIndelStrelkaTumor = 0
#     variantCountIndelMutectTumorPrivate = 0
#     variantCountIndelVarscanTumorPrivate = 0
#     variantCountIndelFreebayesTumorPrivate = 0
#     variantCountIndelStrelkaTumorPrivate = 0
# 
#     def __init__(self, fullstring, mutationType, chrom,pos, id, ref,alt,format,freebayesNormalGenotype,mutectNormalGenotype, \
#                  freebayesTumorGenotype,mutectTumorGenotype,strelkaNormalGenotype,varscanNormalGenotype, \
#                  strelkaTumorGenotype,varscanTumorGenotype):
#         self.tumorCallerCount = 0
#         self.tumorCallerCountIndel = 0
#         self.fullstring = fullstring
#         self.mutationType = mutationType
#         self.chrom = chrom
#         self.pos = pos
#         self.id = id
#         self.ref = ref
#         self.alt = alt
#         self.alleles = alt.count(',') + 1
#         self.format = format
#         self.mutectNormalGenotype = mutectNormalGenotype
#         self.mutectTumorGenotype = mutectTumorGenotype
#         self.freebayesNormalGenotype = freebayesNormalGenotype
#         self.freebayesTumorGenotype = freebayesTumorGenotype
#         self.strelkaNormalGenotype = strelkaNormalGenotype
#         self.strelkaTumorGenotype = strelkaTumorGenotype
#         self.varscanNormalGenotype = varscanNormalGenotype
#         self.varscanTumorGenotype = varscanTumorGenotype
#         self.mutectCall = False
#         self.varscanCall = False
#         self.freebayesCall = False
#         self.strelkaCall = False
#         self.mutectIndelCall = False
#         self.varscanIndelCall = False
#         self.freebayesIndelCall = False
#         self.strelkaIndelCall = False
# 
#         if mutationType == "SNP":
#             if varscanNormalGenotype != './.' and varscanNormalGenotype[1] == '/':
#                 somaticVariant.variantCountVarscanNormal += 1
#             if mutectNormalGenotype != './.' and mutectNormalGenotype[1] == '/':
#                 somaticVariant.variantCountMutectNormal += 1
# 
#             if freebayesNormalGenotype != './.' and freebayesNormalGenotype[1] == '/':
#                 somaticVariant.variantCountFreebayesNormal += 1
#             if strelkaNormalGenotype != './.' and strelkaNormalGenotype[1] == '/':
#                 somaticVariant.variantCountStrelkaNormal += 1
# 
#             if varscanTumorGenotype != './.' and varscanTumorGenotype[1] == '/':
#                 somaticVariant.variantCountVarscanTumor += 1
#                 self.tumorCallerCount += 1
#                 self.varscanCall = True
#             if mutectTumorGenotype != './.' and mutectTumorGenotype[1] == '/':
#                 somaticVariant.variantCountMutectTumor += 1
#                 self.tumorCallerCount += 1
#                 self.mutectCall = True
#             if freebayesTumorGenotype != './.' and freebayesTumorGenotype[1] == '/':
#                 somaticVariant.variantCountFreebayesTumor += 1
#                 self.tumorCallerCount += 1
#                 self.freebayesCall = True
#             if strelkaTumorGenotype != './.' and strelkaTumorGenotype[1] == '/':
#                 somaticVariant.variantCountStrelkaTumor += 1
#                 self.tumorCallerCount += 1
#                 self.strelkaCall = True
# 
#             if self.tumorCallerCount == 1:
#                 somaticVariant.variantCountPrivate += 1
#                 if self.varscanCall:
#                     somaticVariant.variantCountVarscanTumorPrivate += 1
#                 if self.mutectCall:
#                     somaticVariant.variantCountMutectTumorPrivate += 1
#                 if self.strelkaCall:
#                     somaticVariant.variantCountStrelkaTumorPrivate += 1
#                 if self.freebayesCall:
#                     somaticVariant.variantCountFreebayesTumorPrivate += 1
#             elif self.tumorCallerCount == 2:
#                 somaticVariant.variantCountTwoCallers += 1
#             elif self.tumorCallerCount == 3:
#                 somaticVariant.variantCountThreeCallers += 1
#             elif self.tumorCallerCount == 4:
#                 somaticVariant.variantCountFourCallers += 1
# 
#             somaticVariant.variantCountTotal += 1
# 
#         elif mutationType == 'Indel':
# 
#             if varscanNormalGenotype != './.':
#                 somaticVariant.variantCountIndelVarscanNormal += 1
#             if mutectNormalGenotype != './.' and mutectNormalGenotype[1] == '/':  #MutectGivesDifferentInfoforNormalSample
#                 somaticVariant.variantCountIndelMutectNormal += 1
# 
#             if freebayesNormalGenotype != './.':
#                 somaticVariant.variantCountIndelFreebayesNormal += 1
#             if strelkaNormalGenotype != './.':
#                 somaticVariant.variantCountIndelStrelkaNormal += 1
# 
#             if varscanTumorGenotype != './.':
#                 somaticVariant.variantCountIndelVarscanTumor += 1
#                 self.tumorCallerCountIndel += 1
#                 self.varscanIndelCall = True
#             if mutectTumorGenotype != './.':
#                 somaticVariant.variantCountIndelMutectTumor += 1
#                 self.tumorCallerCountIndel += 1
#                 self.mutectIndelCall = True
#                 print fullstring
#             if freebayesTumorGenotype != './.':
#                 somaticVariant.variantCountIndelFreebayesTumor += 1
#                 self.tumorCallerCountIndel += 1
#                 self.freebayesIndelCall = True
#             if strelkaTumorGenotype != './.':
#                 somaticVariant.variantCountIndelStrelkaTumor += 1
#                 self.tumorCallerCountIndel += 1
#                 self.strelkaIndelCall = True
# 
#             if self.tumorCallerCountIndel == 1:
#                 somaticVariant.variantCountIndelPrivate += 1
#                 if self.varscanIndelCall:
#                     somaticVariant.variantCountIndelVarscanTumorPrivate += 1
#                 if self.mutectIndelCall:
#                     somaticVariant.variantCountIndelMutectTumorPrivate += 1
#                 if self.strelkaIndelCall:
#                     somaticVariant.variantCountIndelStrelkaTumorPrivate += 1
#                 if self.freebayesIndelCall:
#                     somaticVariant.variantCountIndelFreebayesTumorPrivate += 1
#             elif self.tumorCallerCountIndel == 2:
#                 somaticVariant.variantCountIndelTwoCallers += 1
#             elif self.tumorCallerCountIndel == 3:
#                 somaticVariant.variantCountIndelThreeCallers += 1
#             elif self.tumorCallerCountIndel == 4:
#                 somaticVariant.variantCountIndelFourCallers += 1
# 
#             somaticVariant.variantCountIndelTotal += 1
# 
# 
#     def displaySNPCount(self):
#         nonPrivateTumors = somaticVariant.variantCountTotal - somaticVariant.variantCountPrivate
#         print "Total Variant count %d" % somaticVariant.variantCountTotal
#         print "Truth Set %d\n" % nonPrivateTumors
#         print "BLOOD\nVarscan %d" % somaticVariant.variantCountVarscanNormal
#         print "Mutect %d" % somaticVariant.variantCountMutectNormal
#         print "Freebayes %d" % somaticVariant.variantCountFreebayesNormal
#         print "Strelka %d\n" % somaticVariant.variantCountStrelkaNormal
#         print "TUMORS\nVarscan %d" % somaticVariant.variantCountVarscanTumor
#         print "Mutect %d" % somaticVariant.variantCountMutectTumor
#         print "Freebayes %d" % somaticVariant.variantCountFreebayesTumor
#         print "Strelka %d\n" % somaticVariant.variantCountStrelkaTumor
#         print "CALLER COUNT\nPrivate %d" % somaticVariant.variantCountPrivate
#         print "Two Callers %d" % somaticVariant.variantCountTwoCallers
#         print "Three Callers %d" % somaticVariant.variantCountThreeCallers
#         print "All Callers %d\n" % somaticVariant.variantCountFourCallers
#         print "PRIVATE COUNT\nVarscan %d" % somaticVariant.variantCountVarscanTumorPrivate
#         print "Freebayes %d" % somaticVariant.variantCountFreebayesTumorPrivate
#         print "Mutect %d" % somaticVariant.variantCountMutectTumorPrivate
#         print "Strelka %d\n" % somaticVariant.variantCountStrelkaTumorPrivate
#         print "SENSITIVITY\nVarscan",round((somaticVariant.variantCountVarscanTumor- \
#                                             somaticVariant.variantCountVarscanTumorPrivate)/float(nonPrivateTumors),ROUNDING_PRECISION)
#         print "Freebayes", round((somaticVariant.variantCountFreebayesTumor - \
#                                   somaticVariant.variantCountFreebayesTumorPrivate) / float(nonPrivateTumors), ROUNDING_PRECISION)
#         print "Mutect", round((somaticVariant.variantCountMutectTumor - \
#                                somaticVariant.variantCountMutectTumorPrivate) / float(nonPrivateTumors), ROUNDING_PRECISION)
#         print "Strelka", round((somaticVariant.variantCountStrelkaTumor - \
#                                somaticVariant.variantCountStrelkaTumorPrivate) / float(nonPrivateTumors), ROUNDING_PRECISION)
#         print "\nPRECISION\nVarscan",1-round(somaticVariant.variantCountVarscanTumorPrivate \
#                                          /float(somaticVariant.variantCountVarscanTumor),ROUNDING_PRECISION)
#         print "Freebayes",1-round(somaticVariant.variantCountFreebayesTumorPrivate \
#                                          /float(somaticVariant.variantCountFreebayesTumor),ROUNDING_PRECISION)
#         print "Mutect",1-round(somaticVariant.variantCountMutectTumorPrivate \
#                                          /float(somaticVariant.variantCountMutectTumor),ROUNDING_PRECISION)
#         print "Strelka",1-round(somaticVariant.variantCountStrelkaTumorPrivate \
#                                          /float(somaticVariant.variantCountStrelkaTumor),ROUNDING_PRECISION)
# 
#     def displayIndelCount(self):
#         nonPrivateTumors = somaticVariant.variantCountIndelTotal - somaticVariant.variantCountIndelPrivate
#         print "Total Variant count %d" % somaticVariant.variantCountIndelTotal
#         print "Truth Set %d\n" % nonPrivateTumors
#         print "BLOOD\nVarscan %d" % somaticVariant.variantCountIndelVarscanNormal
#         print "Mutect %d" % somaticVariant.variantCountIndelMutectNormal
#         print "Freebayes %d" % somaticVariant.variantCountIndelFreebayesNormal
#         print "Strelka %d\n" % somaticVariant.variantCountIndelStrelkaNormal
#         print "TUMORS\nVarscan %d" % somaticVariant.variantCountIndelVarscanTumor
#         print "Mutect %d" % somaticVariant.variantCountIndelMutectTumor
#         print "Freebayes %d" % somaticVariant.variantCountIndelFreebayesTumor
#         print "Strelka %d\n" % somaticVariant.variantCountIndelStrelkaTumor
#         print "CALLER COUNT\nPrivate %d" % somaticVariant.variantCountIndelPrivate
#         print "Two Callers %d" % somaticVariant.variantCountIndelTwoCallers
#         print "Three Callers %d" % somaticVariant.variantCountIndelThreeCallers
#         print "All Callers %d\n" % somaticVariant.variantCountIndelFourCallers
#         print "PRIVATE COUNT\nVarscan %d" % somaticVariant.variantCountIndelVarscanTumorPrivate
#         print "Freebayes %d" % somaticVariant.variantCountIndelFreebayesTumorPrivate
#         print "Mutect %d" % somaticVariant.variantCountIndelMutectTumorPrivate
#         print "Strelka %d\n" % somaticVariant.variantCountIndelStrelkaTumorPrivate
#         print "SENSITIVITY\nVarscan", round((somaticVariant.variantCountIndelVarscanTumor - \
#                                              somaticVariant.variantCountIndelVarscanTumorPrivate) / float(nonPrivateTumors),
#                                             ROUNDING_PRECISION)
#         print "Freebayes", round((somaticVariant.variantCountIndelFreebayesTumor - \
#                                   somaticVariant.variantCountIndelFreebayesTumorPrivate) / float(nonPrivateTumors),
#                                  ROUNDING_PRECISION)
#         print "Mutect", round((somaticVariant.variantCountIndelMutectTumor - \
#                                somaticVariant.variantCountIndelMutectTumorPrivate) / float(nonPrivateTumors),
#                               ROUNDING_PRECISION)
#         print "Strelka", round((somaticVariant.variantCountIndelStrelkaTumor - \
#                                 somaticVariant.variantCountIndelStrelkaTumorPrivate) / float(nonPrivateTumors),
#                                ROUNDING_PRECISION)
#         print "\nPRECISION\nVarscan", 1 - round(somaticVariant.variantCountIndelVarscanTumorPrivate \
#                                                 / float(somaticVariant.variantCountIndelVarscanTumor), ROUNDING_PRECISION)
#         print "Freebayes", 1 - round(somaticVariant.variantCountIndelFreebayesTumorPrivate \
#                                      / float(somaticVariant.variantCountIndelFreebayesTumor), ROUNDING_PRECISION)
#         print "Mutect",0.0000 # Mutect does not have indels
#         print "Strelka", 1 - round(somaticVariant.variantCountIndelStrelkaTumorPrivate \
#                                    / float(somaticVariant.variantCountIndelStrelkaTumor), ROUNDING_PRECISION)
# 
#     def displayVariant(self):
#         print (self.chrom,self.pos,self.id,self.ref,self.alt,self.mutectTumorGenotype,self.freebayesTumorGenotype, \
#                self.strelkaTumorGenotype,self.varscanTumorGenotype)
# 


def loadVaraintsFromVCF(aPath, aVCFFile,aVariants):
    print "reading vcf file. . .\n"
    with open(aPath + aVCFFile, 'r') as f:
        i=0
        for line in f:
            a = [x for x in line.split('\t')]
            if a[0][:1] != '#':
                #print a
                #if len(a[3]) == 1 and (len(a[4]) == 1 or a[4].find(",") == 1):
                aVariants.append(somaticVariant(a,a[0], a[1], a[2], a[3], a[4],a[8],a[9][:3], \
                        a[10][:3],a[11][:3],a[12][:3],a[13][:3],a[14][:3],a[15][:3],a[16][:3]))
                #else:
                #    aVariants.append(somaticVariant(a, "Indel", a[0], a[1], a[2], a[3], a[4],a[8], \
                #            a[9][:3], a[10][:3], a[11][:3], a[12][:3],a[13][:3], a[14][:3], a[15][:3], a[16][:3]))
                i=i+1
                if i > 1000:
                    return 1
    return 1


variants = []
loadVaraintsFromVCF(Path,VCFFile,variants)

#print "\n*********SNP********\n"
#variants[0].displaySNPCount()
#print "\n********INDELS*******\n"
#variants[0].displayIndelCount()


print "Normal Variants: ",variants[0].variantGenotypes['mutect'].variantCountNormal
print "Tumor variants: ",variants[0].variantGenotypes['mutect'].variantCountTumor
print "Private Variants: ",variants[0].variantGenotypes['mutect'].variantCountTumorPrivate
print "Total Indels: ",variants[0].variantCountIndelTotal
print "Total SNP: ",variants[0].variantCountSNPTotal
print "Number of Callers SNP: ",variants[0].variantCountSNPNumberCallers
print "Number of Callers Indel: ", variants[0].variantCountIndelNumberCallers

    #
# j=0
# genotypes = {}
# alleles = {}
# for i in range(0,variants[0].variantCountTotal+variants[0].variantCountIndelTotal):
#
#     myAllele = variants[i].alleles
#     if alleles.has_key(myAllele):
#         alleles[myAllele] +=1
#     else:
#         alleles[myAllele] = 1
#     if variants[i].freebayesIndelCall == True:
#         j += 1
#         print variants[i].freebayesTumorGenotype,variants[i].varscanTumorGenotype,variants[i].strelkaTumorGenotype,variants[i].ref, variants[i].alt
#         myGenotype = variants[i].freebayesTumorGenotype
#         if genotypes.has_key(myGenotype):
#             genotypes[myGenotype] += 1
#         else:
#             genotypes[myGenotype] = 1
#         if j > 1:
#             break
#
# #print genotypes
# print genotypes

#SNP in VCF STAT
# LEN(GT1) = LEN(GT2) = LEN(REF)
# Can be a SNP for 1 and a INDEL for another!
# RTGTOOLS counts this as an SNP
    #REF =CTTTTTTTTTTTTTTA  ; ALT = CTTTTTTTTTTTTTTT ;   GT: 1/1
# If one part of a genotype is a INDEL, then gets counted as an indel
# Missing Genotype = ./.
# Same as Reference = 0/0
    # Mutect does not use
    # Others use for Normal sample only
# Mysteries
    # Freebayes has a different missiing genotype count for R & T????
    # Why do varscan and freebayes have indels in Blood
    # Why does freebayes have SNP in Blood
    # Varscan Tumor count
        #INDEL = 187793+241939 =   429,732  (435,351)  ->
        #SNP = 1049522   (1,053,204)     -> 3682
    # how can the number of SNPs etc change when the file is combined????
    # What is the mutect normal info

#Variant
#INSERTION or DELETION or SNP or Reference or ...
#SNPs
#Freebayes {'0/0|0/1': 199125, './.|./.': 1064230, '0/1|1/1': 53, '0/0|1/1': 1}
#Strelka {'0/0|0/1': 1160332, './.|./.': 103077}
#varscan {'./.|./.': 210232, '0/0|0/1': 1053172, '0/0|1/1': 5}
#mutect {'./.|./.': 136749, './.|0/1': 1126660}

# FORMAT INFO
#GT:AD:BQ:DP:FA:SS     Mutect
#GT:AU:CU:DP:FDP:GU:SDP:SUBDP:TU    Strelka
#GT:AD:DP:DP4:FREQ:RD     Varscan
#GT:AO:DP:GQ:PL:QA:QR:RO   Freebayes
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of 'A' alleles used in tiers 1,2">
##FORMAT=<ID=BQ,Number=A,Type=Float,Description="Average base quality for reads supporting alleles">
##FORMAT=<ID=CU,Number=2,Type=Integer,Description="Number of 'C' alleles used in tiers 1,2">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=DP4,Number=1,Type=String,Description="Strand read counts: ref/fwd, ref/rev, var/fwd, var/rev">
##FORMAT=<ID=FA,Number=A,Type=Float,Description="Allele fraction of the alternate allele with regard to reference">
##FORMAT=<ID=FDP,Number=1,Type=Integer,Description="Number of basecalls filtered from original read depth for tier1">
##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GU,Number=2,Type=Integer,Description="Number of 'G' alleles used in tiers 1,2">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=QA,Number=A,Type=Integer,Description="Sum of quality of the alternate observations">
##FORMAT=<ID=QR,Number=1,Type=Integer,Description="Sum of quality of the reference observations">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">
##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Number of reads with deletions spanning this site at tier1">
##FORMAT=<ID=SS,Number=1,Type=Integer,Description="Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown">
##FORMAT=<ID=SUBDP,Number=1,Type=Integer,Description="Number of reads below tier1 mapping quality threshold aligned across this site">
##FORMAT=<ID=TU,Number=2,Type=Integer,Description="Number of 'T' alleles used in tiers 1,2">


#Freebayes
#SNP  199305    (199392)   87!!
#Indel 14413    (14326)