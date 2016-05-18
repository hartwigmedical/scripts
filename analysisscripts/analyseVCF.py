#!/usr/local/bin/python
#import math
import numpy as np
from array import array

ROUNDING_PRECISION = 4
Path = "/Users/peterpriestley/hmf/70-30sample/"
VCFFile = "combined.vcf"

class genotype:

    variantCountSubTypeNormal = {}
    variantCountTumor = {}
    variantCountTumorPrivate = {}
    variantCountSubTypeTumor = {}
    variantCountSubTypeTumorPrivate = {}

    def __init__(self,caller,ref,alt,inputGenotype):
        self.private = False
        self.caller = caller
        self.genotypeTumor = inputGenotype[1]
        self.genotypeNormal = inputGenotype[0]
        self.altsplit = (ref + ","+ alt).split(',')
        self.tumorVariantSubType = ""
        self.normalVariantSubType = ""

        #TUMOR SAMPLE
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
                if len(self.alleleTumor1) <= len(ref) and len(self.alleleTumor2) <= len(ref):
                    self.tumorVariantSubType = "Delete"
                elif len(self.alleleTumor1) >= len(ref) and len(self.alleleTumor2) >= len(ref):
                    self.tumorVariantSubType = "Insert"

        if genotype.variantCountSubTypeTumor.has_key(self.tumorVariantType + self.tumorVariantSubType + " " + self.caller):
            genotype.variantCountSubTypeTumor[self.tumorVariantType + self.tumorVariantSubType + " " + self.caller] += 1
        else:
            genotype.variantCountSubTypeTumor[self.tumorVariantType + self.tumorVariantSubType + " " + self.caller] = 1

        if genotype.variantCountTumor.has_key(self.tumorVariantType + " " + self.caller):
            genotype.variantCountTumor[self.tumorVariantType + " " + self.caller] += 1
        else:
            genotype.variantCountTumor[self.tumorVariantType + " " + self.caller] = 1


        # NORMAL SAMPLE
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
                if len(self.alleleNormal1) <= len(ref) and len(self.alleleNormal2) <= len(ref):
                    self.normalVariantSubType = "Delete"
                elif len(self.alleleNormal1) >= len(ref) and len(self.alleleNormal2) >= len(ref):
                    self.normalVariantSubType = "Insert"

        if genotype.variantCountSubTypeNormal.has_key(self.normalVariantType + self.normalVariantSubType + " " + self.caller):
            genotype.variantCountSubTypeNormal[self.normalVariantType + self.normalVariantSubType + " " + self.caller] += 1
        else:
            genotype.variantCountSubTypeNormal[self.normalVariantType + self.normalVariantSubType + " " + self.caller] = 1


    def markPrivate(self):
        self.private = True
        if genotype.variantCountSubTypeTumorPrivate.has_key(self.tumorVariantType + self.tumorVariantSubType + " " + self.caller):
            genotype.variantCountSubTypeTumorPrivate[self.tumorVariantType + self.tumorVariantSubType + " " + self.caller] += 1
        else:
            genotype.variantCountSubTypeTumorPrivate[self.tumorVariantType + self.tumorVariantSubType + " " + self.caller] = 1

        if genotype.variantCountTumorPrivate.has_key(self.tumorVariantType + " " + self.caller):
            genotype.variantCountTumorPrivate[self.tumorVariantType + " " + self.caller] += 1
        else:
            genotype.variantCountTumorPrivate[self.tumorVariantType + " " + self.caller] = 1



class somaticVariant:
    variantCountSNPTotal = 0
    variantCountSNPNumberCallers = {}
    variantCountIndelNumberCallers = {}
    variantCountIndelTotal = 0

    def __init__(self, chrom, pos, id, ref, alt, filter, format, inputGenotypes):

        if filter == "PASS" or filter == ".":
            self.tumorCallerCountSNP = 0
            self.tumorCallerCountIndel = 0
            self.chrom = chrom
            self.pos = pos
            self.id = id
            self.format = format
            self.variantGenotypes = {}

            for key in inputGenotypes.iterkeys():
                self.variantGenotypes[key] = genotype(key, ref, alt, inputGenotypes[key])

            #####DEBUG#########
            #if self.chrom == "1" and self.pos == "170993444":
            #    print filter, chrom, pos, id, ref, alt, inputGenotypes
            #####DEBUG END#####

            for key, value in self.variantGenotypes.items():
                if value.tumorVariantType == "SNP":
                    self.tumorCallerCountSNP += 1
                if value.tumorVariantType == "Indel":
                    self.tumorCallerCountIndel += 1

            if self.tumorCallerCountSNP > 0:
                if somaticVariant.variantCountSNPNumberCallers.has_key(self.tumorCallerCountSNP):
                    somaticVariant.variantCountSNPNumberCallers[self.tumorCallerCountSNP] += 1
                else:
                    somaticVariant.variantCountSNPNumberCallers[self.tumorCallerCountSNP] = 1

            if self.tumorCallerCountIndel > 0:
                if somaticVariant.variantCountIndelNumberCallers.has_key(self.tumorCallerCountIndel):
                    somaticVariant.variantCountIndelNumberCallers[self.tumorCallerCountIndel] += 1
                else:
                    somaticVariant.variantCountIndelNumberCallers[self.tumorCallerCountIndel] = 1

            for caller, variantGenotype in self.variantGenotypes.items():
                if variantGenotype.tumorVariantType == "SNP" and self.tumorCallerCountSNP == 1:
                    variantGenotype.markPrivate()
                if variantGenotype.tumorVariantType == "Indel" and self.tumorCallerCountIndel == 1:
                    variantGenotype.markPrivate()

            if self.tumorCallerCountSNP > 0:
                somaticVariant.variantCountSNPTotal += 1

            if self.tumorCallerCountIndel > 0:
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
                aVariants.append(somaticVariant(a[0], a[1], a[2], a[3], a[4], a[6], a[8],myGenotypes))

    return 1


variants = []
loadVaraintsFromVCF(Path,VCFFile,variants)


print "NORMAL VARIANTS"
for key, value in sorted(variants[0].variantGenotypes['varscan'].variantCountSubTypeNormal.items()):
    print "%s: %s" % (key, value)
print "\nTUMOR VARIANTS"
for key, value in sorted(variants[0].variantGenotypes['varscan'].variantCountSubTypeTumor.items()):
    print "%s: %s" % (key, value)

indelTruthSet = variants[0].variantCountIndelTotal - variants[0].variantCountIndelNumberCallers[1]
snpTruthSet = variants[0].variantCountSNPTotal - variants[0].variantCountSNPNumberCallers[1]
print "\nIndel Total: ",variants[0].variantCountIndelTotal
print "SNP Total: ",variants[0].variantCountSNPTotal
print "\nIndel 'Truth Set': ",indelTruthSet
print "SNP 'Truth Set': ",snpTruthSet
print "\nSENSITIVITY"
for key,value in sorted(variants[0].variantGenotypes['varscan'].variantCountTumor.items()):
            if key[:3] == 'SNP':
                print key,":",round(float(value-variants[0].variantGenotypes['varscan'].variantCountTumorPrivate[key])/snpTruthSet,4)
            elif key[:5]  == 'Indel':
                print key,":",round(float(value-variants[0].variantGenotypes['varscan'].variantCountTumorPrivate[key])/indelTruthSet,4)
print "\nPRECISION"
for key,value in sorted(variants[0].variantGenotypes['varscan'].variantCountSubTypeTumorPrivate.items()):
            print key,":",round(1-float(value)/float(variants[0].variantGenotypes['varscan'].variantCountSubTypeTumor[key]),4)
print "\nNumber of Callers SNP: ",variants[0].variantCountSNPNumberCallers
print "Number of Callers Indel: ", variants[0].variantCountIndelNumberCallers



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

