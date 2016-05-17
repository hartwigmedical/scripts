#!/usr/local/bin/python
#import math
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
        self.tumorVariantSubType = ""
        self.normalVariantSubType = ""

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

        if genotype.variantCountTumor.has_key(self.tumorVariantType + self.tumorVariantSubType + " " + self.caller):
            genotype.variantCountTumor[self.tumorVariantType + self.tumorVariantSubType + " " + self.caller] += 1
        else:
            genotype.variantCountTumor[self.tumorVariantType + self.tumorVariantSubType + " " + self.caller] = 1

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

        if genotype.variantCountNormal.has_key(self.normalVariantType + self.normalVariantSubType + " " + self.caller):
            genotype.variantCountNormal[self.normalVariantType + self.normalVariantSubType + " " + self.caller] += 1
        else:
            genotype.variantCountNormal[self.normalVariantType + self.normalVariantSubType + " " + self.caller] = 1

    def markPrivate(self):
        self.private = True
        if genotype.variantCountTumorPrivate.has_key(self.tumorVariantType + " " + self.caller):
            genotype.variantCountTumorPrivate[self.tumorVariantType + " " + self.caller] += 1
        else:
            genotype.variantCountTumorPrivate[self.tumorVariantType + " " + self.caller] = 1

class somaticVariant:
    variantCountSNPTotal = 0
    variantCountSNPNumberCallers = {}
    variantCountIndelNumberCallers = {}
    variantCountIndelTotal = 0

    countSNPIndelOverlap = 0

    def __init__(self, chrom,pos, id, ref,alt,format,freebayesNormalGenotype,mutectNormalGenotype, \
                 freebayesTumorGenotype,mutectTumorGenotype,strelkaNormalGenotype,varscanNormalGenotype, \
                 strelkaTumorGenotype,varscanTumorGenotype):
        self.tumorCallerCountSNP = 0
        self.tumorCallerCountIndel = 0
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

        if self.tumorCallerCountSNP > 0 and self.tumorCallerCountIndel > 0:
            somaticVariant.countSNPIndelOverlap += 1

def loadVaraintsFromVCF(aPath, aVCFFile,aVariants):
    print "reading vcf file. . .\n"
    with open(aPath + aVCFFile, 'r') as f:
        i=0
        for line in f:
            a = [x for x in line.split('\t')]
            if a[0][:1] != '#':
                aVariants.append(somaticVariant(a[0], a[1], a[2], a[3], a[4],a[8],a[9][:3], \
                        a[10][:3],a[11][:3],a[12][:3],a[13][:3],a[14][:3],a[15][:3],a[16][:3]))

                i=i+1
                if i > 10000000:
                    return 1
    return 1


variants = []
loadVaraintsFromVCF(Path,VCFFile,variants)


print "NORMAL VARIANTS"
for key, value in sorted(variants[0].variantGenotypes['mutect'].variantCountNormal.items()):
    print "%s: %s" % (key, value)
print "\nTUMOR VARIANTS"
for key, value in sorted(variants[0].variantGenotypes['mutect'].variantCountTumor.items()):
    print "%s: %s" % (key, value)
print "\nIndel Total: ",variants[0].variantCountIndelTotal
print "SNP Total: ",variants[0].variantCountSNPTotal
print "SNP-INDEL Overlap: ",variants[0].countSNPIndelOverlap
print "\nPRIVATE VARIANTS"
for key,value in sorted(variants[0].variantGenotypes['mutect'].variantCountTumorPrivate.items()):
        print "%s: %s" % (key, value)
print "\nNumber of Callers SNP: ",variants[0].variantCountSNPNumberCallers
print "Number of Callers Indel: ", variants[0].variantCountIndelNumberCallers


#NOTES ON VCF STAT
#SNP in VCF STAT
    # LEN(GT1) = LEN(GT2) = LEN(REF)
    # Could potentially be a SNP for 1 and a INDEL for another but in practice NOT.
    # RTGTOOLS counts this as an SNP
        #REF =CTTTTTTTTTTTTTTA  ; ALT = CTTTTTTTTTTTTTTT ;   GT: 1/1
# INDELS
    # If one part of a genotype is an INDEL, then gets counted as an indel
# MISSING GENOTYPE = ./.
# Same as Reference = 0/0
    # Mutect does not use
    # Others use for Normal sample only
# Mysteries
    # why does vcfstat exclude 9301 samples! (all Varscan!)
    #  Freebayes has a different missiing genotype count (by 13) for R & T????
    # need to look at 'indels'
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