library(IRanges)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
detach("package:purple", unload=TRUE);
library(purple)
library(RMySQL)

outputDir = "~/garvan/RData/"
outputDir = "~/Documents/LKCGP_projects/RData/"

referenceDir = paste0(outputDir, "reference/")
processedDir = paste0(outputDir, "processed/")

load(paste0(processedDir, "highestPurityCohortSummary.RData"))
load(paste0(processedDir, "driverGenes.RData"))
load(paste0(processedDir, "HmfRefCDSCv.RData"))
load(paste0(processedDir, "dndsUnfilteredAnnotatedMutations.RData"))
load(paste0(referenceDir, "cohortExonicSomatics.RData"))


############################################   HIGHEST PURITY COHORT ############################################
hpcSomaticCounts = highestPurityCohortSummary %>% select(sampleId, ends_with("SNV"), ends_with("INDEL"))
hpcSomaticCounts[is.na(hpcSomaticCounts)] <- 0
hpcSomaticCounts = hpcSomaticCounts %>%
    mutate(sample_SNV = TOTAL_SNV, sample_INDEL = TOTAL_INDEL) %>%
    select(starts_with("sample"))

hpcSomatics = cohortExonicSomatics %>%
  filter(sampleId %in% highestPurityCohortSummary$sampleId) %>%
    select(sampleId, chromosome, position, ref, alt, type, worstCodingEffect, canonicalCodingEffect, hotspot, biallelic, clonality, repeatCount) %>%
    mutate(shared = F)
hpcDndsExpectedDriversPerGene = dnds_expected_drivers(HmfRefCDSCv, dndsUnfilteredAnnotatedMutations, hpcSomatics)
save(hpcDndsExpectedDriversPerGene, file=paste0(processedDir, "hpcDndsExpectedDriversPerGene.RData"))

hpcMutations = dnds_annotate_somatics(dndsUnfilteredAnnotatedMutations, hpcSomatics) %>% mutate(pHGVS = "")
hpcMutations = hpcMutations %>% filter(gene %in% tsGenes$gene_name | gene %in% oncoGenes$gene_name , impact != "")

hpcDndsTsgDriversOutput = dnds_tsg_drivers(hpcSomaticCounts, hpcMutations %>% filter(gene %in% tsGenes$gene_name), hpcDndsExpectedDriversPerGene)
hpcDndsTsgMutations = hpcDndsTsgDriversOutput[["tsgMutations"]]; 
hpcDndsTsgDriverRates = hpcDndsTsgDriversOutput[["tsgDriverRates"]]; 
hpcDndsTsgUnknownDriversTotals = hpcDndsTsgDriversOutput[["tsgUnknownDriversTotals"]]; 
hpcDndsTsgDrivers = hpcDndsTsgDriversOutput[["tsgDrivers"]]; 

save(hpcDndsTsgMutations, file=paste0(processedDir, "hpcDndsTsgMutations.RData"))
save(hpcDndsTsgDriverRates, file=paste0(processedDir, "hpcDndsTsgDriverRates.RData"))
save(hpcDndsTsgUnknownDriversTotals, file=paste0(processedDir, "hpcDndsTsgUnknownDriversTotals.RData"))
save(hpcDndsTsgDrivers, file=paste0(processedDir, "hpcDndsTsgDrivers.RData"))


hpcDndsOncoDriversOutput = dnds_onco_drivers(hpcSomaticCounts, hpcMutations %>% filter(gene %in% oncoGenes$gene_name), hpcDndsExpectedDriversPerGene)
hpcDndsOncoMutations = hpcDndsOncoDriversOutput[["oncoMutations"]]; 
hpcDndsOncoDriverRates = hpcDndsOncoDriversOutput[["oncoDriverRates"]];
hpcDndsOncoUnknownDriversTotals = hpcDndsOncoDriversOutput[["oncoUnknownDriversTotals"]]; 
hpcDndsOncoDrivers = hpcDndsOncoDriversOutput[["oncoDrivers"]]; 

save(hpcDndsOncoMutations, file=paste0(processedDir, "hpcDndsOncoMutations.RData"))
save(hpcDndsOncoDriverRates, file=paste0(processedDir, "hpcDndsOncoDriverRates.RData"))
save(hpcDndsOncoUnknownDriversTotals, file=paste0(processedDir, "hpcDndsOncoUnknownDriversTotals.RData"))
save(hpcDndsOncoDrivers, file=paste0(processedDir, "hpcDndsOncoDrivers.RData"))


