library(IRanges)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
detach("package:purple", unload=TRUE);
library(purple)
library(RMySQL)

############################################   HIGHEST PURITY COHORT ############################################
load(file = "~/hmf/RData/processed/highestPurityCohortSummary.RData")
hpcSomaticCounts = highestPurityCohortSummary %>% select(sampleId, ends_with("SNV"), ends_with("INDEL"))
hpcSomaticCounts[is.na(hpcSomaticCounts)] <- 0
hpcSomaticCounts = hpcSomaticCounts %>%
  mutate(sample_SNV = TOTAL_SNV, sample_INDEL = TOTAL_INDEL) %>%
  select(starts_with("sample"))

load(file = "~/hmf/RData/processed/driverGenes.RData")
genePanel = bind_rows(oncoGenes, tsGenes)

load(file = "~/hmf/Rdata/Processed/HmfRefCDSCv.RData")
load(file = "~/hmf/RData/processed/hpcExonicSomaticsDndsAnnotated.RData")
load(file = "~/hmf/RData/reference/hpcExonicSomatics.RData")
hpcSomatics = hpcExonicSomatics %>%
  select(sampleId, chromosome, position, ref, alt, type, worstCodingEffect, canonicalCodingEffect, hotspot, biallelic, subclonalLikelihood, repeatCount) %>%
  mutate(shared = F)
hpcDndsExpectedDriversPerGene = dnds_expected_drivers(HmfRefCDSCv, dndsUnfilteredAnnotatedMutations, hpcSomatics)
save(hpcDndsExpectedDriversPerGene, file = "~/hmf/RData/Processed/hpcDndsExpectedDriversPerGene.RData")

phgvs = read.csv("~/hmf/resources/phgvs.csv", stringsAsFactors = F) %>%
  select(patient = SAMPLEID, chromosome = CHROM, position = POS, ref = REF, alt = ALTS, pHGVS = HGVS_PROTEIN) %>%
  mutate(ref = gsub("\\*","", ref)) %>%
  distinct()

hpcMutations = dnds_annotate_somatics(dndsUnfilteredAnnotatedMutations, hpcSomatics) %>%
  mutate(patient = substr(sampleId, 1,12)) %>%
  left_join(phgvs, by = c("patient","chromosome","position","ref","alt")) %>% 
  select(-patient)

#load(file = "~/hmf/RData/Reference/canonicalTranscripts.RData")
#hpcMutations = hpcMutations %>% left_join(canonicalTranscripts %>% select(gene, transcriptId), by = "gene")
#save(hpcMutations, file = "~/hmf/RData/Processed/hpcMutations.RData")
#load(file = "~/hmf/RData/Processed/hpcMutations.RData")

hpcMutations = hpcMutations %>% filter(gene %in% genePanel$gene_name, impact != "")

hpcDndsTsgDriversOutput = dnds_tsg_drivers(hpcSomaticCounts, hpcMutations %>% filter(gene %in% tsGenes$gene_name), hpcDndsExpectedDriversPerGene)
hpcDndsTsgMutations = hpcDndsTsgDriversOutput[["tsgMutations"]]; save(hpcDndsTsgMutations, file = "~/hmf/RData/Processed/hpcDndsTsgMutations.RData")
hpcDndsTsgDriverRates = hpcDndsTsgDriversOutput[["tsgDriverRates"]]; save(hpcDndsTsgDriverRates, file = "~/hmf/RData/Processed/hpcDndsTsgDriverRates.RData")
hpcDndsTsgUnknownDriversTotals = hpcDndsTsgDriversOutput[["tsgUnknownDriversTotals"]]; save(hpcDndsTsgUnknownDriversTotals, file = "~/hmf/RData/Processed/hpcDndsTsgUnknownDriversTotals.RData")
hpcDndsTsgDrivers = hpcDndsTsgDriversOutput[["tsgDrivers"]]; save(hpcDndsTsgDrivers, file = "~/hmf/RData/Processed/hpcDndsTsgDrivers.RData")

hpcDndsOncoDriversOutput = dnds_onco_drivers(hpcSomaticCounts, hpcMutations %>% filter(gene %in% oncoGenes$gene_name), hpcDndsExpectedDriversPerGene)
hpcDndsOncoMutations = hpcDndsOncoDriversOutput[["oncoMutations"]]; save(hpcDndsOncoMutations, file = "~/hmf/RData/Processed/hpcDndsOncoMutations.RData")
hpcDndsOncoDriverRates = hpcDndsOncoDriversOutput[["oncoDriverRates"]]; save(hpcDndsOncoDriverRates, file = "~/hmf/RData/Processed/hpcDndsOncoDriverRates.RData")
hpcDndsOncoUnknownDriversTotals = hpcDndsOncoDriversOutput[["oncoUnknownDriversTotals"]]; save(hpcDndsOncoUnknownDriversTotals, file = "~/hmf/RData/Processed/hpcDndsOncoUnknownDriversTotals.RData")
hpcDndsOncoDrivers = hpcDndsOncoDriversOutput[["oncoDrivers"]]; save(hpcDndsOncoDrivers, file = "~/hmf/RData/Processed/hpcDndsOncoDrivers.RData")

############################################ MULTIPLE BIOPYSY COHORT ############################################
load(file = "~/hmf/RData/reference/multipleBiopsySomaticsWithScope.Rdata")
mbSomaticCounts = multipleBiopsySomaticsWithScope %>%
  ungroup() %>%
  group_by(sampleId, type) %>%
  summarise(count = n()) %>%
  spread(type, count) %>%
  mutate(sample_SNV = SNP + MNP, sample_INDEL = INDEL)

load(file = "~/hmf/RData/processed/driverGenes.RData")
genePanel = bind_rows(oncoGenes, tsGenes)

load(file = "~/hmf/Rdata/Processed/HmfRefCDSCv.RData")
load(file = "~/hmf/RData/processed/mbExonicSomaticsDndsAnnotated.RData")
load(file = "~/hmf/RData/reference/mbExonicSomatics.RData")
mbSomatics = mbExonicSomatics %>%
  select(sampleId, chromosome, position, ref, alt, type, worstCodingEffect, canonicalCodingEffect, hotspot, biallelic, clonality, scope, repeatCount) %>%
  mutate(shared = scope == "Shared")
mbDndsExpectedDriversPerGene = dnds_expected_drivers(HmfRefCDSCv, dndsUnfilteredMultipleBiopsyMutations, mbSomatics)
save(mbDndsExpectedDriversPerGene, file = "~/hmf/RData/Processed/mbDndsExpectedDriversPerGene.RData")

mbMutations = dnds_annotate_somatics(dndsUnfilteredMultipleBiopsyMutations, mbSomatics)
mbMutations = mbMutations %>% filter(gene %in% genePanel$gene_name, impact != "")

mbDndsTsgDriversOutput = dnds_tsg_drivers(mbSomaticCounts, mbMutations %>% filter(gene %in% tsGenes$gene_name), mbDndsExpectedDriversPerGene)
mbDndsTsgMutations = mbDndsTsgDriversOutput[["tsgMutations"]]; save(mbDndsTsgMutations, file = "~/hmf/RData/Processed/mbDndsTsgMutations.RData")
mbDndsTsgDriverRates = mbDndsTsgDriversOutput[["tsgDriverRates"]]; save(mbDndsTsgDriverRates, file = "~/hmf/RData/Processed/mbDndsTsgDriverRates.RData")
mbDndsTsgUnknownDriversTotals = mbDndsTsgDriversOutput[["tsgUnknownDriversTotals"]]; save(mbDndsTsgUnknownDriversTotals, file = "~/hmf/RData/Processed/mbDndsTsgUnknownDriversTotals.RData")
mbDndsTsgDrivers = mbDndsTsgDriversOutput[["tsgDrivers"]]; save(mbDndsTsgDrivers, file = "~/hmf/RData/Processed/mbDndsTsgDrivers.RData")

mbDndsOncoDriversOutput = dnds_onco_drivers(mbSomaticCounts, mbMutations %>% filter(gene %in% oncoGenes$gene_name), mbDndsExpectedDriversPerGene)
mbDndsOncoMutations = mbDndsOncoDriversOutput[["oncoMutations"]]; save(mbDndsOncoMutations, file = "~/hmf/RData/Processed/mbDndsOncoMutations.RData")
mbDndsOncoDriverRates = mbDndsOncoDriversOutput[["oncoDriverRates"]]; save(mbDndsOncoDriverRates, file = "~/hmf/RData/Processed/mbDndsOncoDriverRates.RData")
mbDndsOncoUnknownDriversTotals = mbDndsOncoDriversOutput[["oncoUnknownDriversTotals"]]; save(mbDndsOncoUnknownDriversTotals, file = "~/hmf/RData/Processed/mbDndsOncoUnknownDriversTotals.RData")
mbDndsOncoDrivers = mbDndsOncoDriversOutput[["oncoDrivers"]]; save(mbDndsOncoDrivers, file = "~/hmf/RData/Processed/mbDndsOncoDrivers.RData")