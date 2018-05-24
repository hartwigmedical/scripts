library(IRanges)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
detach("package:purple", unload=TRUE);
library(purple)
library(RMySQL)

load(file = "~/hmf/RData/processed/highestPurityCohortSummary.RData")
cohortSize = nrow(highestPurityCohortSummary)
sampleSomatics = highestPurityCohortSummary %>% select(sampleId, ends_with("SNP"), ends_with("INDEL"))
sampleSomatics[is.na(sampleSomatics)] <- 0
sampleSomatics = sampleSomatics %>% mutate(sample_SNV = INCONSISTENT_SNP + SUBCLONAL_SNP + CLONAL_SNP, sample_INDEL = INCONSISTENT_INDEL + SUBCLONAL_INDEL + CLONAL_INDEL) %>%
  select(starts_with("sample"))
totalSomatics = sampleSomatics %>% ungroup() %>% summarise(total_SNV = sum(sample_SNV), total_INDEL = sum(sample_INDEL))


load(file = "~/hmf/RData/processed/genePanel.RData")
genePanel = genePanel %>% filter(martincorena | hmf | cosmicCurated)

load(file = "~/hmf/Rdata/Processed/HmfRefCDSCv.RData")
load(file = "~/hmf/RData/processed/dndsUnfilteredAnnotatedMutations.RData")
load(file = "~/hmf/RData/reference/hpcExonicSomatics.RData")
somatics = hpcExonicSomatics %>% 
  select(sampleId, chromosome, position, ref, alt, type, worstCodingEffect, canonicalCodingEffect, hotspot, biallelic, clonality)
expectedDriversPerGene = dnds_expected_drivers(HmfRefCDSCv, dndsUnfilteredAnnotatedMutations, somatics)
save(expectedDriversPerGene, file = "~/hmf/RData/Processed/expectedDriversPerGene.RData")

mutations = dnds_annotate_somatics(dndsUnfilteredAnnotatedMutations, somatics)
mutations = mutations %>% filter(gene %in% genePanel$gene_name, impact != "")

tsgMutations = tsg_mutations(mutations)
save(tsgMutations, file = "~/hmf/RData/Processed/tsgMutations.RData")

tsgDriverRates = dnds_driver_likelihood(tsgMutations, expectedDriversPerGene)
save(tsgDriverRates, file = "~/hmf/RData/Processed/tsgDriverRates.RData")

tsgDrivers = tsgMutations %>%
  filter(redundant == F) %>%
  left_join(tsgDriverRates %>% select(gene, impact, driverLikelihood), by = c("gene","impact")) %>% 
  mutate(driverLikelihood = ifelse(knownDriver, 1, driverLikelihood))
  
tsgUnknownDriversTotals = tsgDrivers %>% 
  group_by(gene, impact) %>%
  filter(!knownDriver) %>% 
  summarise(gene_drivers = sum(driverLikelihood), gene_non_drivers = sum(1 - driverLikelihood))
save(tsgUnknownDriversTotals, file = "~/hmf/RData/Processed/tsgUnknownDriversTotals.RData")

tsgDrivers = tsgDrivers %>%
  left_join(sampleSomatics, by = "sampleId") %>%
  left_join(tsgUnknownDriversTotals, by = c("gene","impact")) %>%
  mutate(
    p_variant_nondriver_multihit = 1 - ppois(1, sample_SNV / totalSomatics$total_SNV  * gene_non_drivers),
    p_variant_nondriver_singlehit = 1 - ppois(0, sample_SNV / totalSomatics$total_SNV  * gene_non_drivers),
    p_variant_nondriver_indel = 1 - ppois(0, sample_INDEL / totalSomatics$total_INDEL  * gene_non_drivers),
    p_variant_nondriver = ifelse(impact == "Indel", p_variant_nondriver_indel, p_variant_nondriver_singlehit),
    p_variant_nondriver = ifelse(driverType == "MultiHit", p_variant_nondriver_multihit, p_variant_nondriver),
    p_driver = gene_drivers / cohortSize,  
    p_driver_variant = p_driver / (p_driver + p_variant_nondriver * (1-p_driver)),
    driverLikelihoodAdjusted = ifelse(driverLikelihood > 0 & driverLikelihood < 1, p_driver_variant, driverLikelihood)
  ) %>%
  group_by(sampleId, gene) %>%
  summarise(
    multihit = n() > 1,
    biallelic = any(biallelic), 
    hotspot = any(hotspot), 
    impact = paste(impact, collapse = ","),
    knownDriver = any(knownDriver),
    driverLikelihood = max(driverLikelihood),
    driverLikelihoodAdjusted = max(driverLikelihoodAdjusted)
  ) %>%
  mutate(driver = ifelse(multihit, "Multihit", impact), type = 'TSG') %>%
  select(sampleId, gene, driver, impact, type, multihit, biallelic, hotspot, knownDriver, driverLikelihood, driverLikelihoodAdjusted) %>%
  ungroup()  
save(tsgDrivers, file = "~/hmf/RData/Processed/tsgDrivers.RData")

oncoMutations = onco_mutations(mutations)
save(oncoMutations, file = "~/hmf/RData/Processed/oncoMutations.RData")

oncoDriverRates = dnds_driver_likelihood(oncoMutations, expectedDriversPerGene)
save(oncoDriverRates, file = "~/hmf/RData/Processed/oncoDriverRates.RData")

oncoDrivers = oncoMutations %>%
  filter(redundant == F) %>%
  left_join(oncoDriverRates %>% select(gene, impact, driverLikelihood), by = c("gene","impact")) %>% 
  mutate(driverLikelihood = ifelse(knownDriver, 1, driverLikelihood)) %>%
  mutate(driverLikelihood = ifelse(impact == "Frameshift", 0, driverLikelihood))

oncoUnknownDriversTotals = oncoDrivers %>% 
  group_by(gene, impact) %>%
  filter(!knownDriver) %>% 
  summarise(gene_drivers = sum(driverLikelihood), gene_non_drivers = sum(1 - driverLikelihood))
save(oncoUnknownDriversTotals, file = "~/hmf/RData/Processed/oncoUnknownDriversTotals.RData")

oncoDrivers = oncoDrivers %>%
  left_join(sampleSomatics, by = "sampleId") %>%
  left_join(oncoUnknownDriversTotals, by = c("gene","impact")) %>%
  mutate(
    p_variant_nondriver = 1 - ppois(0, sample_SNV / totalSomatics$total_SNV  * gene_non_drivers),
    p_driver = gene_drivers / cohortSize,  
    p_driver_variant = p_driver / (p_driver + p_variant_nondriver * (1-p_driver)),
    driverLikelihoodAdjusted = ifelse(driverLikelihood > 0 & driverLikelihood < 1, p_driver_variant, driverLikelihood)
  ) %>%
  group_by(sampleId, gene) %>%
  summarise(
    nearHotspot = any(nearHotspot),
    hotspot = sum(hotspot) > 0, 
    impact = paste(impact, collapse = ","),
    knownDriver = any(knownDriver),
    driverLikelihood = max(driverLikelihood),
    driverLikelihoodAdjusted = max(driverLikelihoodAdjusted),
    sample_SNV = max(sample_SNV)
  ) %>%
  mutate(driver = impact, type = 'ONCO') %>%
  select(sampleId, gene, driver, impact, type, hotspot, nearHotspot, knownDriver, driverLikelihood, driverLikelihoodAdjusted, sample_SNV) %>%
  ungroup()  
save(oncoDrivers, file = "~/hmf/RData/Processed/oncoDrivers.RData")
