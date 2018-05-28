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

load(file = "~/hmf/RData/processed/driverGenes.RData")
genePanel = bind_rows(oncoGenes, tsGenes)

load(file = "~/hmf/Rdata/Processed/HmfRefCDSCv.RData")
load(file = "~/hmf/RData/processed/hpcExonicSomaticsDndsAnnotated.RData")
load(file = "~/hmf/RData/reference/hpcExonicSomatics.RData")
somatics = hpcExonicSomatics %>% 
  select(sampleId, chromosome, position, ref, alt, type, worstCodingEffect, canonicalCodingEffect, hotspot, biallelic, clonality)
expectedDriversPerGene = dnds_expected_drivers(HmfRefCDSCv, dndsUnfilteredAnnotatedMutations, somatics)
save(expectedDriversPerGene, file = "~/hmf/RData/Processed/expectedDriversPerGene.RData")

mutations = dnds_annotate_somatics(dndsUnfilteredAnnotatedMutations, somatics)
mutations = mutations %>% filter(gene %in% genePanel$gene_name, impact != "")

hpcTsgMutations = tsg_mutations(mutations %>% filter(gene %in% tsGenes$gene_name))
save(hpcTsgMutations, file = "~/hmf/RData/Processed/hpcTsgMutations.RData")

hpcTsgDriverRates = dnds_driver_likelihood(hpcTsgMutations, expectedDriversPerGene)
save(hpcTsgDriverRates, file = "~/hmf/RData/Processed/hpcTsgDriverRates.RData")

hpcTsgDrivers = hpcTsgMutations %>%
  filter(redundant == F) %>%
  left_join(hpcTsgDriverRates %>% select(gene, impact, driverLikelihood), by = c("gene","impact")) %>% 
  mutate(driverLikelihood = ifelse(knownDriver, 1, driverLikelihood))
  
hpcTsgUnknownDriversTotals = hpcTsgDrivers %>% 
  group_by(gene, impact) %>%
  filter(!knownDriver) %>% 
  summarise(gene_drivers = sum(driverLikelihood), gene_non_drivers = sum(1 - driverLikelihood))
save(hpcTsgUnknownDriversTotals, file = "~/hmf/RData/Processed/hpcTsgUnknownDriversTotals.RData")

hpcTsgDrivers = hpcTsgDrivers %>%
  left_join(sampleSomatics, by = "sampleId") %>%
  left_join(hpcTsgUnknownDriversTotals, by = c("gene","impact")) %>%
  mutate(
    p_variant_nondriver_snv_single = 1 - ppois(0, sample_SNV / totalSomatics$total_SNV  * gene_non_drivers),
    p_variant_nondriver_snv_multi = 1 - ppois(1, sample_SNV / totalSomatics$total_SNV  * gene_non_drivers),
    p_variant_nondriver_indel_single = 1 - ppois(0, sample_INDEL / totalSomatics$total_INDEL  * gene_non_drivers),
    p_variant_nondriver_indel_multi = 1 - ppois(1, sample_INDEL / totalSomatics$total_INDEL  * gene_non_drivers),
    p_variant_nondriver_single = ifelse(impact == "Indel", p_variant_nondriver_indel_single, p_variant_nondriver_snv_single),
    p_variant_nondriver_multi = ifelse(impact == "Indel", p_variant_nondriver_indel_multi, p_variant_nondriver_snv_multi)
  ) %>%
  group_by(sampleId, gene, impact) %>%
  mutate(
    sameImpact = n() > 1,
    p_variant_nondriver = ifelse(sameImpact, p_variant_nondriver_multi, p_variant_nondriver_single)
  )%>%
  ungroup() %>%
  group_by(sampleId, gene) %>%
  summarise(
    multihit = n() > 1,
    biallelic = any(biallelic), 
    hotspot = any(hotspot), 
    impact = paste(impact, collapse = ","),
    knownDriver = any(knownDriver),
    driverLikelihood = max(driverLikelihood),
    sameImpact = any(sameImpact),
    p_variant_nondriver = ifelse(sameImpact, max(p_variant_nondriver), prod(p_variant_nondriver)),
    gene_drivers = max(gene_drivers),
    sample_SNV = max(sample_SNV),
    clonality = first(clonality)
  ) %>%
  mutate(
    p_driver = gene_drivers / cohortSize,  
    p_driver_variant = p_driver / (p_driver + p_variant_nondriver * (1-p_driver)),
    driverLikelihoodAdjusted = ifelse(driverLikelihood > 0 & driverLikelihood < 1, p_driver_variant, driverLikelihood),
    driver = ifelse(multihit, "Multihit", impact), type = 'TSG') %>%
  select(sampleId, gene, driver, impact, type, multihit, biallelic, hotspot, clonality, knownDriver, driverLikelihood, driverLikelihoodAdjusted, sample_SNV) %>%
  ungroup()  
save(hpcTsgDrivers, file = "~/hmf/RData/Processed/hpcTsgDrivers.RData")



hpcOncoMutations = onco_mutations(mutations %>% filter(gene %in% oncoGenes$gene_name))
save(hpcOncoMutations, file = "~/hmf/RData/Processed/hpcOncoMutations.RData")

hpcOncoDriverRates = dnds_driver_likelihood(hpcOncoMutations, expectedDriversPerGene)
save(hpcOncoDriverRates, file = "~/hmf/RData/Processed/hpcOncoDriverRates.RData")

hpcOncoDrivers = hpcOncoMutations %>%
  filter(redundant == F) %>%
  left_join(hpcOncoDriverRates %>% select(gene, impact, driverLikelihood), by = c("gene","impact")) %>% 
  mutate(driverLikelihood = ifelse(knownDriver, 1, driverLikelihood)) %>%
  mutate(driverLikelihood = ifelse(impact == "Frameshift", 0, driverLikelihood))

hpcOncoUnknownDriversTotals = hpcOncoDrivers %>% 
  group_by(gene, impact) %>%
  filter(!knownDriver) %>% 
  summarise(gene_drivers = sum(driverLikelihood), gene_non_drivers = sum(1 - driverLikelihood))
save(hpcOncoUnknownDriversTotals, file = "~/hmf/RData/Processed/hpcOncoUnknownDriversTotals.RData")

hpcOncoDrivers = hpcOncoDrivers %>%
  left_join(sampleSomatics, by = "sampleId") %>%
  left_join(hpcOncoUnknownDriversTotals, by = c("gene","impact")) %>%
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
    sample_SNV = max(sample_SNV),
    clonality = first(clonality)
  ) %>%
  mutate(driver = impact, type = 'ONCO') %>%
  select(sampleId, gene, driver, impact, type, hotspot, nearHotspot, clonality, knownDriver, driverLikelihood, driverLikelihoodAdjusted, sample_SNV) %>%
  ungroup()  
save(hpcOncoDrivers, file = "~/hmf/RData/Processed/hpcOncoDrivers.RData")
