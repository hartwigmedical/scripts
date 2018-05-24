load(file = "~/hmf/RData/Processed/highestPurityCohortSummary.RData")
load(file = "~/hmf/RData/Processed/tsgDrivers.RData")
load(file = "~/hmf/RData/Processed/oncoDrivers.RData")

sampleSomatics = highestPurityCohortSummary %>% select(sampleId, ends_with("SNP"), ends_with("INDEL"))
sampleSomatics[is.na(sampleSomatics)] <- 0
sampleSomatics = sampleSomatics %>% mutate(sample_SNV = INCONSISTENT_SNP + SUBCLONAL_SNP + CLONAL_SNP, sample_INDEL = INCONSISTENT_INDEL + SUBCLONAL_INDEL + CLONAL_INDEL) %>%
  select(starts_with("sample"))

totalSomatics = sampleSomatics %>% ungroup() %>% summarise(total_SNV = sum(sample_SNV), total_INDEL = sum(sample_INDEL))

cohortSize = nrow(highestPurityCohortSummary)



tsgDriverByGene = tsgDrivers %>%  filter(driver > 0)
tsgDriverByGene$firstImpact = sapply(tsgDriverByGene$impact, function (x) {unlist(strsplit(x, split = ","))[[1]]})

tsgSingleHitsByGene = tsgDriverByGene %>% 
  filter(driverLikelihood < 1, driver != 'Multihit') %>%
  ungroup() %>% 
  group_by(gene, driver) %>% summarise(gene_single_drivers = sum(driverLikelihood), gene_single_non_drivers = sum(1 - driverLikelihood))

tsgDriverByGeneBayesian = tsgDriverByGene %>% 
  left_join(tsgSingleHitsByGene, by = c("gene", "firstImpact" = "driver")) %>% 
  left_join(sampleSomatics, by = "sampleId") %>%
  ungroup()

tsgDriverByGeneBayesian$p_variant_nondriver_snv = 1 - ppois(0, tsgDriverByGeneBayesian$sample_SNV / totalSomatics$total_SNV  * tsgDriverByGeneBayesian$gene_single_non_drivers)
tsgDriverByGeneBayesian$p_variant_nondriver_multi = 1 - ppois(1, tsgDriverByGeneBayesian$sample_SNV / totalSomatics$total_SNV  * tsgDriverByGeneBayesian$gene_single_non_drivers)
tsgDriverByGeneBayesian$p_variant_nondriver_indel = 1 - ppois(0, tsgDriverByGeneBayesian$sample_INDEL / totalSomatics$total_INDEL  * tsgDriverByGeneBayesian$gene_single_non_drivers)
tsgDriverByGeneBayesian$p_variant_nondriver = ifelse(tsgDriverByGeneBayesian$driver %in% c("Indel"), tsgDriverByGeneBayesian$p_variant_nondriver_indel, tsgDriverByGeneBayesian$p_variant_nondriver_snv)
tsgDriverByGeneBayesian$p_variant_nondriver = ifelse(tsgDriverByGeneBayesian$driver %in% c("Multihit"), tsgDriverByGeneBayesian$p_variant_nondriver_multi, tsgDriverByGeneBayesian$p_variant_nondriver)

tsgDriverByGeneBayesian = tsgDriverByGeneBayesian %>%
  mutate(p_driver = gene_single_drivers / cohortSize,  p_driver_variant = p_driver / (p_driver + p_variant_nondriver * (1-p_driver))) 

tsgDriverByGeneBayesian$adjustedDriverLikelihood = ifelse(tsgDriverByGeneBayesian$driverLikelihood < 1 & tsgDriverByGeneBayesian$driverLikelihood > 0, tsgDriverByGeneBayesian$p_driver_variant, tsgDriverByGeneBayesian$driverLikelihood)



oncoDriverByGene = oncoDrivers %>%
  group_by(sampleId, gene) %>%
  summarise(
    hotspot = sum(hotspot) > 0,
    nearHotspot = any(nearHotspot),
    impact = paste(impact, collapse = ","),
    driver = max(driver),
    driverLikelihood = max(driver))  %>%
  mutate(driver = impact, type = 'ONCO') %>%
  select(sampleId, gene, impact, driver, driverLikelihood, type) %>% 
  filter(driver != 'Frameshift') %>%
  mutate(driverLikelihood = ifelse(driver == "Inframe", 1, driverLikelihood)) %>%
  filter(driver > 0)

oncoHitsByGene = oncoDriverByGene %>% 
  filter(driverLikelihood < 1, driver == "Missense") %>%
  ungroup() %>% 
  group_by(gene, driver) %>% summarise(gene_single_drivers = sum(driverLikelihood), gene_single_non_drivers = sum(1 - driverLikelihood))

oncoDriverByGeneBayesian = oncoDriverByGene %>% 
  left_join(oncoHitsByGene, by = c("gene", "driver")) %>% 
  left_join(sampleSomatics, by = "sampleId") %>%
  mutate(p_variant_nondriver = 1 - ppois(0, sample_SNV / totalSomatics$total_SNV  * gene_single_non_drivers)) %>% 
  mutate(p_driver = gene_single_drivers / cohortSize,  p_driver_variant = p_driver / (p_driver + p_variant_nondriver * (1-p_driver))) 
oncoDriverByGeneBayesian$adjustedDriverLikelihood = ifelse(oncoDriverByGeneBayesian$driverLikelihood < 1 & oncoDriverByGeneBayesian$driverLikelihood > 0, oncoDriverByGeneBayesian$p_driver_variant, oncoDriverByGeneBayesian$driverLikelihood)


save(tsgDriverByGeneBayesian, oncoDriverByGeneBayesian, file = "~/hmf/RData/adjustedDriverLikelihood.RData")

