load(file = "~/hmf/RData/Processed/highestPurityCohortSummary.RData")
load(file = "~/hmf/RData/Processed/tsgDrivers.RData")
load(file = "~/hmf/RData/Processed/oncoDrivers.RData")

sampleSomatics = highestPurityCohortSummary %>% select(sampleId, ends_with("SNP"), ends_with("INDEL"))
sampleSomatics[is.na(sampleSomatics)] <- 0
sampleSomatics = sampleSomatics %>% mutate(sample_SNV = INCONSISTENT_SNP + SUBCLONAL_SNP + CLONAL_SNP, sample_INDEL = INCONSISTENT_INDEL + SUBCLONAL_INDEL + CLONAL_INDEL) %>%
  select(starts_with("sample"))

totalSomatics = sampleSomatics %>% ungroup() %>% summarise(total_SNV = sum(sample_SNV), total_INDEL = sum(sample_INDEL))

cohortSize = nrow(highestPurityCohortSummary)

tsgDriverByGene = tsgDrivers %>%
  group_by(sampleId, gene) %>%
  summarise(
    multihit = n() > 1,
    biallelic = sum(biallelic) > 0, 
    hotspot = sum(hotspot) > 0, 
    impact = paste(impact, collapse = ","),
    driver = max(driver),
    driverLikelihood = max(driver)) %>%
  mutate(driver = ifelse(multihit, "Multihit", impact), type = 'TSG') %>%
  select(sampleId, gene, impact, driver, driverLikelihood, type, biallelic)

tsgSingleHitsByGene = tsgDriverByGene %>% 
  filter(driverLikelihood < 1, driver != 'Multihit') %>%
  ungroup() %>% 
  group_by(gene, driver) %>% summarise(gene_single_drivers = sum(driverLikelihood), gene_single_non_drivers = sum(1 - driverLikelihood))

tsgMultiHitsByGene = tsgDriverByGene %>% 
  filter(driverLikelihood < 1, driver == 'Multihit') %>%
  ungroup() %>% 
  group_by(gene, driver) %>% summarise(gene_multihit_drivers = sum(driverLikelihood), gene_multihit_non_drivers = sum(1 - driverLikelihood))

tsgDriverByGeneBayesian = tsgDriverByGene %>% left_join(tsgSingleHitsByGene, by = c("gene", "driver")) %>% left_join(tsgMultiHitsByGene, by = c("gene", "driver")) %>% left_join(sampleSomatics, by = "sampleId")
tsgDriverByGeneBayesian$p_variant_nondriver = 1 - ppois(0, tsgDriverByGeneBayesian$sample_SNV / totalSomatics$total_SNV  * tsgDriverByGeneBayesian$gene_single_non_drivers)
tsgDriverByGeneBayesian$p_variant_nondriver = ifelse(tsgDriverByGeneBayesian$driver %in% c("Frameshift","Inframe"), 1 - ppois(0, tsgDriverByGeneBayesian$sample_INDEL / totalSomatics$total_INDEL  * tsgDriverByGeneBayesian$gene_single_non_drivers), tsgDriverByGeneBayesian$p_variant_nondriver)
tsgDriverByGeneBayesian$p_variant_nondriver = ifelse(tsgDriverByGeneBayesian$driver %in% c("Multihit"), 1 - ppois(1, tsgDriverByGeneBayesian$sample_SNV / totalSomatics$total_SNV  * (tsgDriverByGeneBayesian$gene_single_non_drivers + tsgDriverByGeneBayesian$gene_multihit_non_drivers)), tsgDriverByGeneBayesian$p_variant_nondriver)
tsgDriverByGeneBayesian = tsgDriverByGeneBayesian %>%
  mutate(p_driver = ifelse(driver == "Multihit", gene_multihit_drivers / cohortSize, gene_single_drivers / cohortSize),  p_driver_variant = p_driver / (p_driver + p_variant_nondriver * (1-p_driver))) 

tsgDriverByGeneBayesian$adjustedDriverLikelihood = ifelse(tsgDriverByGeneBayesian$driverLikelihood < 1, tsgDriverByGeneBayesian$p_driver_variant, tsgDriverByGeneBayesian$driverLikelihood)


oncoDriverByGene = oncoDrivers %>%
  group_by(sampleId, gene) %>%
  summarise(
    hotspot = sum(hotspot) > 0,
    nearHotspot = any(nearHotspot),
    impact = paste(impact, collapse = ","),
    driver = max(driver),
    driverLikelihood = max(driver))  %>%
  mutate(driver = impact, type = 'ONCO') %>%
  select(sampleId, gene, impact, driver, driverLikelihood, type)

oncoHitsByGene = oncoDriverByGene %>% 
  filter(driverLikelihood < 1, driver == "Missense") %>%
  ungroup() %>% 
  group_by(gene, driver) %>% summarise(gene_single_drivers = sum(driverLikelihood), gene_single_non_drivers = sum(1 - driverLikelihood))

oncoDriverByGeneBayesian = oncoDriverByGene %>% left_join(oncoHitsByGene, by = c("gene", "driver")) %>% left_join(sampleSomatics, by = "sampleId")
oncoDriverByGeneBayesian$p_variant_nondriver = 1 - ppois(0, oncoDriverByGeneBayesian$sample_SNV / totalSomatics$total_SNV  * oncoDriverByGeneBayesian$gene_single_non_drivers)
oncoDriverByGeneBayesian = oncoDriverByGeneBayesian %>% mutate(p_driver = gene_single_drivers / cohortSize,  p_driver_variant = p_driver / (p_driver + p_variant_nondriver * (1-p_driver))) 
oncoDriverByGeneBayesian$adjustedDriverLikelihood = ifelse(!is.na(oncoDriverByGeneBayesian$p_driver_variant) & oncoDriverByGeneBayesian$driverLikelihood < 1, oncoDriverByGeneBayesian$p_driver_variant, oncoDriverByGeneBayesian$driverLikelihood)



