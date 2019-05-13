detach("package:purple", unload=TRUE);
library(purple)
library(RMySQL)
library(dplyr)
library(tidyr)
library(dndscv)
library(GenomicRanges)


############################################ LOAD NEW DATA
load(file = "~/hmf/analysis/genepanel/HmfRefCDSCv.RData")
load(file = "~/hmf/analysis/genepanel/dndsUnfilteredAnnotatedMutations.RData")
load(file = "~/hmf/analysis/genepanel/highestPurityCohort.RData")
load(file = "~/hmf/analysis/genepanel/somaticCounts.RData")
load(file = "~/hmf/analysis/genepanel/highestPurityCohort.RData")
load(file = "~/hmf/analysis/genepanel/somaticCounts.RData")
load(file = "~/hmf/analysis/genepanel/driverGenes.RData")
load(file = "~/hmf/analysis/genepanel/hpcExonicSomatics.RData")

somaticCounts[is.na(somaticCounts)] <- 0
somaticCounts = somaticCounts %>%
  mutate(sample_SNV = TOTAL_SNP, sample_INDEL = TOTAL_INDEL) %>%
  select(starts_with("sample"))
genePanel = bind_rows(oncoGenes, tsGenes)



############################################CODE RUNS FROM HERE  ############################################
somatics = hpcExonicSomatics %>%
  mutate(subclonalLikelihood = 0) %>%
  select(sampleId, chromosome, position, ref, alt, type, worstCodingEffect, canonicalCodingEffect, hotspot, biallelic, subclonalLikelihood, repeatCount) %>%
  mutate(shared = F)

mutations = dnds_annotate_somatics(dndsUnfilteredAnnotatedMutations %>% filter(gene %in% genePanel$gene_name), somatics) %>%
  mutate(patient = substr(sampleId, 1,12), pHGVS = "") %>%
  select(-patient) %>% filter(impact != "")

tsgMutations = tsg_mutations(mutations %>% filter(gene %in% tsGenes$gene_name))
tsgDriverRates = dnds_driver_rates(HmfRefCDSCv, tsgMutations)

oncoMutations = onco_mutations(mutations %>% filter(gene %in% oncoGenes$gene_name))
oncoDriverRates = dnds_driver_rates(HmfRefCDSCv, oncoMutations)

save(tsgDriverRates, oncoDriverRates, file = "~/hmf/analysis/genepanel/driverRates.RData")

############################################  TSG ############################################
cohortSize = nrow(highestPurityCohort)
totalSomatics = somaticCounts %>% ungroup() %>% summarise(total_SNV = sum(sample_SNV), total_INDEL = sum(sample_INDEL))

tsgDriverLikelihood = tsgDriverRates %>% select(gene, impact, driverLikelihood) %>% 
  spread(impact, driverLikelihood, fill = 0) %>%
  select(gene, dndsDriverLikelihoodMissense = Missense, dndsDriverLikelihoodNonsense = Nonsense, dndsDriverLikelihoodSplice = Splice, dndsDriverLikelihoodIndel = Indel)

tsgPDrivers = tsgDriverRates %>% select(gene, impact, gene_drivers) %>%
  mutate(p_driver = gene_drivers / cohortSize) %>%
  select(-gene_drivers) %>%
  spread(impact, p_driver, fill = 0) %>%
  select(gene, pDriverMissense = Missense, pDriverNonsense = Nonsense, pDriverSplice = Splice, pDriverIndel = Indel)

tsgPVariantNonDriverFactor = tsgDriverRates %>% select(gene, impact, gene_non_drivers) %>%
  mutate(p_variant_nondriver_factor = ifelse(impact == 'Indel',  gene_non_drivers / totalSomatics$total_INDEL, gene_non_drivers / totalSomatics$total_SNV)) %>%
  select(-gene_non_drivers) %>%
  spread(impact, p_variant_nondriver_factor, fill = 0) %>%
  select(gene, pVariantNonDriverFactorMissense = Missense, pVariantNonDriverFactorNonsense = Nonsense, pVariantNonDriverFactorSplice = Splice, pVariantNonDriverFactorIndel = Indel)

tsg = tsgDriverLikelihood %>% left_join(tsgPDrivers, by = "gene") %>% left_join(tsgPVariantNonDriverFactor, by = "gene") %>%
  select(gene,
         dndsDriverLikelihoodMissense, pDriverMissense, pVariantNonDriverFactorMissense,
         dndsDriverLikelihoodNonsense, pDriverNonsense, pVariantNonDriverFactorNonsense,
         dndsDriverLikelihoodSplice, pDriverSplice, pVariantNonDriverFactorSplice,
         dndsDriverLikelihoodIndel, pDriverIndel, pVariantNonDriverFactorIndel
  )
write.table(tsg, "/Users/jon/hmf/repos/hmftools/hmf-common/src/main/resources/dnds/DndsDriverLikelihoodTsg.tsv", quote = F, row.names = F, sep = "\t")

############################################  ONCO ############################################
oncoDriverLikelihood = oncoDriverRates %>% select(gene, impact, driverLikelihood) %>% 
  spread(impact, driverLikelihood, fill = 0) %>% 
  select(gene, dndsDriverLikelihoodMissense = Missense)

oncoDriverLikelihoodAdjustment = oncoDriverRates %>%
  filter(impact == 'Missense') %>%
  mutate(
    p_driver = gene_drivers / cohortSize,
    p_variant_nondriver_factor = gene_non_drivers / totalSomatics$total_SNV 
  ) %>% 
  select(gene, pDriverMissense = p_driver, pVariantNonDriverFactorMissense = p_variant_nondriver_factor) 

onco = oncoDriverLikelihood %>% left_join(oncoDriverLikelihoodAdjustment, by = "gene") %>%
  select(gene, dndsDriverLikelihoodMissense, pDriverMissense, pVariantNonDriverFactorMissense) %>%
  mutate(dndsDriverLikelihoodNonsense = 0, pDriverNonsense = 0, pVariantNonDriverFactorNonsense = 0,
         dndsDriverLikelihoodSplice = 0, pDriverSplice = 0, pVariantNonDriverFactorSplice = 0,
         dndsDriverLikelihoodIndel = 0, pDriverIndel = 0, pVariantNonDriverFactorIndel = 0) 
write.table(onco, "/Users/jon/hmf/repos/hmftools/hmf-common/src/main/resources/dnds/DndsDriverLikelihoodOnco.tsv", quote = F, row.names = F, sep = "\t")





############################################ LOAD OLD DATA
load(file = "~/hmf/RData/processed/highestPurityCohortSummary.RData")
load(file = "~/hmf/RData/Processed/hpcDndsTsgDriverRates.RData")
load(file = "~/hmf/RData/Processed/hpcDndsTsgUnknownDriversTotals.RData")
load(file = "~/hmf/RData/reference/hpcExonicSomatics.RData")
load(file = "~/hmf/RData/processed/hpcExonicSomaticsDndsAnnotated.RData")
load(file = "~/hmf/RData/processed/driverGenes.RData")
load(file = "~/hmf/Rdata/Processed/HmfRefCDSCv.RData")

highestPurityCohort = highestPurityCohortSummary
somaticCounts = highestPurityCohortSummary %>% select(sampleId, ends_with("SNV"), ends_with("INDEL"))
somaticCounts = somaticCounts %>%
  mutate(sample_SNV = TOTAL_SNV, sample_INDEL = TOTAL_INDEL) %>%
  select(starts_with("sample"))
genePanel = bind_rows(oncoGenes, tsGenes)

