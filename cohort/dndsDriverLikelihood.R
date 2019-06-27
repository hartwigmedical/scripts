detach("package:purple", unload=TRUE);
library(purple)
library(GenomicRanges)
library(RMySQL)
library(dplyr)
library(tidyr)

############################################ LOAD NEW DATA
load(file = "~/hmf/analysis/cohort/processed/HmfRefCDSCv.RData")
load(file = "~/hmf/analysis/cohort/processed/dndsUnfilteredAnnotatedMutations.RData")
load(file = "~/hmf/analysis/cohort/reference/highestPurityCohort.RData")
load(file = "~/hmf/analysis/cohort/reference/somaticCounts.RData")
load(file = "~/hmf/analysis/cohort/reference/hpcExonicSomatics.RData")
load(file = "~/hmf/analysis/cohort/processed/genePanel.RData")

genePanel = genePanel %>% filter(!is.na(reportablePointMutation))
tsGenes = genePanel %>% filter(reportablePointMutation == "tsg")
oncoGenes = genePanel %>% filter(reportablePointMutation == "onco")

load(file = "~/hmf/analysis/cohort/processed/HmfRefCDSCvBiallelic.RData")
load(file = "~/hmf/analysis/cohort/processed/HmfRefCDSCvNonBiallelic.RData")
load(file = "~/hmf/analysis/cohort/processed/tsGenesToSeparateBiallelic.RData")

tsGenesBiallelic = genePanel %>% filter(reportablePointMutation == "tsg", gene %in% tsGenesToSeparateBiallelic$gene)
HmfRefCDSCvBiallelic = HmfRefCDSCvBiallelic %>% mutate(cancerType = "All") %>% filter(gene_name %in% tsGenesBiallelic$gene)
HmfRefCDSCvNonBiallelic = HmfRefCDSCvNonBiallelic %>% mutate(cancerType = "All") %>% filter(gene_name %in% tsGenesBiallelic$gene)

############################################ Driver Rates  ############################################
somatics = hpcExonicSomatics %>%
  mutate(subclonalLikelihood = 0) %>%
  select(sampleId, chromosome, position, ref, alt, type, worstCodingEffect, canonicalCodingEffect, hotspot, biallelic, subclonalLikelihood, repeatCount) %>%
  mutate(shared = F)

allMutations = dnds_annotate_somatics(dndsUnfilteredAnnotatedMutations %>% filter(gene %in% genePanel$gene), somatics) %>%
  mutate(patient = substr(sampleId, 1,12), pHGVS = "") %>%
  select(-patient) %>% filter(impact != "")

oncoMutations = onco_mutations(allMutations %>% filter(gene %in% oncoGenes$gene))
oncoDriverRates = dnds_driver_rates(HmfRefCDSCv, oncoMutations)

### OLD METHOD
#tsgMutations = tsg_mutations(allMutations %>% filter(gene %in% tsGenes$gene), T)
#tsgDriverRates = dnds_driver_rates(HmfRefCDSCv, tsgMutations)
##save(tsgDriverRates, oncoDriverRates, file = "~/hmf/analysis/cohort/processed/driverRates.RData")


############################################  Biallelic Driver Rates ############################################
biallelicMutations = allMutations %>% filter(gene %in% tsGenesBiallelic$gene, biallelic)
nonBiallelicMutations = allMutations %>% filter(gene %in% tsGenesBiallelic$gene, !biallelic)

tsgMutations = tsg_mutations(allMutations %>% filter(gene %in% tsGenes$gene), F)
tsgBiallelicMutations = tsg_mutations(biallelicMutations, F)
tsgNonBiallelicMutations = tsg_mutations(nonBiallelicMutations, F)

tsgDriverRates = dnds_driver_rates(HmfRefCDSCv, tsgMutations)
tsgBiallelicDriverRates = dnds_driver_rates(HmfRefCDSCvBiallelic, tsgBiallelicMutations)
tsgNonBiallelicDriverRates = dnds_driver_rates(HmfRefCDSCvNonBiallelic, tsgNonBiallelicMutations)

save(tsgDriverRates, oncoDriverRates, tsgBiallelicDriverRates, tsgNonBiallelicDriverRates, file = "~/hmf/analysis/cohort/processed/driverRates.RData")

############################################  TSG ############################################
load(file = "~/hmf/analysis/cohort/processed/driverRates.RData")
load(file = "~/hmf/analysis/cohort/reference/highestPurityCohort.RData")
load(file = "~/hmf/analysis/cohort/reference/somaticCounts.RData")

tsg_dnds_driver_likelihood <- function(cohortSize, somaticCounts, driverRates, allFields = T) {

  totalSomatics = somaticCounts %>% ungroup() %>% summarise(total_SNV = sum(TOTAL_SNP), total_INDEL = sum(TOTAL_INDEL))
  
  tsgDriverLikelihood = driverRates %>% select(gene, impact, driverLikelihood) %>% 
    spread(impact, driverLikelihood, fill = 0) %>%
    select(gene, dndsDriverLikelihoodMissense = Missense, dndsDriverLikelihoodNonsense = Nonsense, dndsDriverLikelihoodSplice = Splice, dndsDriverLikelihoodIndel = Indel)
  
  tsgPDrivers = driverRates %>% select(gene, impact, gene_drivers) %>%
    mutate(p_driver = gene_drivers / cohortSize) %>%
    select(-gene_drivers) %>%
    spread(impact, p_driver, fill = 0) %>%
    select(gene, pDriverMissense = Missense, pDriverNonsense = Nonsense, pDriverSplice = Splice, pDriverIndel = Indel)
  
  tsgPVariantNonDriverFactor = driverRates %>% select(gene, impact, gene_non_drivers) %>%
    mutate(p_variant_nondriver_factor = ifelse(impact == 'Indel',  gene_non_drivers / totalSomatics$total_INDEL, gene_non_drivers / totalSomatics$total_SNV)) %>%
    select(-gene_non_drivers) %>%
    spread(impact, p_variant_nondriver_factor, fill = 0) %>%
    select(gene, pVariantNonDriverFactorMissense = Missense, pVariantNonDriverFactorNonsense = Nonsense, pVariantNonDriverFactorSplice = Splice, pVariantNonDriverFactorIndel = Indel)
  
  tsg = tsgDriverLikelihood %>% left_join(tsgPDrivers, by = "gene") %>% left_join(tsgPVariantNonDriverFactor, by = "gene") %>%
    select(gene,
           dndsDriverLikelihoodMissense, pDriverMissense, pVariantNonDriverFactorMissense,
           dndsDriverLikelihoodNonsense, pDriverNonsense, pVariantNonDriverFactorNonsense,
           dndsDriverLikelihoodSplice, pDriverSplice, pVariantNonDriverFactorSplice,
           dndsDriverLikelihoodIndel, pDriverIndel, pVariantNonDriverFactorIndel)
  
  if (!allFields) {
    tsg = tsg %>% select(gene, dndsDriverLikelihoodMissense, pDriverMissense, pVariantNonDriverFactorMissense)
  }

  return (tsg)
}

cohortSize = nrow(highestPurityCohort)
tsg = tsg_dnds_driver_likelihood(cohortSize, somaticCounts, tsgDriverRates, T)
tsgBiallelic = tsg_dnds_driver_likelihood(cohortSize, biallelicSomaticCounts, tsgBiallelicDriverRates, F)
tsgNonBiallelic = tsg_dnds_driver_likelihood(cohortSize, nonBiallelicSomaticCounts, tsgNonBiallelicDriverRates, F)

write.table(tsg, "/Users/jon/hmf/repos/hmftools/hmf-common/src/main/resources/dnds/DndsDriverLikelihoodTsg.tsv", quote = F, row.names = F, sep = "\t")
write.table(tsgBiallelic, "/Users/jon/hmf/repos/hmftools/hmf-common/src/main/resources/dnds/DndsDriverLikelihoodTsgBiallelic.tsv", quote = F, row.names = F, sep = "\t")
write.table(tsgNonBiallelic, "/Users/jon/hmf/repos/hmftools/hmf-common/src/main/resources/dnds/DndsDriverLikelihoodTsgNonBiallelic.tsv", quote = F, row.names = F, sep = "\t")

############################################  ONCO ############################################
totalSomatics = somaticCounts %>% ungroup() %>% summarise(total_SNV = sum(TOTAL_SNP), total_INDEL = sum(TOTAL_INDEL))

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
  select(gene, dndsDriverLikelihoodMissense, pDriverMissense, pVariantNonDriverFactorMissense) 

write.table(onco, "/Users/jon/hmf/repos/hmftools/hmf-common/src/main/resources/dnds/DndsDriverLikelihoodOnco.tsv", quote = F, row.names = F, sep = "\t")



############################################  TESTING ############################################
current = tsgMutations
old = tsg_mutations(mutations %>% filter(gene %in% tsGenes$gene), T)
new = tsg_mutations(mutations %>% filter(gene %in% tsGenes$gene), F)

intersting = new %>% filter(biallelic, !knownDriver, driverType !="Redundant")
intersting = old %>% filter(biallelic, !knownDriver, driverType !="Redundant")
intersting = current %>% filter(biallelic, !knownDriver, driverType !="Redundant")





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

