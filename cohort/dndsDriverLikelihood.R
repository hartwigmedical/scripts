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

somaticCounts[is.na(somaticCounts)] <- 0
somaticCounts = somaticCounts %>%
  mutate(sample_SNV = TOTAL_SNP, sample_INDEL = TOTAL_INDEL) %>%
  select(starts_with("sample"))


load(file = "~/hmf/analysis/cohort/processed/HmfRefCDSCvBiallelic.RData")
load(file = "~/hmf/analysis/cohort/processed/HmfRefCDSCvNonBiallelic.RData")
load(file = "~/hmf/analysis/cohort/processed/tsGenesToSeparateBiallelic.RData")
#load(file = "~/hmf/analysis/cohort/processed/dndsUnfilteredBiallelicMutations.RData")
#load(file = "~/hmf/analysis/cohort/processed/dndsUnfilteredNonBiallelicMutations.RData")

genePanel = genePanel %>% filter(!is.na(reportablePointMutation))
tsGenes = genePanel %>% filter(reportablePointMutation == "tsg")
oncoGenes = genePanel %>% filter(reportablePointMutation == "onco")

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

### CAREFUL ABOUT WHICH ONE OF THESE YOU CHOOSE
tsgMutations = tsg_mutations(allMutations %>% filter(gene %in% tsGenes$gene), F)
tsgMutations = tsg_mutations(allMutations %>% filter(gene %in% tsGenes$gene), T)
tsgDriverRates = dnds_driver_rates(HmfRefCDSCv, tsgMutations)

oncoMutations = onco_mutations(allMutations %>% filter(gene %in% oncoGenes$gene))
oncoDriverRates = dnds_driver_rates(HmfRefCDSCv, oncoMutations)

save(tsgDriverRates, oncoDriverRates, file = "~/hmf/analysis/cohort/processed/driverRates.RData")
load(file = "~/hmf/analysis/cohort/processed/driverRates.RData")

############################################  Biallelic Driver Rates ############################################
biallelicMutations = allMutations %>% filter(gene %in% tsGenesBiallelic$gene, biallelic)
nonBiallelicMutations = allMutations %>% filter(gene %in% tsGenesBiallelic$gene, !biallelic)
tsgBiallelicMutations = tsg_mutations(biallelicMutations, F)
tsgNonBiallelicMutations = tsg_mutations(nonBiallelicMutations, F)
tsgBiallelicDriverRates = dnds_driver_rates(HmfRefCDSCvBiallelic, tsgBiallelicMutations) %>% filter(impact == "Missense")
tsgNonBiallelicDriverRates = dnds_driver_rates(HmfRefCDSCvNonBiallelic, tsgNonBiallelicMutations) %>% filter(impact == "Missense")
save(tsgDriverRates, tsgBiallelicDriverRates, tsgNonBiallelicDriverRates, file = "~/hmf/analysis/cohort/processed/newDriverRates.RData")


############################################  TESTING ############################################
load(file = "~/hmf/analysis/cohort/processed/newDriverRates.RData")

normal = tsgDriverRates %>% filter(impact == "Missense")

combined =  left_join(tsgBiallelicDriverRates, tsgNonBiallelicDriverRates, by = "gene", suffix = c(".bi", ".non")) %>%
  left_join(normal, by = "gene")
names(combined)

uhoh = combined %>% filter(driverLikelihood.bi < driverLikelihood.non)


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

