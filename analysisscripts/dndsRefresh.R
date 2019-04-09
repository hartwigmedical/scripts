detach("package:purple", unload=TRUE);
library(purple)
library(RMySQL)
library(dplyr)
library(tidyr)
library(dndscv)


####### DB
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis", host = "127.0.0.1")
dbDisconnect(dbProd)
rm(dbProd)

####### Cohort
entireCohort = purple::query_cohort(dbProd)
highestPurityCohort = purple::highest_purity_cohort(entireCohort)
save(highestPurityCohort, file = "~/hmf/analysis/dnds/highestPurityCohort.RData")

sampleIdString = paste("'", highestPurityCohort$sampleId, "'", collapse = ",", sep = "")
somaticQuery = paste0("select * from somaticVariant where filter = 'PASS' and gene <> '' and worstCodingEffect != 'NONE' AND sampleId in (",sampleIdString, ")")
exonicSomatics = dbGetQuery(dbProd, somaticQuery)
save(exonicSomatics, file = "~/hmf/analysis/dnds/exonicSomatics.RData")

somaticCountsQuery = paste0("select sampleId, type, count(*) as total from somaticVariant where filter = 'PASS' and sampleId in (",sampleIdString, ") group by sampleId, type")
somaticCounts = dbGetQuery(dbProd, somaticCountsQuery)
somaticCounts = somaticCounts %>% mutate(type = paste0("TOTAL_", type)) %>% spread(type, total)
save(somaticCounts, file = "~/hmf/analysis/dnds/somaticCounts.RData")

####### DNDS HmfRefCDSCv 
load(file = "~/hmf/analysis/dnds/highestPurityCohort.RData")
load(file = "~/hmf/analysis/dnds/exonicSomatics_question.RData")

somatics = exonicSomatics %>% 
  filter(type == "INDEL" | type == "SNP", repeatCount <= 8) %>%
  select(sampleId, chr = chromosome, pos = position, ref, alt)

data("cancergenes_cgc81", package="dndscv")
data("covariates_hg19", package="dndscv")
refdb="~/hmf/RData/reference/HmfRefCDS.RData"
load(refdb)
newgenes = sapply(RefCDS, function(x) x$gene_name) # New list of genes
oldgenes = sapply(strsplit(newgenes,split=":"), function(x) x[1])
cv = covs[oldgenes,]
rownames(cv) = newgenes
kc = newgenes[oldgenes %in% known_cancergenes]

output = dndscv(somatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)
#save(output, file = "~/hmf/analysis/dnds/dndsOutput.RData")

HmfRefCDSCv  <- output$sel_cv
save(HmfRefCDSCv, file = "~/hmf/analysis/dnds/HmfRefCDSCv.RData")


####### DNDS UNFILTERED ANNOTATED MUTATIONS
unfilteredSomatics = exonicSomatics %>% select(sampleId, chr = chromosome, pos = position, ref, alt)
unfilteredOutput = dndscv(unfilteredSomatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE, max_muts_per_gene_per_sample = 30000000, max_coding_muts_per_sample = 30000000)
dndsUnfilteredAnnotatedMutations = unfilteredOutput$annotmuts
save(dndsUnfilteredAnnotatedMutations, file = "~/hmf/analysis/dnds/exonicSomaticsDndsAnnotated.RData")



############################################   HIGHEST PURITY COHORT ############################################
load(file = "~/hmf/analysis/dnds/somaticCounts.RData")
somaticCounts[is.na(somaticCounts)] <- 0
somaticCounts = somaticCounts %>%
  mutate(sample_SNV = TOTAL_SNP, sample_INDEL = TOTAL_INDEL) %>%
  select(starts_with("sample"))

load(file = "~/hmf/RData/processed/driverGenes.RData")
genePanel = bind_rows(oncoGenes, tsGenes)

load(file = "~/hmf/analysis/dnds/HmfRefCDSCv.RData")
load(file = "~/hmf/analysis/dnds/exonicSomaticsDndsAnnotated.RData")
load(file = "~/hmf/analysis/dnds/exonicSomatics.RData")
HmfRefCDSCv$cancerType <- "All"
somatics = exonicSomatics %>%
  mutate(subclonalLikelihood = 0) %>%
  select(sampleId, chromosome, position, ref, alt, type, worstCodingEffect, canonicalCodingEffect, hotspot, biallelic, subclonalLikelihood, repeatCount) %>%
  mutate(shared = F)
dndsExpectedDriversPerGene = dnds_expected_drivers(HmfRefCDSCv, dndsUnfilteredAnnotatedMutations, somatics)

save(dndsExpectedDriversPerGene, file = "~/hmf/analysis/dnds/dndsExpectedDriversPerGene.RData")

mutations = dnds_annotate_somatics(dndsUnfilteredAnnotatedMutations, somatics) %>%
  mutate(patient = substr(sampleId, 1,12), pHGVS = "") %>%
  select(-patient)
save(mutations, file = "~/hmf/analysis/dnds/dndsMutations.RData")

mutations = mutations %>% filter(gene %in% genePanel$gene_name, impact != "")

dndsTsgDriversOutput = dnds_tsg_drivers(somaticCounts, mutations %>% filter(gene %in% tsGenes$gene_name), dndsExpectedDriversPerGene)
dndsTsgDriverRates = dndsTsgDriversOutput[["tsgDriverRates"]]
dndsTsgUnknownDriversTotals = dndsTsgDriversOutput[["tsgUnknownDriversTotals"]]

dndsOncoDriversOutput = dnds_onco_drivers(somaticCounts, mutations %>% filter(gene %in% oncoGenes$gene_name), dndsExpectedDriversPerGene)
dndsOncoDriverRates = dndsOncoDriversOutput[["oncoDriverRates"]]
dndsOncoUnknownDriversTotals = dndsOncoDriversOutput[["oncoUnknownDriversTotals"]]

############################################  TSG ############################################

cohortSize = nrow(highestPurityCohort)
totalSomatics = somaticCounts %>% ungroup() %>% summarise(total_SNV = sum(sample_SNV), total_INDEL = sum(sample_INDEL))

tsgDriverLikelihood = dndsTsgDriverRates %>% select(gene, impact, driverLikelihood) %>% spread(impact, driverLikelihood, fill = 0) %>%
  select(gene, dndsDriverLikelihoodMissense = Missense, dndsDriverLikelihoodNonsense = Nonsense, dndsDriverLikelihoodSplice = Splice, dndsDriverLikelihoodIndel = Indel)

tsgPDrivers = dndsTsgUnknownDriversTotals %>% select(gene, impact, gene_drivers) %>%
  mutate(p_driver = gene_drivers / cohortSize) %>%
  select(-gene_drivers) %>%
  spread(impact, p_driver, fill = 0) %>%
  select(gene, pDriverMissense = Missense, pDriverNonsense = Nonsense, pDriverSplice = Splice, pDriverIndel = Indel)

tsgPVariantNonDriverFactor = dndsTsgUnknownDriversTotals %>% select(gene, impact, gene_non_drivers) %>%
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
write.table(tsg, "~/hmf/analysis/dnds/DndsDriverLikelihoodTsg.tsv", quote = F, row.names = F, sep = "\t")

############################################  ONCO ############################################

oncoDriverLikelihoodAdjustment = dndsOncoUnknownDriversTotals %>%
  filter(impact == 'Missense') %>%
  mutate(
    p_driver = gene_drivers / cohortSize,
    p_variant_nondriver_factor = gene_non_drivers / totalSomatics$total_SNV 
  ) %>% 
  select(gene, pDriverMissense = p_driver, pVariantNonDriverFactorMissense = p_variant_nondriver_factor) 


oncoDriverLikelihood = dndsOncoDriverRates %>% select(gene, impact, driverLikelihood) %>% spread(impact, driverLikelihood, fill = 0) %>% 
  select(gene, dndsDriverLikelihoodMissense = Missense)

onco = oncoDriverLikelihood %>% left_join(oncoDriverLikelihoodAdjustment, by = "gene") %>%
  select(gene, dndsDriverLikelihoodMissense, pDriverMissense, pVariantNonDriverFactorMissense) %>%
  mutate(dndsDriverLikelihoodNonsense = 0, pDriverNonsense = 0, pVariantNonDriverFactorNonsense = 0,
         dndsDriverLikelihoodSplice = 0, pDriverSplice = 0, pVariantNonDriverFactorSplice = 0,
         dndsDriverLikelihoodIndel = 0, pDriverIndel = 0, pVariantNonDriverFactorIndel = 0) 
write.table(onco, "~/hmf/analysis/dnds/DndsDriverLikelihoodOnco.tsv", quote = F, row.names = F, sep = "\t")



