detach("package:purple", unload=TRUE);
library(purple)
library(RMySQL)
library(dplyr)
library(tidyr)
library(dndscv)
library(GenomicRanges)

create_data_frame <- function(cvList) {
  result = data.frame(stringsAsFactors = F)
  for (cancerType in names(cvList)) {
    df = cvList[[cancerType]]
    df$cancerType <- cancerType
    result = rbind(result, df)  
  }
  return (result)
}

####### DNDS HmfRefCDSCv 
load(file = "~/hmf/analysis/cohort/reference/highestPurityCohort.RData")
load(file = "~/hmf/analysis/cohort/reference/hpcExonicSomatics.RData")

#Verify we have all the somatics
length(unique(hpcExonicSomatics$sampleId)) == nrow(highestPurityCohort)

somatics = hpcExonicSomatics %>% 
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
HmfRefCDSCvList = list()
HmfRefCDSCvList[["All"]]  <- output$sel_cv
cancerTypes = unique(highestPurityCohort$cancerType)
cancerTypes = cancerTypes[!is.na(cancerTypes)]
for (selectedCancerType in cancerTypes) {
  cat("Processing", selectedCancerType)
  cancerTypeSampleIds =  highestPurityCohort %>% filter(!is.na(cancerType), cancerType == selectedCancerType) %>% select(sampleId)
  input = somatics %>% filter(sampleId %in% cancerTypeSampleIds$sampleId)
  output = dndscv(input, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)
  
  HmfRefCDSCvList[[selectedCancerType]] <- output$sel_cv
}

HmfRefCDSCv  <- create_data_frame(HmfRefCDSCvList)
save(HmfRefCDSCv, file = "~/hmf/analysis/cohort/processed/HmfRefCDSCv.RData")


####### DNDS UNFILTERED ANNOTATED MUTATIONS
unfilteredSomatics = hpcExonicSomatics %>% select(sampleId, chr = chromosome, pos = position, ref, alt)
unfilteredOutput = dndscv(unfilteredSomatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE, max_muts_per_gene_per_sample = 30000000, max_coding_muts_per_sample = 30000000)
dndsUnfilteredAnnotatedMutations = unfilteredOutput$annotmuts
save(dndsUnfilteredAnnotatedMutations, file = "~/hmf/analysis/cohort/processed/dndsUnfilteredAnnotatedMutations.RData")




############################################   TESTING ############################################

load(file = "~/hmf/analysis/genepanel/HmfRefCDSCv.RData")



biallelicSomatics = exonicSomatics %>% 
  filter(type == "SNP", repeatCount <= 8, biallelic) %>%
  select(sampleId, chr = chromosome, pos = position, ref, alt)

biallelicOutput = dndscv(biallelicSomatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)
HmfRefCDSCvBiallelic  <- biallelicOutput$sel_cv
save(HmfRefCDSCvBiallelic, file = "~/hmf/analysis/dnds/HmfRefCDSCvBiallelic.RData")

nonBiallelicSomatics = exonicSomatics %>% 
  filter(type == "SNP", repeatCount <= 8, !biallelic) %>%
  select(sampleId, chr = chromosome, pos = position, ref, alt)

nonBiallelicOutput = dndscv(nonBiallelicSomatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)
HmfRefCDSCvNonBiallelic  <- nonBiallelicOutput$sel_cv
save(HmfRefCDSCvNonBiallelic, file = "~/hmf/analysis/dnds/HmfRefCDSCvNonBiallelic.RData")




#load('~/hmf/RData/Processed/genePanel.RData')
#load('Downloads/HmfRefCDSCv.RData')
#Sig = HmfRefCDSCv %>% filter(qglobal_cv<0.01,!gene_name %in% c('POM121L12','TRIM49B'))
#View(merge(genePanel,Sig,by='gene_name',all.y = T) %>% filter(is.na(hmf),is.na(martincorena)) %>% select(cancerType,qglobal_cv,everything()))



############################################   HIGHEST PURITY COHORT ############################################
load(file = "~/hmf/analysis/genepanel/highestPurityCohort.RData")
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

####### JON JON
#dndsExpectedDriversPerGene = dnds_expected_drivers(HmfRefCDSCv, dndsUnfilteredAnnotatedMutations, somatics)
#mutations = dnds_annotate_somatics(dndsUnfilteredAnnotatedMutations, somatics) %>%
#  mutate(patient = substr(sampleId, 1,12), pHGVS = "") %>%
#  select(-patient) 
#mutations = mutations %>% filter(gene %in% genePanel$gene_name, impact != "")

############################################CODE RUNS FROM HERE  ############################################
somatics = exonicSomatics %>%
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

#### CHANGE IN SIGNIFICANT GENES

load(file="~/hmf/RData/processed/HmfRefCDSCv.RData")
sig = 0.01
oldHmfSignificant =  HmfRefCDSCv %>% filter(qglobal_cv < sig, cancerType == 'All') %>% distinct(gene_name)
oldHmfSignificant$InOld = T

load(file = "~/hmf/analysis/dnds/HmfRefCDSCv.RData")
newHmfSignificant =  HmfRefCDSCv %>% filter(qglobal_cv < sig) %>% distinct(gene_name)
newHmfSignificant$InNew = T


combinedSignificant = full_join(newHmfSignificant %>% select(gene_name, InNew), oldHmfSignificant %>% select(gene_name, InOld), by = "gene_name")
combinedSignificant[is.na(combinedSignificant)] <- F

combinedSignificant %>% filter(!InNew | !InOld)



##### IMPACT
load(file = "~/hmf/RData/processed/driverGenes.RData")
genePanel = bind_rows(oncoGenes, tsGenes)

load(file = "~/hmf/Rdata/Processed/HmfRefCDSCv.RData")
jon = HmfRefCDSCv %>% filter(gene_name == "AR")

Old = HmfRefCDSCv %>% filter(cancerType == 'All', gene_name %in% genePanel$gene_name)

load(file = "~/hmf/analysis/dnds/HmfRefCDSCv.RData")
New = HmfRefCDSCv %>% filter(gene_name %in% genePanel$gene_name)
New$change = round((New$wmis_cv - Old$wmis_cv) * New$n_mis, 2)
save(New, Old, file = "~/hmf/Rdata/Processed/ImpactOfNewData.RData")

merged=(merge(Old,New,by='gene_name',suffixes=c('.old','.new')))
View(merged %>% mutate(chMisDrivers=n_mis.new*(pmax(0,(wmis_cv.new-1)/wmis_cv.new)-pmax(0,(wmis_cv.old-1)/wmis_cv.old)),
                       chNonDrivers=n_non.new*(pmax(0,(wnon_cv.new-1)/wnon_cv.new)-pmax((wnon_cv.old-1)/wnon_cv.old,0))) %>% 
       select(gene_name,chMisDrivers,wmis_cv.old,wmis_cv.new,chNonDrivers,everything()))

#### KRAS

load(file = "~/hmf/analysis/dnds/exonicSomaticsDndsAnnotated.RData")
dndsUnfilteredAnnotatedMutations %>% filter(gene == 'KRAS', impact == 'Synonymous')

load(file = "~/hmf/RData/processed/hpcExonicSomaticsDndsAnnotated.RData")
dndsUnfilteredAnnotatedMutations %>% filter(gene == 'KRAS', impact == 'Synonymous')



load(file = "~/hmf/RData/Processed/hpcDndsTsgDriverRates.RData")


load(file = "~/hmf/analysis/dnds/OldDataNewLogicHmfRefCDSCv.RData")
OldDataNewLogic = HmfRefCDSCv

load(file = "~/hmf/Rdata/Processed/HmfRefCDSCv.RData")
Old = HmfRefCDSCv

load(file = "~/hmf/analysis/dnds3/HmfRefCDSCv.RData")
NewDataNewLogic = HmfRefCDSCv

load(file = "~/hmf/analysis/dnds2/HmfRefCDSCv.RData")
Updated = HmfRefCDSCv

Old %>% filter(gene_name == "RB1", cancerType == 'All')
OldDataNewLogic %>% filter(gene_name == "RB1")
NewDataNewLogic %>% filter(gene_name == "RB1")
Updated %>% filter(gene_name == "RB1")


load(file = "~/hmf/RData/Processed/hpcDndsTsgDrivers.RData")

####OLD DATA
load(file = "~/hmf/RData/reference/highestPurityCohort.RData")
load(file = "~/hmf/RData/reference/hpcExonicSomatics.RData")
exonicSomatics = hpcExonicSomatics


load(file = "~/hmf/Rdata/Processed/HmfRefCDSCv.RData")
load(file = "~/hmf/RData/processed/hpcExonicSomaticsDndsAnnotated.RData")
load(file = "~/hmf/RData/processed/driverGenes.RData")



load(file = "~/hmf/RData/processed/highestPurityCohortSummary.RData")
somaticCounts = highestPurityCohortSummary %>% select(sampleId, ends_with("SNV"), ends_with("INDEL"))
somaticCounts = somaticCounts %>%
  mutate(sample_SNV = TOTAL_SNV, sample_INDEL = TOTAL_INDEL) %>%
  select(starts_with("sample"))

