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
dndsFilteredAnnotatedMutations = output$annotmuts
save(dndsFilteredAnnotatedMutations, file = "~/hmf/analysis/cohort/processed/dndsFilteredAnnotatedMutations.RData")

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

############################################   Biallelic V Non-Biallelic ############################################
biallelicSomatics = hpcExonicSomatics %>% 
  filter(type == "SNP", repeatCount <= 8, biallelic) %>%
  select(sampleId, chr = chromosome, pos = position, ref, alt)

biallelicOutput = dndscv(biallelicSomatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)
HmfRefCDSCvBiallelic  <- biallelicOutput$sel_cv
save(HmfRefCDSCvBiallelic, file = "~/hmf/analysis/cohort/processed/HmfRefCDSCvBiallelic.RData")

nonBiallelicSomatics = hpcExonicSomatics %>% 
  filter(type == "SNP", repeatCount <= 8, !biallelic) %>%
  select(sampleId, chr = chromosome, pos = position, ref, alt)

nonBiallelicOutput = dndscv(nonBiallelicSomatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)
HmfRefCDSCvNonBiallelic  <- nonBiallelicOutput$sel_cv
save(HmfRefCDSCvNonBiallelic, file = "~/hmf/analysis/cohort/processed/HmfRefCDSCvNonBiallelic.RData")

############################################   Choose genes to use separate figures for Biallelic / NonBiallelic ############################################
load(file = "~/hmf/analysis/cohort/processed/genePanel.RData")
load(file = "~/hmf/analysis/cohort/processed/HmfRefCDSCv.RData")
HmfRefCDSCv = HmfRefCDSCv %>% filter(cancerType == 'All')

load(file = "~/hmf/analysis/cohort/processed/HmfRefCDSCvBiallelic.RData")
load(file = "~/hmf/analysis/cohort/processed/HmfRefCDSCvNonBiallelic.RData")

tsGenes = genePanel %>% filter(reportablePointMutation == "tsg") %>% select(gene, reportablePointMutation)
allCv = HmfRefCDSCv %>% filter(cancerType == 'All', gene_name %in% tsGenes$gene) %>% select(gene = gene_name, starts_with("n_"), starts_with("w"))
biallelicCv = HmfRefCDSCvBiallelic %>% filter(gene_name %in% tsGenes$gene) %>% select(gene = gene_name, starts_with("n_"), starts_with("w")) 
nonbiallelicCv = HmfRefCDSCvNonBiallelic %>% filter(gene_name %in% tsGenes$gene) %>% select(gene = gene_name, starts_with("n_"), starts_with("w"))

tsGenesToSeparateBiallelic = left_join(biallelicCv, nonbiallelicCv, by = "gene", suffix = c(".bi", ".nonbi")) %>% 
  left_join(allCv, by = "gene") %>% 
  right_join(tsGenes, by = "gene") %>%
  filter(wmis_cv.bi > wmis_cv.nonbi)

save(tsGenesToSeparateBiallelic, file = "~/hmf/analysis/cohort/processed/tsGenesToSeparateBiallelic.RData")
