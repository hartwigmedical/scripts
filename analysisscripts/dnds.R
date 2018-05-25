library(GenomicRanges)
library(dplyr)
detach("package:dndscv", unload=TRUE); 
library(dndscv)

create_data_frame <- function(cvList) {
  result = data.frame(stringsAsFactors = F)
  for (cancerType in names(cvList)) {
    df = cvList[[cancerType]]
    df$cancerType <- cancerType
    result = rbind(result, df)  
  }
  return (result)
}

load(file = "~/hmf/RData/reference/highestPurityCohort.RData")
load(file = "~/hmf/RData/reference/hpcExonicSomatics.RData")

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
#dndsFilteredAnnotatedMutations = output$annotmuts
#save(dndsFilteredAnnotatedMutations, file = "~/hmf/RData/processed/dndsFilteredAnnotatedMutations.RData")
#rm(dndsFilteredAnnotatedMutations)

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

HmfRefCDSCv = create_data_frame(HmfRefCDSCvList)
save(HmfRefCDSCv, file="~/hmf/RData/processed/HmfRefCDSCv.RData")


#################### Unfiltered Annotations #########################
unfilteredSomatics = hpcExonicSomatics %>% select(sampleId, chr = chromosome, pos = position, ref, alt)
unfilteredOutput = dndscv(unfilteredSomatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE, max_muts_per_gene_per_sample = 30000000, max_coding_muts_per_sample = 30000000)
dndsUnfilteredAnnotatedMutations = unfilteredOutput$annotmuts
save(dndsUnfilteredAnnotatedMutations, file = "~/hmf/RData/processed/hpcExonicSomaticsDndsAnnotated.RData")
rm(unfilteredSomatics, unfilteredOutput, dndsUnfilteredAnnotatedMutations)

#################### Unfiltered Multiple Biopsy Annotations #########################
load(file = "~/hmf/RData/reference/mbExonicSomatics.RData")
mbUnfilteredSomatics = mbExonicSomatics  %>% select(sampleId, chr = chromosome, pos = position, ref, alt)
mbUnfilteredOutput = dndscv(mbUnfilteredSomatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE, max_muts_per_gene_per_sample = 30000000, max_coding_muts_per_sample = 30000000)
dndsUnfilteredMultipleBiopsyMutations = mbUnfilteredOutput$annotmuts
save(dndsUnfilteredMultipleBiopsyMutations, file = "~/hmf/RData/processed/mbExonicSomaticsDndsAnnotated.RData")

#################### PCAWG dNdS #########################
PcawgRefCDSCv = read.delim("~/hmf/dnds/dNdScv_output_PANCANCER.txt", stringsAsFactors = FALSE)
PcawgRefCDSCv$cancerType <- "All"
for (i in 2:nrow(CancerMap)) {
  cancerType = CancerMap[i, 1]
  file = paste("~/hmf/dnds/dNdScv_output_", CancerMap[i, 2],".txt", sep = "")
  df = read.delim(file, stringsAsFactors = FALSE)
  df$cancerType <- cancerType
  PcawgRefCDSCv = rbind(PcawgRefCDSCv, df)
}
save(PcawgRefCDSCv, file = "~/hmf/RData/output/PcawgRefCDSCv.RData")



#################### PCAWG Null Hypothosis #########################
loadPcawgNullHypothesis<-function(file, RefCDS) {
  NullHypothesis = read.delim(file, stringsAsFactors = FALSE)[, c("gene_name", "wmis3", "wnon3", "wspl3", "wind")]
  
  RefCDSNames = sapply(RefCDS, function(x) {x$gene_name})
  NullHypothesisNames = NullHypothesis$gene_name
  NullHypothesInd = setNames(1:length(NullHypothesisNames), NullHypothesisNames)
  
  NullHypothesis = NullHypothesis[NullHypothesInd[RefCDSNames], ]
  colnames(NullHypothesis) <- c("gene_name", "wmis","wnon","wspl","wind")
  
  return (NullHypothesis)
}


HmfCancerTypes = 
  c("All", "Urinary tract", "CNS", "Breast", "Colon/Rectum", "Esophagus", "Head and neck", "Kidney", "Liver", "Lung", "Skin", "Mesothelioma", "Ovary", "Pancreas",
    "Prostate", "Bone/Soft tissue", "Stomach", "Uterus")
PcawgCancerTypes = 
  c("PANCANCER","BLCA", "GBM","BRCA", "COREAD", "ESCA", "HNSC","KIRC","LIHC","LUAD","SKCM","MESO","OV","PAAD","PRAD","SARC","STAD","UCEC")
CancerMap = data.frame(hmf = HmfCancerTypes, pcawg = PcawgCancerTypes, stringsAsFactors = F)
CancerMap[!CancerMap$hmf %in% cancerTypes, ]

HmfRefCDSCvNullPcawgList = list()
for (i in 1:nrow(CancerMap)) {
  cancerType = CancerMap[i, 1]
  file = paste("~/hmf/dnds/dNdScv_output_", CancerMap[i, 2],".txt", sep = "")
  null_hypothesis = loadPcawgNullHypothesis(file, RefCDS)
  
  cat(file, nrow(null_hypothesis), "\n")
  cat("Processing", cancerType)
  cancerTypeSampleIds = highestPurityCohort[!is.na(highestPurityCohort$primaryTumorLocation) & highestPurityCohort$primaryTumorLocation == cancerType, c("sampleId")]
  if (cancerType == "All") {
    input = somatics
  } else {
    input = somatics[somatics$sampleId %in% cancerTypeSampleIds, c("sampleId", "chr", "pos", "ref", "alt")]
  }
  output = dndscv(input, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE, null_hypothesis = null_hypothesis)
  
  HmfRefCDSCvNullPcawgList[[cancerType]] <- output$sel_cv
}

HmfRefCDSCvNullPcawg = create_data_frame(HmfRefCDSCvNullPcawgList)
save(HmfRefCDSCvNullPcawg, file = "~/hmf/RData/output/HmfRefCDSCvNullPcawg.RData")


#################### HmfRefCDSCv Excluding Cancer Type #########################
HmfRefCDSCvExcludingCancerTypeList = list()
for (cancerType in cancerTypes) {
  cat("Processing", cancerType)
  cancerTypeSampleIds = highestPurityCohort[!is.na(highestPurityCohort$primaryTumorLocation) & highestPurityCohort$primaryTumorLocation != cancerType, c("sampleId")]
  input = somatics[somatics$sampleId %in% cancerTypeSampleIds, c("sampleId", "chr", "pos", "ref", "alt")]
  output = dndscv(input, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)
  
  HmfRefCDSCvExcludingCancerTypeList[[cancerType]] <- output$sel_cv
}

HmfRefCDSCvExcludingCancerType = create_data_frame(HmfRefCDSCvExcludingCancerTypeList)
save(HmfRefCDSCvExcludingCancerType, file = "~/hmf/RData/output/HmfRefCDSCvExcludingCancerType.RData")


#################### Global sans cancertype Null Hypothesis #########################
createNullHypothesisFromSelCV <-function(sel_cv, RefCDS) {
  NullHypothesis = sel_cv[, c("gene_name", "wmis_cv","wnon_cv","wspl_cv","wind_cv")]
  
  RefCDSNames = sapply(RefCDS, function(x) {x$gene_name})
  NullHypothesisNames = NullHypothesis$gene_name
  NullHypothesInd = setNames(1:length(NullHypothesisNames), NullHypothesisNames)
  
  NullHypothesis = NullHypothesis[NullHypothesInd[RefCDSNames], ]
  colnames(NullHypothesis) <- c("gene_name", "wmis","wnon","wspl","wind")
  
  return (NullHypothesis)
}


#load("~/hmf/RData/output/HmfRefCDSCv.RData")
HmfRefCDSCvNullGlobalSansCancerTypeList = list()
HmfRefCDSCvNullGlobalSansCancerType = list()
for (cancerType in cancerTypes) {
  cat("Processing", cancerType)
  null_hypothesis = createNullHypothesisFromSelCV(HmfRefCDSCvExcludingCancerType[HmfRefCDSCvExcludingCancerType$cancerType == cancerType, ], RefCDS)
  cancerTypeSampleIds = highestPurityCohort[!is.na(highestPurityCohort$primaryTumorLocation) & highestPurityCohort$primaryTumorLocation == cancerType, c("sampleId")]
  input = somatics[somatics$sampleId %in% cancerTypeSampleIds, c("sampleId", "chr", "pos", "ref", "alt")]
  output = dndscv(input, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE, null_hypothesis = null_hypothesis)
  
  HmfRefCDSCvNullGlobalSansCancerTypeList[[cancerType]] <- output$sel_cv
}

HmfRefCDSCvNullGlobalSansCancerType = create_data_frame(HmfRefCDSCvNullGlobalSansCancerTypeList)
save(HmfRefCDSCvNullGlobalSansCancerType, file = "~/hmf/RData/output/HmfRefCDSCvNullGlobalSansCancerType.RData")






