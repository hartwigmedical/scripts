library(RMySQL)

library(IRanges)
detach("package:purple", unload=TRUE); 
library(purple);
library(Biostrings)
library(GenomicRanges)
library(MASS)
library(seqinr)
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


createNullHypothesisFromSelCV <-function(sel_cv, RefCDS) {
  NullHypothesis = sel_cv[, c("gene_name", "wmis_cv","wnon_cv","wspl_cv","wind_cv")]
  
  RefCDSNames = sapply(RefCDS, function(x) {x$gene_name})
  NullHypothesisNames = NullHypothesis$gene_name
  NullHypothesInd = setNames(1:length(NullHypothesisNames), NullHypothesisNames)
  
  NullHypothesis = NullHypothesis[NullHypothesInd[RefCDSNames], ]
  colnames(NullHypothesis) <- c("gene_name", "wmis","wnon","wspl","wind")
  
  return (NullHypothesis)
}

loadPcawgNullHypothesis<-function(file, RefCDS) {
  NullHypothesis = read.delim(file, stringsAsFactors = FALSE)[, c("gene_name", "wmis3", "wnon3", "wspl3", "wind")]
  
  RefCDSNames = sapply(RefCDS, function(x) {x$gene_name})
  NullHypothesisNames = NullHypothesis$gene_name
  NullHypothesInd = setNames(1:length(NullHypothesisNames), NullHypothesisNames)
  
  NullHypothesis = NullHypothesis[NullHypothesInd[RefCDSNames], ]
  colnames(NullHypothesis) <- c("gene_name", "wmis","wnon","wspl","wind")
  
  return (NullHypothesis)
}

load(file = "~/hmf/RData/highestPurityCohort.RData")
load(file = "~/hmf/RData/highestPurityExonicSomatics.RData")



somatics = highestPurityExonicSomatics %>% 
  filter(type == "INDEL" | type == "SNP") %>%
  select(sampleId, chr = chromosome, pos = position, ref = ref, alt = alt)
colnames(somatics) <- c("sampleId", "chr", "pos", "ref", "alt")

cancerTypes = unique(highestPurityCohort$primaryTumorLocation)
cancerTypes = cancerTypes[!is.na(cancerTypes)]
cancerTypes

data("cancergenes_cgc81", package="dndscv")
data("covariates_hg19", package="dndscv")
refdb="~/hmf/RData/HmfRefCDS.RData"
load(refdb)
newgenes = sapply(RefCDS, function(x) x$gene_name) # New list of genes
oldgenes = sapply(strsplit(newgenes,split=":"), function(x) x[1])
cv = covs[oldgenes,]
rownames(cv) = newgenes
kc = newgenes[oldgenes %in% known_cancergenes]

#RefCDSNames = sapply(RefCDS, function (x) {x$gene_name})
#RefCDSInd = setNames(1:length(RefCDSNames), RefCDSNames)
#sum(RefCDS[[RefCDSInd["KRAS"]]]$L[1:192,3])

HmfRefCDSGlobalList = list()
HmfRefCDSCvList = list()
for (cancerType in cancerTypes) {
  cat("Processing", cancerType)
  cancerTypeSampleIds = highestPurityCohort[!is.na(highestPurityCohort$primaryTumorLocation) & highestPurityCohort$primaryTumorLocation == cancerType, c("sampleId")]
  input = somatics[somatics$sampleId %in% cancerTypeSampleIds, c("sampleId", "chr", "pos", "ref", "alt")]
  output = dndscv(input, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)
  
  HmfRefCDSCvList[[cancerType]] <- output$sel_cv
  HmfRefCDSGlobalList[[cancerType]] <- output$globaldnds
  #save(HmfRefCDSCvList, file="~/hmf/RData/HmfRefCDSCvList")
}

output = dndscv(somatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)
names(output)
globaldnds = output$globaldnds

HmfRefCDSGlobalList[["All"]] <- output$globaldnds
HmfRefCDSCvList[["All"]]  <- output$sel_cv

HmfRefCDSCv = create_data_frame(HmfRefCDSCvList)
save(HmfRefCDSCv, file="~/hmf/RData/HmfRefCDSCv.RData")











# Create HmfAllNullHypothesis
NullHypothesisAll = createNullHypothesisFromSelCV(HmfRefCDSCvList[["All"]], RefCDS)


##### Create Cancer Map
HmfCancerTypes = 
  c("All", "Bladder", "Brain", "Breast", "Colorectal", "Esophagus", "Head and neck", "Kidney", "Liver", "Lung", "Melanoma", "Mesothelioma",
    "Neuroendocrine", "Ovary", "Pancreas", "Prostate", "Sarcoma", "Stomach", "Testis", "Uterus")
PcawgCancerTypes = 
  c("PANCANCER","BLCA", "GBM","BRCA", "COREAD", "ESCA", "HNSC","KIRC","LIHC","LUAD","SKCM","MESO","THCA","OV","PAAD","PRAD","SARC","STAD","TGCT","UCEC")
CancerMap = data.frame(hmf = HmfCancerTypes, pcawg = PcawgCancerTypes, stringsAsFactors = F)
CancerMap[CancerMap$hmf %in% cancerTypes, ]
CancerMap[1, 2]

##### PCAWG NULL HYP
HmfRefCDSCvNullPcawgList = list()
for (i in 2:nrow(CancerMap)) {

  cancerType = CancerMap[i, 1]
  file = paste("~/hmf/dnds/dNdScv_output_", CancerMap[i, 2],".txt", sep = "")
  null_hypothesis = loadPcawgNullHypothesis(file, RefCDS)
  
  cat(file, nrow(null_hypothesis), "\n")
  cat("Processing", cancerType)
  cancerTypeSampleIds = highestPurityCohort[!is.na(highestPurityCohort$cancerType) & highestPurityCohort$cancerType == cancerType, c("sampleId")]
  input = somatics[somatics$sampleId %in% cancerTypeSampleIds, c("sampleId", "chr", "pos", "ref", "alt")]
  output = jondndscv(input, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE, null_hypothesis = null_hypothesis)
  
  HmfRefCDSCvNullPcawgList[[cancerType]] <- output$sel_cv
}

i = 1
file = paste("~/hmf/dnds/dNdScv_output_", CancerMap[i, 2],".txt", sep = "")
null_hypothesis = loadPcawgNullHypothesis(file, RefCDS)
output = jondndscv(somatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE, null_hypothesis = null_hypothesis)
HmfRefCDSCvNullPcawgList[["All"]] <- output$sel_cv

HmfRefCDSCvNullPcawg = HmfRefCDSCvNullPcawgList[["All"]]
HmfRefCDSCvNullPcawg$cancerType <- "All"
for (cancerType in CancerMap[2:nrow(CancerMap), c("hmf")]) {
  df = HmfRefCDSCvNullPcawgList[[cancerType]]
  df$cancerType <- cancerType
  HmfRefCDSCvNullPcawg = rbind(HmfRefCDSCvNullPcawg, df)
}
save(HmfRefCDSCvNullPcawg, file = "~/hmf/RData/HmfRefCDSCvNullPcawg.RData")

##### PCAWG DNDS
PcawgRefCDSCv = read.delim("~/hmf/dnds/dNdScv_output_PANCANCER.txt", stringsAsFactors = FALSE)
PcawgRefCDSCv$cancerType <- "All"
for (i in 2:nrow(CancerMap)) {
  cancerType = CancerMap[i, 1]
  file = paste("~/hmf/dnds/dNdScv_output_", CancerMap[i, 2],".txt", sep = "")
  df = read.delim(file, stringsAsFactors = FALSE)
  df$cancerType <- cancerType
  PcawgRefCDSCv = rbind(PcawgRefCDSCv, df)
}
save(PcawgRefCDSCv, file = "~/hmf/RData/PcawgRefCDSCv.RData")

#################### dNdS EXCLUDING Cancer Types #########################

HmfRefCDSCvExcludingCancerTypeList = list()
for (cancerType in cancerTypes) {
  cat("Processing", cancerType)
  cancerTypeSampleIds = highestPurityCohort[!is.na(highestPurityCohort$primaryTumorLocation) & highestPurityCohort$primaryTumorLocation != cancerType, c("sampleId")]
  input = somatics[somatics$sampleId %in% cancerTypeSampleIds, c("sampleId", "chr", "pos", "ref", "alt")]
  output = jondndscv(input, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)
  
  HmfRefCDSCvExcludingCancerTypeList[[cancerType]] <- output$sel_cv
}

output = jondndscv(somatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)

HmfRefCDSExcludingCancerTypeCV = output$sel_cv
HmfRefCDSExcludingCancerTypeCV$cancerType <- "All"
for (cancerType in cancerTypes[!is.na(cancerTypes)]) {
  df = HmfRefCDSCvExcludingCancerTypeList[[cancerType]]
  df$cancerType <- cancerType
  HmfRefCDSExcludingCancerTypeCV = rbind(HmfRefCDSExcludingCancerTypeCV, df)
}
save(HmfRefCDSExcludingCancerTypeCV, file = "~/hmf/RData/HmfRefCDSExcludingCancerTypeCV.RData")
rm(HmfRefCDSCvExcludingCancerTypeList)

#################### dNdS Global SANS cancertype  #########################
load("~/hmf/RData/HmfRefCDSCv.RData")
HmfRefCDSCvNullGlobalSansCancerType = list()
for (cancerType in cancerTypes[!is.na(cancerTypes)]) {
  cat("Processing", cancerType)
  null_hypothesis = createNullHypothesisFromSelCV(HmfRefCDSExcludingCancerTypeCV[HmfRefCDSExcludingCancerTypeCV$cancerType == cancerType, ], RefCDS)
  cancerTypeSampleIds = highestPurityCohort[!is.na(highestPurityCohort$cancerType) & highestPurityCohort$cancerType == cancerType, c("sampleId")]
  input = somatics[somatics$sampleId %in% cancerTypeSampleIds, c("sampleId", "chr", "pos", "ref", "alt")]
  output = jondndscv(input, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE, null_hypothesis = null_hypothesis)
  
  HmfRefCDSCvNullGlobalSansCancerType[[cancerType]] <- output$sel_cv
}

null_hypothesis = createNullHypothesisFromSelCV(HmfRefCDSExcludingCancerTypeCV[HmfRefCDSExcludingCancerTypeCV$cancerType == "All", ], RefCDS)
output = jondndscv(somatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE, null_hypothesis = null_hypothesis)


HmfRefCDSNullGlobalSansCancerTypeCV = output$sel_cv
HmfRefCDSNullGlobalSansCancerTypeCV$cancerType <- "All"
for (cancerType in cancerTypes[!is.na(cancerTypes)]) {
  df = HmfRefCDSCvNullGlobalSansCancerType[[cancerType]]
  df$cancerType <- cancerType
  HmfRefCDSNullGlobalSansCancerTypeCV = rbind(HmfRefCDSNullGlobalSansCancerTypeCV, df)
}
save(HmfRefCDSNullGlobalSansCancerTypeCV, file = "~/hmf/RData/HmfRefCDSNullGlobalSansCancerTypeCV.RData")

HmfRefCDSNullGlobalSansCancerTypeCV[HmfRefCDSNullGlobalSansCancerTypeCV$qglobal_cv < 0.1, ]

jon = HmfRefCDSNullGlobalSansCancerTypeCV %>% filter(qglobal_cv < 0.1)

#################### MUCK AROUND #########################
library(ggplot2)
load("~/hmf/RData/HmfRefCDSCv.RData")
load("~/hmf/RData/HmfRefCDSCvNullGlobal.RData")
load("~/hmf/RData/HmfRefCDSCvNullPcawg.RData")
load("~/hmf/RData/PcawgRefCDSCv.RData")

#Attach pcawg score to hmf
PcawgScore = PcawgRefCDSCv[, c("gene_name", "cancerType", "qglobal")]
colnames(PcawgScore) <- c("gene_name", "cancerType", "pcawg_qglobal_cv")

jon = left_join(HmfRefCDSCv, PcawgScore, by = c("gene_name", "cancerType"))
jon$status = "Unchanged"
jon$status = ifelse(jon$qglobal_cv < 0.1 & jon$pcawg_qglobal_cv > 0.1, "Added", jon$status)
jon$status = ifelse(jon$qglobal_cv > 0.1 & jon$pcawg_qglobal_cv < 0.1, "Removed", jon$status)

significant = jon %>% filter(qglobal_cv < 0.1 | pcawg_qglobal_cv < 0.1) %>%select(gene_name, cancerType, status) %>% spread(cancerType, status)
significant2 = jon %>% filter(qglobal_cv < 0.001 | pcawg_qglobal_cv < 0.001)
ggplot(significant2, aes(x = cancerType, y = gene_name))+ geom_tile(aes(fill=status)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
