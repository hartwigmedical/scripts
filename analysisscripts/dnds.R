library(RMySQL)
library(dndscv)
library(IRanges)
detach("package:purple", unload=TRUE); 
library(purple);
library(Biostrings)
library(GenomicRanges)
library(MASS)
library(seqinr)

createNullHypothesisFromSelCV <-function(sel_cv, RefCDS) {
  NullHypothesis = sel_cv[, c("gene_name", "wmis_cv","wnon_cv","wspl_cv","wind_cv")]
  
  RefCDSNames = sapply(RefCDS, function(x) {x$gene_name})
  NullHypothesisNames = NullHypothesis$gene_name
  NullHypothesInd = setNames(1:length(NullHypothesisNames), NullHypothesisNames)
  
  NullHypothesis = NullHypothesis[NullHypothesInd[RefCDSNames], ]
  colnames(NullHypothesis) <- c("gene_name", "wmis","wnon","wspl","wind")
  
  return (NullHypothesis)
}


dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
cat("Querying purple")
rawCohort = purple::query_purity(dbProd)
cohort = rawCohort

# PatientIds
cat("Mapping samples to patients")
patientIdLookups = query_patient_id_lookup(dbProd)
patientIds = purple::apply_to_cohort(cohort, function(x) {purple::sample_to_patient_id(x$sampleId, patientIdLookups)})
cohort$patientId <- patientIds$V1

#Clinical Data
cat("Querying clinical data")
clinicalData = purple::query_clinical_data(dbProd)
cohort = left_join(cohort, clinicalData[, c("sampleId", "cancerType")])

# Cohort
highestPurityCohort = purple::highest_purity_patients(cohort)
save(highestPurityCohort, file="~/hmf/RData/highestPurityCohort.RData")
rm(cohort)
rm(rawCohort)

# Somatics 
#cat("Querying somatics")
#rawSomatics = purple::query_somatic_variants(dbProd, cohort, filterEmptyGenes = TRUE)
#save(rawSomatics, file="/Users/jon/hmf/RData/dnds/dndsSomatics.RData")
#somatics = rawSomatics[rawSomatics$type == "INDEL" | rawSomatics$type == "SNP", c("sampleId", "chromosome", "position", "ref", "alt")]
#colnames(somatics) <- c("sampleId", "chr", "pos", "ref", "alt")

load("~/hmf/RData/highestPurityCohort.RData")
load("~/hmf/RData/allHighestPuritySomaticsProd.RData")
somatics = highestPuritySomaticsProd[highestPuritySomaticsProd$type == "INDEL" | highestPuritySomaticsProd$type == "SNP", c("sampleId", "chromosome", "position", "ref", "alt")]
rm(highestPuritySomaticsProd)
colnames(somatics) <- c("sampleId", "chr", "pos", "ref", "alt")

# Takes AGES to annotate with cancer type... don't bother!
#somatics$cancerType <- "UNKNOWN"
#for (sampleId in unique(somatics$sampleId)) {
#  matchedCancerType = cohort[match(sampleId, cohort$sampleId), c("cancerType")]
#  cat("SampleId:", sampleId, ", cancerType:", matchedCancerType, "\n")
#  somatics[somatics$sampleId == sampleId, ]$cancerType <- matchedCancerType
#}

# Clean up DB Connection
dbDisconnect(dbProd)
rm(dbProd)
rm(patientIdLookups)
rm(patientIds)
rm(clinicalData)

cancerTypes = unique(highestPurityCohort$cancerType)
cancerTypes = cancerTypes[!is.na(cancerTypes)]
cancerTypes = cancerTypes[cancerTypes != "Unknown primary"]
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

RefCDSNames = sapply(RefCDS, function (x) {x$gene_name})
RefCDSInd = setNames(1:length(RefCDSNames), RefCDSNames)
sum(RefCDS[[RefCDSInd["KRAS"]]]$L[1:192,3])

HmfRefCDSCvList = list()
for (cancerType in cancerTypes[!is.na(cancerTypes)]) {
  cat("Processing", cancerType)
  cancerTypeSampleIds = highestPurityCohort[!is.na(highestPurityCohort$cancerType) & highestPurityCohort$cancerType == cancerType, c("sampleId")]
  input = somatics[somatics$sampleId %in% cancerTypeSampleIds, c("sampleId", "chr", "pos", "ref", "alt")]
  output = jondndscv(input, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)
  
  HmfRefCDSCvList[[cancerType]] <- output$sel_cv
  #save(HmfRefCDSCvList, file="~/hmf/RData/HmfRefCDSCvList")
}

output = jondndscv(somatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)
HmfRefCDSCvList[["All"]]  <- output$sel_cv
save(HmfRefCDSCvList, file="~/hmf/RData/HmfRefCDSCvList.Rdat")

#Combine together into one big happy data frame
HmfRefCDSCv = HmfRefCDSCvList[["All"]]
HmfRefCDSCv$cancerType <- "All"
for (cancerType in cancerTypes[!is.na(cancerTypes)]) {
  df = HmfRefCDSCvList[[cancerType]]
  df$cancerType <- cancerType
  HmfRefCDSCv = rbind(HmfRefCDSCv, df)
}
save(HmfRefCDSCv, file="~/hmf/RData/HmfRefCDSCv.RData")

# Create HmfAllNullHypothesis
NullHypothesisAll = createNullHypothesisFromSelCV(HmfRefCDSCvList[["All"]], RefCDS)




#################### MUCK AROUND #########################
library(dplyr)
library(tidyr)
SignificanceMatrix = HmfRefCDSCv %>% filter(qglobal_cv < 0.1) %>% select(gene_name, cancerType, qglobal_cv) %>% spread(cancerType, qglobal_cv)

jon2 = dcast(HmfRefCDSCv[HmfRefCDSCv$qglobal_cv < 0.1, ], gene_name ~ cancerType, value.var = "qglobal_cv")


library(data.table)
jon = HmfRefCDSCv[HmfRefCDSCv$qglobal_cv < 0.1, ]
jon = dcast(jon, gene_name ~ cancerType, value.var = "qglobal_cv")




loadPcawgNullHypothesis<-function(file = "~/hmf/dnds/dNdScv_output_PANCANCER.txt", RefCDS) {
  NullHypothesis = read.delim(file, stringsAsFactors = FALSE)[, c("gene_name", "wmis3", "wnon3", "wspl3", "wind")]
  
  RefCDSNames = sapply(RefCDS, function(x) {x$gene_name})
  NullHypothesisNames = NullHypothesis$gene_name
  NullHypothesInd = setNames(1:length(NullHypothesisNames), NullHypothesisNames)
  
  NullHypothesis = NullHypothesis[NullHypothesInd[RefCDSNames], ]
  colnames(NullHypothesis) <- c("gene_name", "wmis","wnon","wspl","wind")
  
  return (NullHypothesis)
}
