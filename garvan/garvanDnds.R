library(GenomicRanges)
library(dplyr)
detach("package:dndscv", unload=TRUE);
library(dndscv)

outputDir = "~/garvan/RData/"
resourceDir = "~/hmf/RData/reference/"

outputDir = "~/Documents/LKCGP_projects/RData/"
referenceDir = paste0(outputDir, "reference/")
processedDir = paste0(outputDir, "processed/")
resourceDir = "/Users/marwo2/Documents/LKCGP_projects/RData/Resources/"

load(paste0(referenceDir, "cohortExonicSomatics.RData"))
load(paste0(processedDir, "highestPurityCohortSummary.RData"))

create_data_frame <- function(cvList) {
  result = data.frame(stringsAsFactors = F)
  for (cancerType in names(cvList)) {
    df = cvList[[cancerType]]
    if (nrow(df) > 0 & ncol(df) == 24) {
      df$cancerType <- cancerType
      result = rbind(result, df)
    }
  }
  return (result)
}

somatics = cohortExonicSomatics %>%
  filter(type == "INDEL" | type == "SNP", repeatCount <= 8) %>%
  select(sampleId, chr = chromosome, pos = position, ref, alt)

data("cancergenes_cgc81", package="dndscv")
data("covariates_hg19", package="dndscv")
refdb=paste0(resourceDir, "HmfRefCDS.RData");
load(refdb)
newgenes = sapply(RefCDS, function(x) x$gene_name) # New list of genes
oldgenes = sapply(strsplit(newgenes,split=":"), function(x) x[1])
cv = covs[oldgenes,]
rownames(cv) = newgenes
kc = newgenes[oldgenes %in% known_cancergenes]


#### Use dNdS to annotate all somatics
unfilteredSomatics = cohortExonicSomatics %>% select(sampleId, chr = chromosome, pos = position, ref, alt)
unfilteredOutput = dndscv(unfilteredSomatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE, max_muts_per_gene_per_sample = 30000000, max_coding_muts_per_sample = 30000000)
dndsUnfilteredAnnotatedMutations = unfilteredOutput$annotmuts
save(dndsUnfilteredAnnotatedMutations, file=paste0(processedDir, "dndsUnfilteredAnnotatedMutations.RData"))

#### Standard dNdS for PANCANCER
output = dndscv(somatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)
HmfRefCDSCvList = list()
HmfRefCDSCvList[["PanCancer"]]  <- output$sel_cv
HmfRefCDSCv = create_data_frame(HmfRefCDSCvList)
save(HmfRefCDSCv, file=paste0(processedDir, "HmfRefCDSCv.RData"))

#### IF YOU WANT TO RUN dNdS for ALL CANCER TYPES YOU CAN RUN THE FOLLOWING AS WELL
cancerTypes = unique(highestPurityCohortSummary$cancerType)
cancerTypes = cancerTypes[!is.na(cancerTypes)]
for (selectedCancerType in cancerTypes) {
  cat("Processing", selectedCancerType)
  cancerTypeSampleIds =  highestPurityCohortSummary %>% filter(!is.na(cancerType), cancerType == selectedCancerType) %>% select(sampleId)
  input = somatics %>% filter(sampleId %in% cancerTypeSampleIds$sampleId)
  output = dndscv(input, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)
  HmfRefCDSCvList[[selectedCancerType]] <- output$sel_cv
}