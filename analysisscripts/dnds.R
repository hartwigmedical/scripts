library(RMySQL)
library(dndscv)
library(IRanges)
detach("package:purple", unload=TRUE)
library(purple)


load(file="~/hmf/cohort.RData")
distinctCohort = purple::highest_purity_patients(cohort)

#Select snps
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
snps = query_snps(dbProd, distinctCohort[1:3,])
dbDisconnect(dbProd)
rm(dbProd)

#Enrich snps with cancer type
snps = snps[,c("sampleId", "chromosome", "position", "ref", "alt")]
colnames(snps) <- c("sampleId", "chr", "pos", "ref", "alt")
snps$cancerType = slookup(jon, distinctCohort[, c("sampleId", "cancerType")])


panResults = dndscv(snps[, c("sampleId", "chr", "pos", "ref", "alt")])
results = list()

cancerTypes = unique(snps$cancerType)
for (cancerType in cancerTypes) {
  cat("Processing ", cancerType)
  input = snps[snps$cancerType == cancerType, c("sampleId", "chr", "pos", "ref", "alt")]
  output = dndscv(input)
  results[[cancerType]] <- output
}

# Example of writing globadnds to file...
for (cancerType in cancerTypes) {
  write.csv(results[[cancerType]]$globaldnds, paste('~/hmf/dnds_', cancerType, ".RData", sep=""))
}

# Accessing individual ones
breastResults = results$Breast
lungResults = results$Lung
