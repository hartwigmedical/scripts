library(RMySQL)
library(dndscv)
library(IRanges)
detach("package:purple", unload=TRUE)
library(purple)


load(file="~/hmf/cohort.RData")
distinctCohort = purple::highest_purity_patients(cohort)

#Select snps
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
snps = purple::query_snps_cohort(dbProd, distinctCohort[1:3,])
### ALTERNATIVE WAY OF DOING IT ONE BY ONE: snps = purple::apply_to_cohort(distinctCohort[1:3,], function(x) {purple::query_snps_sample(dbProd, x$sampleId)})
dbDisconnect(dbProd)
rm(dbProd)

#Select gene panel - NOTE THIS COMES FROM PILOT
dbPilot = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
genePanel = purple::query_gene_panel(dbPilot)
dbDisconnect(dbPilot)
rm(dbPilot)

#Enrich snps with cancer type
snps$cancerType = slookup(snps, distinctCohort[, c("sampleId", "cancerType")])

dndsSelCV = list()
dndsResults = list()
dndsResults$pan = dndscv(snps[, c("sampleId", "chr", "pos", "ref", "alt")])
dndsSelCV$pan = dndsResults$pan$sel_cv

cancerTypes = unique(snps$cancerType)
for (cancerType in cancerTypes) {
  cat("Processing", cancerType)
  input = snps[snps$cancerType == cancerType, c("sampleId", "chr", "pos", "ref", "alt")]
  output = dndscv(input)
  dndsResults[[cancerType]] <- output
}

# Example of selecting only things we are interested in:
for (cancerType in c("pan", cancerTypes)) {
  output = dndsResults[[cancerType]]   
  sel_cv = output$sel_cv
  dndsSelCV[[cancerType]] = sel_cv[sel_cv$gene_name %in% genePanel$gene, ]
}

# Accessing individual ones
breastResults = dndsResults$Breast
lungResults = dndsResults$Lung

save(dndsSelCV, file="~/hmf/dndsSelCV.RData")

