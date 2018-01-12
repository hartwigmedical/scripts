library(RMySQL)
detach("package:purple", unload=TRUE)
library(purple)

load(file="~/hmf/dndsSelCV.RData")
load(file="~/hmf/cohort.RData")
distinctCohort = purple::highest_purity_patients(cohort)

#Select gene panel - NOTE THIS COMES FROM PILOT
dbPilot = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
genePanel = purple::query_gene_panel(dbPilot)
dbDisconnect(dbPilot)
rm(dbPilot)

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
somaticDrivers = purple::query_somatic_drivers(dbProd, distinctCohort[1:3, ], genePanel)
cancerTypes = unique(distinctCohort$cancerType)

#Enrich snps with cancer type
somaticDrivers$cancerType = slookup2(somaticDrivers$sampleId, distinctCohort[, c("sampleId", "cancerType")])
cancerTypes = unique(somaticDrivers$cancerType)

# Enrich with pan cancer dndsSelCV
somaticDrivers$pan_wmis_cv = slookup2(somaticDrivers$gene, dndsSelCV$pan[, c("gene_name", "wmis_cv")], on = "gene_name")

# Enriche with cancer type dnds data
somaticDrivers$wmis_cv <- NA
for (cancerType in cancerTypes) {
  dndsData = dndsSelCV[[cancerType]]
  lookupCancerType = slookup2(somaticDrivers$gene, dndsData[, c("gene_name", "wmis_cv")], on = "gene_name")
  somaticDrivers$wmis_cv = ifelse(somaticDrivers$cancerType == cancerType, lookupCancerType, somaticDrivers$wmis_cv)
}

