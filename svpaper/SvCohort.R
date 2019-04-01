detach("package:purple", unload=TRUE)
library(purple)
library(tidyr)
library(dplyr)
library(RMySQL)

####### DB
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis", host = "127.0.0.1")
dbDisconnect(dbProd)
rm(dbProd)

####### Cohort
entireCohort = purple::query_cohort(dbProd)
highestPurityCohort = purple::highest_purity_cohort(entireCohort)
multipleBiopsyCohort = multiple_biopsy_cohort(entireCohort)
save(entireCohort, highestPurityCohort, multipleBiopsyCohort, file = "~/hmf/analysis/svPaper/cohort.RData")

####### Hpc Copy Numbers
load(file = "~/hmf/analysis/svPaper/cohort.RData")
hpcCopyNumbers = purple::query_copy_number(dbProd, highestPurityCohort)
save(hpcCopyNumbers,file = "~/hmf/analysis/svPaper/hpcCopyNumbers.RData")
