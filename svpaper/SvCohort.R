detach("package:purple", unload=TRUE)
library(purple)
library(tidyr)
library(dplyr)
library(RMySQL)

load(file = "~/hmf/analysis/genepanel/cohort.RData")

save(highestPurityCohort, file = "~/hmf/analysis/svPaper/highestPurityCohort.RData")
save(cohort, highestPurityCohort, multipleBiopsyMapping, file = "~/hmf/analysis/svPaper/cohort.RData")

load(file = "~/hmf/analysis/svPaper/highestPurityCohort.RData")
