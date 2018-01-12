detach("package:purple", unload=TRUE)
library(purple)
library(RMySQL)
library(data.table)
LOAD_FROM_FILE = TRUE

# From DB
if (LOAD_FROM_FILE) {
  load("~/hmf/cohortRawData.RData")
} else {
  dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
  cohortRawData=purple::cohort_raw_data(dbProd, limit = 5)
  save(cohortRawData, file="~/hmf/cohortRawData.RData")
  dbDisconnect(dbProd)
  rm(dbProd)
}

cohort = purple::cohort(cohortRawData)
save(cohort, file="~/hmf/cohort.RData")
#write.csv(cohort, file = '~/hmf/cohortOverview.csv')
