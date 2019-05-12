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


### Primary Tumor Location
load(file = "~/hmf/RData/reference/allClinicalData.RData")

allClinicalDataDatabase = purple::query_clinical_data(dbProd)

cohortCancerTypes = allClinicalDataDatabase %>% 
  select(sampleId, primaryTumorLocation) %>% 
  filter(!sampleId %in% allClinicalData$sampleId) %>%
  bind_rows(allClinicalData %>% select(sampleId, primaryTumorLocation)) %>%
  left_join(entireCohort %>% select(sampleId, patientId), by = "sampleId") %>%
  filter(!is.na(patientId)) %>% 
  group_by(patientId) %>% 
  mutate(valid = length(unique(primaryTumorLocation)) ==1) 

save(cohortCancerTypes,file = "~/hmf/analysis/svPaper/cohortCancerTypes.RData")


jon = entireCohort %>% filter(qcStatus !=  'PASS')
load(file = "~/hmf/analysis/svPaper/cohort.RData")


load( file = "~/hmf/RData/reference/highestPurityCohort.RData")
highestPurityCohort %>% filter(qcStatus != 'PASS' | purity < 0.20)
highestPurityCohort %>% filter(sampleId == 'CPCT02030224T')
