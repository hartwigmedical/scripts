library(dplyr)
library(tibble)
library(RMySQL)
library(ggplot2)

rm(list=ls())

# RETRIEVE DATA ------------------------------------------------------------------
dbActin <- dbConnect(MySQL(), dbname='actin_paper', groups="RAnalysis")

query_sample <-"select * from paperSamples;"
query_patient <-"select * from patient;"
query_clinical_status <-"select * from clinicalStatus;"
query_tumor <-"select * from tumor;"
query_eligibleCohorts <-"select * from eligibleCohorts union select * from eligibleCohorts_addition;"
query_molecular <-"select * from molecular;"
query_drivers <- "select * from molecularDrivers;"
query_drivers_gene_count <- "select gene, count(distinct sampleId) as sampleCount from molecularDrivers where driverlikelihood='High' and sampleId in (select sampleId from molecular where containsTumorCells=1) group by 1 order by 2 desc;"
query_patients_second_eval_phase <- "select * from patientsEvaluationPhase;"
query_evaluated_trials_second_eval_phase <- "select distinct patientId, code from trialMatch t LEFT JOIN treatmentMatch tm ON tm.id = t.treatmentMatchId;"

sample <- dbGetQuery(dbActin, query_sample)
patient <- dbGetQuery(dbActin, query_patient)
clinicalStatus <- dbGetQuery(dbActin, query_clinical_status)
tumor <- dbGetQuery(dbActin, query_tumor)
eligibleCohorts <- dbGetQuery(dbActin, query_eligibleCohorts)
molecular <- dbGetQuery(dbActin, query_molecular)
drivers <- dbGetQuery(dbActin, query_drivers)
driversGeneCount <- dbGetQuery(dbActin, query_drivers_gene_count)
patientsSecondEvalPhase <- dbGetQuery(dbActin, query_patients_second_eval_phase)
evaluatedTrialsSecondEvalPhase <- dbGetQuery(dbActin, query_evaluated_trials_second_eval_phase)

dbDisconnect(dbActin)

actionableEvents <- read.csv(paste0(Sys.getenv("HOME"), "/hmf/tmp/TAT and other results overview - WGS relevant drivers+characteristics.csv"), sep = ",")

## CLINICAL ------------------------------------------------------------------
# Tumor type
tumorType <- tumor %>% 
  mutate(primaryTumorLocation = ifelse(primaryTumorType == "Neuroendocrine tumor", "Neuroendocrine (NET)", primaryTumorLocation)) %>%
  group_by(primaryTumorLocation) %>% 
  summarise(count=n()) %>%
  mutate(primaryTumorLocation = ifelse(count == 1 | count == 2, "Other (n<=2)", primaryTumorLocation)) %>%
  group_by(primaryTumorLocation) %>%
  summarize_all(sum) %>% 
  arrange(count)

tumorType %>% filter(is.na(tumorType$primaryTumorLocation)) %>% nrow()
row_to_move_1 <- tumorType[tumorType$primaryTumorLocation=="Other (n<=2)", ]
row_to_move_2 <- tumorType[!tumorType$primaryTumorLocation == "Other (n<=2)", ]

tumorTypeRearranged <- rbind(row_to_move_1, row_to_move_2)
tumorTypeRearranged$primaryTumorLocation <- factor(tumorTypeRearranged$primaryTumorLocation, levels = tumorTypeRearranged$primaryTumorLocation)
tumorTypeRearranged %>% ggplot(aes(x=count, y=primaryTumorLocation)) + geom_bar(stat="identity") + theme_light()

# Biopsy location
biopsyCuration <- read.csv(paste0(Sys.getenv("HOME"), "/hmf/tmp/Paper Curation - Biopsy location.csv"), sep = ",")
biopsy <- left_join(tumor, biopsyCuration, by=c('biopsyLocation'))
biopsyLocation <- biopsy %>% 
  group_by(biopsyLocationCurated) %>% 
  summarise(count=n()) %>%
  mutate(biopsyLocationCurated = ifelse(count == 1 | count == 2, "Other (n<=2)", biopsyLocationCurated)) %>%
  group_by(biopsyLocationCurated) %>%
  summarize_all(sum) %>%
  arrange(count)

biopsyLocation %>% filter(is.na(biopsyLocation$biopsyLocationCurated)) %>% nrow()
row_to_move_1 <- biopsyLocation[biopsyLocation$biopsyLocationCurated=="Other (n<=2)", ]
row_to_move_2 <- biopsyLocation[!biopsyLocation$biopsyLocationCurated == "Other (n<=2)", ]

biopsyLocationRearranged <- rbind(row_to_move_1, row_to_move_2)
biopsyLocationRearranged$biopsyLocationCurated <- factor(biopsyLocationRearranged$biopsyLocationCurated, levels = biopsyLocationRearranged$biopsyLocationCurated)
biopsyLocationRearranged %>% ggplot(aes(x=count, y=biopsyLocationCurated)) + geom_bar(stat="identity") + theme_light()

# Gender, age & WHO
patientCount <- nrow(patient)
patient %>% filter(is.na(patient$gender)) %>% nrow()
patient %>% group_by(gender) %>% summarise(count=n(), percentage=round(n()/patientCount*100,1))

patient <- patient %>% add_column(registrationYear = substr(patient$registrationDate,1,4))
patient$registrationYear <- as.integer(patient$registrationYear)
patient$birthYear <- as.integer(patient$birthYear)

patient %>% filter(is.na(patient$registrationYear)) %>% nrow()
patient %>% filter(is.na(patient$birthYear)) %>% nrow()
patient <- patient %>% add_column(ageAtRegistration = patient$registrationYear-patient$birthYear)

median(patient$ageAtRegistration)
min(patient$ageAtRegistration)
max(patient$ageAtRegistration)

clinicalStatus %>% filter(is.na(clinicalStatus$who)) %>% nrow()
clinicalStatus %>% group_by(who) %>% summarise(count=n(), percentage=round(n()/patientCount*100,1))

## MOLECULAR ------------------------------------------------------------------
# Successful WGS reports
molecularSuccessful <- molecular %>% filter(containsTumorCells==1)
n_successful_reports <- molecularSuccessful %>% nrow()

# Filter for relevant drivers
driversFiltered <- drivers %>% filter(category=='c_fusions' | driverLikelihood=='High') %>% filter(sampleId %in% molecularSuccessful$sampleId)

# Drivers per successful sample
driversPerSample <- driversFiltered %>% group_by(sampleId) %>% summarise(count = n())
driversPerSampleOverview <- molecularSuccessful %>% left_join(driversPerSample, by='sampleId') %>% mutate(count = ifelse(is.na(count), 0, count)) %>% select(sampleId, count)

median(driversPerSampleOverview$count)
min(driversPerSampleOverview$count)
max(driversPerSampleOverview$count)
driversPerSampleOverview %>% filter(count>0) %>% nrow()

# High TML, TMB, HRD and MSI
molecularSuccessful %>% filter(isMicrosatelliteUnstable==1) %>% nrow()
molecularSuccessful %>% filter(isHomologousRepairDeficient==1) %>% nrow()
molecularSuccessful %>% filter(tumorMutationalBurden>=10) %>% nrow()
molecularSuccessful %>% filter(tumorMutationalLoad>=140) %>% nrow()

# Recurring genes in driver events (excluding fusions and viruses)
head(driversGeneCount %>% add_column(percentage = round(driversGeneCount$sampleCount/n_successful_reports*100,1)),5)

# High-driver likelihood WGS aberrations per gene (excluding fusions and viruses) + plot
driversFilteredPerGene <- driversFiltered %>%
  group_by(gene) %>%
  summarise(count=n_distinct(sampleId)) %>%
  mutate(gene = ifelse(count <= 5, "Other (n<=5)", gene)) %>%
  group_by(gene) %>%
  summarize_all(sum) %>% 
  arrange(count)

rows_to_keep <- driversFilteredPerGene[!is.na(driversFilteredPerGene$gene) & !driversFilteredPerGene$gene=="Other (n<=5)" & !driversFilteredPerGene$gene =="NA", ]

driversFilteredPerGene <- rbind(rows_to_keep)
driversFilteredPerGene$gene <- factor(driversFilteredPerGene$gene, levels = driversFilteredPerGene$gene)
driversFilteredPerGene %>% ggplot(aes(x=count, y=gene)) + geom_bar(stat="identity") + theme_light()

# Fusions and viruses
driversFiltered %>% 
  filter(grepl('fusion', category)) %>% nrow()

driversFiltered %>% 
  filter(grepl('fusion', category)) %>% select(sampleId) %>% n_distinct()

driversFiltered %>% 
  filter(grepl('fusion', category) & driverLikelihood == 'High') %>% nrow()

driversFiltered %>% 
  filter(grepl('virus', category)) %>% nrow()

driversFiltered %>% 
  filter(grepl('virus', category)) %>% select(sampleId) %>% n_distinct()

# Actionable events plot
actionableEventsFiltered <- actionableEvents %>% filter(sampleId %in% molecularSuccessful$sampleId)

actionableEventsFiltered %>% 
  filter(potentially.actionable. == 'Yes') %>% select(sampleId) %>% n_distinct()

actionableEventsFilteredPerEvent <- actionableEventsFiltered %>%
  filter(potentially.actionable. == 'Yes') %>%
  group_by(curated.event) %>%
  summarise(count=n_distinct(sampleId)) %>%
  mutate(curated.event = ifelse(count <= 5, "Other (n<=5)", curated.event)) %>%
  group_by(curated.event) %>%
  summarize_all(sum) %>% 
  arrange(count)

actionableEventsFilteredPerEvent %>% tail(6)

rows_to_keep <- actionableEventsFilteredPerEvent[!is.na(actionableEventsFilteredPerEvent$curated.event) & !actionableEventsFilteredPerEvent$curated.event=="Other (n<=5)", ]

actionableEventsFilteredPerEventRearranged <- rbind(rows_to_keep)
actionableEventsFilteredPerEventRearranged$curated.event <- factor(actionableEventsFilteredPerEventRearranged$curated.event, levels = actionableEventsFilteredPerEventRearranged$curated.event)
actionableEventsFilteredPerEventRearranged %>% ggplot(aes(x=count, y=curated.event)) + geom_bar(stat="identity") + theme_light()

## COHORTS & TRIALS NR ------------------------------------------------------------------
eligibleCohortsJoin <- left_join(patientsSecondEvalPhase, eligibleCohorts, by="patientId") %>% group_by(patientId) %>% summarise(numberOfCohorts=sum(!is.na(trialId))) 
median(eligibleCohortsJoin$numberOfCohorts)
min(eligibleCohortsJoin$numberOfCohorts)
max(eligibleCohortsJoin$numberOfCohorts)

evaluatedTrialsJoinAll <- left_join(patientsSecondEvalPhase, evaluatedTrialsSecondEvalPhase, by="patientId") %>% group_by(patientId) %>% summarise(numberOfTrials=sum(!is.na(code)))
evaluatedTrialsJoin <- evaluatedTrialsJoinAll %>% filter(numberOfTrials != 0)
median(evaluatedTrialsJoin$numberOfTrials)
min(evaluatedTrialsJoin$numberOfTrials)
max(evaluatedTrialsJoin$numberOfTrials)
