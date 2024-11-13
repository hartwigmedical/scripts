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
query_drivers <- "select * from molecularDrivers where (driverlikelihood='High' or category='c_fusions');"
query_drivers_per_sample <- "select m.sampleId, count(md.driverLikelihood) as driverCount from molecular m left join molecularDrivers md on m.sampleId=md.sampleId and (driverlikelihood='High' or category='c_fusions') group by 1;"
query_drivers_gene_count <- "select gene, count(distinct sampleId) as sampleCount from molecularDrivers where (driverlikelihood='High' or category='c_fusions') and sampleId in (select sampleId from molecular where containsTumorCells=1) group by 1 order by 2 desc;"
query_patients_second_eval_phase <- "select * from patientsEvaluationPhase;"
query_evaluated_trials_second_eval_phase <- "select distinct patientId, code from trialMatch t LEFT JOIN treatmentMatch tm ON tm.id = t.treatmentMatchId;"

sample <- dbGetQuery(dbActin, query_sample)
patient <- dbGetQuery(dbActin, query_patient)
clinicalStatus <- dbGetQuery(dbActin, query_clinical_status)
tumor <- dbGetQuery(dbActin, query_tumor)
eligibleCohorts <- dbGetQuery(dbActin, query_eligibleCohorts)
molecular <- dbGetQuery(dbActin, query_molecular)
drivers <- dbGetQuery(dbActin, query_drivers)
driversPerSample <- dbGetQuery(dbActin, query_drivers_per_sample)
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

#row_to_move_1 <- tumorType[!is.na(tumorType$primaryTumorLocation) & tumorType$primaryTumorLocation=="Other (n<=2)", ]
#row_to_move_2 <- tumorType[!is.na(tumorType$primaryTumorLocation) & !tumorType$primaryTumorLocation == "Other (n<=2)", ]
#row_to_drop <- tumorType[is.na(tumorType$primaryTumorLocation), ]

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

#row_to_move_1 <- biopsyLocation[!is.na(biopsyLocation$biopsyLocationCurated) & biopsyLocation$biopsyLocationCurated=="Other (n<=2)", ]
#row_to_move_2 <- biopsyLocation[!is.na(biopsyLocation$biopsyLocationCurated) & !biopsyLocation$biopsyLocationCurated == "Other (n<=2)", ]
#row_to_drop <- biopsyLocation[is.na(biopsyLocation$biopsyLocationCurated), ]

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

# Driver events per sample (only selecting successful samples)
driversPerSuccessfulSample <- driversPerSample %>% filter(sampleId %in% molecularSuccessful$sampleId)

median(driversPerSuccessfulSample$driverCount)
min(driversPerSuccessfulSample$driverCount)
max(driversPerSuccessfulSample$driverCount)
driversPerSuccessfulSample %>% filter(driverCount>0) %>% nrow()

# High TML, TMB, HRD and MSI
molecularSuccessful %>% filter(isMicrosatelliteUnstable==1) %>% nrow()
molecularSuccessful %>% filter(isHomologousRepairDeficient==1) %>% nrow()
molecularSuccessful %>% filter(tumorMutationalBurden>=10) %>% nrow()
molecularSuccessful %>% filter(tumorMutationalLoad>=140) %>% nrow()

# Recurring genes in driver events (excluding fusions and viruses)
head(driversGeneCount %>% add_column(percentage = round(driversGeneCount$sampleCount/n_successful_reports*100,1)),5)

# High-driver likelihood WGS aberrations per gene (excluding fusions and viruses) + plot
driversSuccessfulRelevant <- drivers %>% filter(sampleId %in% molecularSuccessful$sampleId)

driversSuccessfulRelevantPerGene <- driversSuccessfulRelevant %>%
  group_by(gene) %>%
  summarise(count=n_distinct(sampleId)) %>%
  mutate(gene = ifelse(count <= 5, "Other (n<=5)", gene)) %>%
  group_by(gene) %>%
  summarize_all(sum) %>% 
  arrange(count)

rows_to_keep <- driversSuccessfulRelevantPerGene[!is.na(driversSuccessfulRelevantPerGene$gene) & !driversSuccessfulRelevantPerGene$gene=="Other (n<=5)" & !driversSuccessfulRelevantPerGene$gene =="NA", ]

driversSuccessfulRelevantPerGeneRearranged <- rbind(rows_to_keep)
driversSuccessfulRelevantPerGeneRearranged$gene <- factor(driversSuccessfulRelevantPerGeneRearranged$gene, levels = driversSuccessfulRelevantPerGeneRearranged$gene)
driversSuccessfulRelevantPerGeneRearranged %>% ggplot(aes(x=count, y=gene)) + geom_bar(stat="identity") + theme_light()

# Fusions and viruses
driversSuccessfulRelevant %>% 
  filter(grepl('fusion', category)) %>% nrow()

driversSuccessfulRelevant %>% 
  filter(grepl('fusion', category)) %>% select(sampleId) %>% n_distinct()

driversSuccessfulRelevant %>% 
  filter(grepl('fusion', category) & driverLikelihood == 'High') 

driversSuccessfulRelevant %>% 
  filter(grepl('virus', category)) %>% nrow()

driversSuccessfulRelevant %>% 
  filter(grepl('virus', category)) %>% select(sampleId) %>% n_distinct()

# Actionable events plot
actionableEventsSuccessfulRelevant <- actionableEvents %>% filter(sampleId %in% molecularSuccessful$sampleId)

actionableEventsSuccessfulRelevant %>% 
  filter(potentially.actionable. == 'Yes') %>% select(sampleId) %>% n_distinct()

actionableEventsSuccessfulRelevantPerEvent <- actionableEventsSuccessfulRelevant %>%
  filter(potentially.actionable. == 'Yes') %>%
  group_by(curated.event) %>%
  summarise(count=n_distinct(sampleId)) %>%
  mutate(curated.event = ifelse(count <= 5, "Other (n<=5)", curated.event)) %>%
  group_by(curated.event) %>%
  summarize_all(sum) %>% 
  arrange(count)

actionableEventsSuccessfulRelevantPerEvent %>% tail(6)

rows_to_keep <- actionableEventsSuccessfulRelevantPerEvent[!is.na(actionableEventsSuccessfulRelevantPerEvent$curated.event) & !actionableEventsSuccessfulRelevantPerEvent$curated.event=="Other (n<=5)", ]

actionableEventsSuccessfulRelevantPerEventRearranged <- rbind(rows_to_keep)
actionableEventsSuccessfulRelevantPerEventRearranged$curated.event <- factor(actionableEventsSuccessfulRelevantPerEventRearranged$curated.event, levels = actionableEventsSuccessfulRelevantPerEventRearranged$curated.event)
actionableEventsSuccessfulRelevantPerEventRearranged %>% ggplot(aes(x=count, y=curated.event)) + geom_bar(stat="identity") + theme_light()

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
