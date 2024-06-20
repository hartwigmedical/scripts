library(dplyr)
library(tibble)

rm(list=ls())

# Retrieve data ------------------------------------------------------------------
dbActin <- dbConnect(MySQL(), dbname='actin_pilot', groups="RAnalysis")

query_sample <-"select * from actinLastResortPaperSamples;"
query_patient <-"select * from patient where patientId in (select patientId from actinLastResortPaperSamples);"
query_clinical_status <-"select * from clinicalStatus where patientId in (select patientId from actinLastResortPaperSamples);"
query_molecular <-"select * from molecular where sampleId in (select sampleId from actinLastResortPaperSamples);"
query_tumor <-"select * from tumor where patientId in (select patientId from actinLastResortPaperSamples);"
query_drivers <- "select distinct patientId, trialId, inclusionMolecularEvents from trialEvaluation where inclusionMolecularEvents not in ('') and isEligibleTrial and isEligibleCohort and patientId in (select patientId from actinLastResortPaperSamples);"
query_drivers_2 <- "select gene, count(distinct sampleId) from hmfpatients.driverCatalog where driverLikelihood>0.8 and sampleId in (select sampleId from actinLastResortPaperSamples) and sampleId in (select sampleId from hmfpatients.purity where qcStatus not like '%FAIL%') group by 1;"
query_drivers_3 <- "select sampleId, count(*) from hmfpatients.driverCatalog where driverLikelihood>0.8 and sampleId in (select sampleId from actinLastResortPaperSamples) and sampleId in (select sampleId from hmfpatients.purity where qcStatus not like '%FAIL%') group by 1;"

sample <- dbGetQuery(dbActin, query_sample)
patient <- dbGetQuery(dbActin, query_patient)
clinicalStatus <- dbGetQuery(dbActin, query_clinical_status)
molecular <- dbGetQuery(dbActin, query_molecular)
tumor <- dbGetQuery(dbActin, query_tumor)
drivers <- dbGetQuery(dbActin, query_drivers)
drivers_2 <- dbGetQuery(dbActin, query_drivers_2)
drivers_3 <- dbGetQuery(dbActin, query_drivers_3)

dbDisconnect(dbActin)

# Tumor type WIP
tumorType <- tumor %>% group_by(primaryTumorLocation) %>% summarise(count=n()) %>%
  mutate(primaryTumorLocation = ifelse(count == 1 | count == 2, "Other (n<=2)", primaryTumorLocation)) %>%
  group_by(primaryTumorLocation) %>%
  summarize_all(sum) 

tumorType$primaryTumorLocation <- factor(tumorType$primaryTumorLocation, levels = tumorType$primaryTumorLocation[order(tumorType$count)])
tumorType %>% ggplot(aes(x=primaryTumorLocation, y=count)) + geom_bar(stat="identity") + coord_flip()

biopsyCuration <- read.csv(paste0(Sys.getenv("HOME"), "/hmf/tmp/Curation - Biopsy location"), sep = ",")
biopsy <- inner_join(tumor, biopsyCuration, by=c('biopsyLocation'))
biopsyLocation <- biopsy %>% group_by(biopsyLocationCurated) %>% summarise(count=n()) %>%
  mutate(biopsyLocationCurated = ifelse(count == 1 | count == 2, "Other (n<=2)", biopsyLocationCurated)) %>%
  group_by(biopsyLocationCurated) %>%
  summarize_all(sum) 

biopsyLocation$biopsyLocationCurated <- factor(biopsyLocation$biopsyLocationCurated, levels = biopsyLocation$biopsyLocationCurated[order(biopsyLocation$count)])
biopsyLocation %>% ggplot(aes(x=biopsyLocationCurated, y=count)) + geom_bar(stat="identity") + coord_flip()

# clinical WIP
patient %>% group_by(gender) %>% summarise(count=n())
patient %>% group_by(birthYear) %>% summarise(count=n())
clinicalStatus %>% group_by(who) %>% summarise(count=n())

# drivers WIP
driverCuration <- read.csv(paste0(Sys.getenv("HOME"), "/hmf/tmp/Curation - Drivers.csv"), sep = ",")
driver <- inner_join(drivers, driverCuration, by=c('inclusionMolecularEvents'))
arranged <- driver %>% group_by(inclusionMolecularEventsCurated) %>% summarise(count=n_distinct(patientId))
