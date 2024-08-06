library(dplyr)
library(tibble)
library(RMySQL)

rm(list=ls())

# Retrieve data ------------------------------------------------------------------
dbActin <- dbConnect(MySQL(), dbname='actin_paper', groups="RAnalysis")

query_sample <-"select * from paperSamples;"
query_patient <-"select * from patient;"
query_clinical_status <-"select * from clinicalStatus;"
query_tumor <-"select * from tumor;"
query_eligibleCohorts <-"select * from eligibleCohorts union select * from eligibleCohorts_addition;"
query_molecular <-"select * from molecular;"
query_drivers <- "select * from molecularDrivers;"
query_drivers_per_sample <- "select m.sampleId, count(md.driverLikelihood) as highDriverCount from molecular m left join molecularDrivers md on m.sampleId=md.sampleId and driverlikelihood='high' where containsTumorCells=1 group by 1;"
query_drivers_gene_count <- "select gene, count(distinct sampleId) as sampleCount from molecularDrivers where driverLikelihood='High' and sampleId in (select sampleId from molecular where containsTumorCells=1) group by 1 order by 2 desc;"

sample <- dbGetQuery(dbActin, query_sample)
patient <- dbGetQuery(dbActin, query_patient)
clinicalStatus <- dbGetQuery(dbActin, query_clinical_status)
tumor <- dbGetQuery(dbActin, query_tumor)
eligibleCohorts <- dbGetQuery(dbActin, query_eligibleCohorts)
molecular <- dbGetQuery(dbActin, query_molecular)
drivers <- dbGetQuery(dbActin, query_drivers)
drivers_per_sample <- dbGetQuery(dbActin, query_drivers_per_sample)
drivers_gene_count <- dbGetQuery(dbActin, query_drivers_gene_count)

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

## DRIVERS - WIP
#driverCuration <- read.csv(paste0(Sys.getenv("HOME"), "/hmf/tmp/Curation - Drivers.csv"), sep = ",")
#driver <- inner_join(drivers, driverCuration, by=c('inclusionMolecularEvents'))
#arranged <- driver %>% group_by(inclusionMolecularEventsCurated) %>% summarise(count=n_distinct(patientId))

# Driver events per sample
median(drivers_per_sample$highDriverCount)
min(drivers_per_sample$highDriverCount)
max(drivers_per_sample$highDriverCount)
drivers_per_sample %>% dplyr::filter(highDriverCount>0) %>% nrow()

# Recurring genes in driver events (excluding fusions and viruses)
sampleCount <- nrow(drivers_per_sample)
head(drivers_gene_count %>% add_column(percentage = round(drivers_gene_count$sampleCount/sampleCount*100,1)),3)
