library(dplyr)
library(tibble)
library(RMySQL)
library(ggplot2)

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
query_drivers_per_sample <- "select m.sampleId, count(md.driverLikelihood) as highDriverCount from molecular m left join molecularDrivers md on m.sampleId=md.sampleId and driverlikelihood='high' group by 1;"
query_drivers_gene_count <- "select gene, count(distinct sampleId) as sampleCount from molecularDrivers where driverLikelihood='High' and sampleId in (select sampleId from molecular where containsTumorCells=1) group by 1 order by 2 desc;"

sample <- dbGetQuery(dbActin, query_sample)
patient <- dbGetQuery(dbActin, query_patient)
clinicalStatus <- dbGetQuery(dbActin, query_clinical_status)
tumor <- dbGetQuery(dbActin, query_tumor)
eligibleCohorts <- dbGetQuery(dbActin, query_eligibleCohorts)
molecular <- dbGetQuery(dbActin, query_molecular)
drivers <- dbGetQuery(dbActin, query_drivers)
driversPerSample <- dbGetQuery(dbActin, query_drivers_per_sample)
driversGeneCount <- dbGetQuery(dbActin, query_drivers_gene_count)

dbDisconnect(dbActin)

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

# clinical
patientCount <- nrow(patient)
patient %>% group_by(gender) %>% summarise(count=n(), percentage=round(n()/patientCount*100,1))

patient <- patient %>% add_column(registrationYear = substr(patient$registrationDate,1,4))
patient$registrationYear <- as.integer(patient$registrationYear)
patient$birthYear <- as.integer(patient$birthYear)
patient <- patient %>% add_column(ageAtRegistration = patient$registrationYear-patient$birthYear)

median(patient$ageAtRegistration)
min(patient$ageAtRegistration)
max(patient$ageAtRegistration)

clinicalStatus %>% group_by(who) %>% summarise(count=n(), percentage=round(n()/patientCount*100,1))

## DRIVERS
#driverCuration <- read.csv(paste0(Sys.getenv("HOME"), "/hmf/tmp/Curation - Drivers.csv"), sep = ",")
#driver <- inner_join(drivers, driverCuration, by=c('inclusionMolecularEvents'))
#arranged <- driver %>% group_by(inclusionMolecularEventsCurated) %>% summarise(count=n_distinct(patientId))

# Successful WGS reports
molecularSuccessful <- molecular %>% filter(containsTumorCells==1)
n_successful_reports <- molecularSuccessful %>% nrow()

# Driver events per sample
driversPerSuccessfulSample <- driversPerSample %>% filter(sampleId %in% molecularSuccessful$sampleId)

median(driversPerSuccessfulSample$highDriverCount)
min(driversPerSuccessfulSample$highDriverCount)
max(driversPerSuccessfulSample$highDriverCount)
driversPerSuccessfulSample %>% filter(highDriverCount>0) %>% nrow()

# High TML, TMB, HRD and MSI
molecularSuccessful %>% filter(isMicrosatelliteUnstable==1) %>% nrow()
molecularSuccessful %>% filter(isHomologousRepairDeficient==1) %>% nrow()
molecularSuccessful %>% filter(tumorMutationalBurden>=10) %>% nrow()
molecularSuccessful %>% filter(tumorMutationalLoad>=140) %>% nrow()

# Recurring genes in driver events (excluding fusions and viruses)
head(driversGeneCount %>% add_column(percentage = round(driversGeneCount$sampleCount/n_successful_reports*100,1)),5)

# High-driver likelihood WGS aberrations per gene (excluding fusions and viruses)
driversSuccessfulHigh <- drivers %>% filter(driverLikelihood=='High' & sampleId %in% molecularSuccessful$sampleId)

driversSuccessfulHighPerGene <- driversSuccessfulHigh %>%
  group_by(gene) %>%
  summarise(count=n_distinct(sampleId)) %>%
  mutate(gene = ifelse(count <= 5, "Other (n<=5)", gene)) %>%
  group_by(gene) %>%
  summarize_all(sum) %>% 
  arrange(count)

rows_to_keep <- driversSuccessfulHighPerGene[!is.na(driversSuccessfulHighPerGene$gene) & !driversSuccessfulHighPerGene$gene=="Other (n<=5)" & !driversSuccessfulHighPerGene$gene =="NA", ]
rows_to_drop <- driversSuccessfulHighPerGene[is.na(driversSuccessfulHighPerGene$gene) | driversSuccessfulHighPerGene$gene=="Other (n<=5)" | driversSuccessfulHighPerGene$gene =="NA", ]

driversSuccessfulHighPerGeneRearranged <- rbind(rows_to_keep)
driversSuccessfulHighPerGeneRearranged$gene <- factor(driversSuccessfulHighPerGeneRearranged$gene, levels = driversSuccessfulHighPerGeneRearranged$gene)
driversSuccessfulHighPerGeneRearranged %>% ggplot(aes(x=count, y=gene)) + geom_bar(stat="identity") + theme_light()

# Fusions and viruses
driversSuccessfulHigh %>% 
  filter(grepl('fusion', category)) %>% nrow()

driversSuccessfulHigh %>% 
  filter(grepl('virus', category)) %>% nrow()
