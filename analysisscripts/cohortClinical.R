library(RMySQL)
library(dplyr)
library(tidyr)

#### Curated Data
allClinicalDataCurated = read.csv(file = "~/hmf/resources/ClinicalData20180516.csv", stringsAsFactors = F, header = T)

### Database Data
dbProd = dbConnect(MySQL(), dbname='hmfpatients_20180418', groups="RAnalysis")
allClinicalDataDatabase = purple::query_clinical_data(dbProd) %>%
  select(-starts_with("biopsy"), -primaryTumorLocation, -cancerSubtype, -birthYear)
dbDisconnect(dbProd)
rm(dbProd)

### Combined
allClinicalData = left_join(allClinicalDataCurated, allClinicalDataDatabase, by = "sampleId") %>%
  group_by(primaryTumorLocation) %>% 
  mutate(primaryTumorLocationSamples = n()) %>%
  mutate(cancerType = ifelse(primaryTumorLocationSamples > 15, primaryTumorLocation, "Other")) %>%
  ungroup() %>%
  select(-primaryTumorLocationSamples)
save(allClinicalData, file = "~/hmf/RData/Reference/allClinicalData.RData")
rm(allClinicalDataCurated, allClinicalDataDatabase)


############################## MANUAL MODIFICATIONS ##############################
allClinicalDataCurated = read.csv(file = "~/hmf/resources/ClinicalData20180512.csv", stringsAsFactors = F, header = T) %>%
  select(-treatment, -treatmentType) 
allClinicalDataCurated[allClinicalDataCurated$sampleId == "DRUP01020004T", "birthYear"] <- 1985

allClinicalDataCurated[allClinicalDataCurated$sampleId == "CPCT02070029T", "primaryTumorLocation"] <- "Colon/Rectum"
allClinicalDataCurated[allClinicalDataCurated$sampleId == "CPCT02070029T", "cancerSubtype"] <- ""
allClinicalDataCurated[allClinicalDataCurated$sampleId == "DRUP01070039T", "primaryTumorLocation"] <- "Colon/Rectum"
allClinicalDataCurated[allClinicalDataCurated$sampleId == "DRUP01070039T", "cancerSubtype"] <- ""

allClinicalDataCurated[allClinicalDataCurated$sampleId == "CPCT02080086T", "primaryTumorLocation"] <- "NET"
allClinicalDataCurated[allClinicalDataCurated$sampleId == "CPCT02080086T", "cancerSubtype"] <- ""
allClinicalDataCurated[allClinicalDataCurated$sampleId == "DRUP01080005T", "primaryTumorLocation"] <- "NET"
allClinicalDataCurated[allClinicalDataCurated$sampleId == "DRUP01080005T", "cancerSubtype"] <- ""

allClinicalDataCurated[allClinicalDataCurated$sampleId == "CPCT02050116T", "primaryTumorLocation"] <- "CNS"
allClinicalDataCurated[allClinicalDataCurated$sampleId == "CPCT02050116T", "cancerSubtype"] <- "Neuroblastoma"
allClinicalDataCurated[allClinicalDataCurated$sampleId == "DRUP01050008T", "primaryTumorLocation"] <- "CNS"
allClinicalDataCurated[allClinicalDataCurated$sampleId == "DRUP01050008T", "cancerSubtype"] <- "Neuroblastoma"

allClinicalDataCurated[allClinicalDataCurated$sampleId == "CPCT02030349T", "primaryTumorLocation"] <- "CUP"
allClinicalDataCurated[allClinicalDataCurated$sampleId == "CPCT02030349T", "cancerSubtype"] <- "Endometrial or Ovarium"
allClinicalDataCurated[allClinicalDataCurated$sampleId == "DRUP01030013T", "primaryTumorLocation"] <- "CUP"
allClinicalDataCurated[allClinicalDataCurated$sampleId == "DRUP01030013T", "cancerSubtype"] <- "Endometrial or Ovarium"

allClinicalDataCurated[allClinicalDataCurated$sampleId == "CPCT02050240T", "primaryTumorLocation"] <- "Skin"
allClinicalDataCurated[allClinicalDataCurated$sampleId == "CPCT02050240T", "cancerSubtype"] <- "Atypical fibroxanthoma"
allClinicalDataCurated[allClinicalDataCurated$sampleId == "DRUP01050021T", "primaryTumorLocation"] <- "Skin"
allClinicalDataCurated[allClinicalDataCurated$sampleId == "DRUP01050021T", "cancerSubtype"] <- "Atypical fibroxanthoma"

allClinicalDataCurated[allClinicalDataCurated$sampleId == "CPCT02230072T", "primaryTumorLocation"] <- "Skin"
allClinicalDataCurated[allClinicalDataCurated$sampleId == "CPCT02230072T", "cancerSubtype"] <- "Skin squamous cell carcinoma"
allClinicalDataCurated[allClinicalDataCurated$sampleId == "DRUP01230001T", "primaryTumorLocation"] <- "Skin"
allClinicalDataCurated[allClinicalDataCurated$sampleId == "DRUP01230001T", "cancerSubtype"] <- "Skin squamous cell carcinoma"

write.csv(allClinicalDataCurated, file = "~/hmf/resources/ClinicalData20180516.csv", row.names = F)

############################## OLD ##############################
clean_data <- function(vector) {
  result <- ifelse(vector == '#N/A', NA, vector)
  result <- ifelse(substr(result,1, 1) == '?', NA, result)
  result <- ifelse(result == 'NULL', NA, result)
  result <- ifelse(result == 'NA', NA, result)
  return (result)
}

disagree <- function(database, curated) {
  !is.na(database) & (is.na(curated) | database != curated)
}

improve <- function(database, curated) {
  is.na(database) & !is.na(curated)
}

capitalise <- function(vector) {
  ifelse(is.na(vector), NA, paste0(toupper(substr(vector, 1, 1)), substring(vector, 2)))
}

### Database Data
databaseClinicalData = purple::query_clinical_data(dbProd) %>% 
  mutate(biopsyLocation = clean_data(biopsyLocation))

databaseClinicalData = databaseClinicalData %>% 
  select(sampleId, primaryTumorLocation, cancerSubtype, biopsyDate, biopsySite, biopsyType, biopsyLocation, treatment, treatmentType, birthYear)

curatedClinicalDatar
curatedClinicalData = read.csv(file = "~/hmf/RData/ClinicalData20180412.csv", stringsAsFactors = F, header = T) 

### Sanity Checks (should be false)
any(is.na(curatedClinicalData$primaryTumorLocation))
any(!curatedClinicalData$primaryTumorLocation %in% databaseClinicalData$primaryTumorLocation)

### Join data
clinicalData = left_join(databaseClinicalData, curatedClinicalData, suffix = c('.database','.curated'), by = "sampleId") %>%
  mutate(birthYear = birthYear.curated,
         primaryTumorLocation = primaryTumorLocation.curated,
         biopsySite = biopsySite.curated,
         biopsyDate = biopsyDate.curated,
         cancerSubtype = cancerSubtype.curated,
         biopsyType = biopsyType.curated,
         biopsyLocation = ifelse(disagree(biopsyLocation.database, biopsyLocation.curated), biopsyLocation.database, biopsyLocation.curated)
         )

#### Check fields that disagree
jon = clinicalData[disagree(clinicalData$birthYear.database, clinicalData$birthYear.curated), ]
jon = clinicalData[disagree(clinicalData$primaryTumorLocation.database, clinicalData$primaryTumorLocation.curated), ]
jon = clinicalData[disagree(clinicalData$biopsySite.database, clinicalData$biopsySite.curated), ]
jon = clinicalData[disagree(clinicalData$biopsyDate.database, clinicalData$biopsyDate.curated), ]
jon = clinicalData[disagree(clinicalData$cancerSubtype.database, clinicalData$cancerSubtype.curated), ]
jon = clinicalData[disagree(clinicalData$biopsyType.database, clinicalData$biopsyType.curated), ]

clinicalData[disagree(clinicalData$biopsyLocation.database, clinicalData$biopsyLocation.curated), ]
clinicalData[disagree(clinicalData$biopsyLocation.database, clinicalData$biopsyLocation), ]

clinicalData = clinicalData %>% select(-ends_with(".curated"), -ends_with(".database")) 
save(clinicalData, file = "~/hmf/RData/reference/clinicalData.RData")

jon = read.csv(file = "~/hmf/RData/JonClinicalData.csv", stringsAsFactors = F)

edwin = read.csv(file = "~/hmf/RData/ClinicalDataEdwin.csv", stringsAsFactors = F)
sort(unique(edwin$cancerType))

edwin = edwin %>% group_by(primaryTumorLocation) %>% mutate(primaryTumorLocationSamples = n()) %>%
  mutate(cancerType = ifelse(primaryTumorLocationSamples> 15, primaryTumorLocation, "Other"))


jon2 = originalCurated %>% filter(sampleId %in% jon$sampleId)
jon2$sampleId


############################################ PREVIOUS CURATED DATA

### Curated Data
originalCurated = read.csv(file = "~/hmf/resources/CuratedClinicalData.csv", stringsAsFactors = F, header = T) %>%
  mutate(
    gender = tolower(clean_data(gender)),
    biopsyDate = clean_data(biopsyDate), 
    biopsySite = clean_data(biopsySite), 
    biopsyType = clean_data(biopsyType), 
    biopsyLocation = clean_data(biopsyLocation), 
    birthYear = clean_data(birthYear), 
    primaryTumorLocation = capitalise(clean_data(primaryTumorLocation))) %>%
  select(-gender, -patientId)

### Curated Data Corrections
originalCurated[originalCurated$primaryTumorLocation == "GIST",  "cancerSubtype"] <- "Gastrointestinal stromal tumor (GIST)"
originalCurated[originalCurated$primaryTumorLocation == "GIST",  "primaryTumorLocation"] <- "Bone/Soft tissue"
originalCurated[originalCurated$primaryTumorLocation == "Sarcoma",  "cancerSubtype"] <- "Sarcoma"
originalCurated[originalCurated$primaryTumorLocation == "Sarcoma",  "primaryTumorLocation"] <- "Bone/Soft tissue"



