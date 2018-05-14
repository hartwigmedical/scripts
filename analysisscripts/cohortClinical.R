library(RMySQL)
library(dplyr)
library(tidyr)

#### Curated Data
allClinicalDataCurated = read.csv(file = "~/hmf/resources/ClinicalData20180512.csv", stringsAsFactors = F, header = T) %>%
  select(-treatment, -treatmentType) 
allClinicalDataCurated[allClinicalDataCurated$sampleId == "DRUP01020004T", "birthYear"] <- 1985

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

### Cancer Type Colours
load(file = "~/hmf/RData/reference/allClinicalData.RData")
cancerTypes = sort(unique(allClinicalData$cancerType))
cosmicSignatureColours = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
                           "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a",
                           "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
                           "#dea185","#a0729d","#8a392f")
cancerTypeColours = setNames(cosmicSignatureColours[1:length(cancerTypes)], cancerTypes)
save(cancerTypeColours, file = "~/hmf/RData/reference/cancerTypeColours.RData")



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



