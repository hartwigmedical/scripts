library(RMySQL)
library(dplyr)
library(tidyr)

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
dbProd = dbConnect(MySQL(), dbname='hmfpatients_20180418', groups="RAnalysis")
databaseClinicalData = purple::query_clinical_data(dbProd) %>% 
  mutate(biopsyLocation = clean_data(biopsyLocation))

databaseClinicalData = databaseClinicalData %>% 
  select(sampleId, primaryTumorLocation, cancerSubtype, biopsyDate, biopsySite, biopsyType, biopsyLocation, treatment, treatmentType, birthYear)

dbDisconnect(dbProd)
rm(dbProd)

### Curated Data
curatedClinicalData = read.csv(file = "~/hmf/resources/CuratedClinicalData.csv", stringsAsFactors = F, header = T) %>%
  mutate(
    gender = tolower(clean_data(gender)),
    biopsyDate = clean_data(biopsyDate), 
    biopsySite = clean_data(biopsySite), 
    biopsyType = clean_data(biopsyType), 
    biopsyLocation = clean_data(biopsyLocation), 
    birthYear = clean_data(birthYear), 
    primaryTumorLocation = capitalise(clean_data(primaryTumorLocation))) %>%
  filter(!is.na(primaryTumorLocation)) %>%
  select(-gender, -patientId)


### Curated Data Corrections
curatedClinicalData[curatedClinicalData$primaryTumorLocation == "GIST",  "cancerSubtype"] <- "Gastrointestinal stromal tumor (GIST)"
curatedClinicalData[curatedClinicalData$primaryTumorLocation == "GIST",  "primaryTumorLocation"] <- "Bone/Soft tissue"
curatedClinicalData[curatedClinicalData$primaryTumorLocation == "Sarcoma",  "cancerSubtype"] <- "Sarcoma"
curatedClinicalData[curatedClinicalData$primaryTumorLocation == "Sarcoma",  "primaryTumorLocation"] <- "Bone/Soft tissue"

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
clinicalData[disagree(clinicalData$birthYear.database, clinicalData$birthYear.curated), ]
clinicalData[disagree(clinicalData$primaryTumorLocation.database, clinicalData$primaryTumorLocation.curated), ]
clinicalData[disagree(clinicalData$biopsySite.database, clinicalData$biopsySite.curated), ]
clinicalData[disagree(clinicalData$biopsyDate.database, clinicalData$biopsyDate.curated), ]
clinicalData[disagree(clinicalData$cancerSubtype.database, clinicalData$cancerSubtype.curated), ]
clinicalData[disagree(clinicalData$biopsyType.database, clinicalData$biopsyType.curated), ]

clinicalData[disagree(clinicalData$biopsyLocation.database, clinicalData$biopsyLocation.curated), ]
clinicalData[disagree(clinicalData$biopsyLocation.database, clinicalData$biopsyLocation), ]

clinicalData = clinicalData %>% select(-ends_with(".curated"), -ends_with(".database")) 
save(clinicalData, file = "~/hmf/RData/reference/clinicalData.RData")


#### Primary Tumor Location Colours
load(file = "~/hmf/RData/reference/clinicalData.RData")
primaryTumorLocations = sort(unique(clinicalData$primaryTumorLocation))
primaryTumorLocations= primaryTumorLocations[!is.na(primaryTumorLocations)]

cosmicSignatureColours = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
                           "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a",
                           "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
                           "#dea185","#a0729d","#8a392f")
primaryTumorLocationColours = setNames(cosmicSignatureColours[1:length(primaryTumorLocations)], primaryTumorLocations)
save(primaryTumorLocationColours, file = "~/hmf/RData/reference/primaryTumorLocationColours.RData")
