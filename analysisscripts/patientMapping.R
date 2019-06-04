library(dplyr)
library(tidyr)
library(RMySQL)

sample_mapping <- function(amberBafs, percentCutoff = 0.6) {
  simpleBafs = amberBafs %>% mutate(loci = paste0(chromosome ,":", position)) %>% select(sampleId, loci) %>% mutate(het = T)
  spreadBafs = simpleBafs %>% spread(loci, het, fill = F)
  bafMatrix = as.matrix(spreadBafs %>% select(-sampleId))
  rownames(bafMatrix) <- spreadBafs$sampleId
  matrixResult = bafMatrix %*% t(bafMatrix)
  dfResult = data.frame(matrixResult)
  dfResult$sample1 <- spreadBafs$sampleId
  dfResult = dfResult %>% gather(sample2, match,  spreadBafs$sampleId)
  
  sampleCounts = dfResult %>% filter(sample1 == sample2) %>% select(sample1, count = match)
  
  sampleMapping = dfResult %>%
    filter(sample1 != sample2) %>%
    left_join(sampleCounts, by = "sample1") %>%
    mutate(matchPercent = match / count) %>%
    filter(matchPercent > percentCutoff) %>%
    mutate(tmp = sample1, sample1 = ifelse(sample1 < sample2, sample1, sample2), sample2 = ifelse(tmp < sample2, sample2, tmp)) %>%
    group_by(sample1, sample2) %>% summarise(matchPercent = max(matchPercent))
  
  return (sampleMapping)
}

patient_mapping <- function(sampleMapping, existingMappings) {
  existingSamples = existingMappings %>% filter(type == 'Sample') %>% select(sampleId = sourceId) %>%
    mutate(match = T)
  
  existingPatientMappings = existingMappings %>% filter(type == 'Patient') %>% select(sourceId, targetId) %>% mutate(new = F)
  
  patientMapping = sampleMapping %>%
    left_join(existingSamples %>% select(sample1 = sampleId, sample1Prior = match), by = "sample1") %>%
    left_join(existingSamples %>% select(sample2 = sampleId, sample2Prior = match), by = "sample2") %>%
    mutate(sample1Prior = ifelse(is.na(sample1Prior), F, T), sample2Prior = ifelse(is.na(sample2Prior), F, T)) %>%
    mutate(patient1 = substr(sample1, 1, 12), patient2 = substr(sample2, 1, 12)) %>%
    filter(patient1 != patient2) %>%
    mutate(
      sourceId = ifelse(patient1 > patient2, patient1, patient2),
      targetId = ifelse(patient1 > patient2, patient2, patient1),
      xorSamples = xor(sample1Prior, sample2Prior),
      sourceId = ifelse(xorSamples, ifelse(sample1Prior, patient2, patient1), sourceId),
      targetId = ifelse(xorSamples, ifelse(sample1Prior, patient1, patient2), targetId)
    ) %>%
    group_by(sourceId, targetId) %>% summarise(matchPercent = max(matchPercent)) %>%
    group_by(targetId) %>% mutate(target_n = n()) %>%
    group_by(sourceId) %>% mutate(source_n = n()) %>%
    arrange(-target_n) %>%
    filter(row_number() == 1) %>%
    select(sourceId, targetId, matchPercent)
  
  patientMapping = patientMapping %>% left_join(existingPatientMappings, by = c("sourceId", "targetId"))
  patientMapping[is.na(patientMapping)] <- T
  
  return (patientMapping %>% arrange(-new, matchPercent))
}

######################## Step 1 - Regenerate patient mappings
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis", host = "127.0.0.1")
amberBafs = dbGetQuery(dbProd, "SELECT * FROM amber")
existingMappings = dbGetQuery(dbProd, "SELECT * FROM idMapping")
dbDisconnect(dbProd)
rm(dbProd)

sampleMapping = sample_mapping(amberBafs)
patientMapping = patient_mapping(sampleMapping, existingMappings)

amberSamples = unique(amberBafs$sampleId)
write.table(amberSamples, file = "~/hmf/resources/idgenerator/samples_20190604.csv", quote = F, row.names = F, col.names = F, sep = ",")
write.table(patientMapping %>% select(sourceId, targetId), file = "~/hmf/resources/idgenerator/patient_mapping_20190604.csv", quote = F, row.names = F, col.names = F, sep = ",")


######################## Step 2 - Update hashes
# com.hartwig.hmftools.idgenerator.HmfIdApplicationKt
# -update_ids
# -password xxxxx 
# -sample_ids_file ~/hmf/resources/idgenerator/samples_20190604.csv
# -patient_mapping_file ~/hmf/resources/idgenerator/patient_mapping_20190604.csv
# -out ~/hmf/resources/idgenerator/sample_hashes_20190604.csv
# -mapping_out ~/hmf/resources/idgenerator/remapping_20190604.csv

######################## Step 3 - Include hash in build
# cp ~/hmf/resources/idgenerator/sample_hashes_20190604.csv ~/hmf/repos/hmftools/hmf-id-generator/src/main/resources/sample_hashes.csv

######################## Step 4 - Run Anonymise Ids
# com.hartwig.hmftools.idgenerator.HmfIdApplicationKt
# -anonymize_ids
# -password xxxxx
# -sample_ids_file ~/hmf/resources/idgenerator/samples_20190604.csv
# -patient_mapping_file ~hmf/resources/idgenerator/patient_mapping_20190604.csv

######################## Step 5 - Copy output
# Need to copy the relevant run output without header to ~/hmf/resources/idgenerator/anonymized_20190604.csv

######################## Step 6 - Load into database
# Note we do not delete type = 'Sample' so we can have a history. 

#delete from idMapping where type in ('Patient', 'Remapping');
#LOAD DATA LOCAL INFILE '~/hmf/resources/idgenerator/patient_mapping_20190604.csv' IGNORE INTO TABLE idMapping  FIELDS TERMINATED BY ',' IGNORE 0 LINES (sourceId, targetId) SET modified = CURRENT_TIMESTAMP, type = 'Patient';
#LOAD DATA LOCAL INFILE '~/hmf/resources/idgenerator/anonymized_20190604.csv' IGNORE INTO TABLE idMapping FIELDS TERMINATED BY ',' IGNORE 0 LINES  (sourceId, targetId) SET modified = CURRENT_TIMESTAMP, type = 'Sample';
#LOAD DATA LOCAL INFILE '~/hmf/resources/idgenerator/remapping_20190604.csv' IGNORE INTO TABLE idMapping  FIELDS TERMINATED BY ',' LINES TERMINATED BY '\r\n' IGNORE 0 LINES (sourceId, targetId) SET modified = CURRENT_TIMESTAMP, type = 'Remapping';

######################## Step 7 - Store files
# copy all new files to datastore, eg
# scp *20190604.csv hmf-datastore:/data/common/dbs/idgenerator/

