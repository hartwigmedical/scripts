library(dplyr)
library(tidyr)
library(RMySQL)

######################## Generate BED file of AMBER heterozygous locations to use
oldHetLocations = read.table(file = "/Users/jon/hmf/resources/CytoScanHD_hg19_SNPs_sorted.bed", sep = "\t", header = F, stringsAsFactors = F) %>% select(chr = V1, start = V2, end = V3) %>% mutate(old = T)
newHetLocations = read.table(file = "/Users/jon/hmf/resources/GermlineHetPon.hg19.bed", sep = "\t", header = F, stringsAsFactors = F)  %>% select(chr = V1, start = V2) %>% mutate(new = T)
combinedHetLocations = inner_join(oldHetLocations, newHetLocations, by = c("chr", "start"))
patientMappingBed = combinedHetLocations %>% filter(chr != 'X') %>%
  mutate(chr = as.numeric(chr)) %>%
  arrange(chr, start) %>%
  select(-old, -new) %>%
  mutate(distance = lead(start) - start) %>%
  filter(distance < 0 | distance > 200000) %>%
  select(-distance)
write.table(patientMappingBed, file = "~/hmf/resources/patientMapping.bed", sep = "\t", row.names = F, col.names = F)

######################## Query AMBER Het Locations
dbProd = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis", host = "127.0.0.1")
amberBafs = dbGetQuery(dbProd, "SELECT * FROM amber")
dbDisconnect(dbProd)
rm(dbProd)
save(amberBafs, file = "~/hmf/analysis/amberancestry/amberBafs.RData")

######################## Sample / Patient Mapping
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

patient_mapping <- function(sampleMapping) {
  patientMapping = sampleMapping %>% 
    mutate(patient1 = substr(sample1, 1, 12), patient2 = substr(sample2, 1, 12)) %>%
    filter(patient1 != patient2) %>% 
    mutate(tmp = patient1, patient1 = ifelse(patient1 > patient2, patient1, patient2), patient2 = ifelse(tmp > patient2, patient2, tmp)) %>%
    group_by(patient1, patient2) %>% summarise(matchPercent = max(matchPercent)) %>% 
    group_by(patient2) %>% mutate(p2n = n()) %>% 
    group_by(patient1) %>% mutate(p1n = n()) %>%
    arrange(-p2n) %>%
    filter(row_number() == 1) %>% 
    select(-p1n, -p2n)
  
  return (patientMapping)
}



#### EXPERIMENTAL
uhoh = read.table(file = "~/hmf/resources/idgenerator/anonymized.csv", header = F, stringsAsFactors = F, sep = ",")

existingMapping = read.table(file = "~/hmf/resources/idgenerator/patient_mapping.csv", header = F, stringsAsFactors = F, sep = ",") %>%
  select(From = V1, To = V2)

previouslyMappedSamples = read.table(file = "~/hmf/resources/idgenerator/samples.csv", header = F, stringsAsFactors = F, sep = ",") %>%
  select(sampleId = V1) %>% mutate(match = T)


patientMapping = sampleMapping %>% 
  left_join(previouslyMappedSamples %>% select(sample1 = sampleId, sample1Prior = match), by = "sample1") %>%
  left_join(previouslyMappedSamples %>% select(sample2 = sampleId, sample2Prior = match), by = "sample2") %>%
  mutate(sample1Prior = ifelse(is.na(sample1Prior), F, T), sample2Prior = ifelse(is.na(sample2Prior), F, T)) %>%
  mutate(patient1 = substr(sample1, 1, 12), patient2 = substr(sample2, 1, 12)) %>%
  filter(patient1 != patient2) %>% 
  mutate(from = ifelse(patient1 > patient2, patient1, patient2), to = ifelse(patient1 > patient2, patient2, patient1))
  group_by(from, to) %>% summarise(matchPercent = max(matchPercent)) %>%
  group_by(to) %>% mutate(to_n = n()) %>% 
  group_by(from) %>% mutate(from_n = n())
#### EXPERIMENTAL




load(file = "~/hmf/analysis/amberancestry/amberBafs.RData")
amberSamples = unique(amberBafs$sampleId)
date()
sampleMapping = sample_mapping(amberBafs)
patientMapping = patient_mapping(sampleMapping)
date()

write.table(amberSamples, file = "~/hmf/resources/idgenerator/amber_samples.csv", quote = F, row.names = F, col.names = F, sep = ",")
write.table(patientMapping, file = "~/hmf/resources/idgenerator/amber_patient_mapping.csv", quote = F, row.names = F, col.names = F, sep = ",")

######################## Verfication of single mapping
patientMapping %>% group_by(patient1) %>% count() %>% filter(n > 1)

######################## Verfication of known duplicates
knownDuplicates = data.frame(sampleId = amberSamples) %>% filter(grepl('TII', sampleId)) %>% mutate(earlierSample = substr(sampleId, 1, 13))
missingKnownDuplicates = knownDuplicates %>% filter(!sampleId %in% sampleMapping$sample1 & !sampleId %in% sampleMapping$sample2 & earlierSample %in% amberBafs$sampleId)
knownDuplicates2 = data.frame(sampleId = amberSamples) %>% filter(grepl('TIII', sampleId)) %>% mutate(earlierSample = substr(sampleId, 1, 14))
missingKnownDuplicates2 = knownDuplicates2 %>% filter(!sampleId %in% sampleMapping$sample1 & !sampleId %in% sampleMapping$sample2 & earlierSample %in% amberBafs$sampleId)

####################### Verfification of JAN14
jan14Samples = read.table("~/hmf/resources/idgenerator/samples.csv", header = F, sep  =",")
jan14AmberBafs = amberBafs %>% filter(sampleId %in% jan14Samples$V1)
oldJan14PatientMapping = read.table("~/hmf/resources/idgenerator/patient_mapping.csv", header = F, sep  =",") %>% mutate(old = T)
newJan14SampleMapping = sample_mapping(jan14AmberBafs)
newJan14PatientMapping = patient_mapping(newJan14SampleMapping) %>% mutate(new = T)
compare = full_join(newJan14PatientMapping, oldJan14PatientMapping, by = c("patient1" = "V1", "patient2" = "V2"))
