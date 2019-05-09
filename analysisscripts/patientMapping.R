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
    group_by(patient1, patient2) %>% summarise(matchPercent = max(matchPercent))
  
  return (patientMapping)
}

load(file = "~/hmf/analysis/amberancestry/amberBafs.RData")
amberSamples = unique(amberBafs$sampleId)
date()
sampleMapping = sample_mapping(amberBafs)
patientMapping = patient_mapping(sampleMapping)
date()

write.table(amberSamples, file = "~/hmf/resources/idgenerator/amber_samples.csv", quote = F, row.names = F, col.names = F, sep = ",")
write.table(patientMapping, file = "~/hmf/resources/idgenerator/amber_patient_mapping.csv", quote = F, row.names = F, col.names = F, sep = ",")

######################## Verfication of known duplicates
knownDuplicates = data.frame(sampleId = amberSamples) %>% filter(grepl('TII', sampleId)) %>% mutate(earlierSample = substr(sampleId, 1, 13))
missingKnownDuplicates = knownDuplicates %>% filter(!sampleId %in% sampleMapping$sample1 & !sampleId %in% sampleMapping$sample2 & earlierSample %in% amberBafs$sampleId)
knownDuplicates2 = data.frame(sampleId = amberSamples) %>% filter(grepl('TIII', sampleId)) %>% mutate(earlierSample = substr(sampleId, 1, 14))
missingKnownDuplicates2 = knownDuplicates2 %>% filter(!sampleId %in% sampleMapping$sample1 & !sampleId %in% sampleMapping$sample2 & earlierSample %in% amberBafs$sampleId)



