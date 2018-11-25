library(dplyr)
library(tidyr)
library(RMySQL)
library(GenomicRanges)

###### COLLECT DATA ########

dbNovaseq = dbConnect(MySQL(), dbname='novaseq_validation', groups="RAnalysisWrite")
dbNovaseqSomatics = dbGetQuery(dbNovaseq, "SELECT * from somaticVariant where sampleId like '180913novaseq%' and filter = 'PASS'")
dbNovaseqSVs= dbGetQuery(dbNovaseq, "SELECT * from structuralVariant where sampleId like '180913novaseq%' and filter = 'PASS'")
dbDisconnect(dbNovaseq)
rm(dbNovaseq)

samples = unique(dbNovaseqSomatics$sampleId)
samples = substr(samples, 15, 27)
sampleIdString = paste("'", samples, "'", collapse = ",", sep = "")

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
dbProdSomatics = dbGetQuery(dbProd, paste0("select * from somaticVariant WHERE sampleId in (",sampleIdString,") and filter = 'PASS'"))
dbProdSVs = dbGetQuery(dbProd, paste0("select * from structuralVariant WHERE sampleId in (",sampleIdString,") and filter = 'PASS'"))

save(dbNovaseqSomatics, dbNovaseqSVs, dbProdSomatics, dbProdSVs, file = "~/hmf/RData/NovaseqValidation.RData")

###### EXECUTE ALGO ########
sv_overlaps <- function(query, subject, maxgap = -1) {
  require(tidyr)
  require(dplyr)
  require(GenomicRanges)
  
  queryStartRange <- GRanges(query$startChromosome, IRanges(query$startPosition, query$startPosition))
  subjectStartRange <- GRanges(subject$startChromosome, IRanges(subject$startPosition, subject$startPosition))
  startOverlaps = data.frame(findOverlaps(queryStartRange, subjectStartRange, type="any", select="all", maxgap = maxgap))
  
  queryEndRange <- GRanges(query$endChromosome, IRanges(query$endPosition, query$endPosition))
  subjectEndRange <- GRanges(subject$endChromosome, IRanges(subject$endPosition, subject$endPosition))
  endOverlaps = data.frame(findOverlaps(queryEndRange, subjectEndRange, type="any", select="all", maxgap = maxgap))
  
  overlaps = inner_join(startOverlaps, endOverlaps, by = c("queryHits", "subjectHits"))
  
  overlapQueryData = query[overlaps$queryHits, ] %>%
    mutate(queryHits = overlaps$queryHits) %>%
    select(queryHits, sampleId, startChromosome, endChromosome, startPosition, endPosition, startOrientation, endOrientation)
  
  overlapSubjectData = subject[overlaps$subjectHits, ] %>%
    mutate(subjectHits = overlaps$subjectHits) %>%
    select(subjectHits, subjectStartPosition = startPosition, subjectEndPosition = endPosition, subjectStartOrientation = startOrientation, subjectEndOrientation = endOrientation)
  
  overlapsData = bind_cols(overlapQueryData, overlapSubjectData) %>%
    filter(startOrientation == subjectStartOrientation, endOrientation == subjectEndOrientation) %>%
    select(-subjectStartOrientation, -subjectEndOrientation) %>%
    mutate(startPositionDiff = abs(startPosition - subjectStartPosition), endPositionDiff = abs(endPosition - subjectEndPosition), positionDiff = startPositionDiff + endPositionDiff) %>%
    group_by(startChromosome, endChromosome, startPosition, endPosition,startOrientation,endOrientation) %>%
    top_n(1, -positionDiff) %>%
    group_by(queryHits) %>%
    top_n(1, -subjectHits)
  
  return (overlapsData %>% select(queryHits, subjectHits))
}

load(file = "~/hmf/RData/NovaseqValidation.RData")
dbNovaseqSomatics$sampleId = substr(dbNovaseqSomatics$sampleId, 15, 27)
dbNovaseqSVs$sampleId = substr(dbNovaseqSVs$sampleId, 15, 27)
dbNovaseqSomatics$source = "NovaSeq"
dbProdSomatics$source  = "HiSeq"
dbNovaseqSVs$source = "NovaSeq"
dbProdSVs$source  = "HiSeq"

combinedSomatics = bind_rows(dbNovaseqSomatics, dbProdSomatics) %>%
  group_by(sampleId, chromosome, position, ref, alt) %>%
  mutate(scope = ifelse(n() > 1, "Shared", source))

combinedSV = data.frame(stringsAsFactors = F)
for (sample in samples) {
  prodSampleSV = dbProdSVs %>% filter(sampleId == sample) %>% mutate(scope = source)
  novaSeqSampleSV = dbNovaseqSVs %>% filter(sampleId == sample)  %>% mutate(scope = source)
  
  overlaps = sv_overlaps(prodSampleSV, novaSeqSampleSV, maxgap = 100)
  prodSampleSV[overlaps$queryHits, "scope"] <- "Shared"
  novaSeqSampleSV[overlaps$subjectHits, "scope"] <- "Shared"
  combinedSV = bind_rows(combinedSV, prodSampleSV)
  combinedSV = bind_rows(combinedSV, novaSeqSampleSV)
}
rm(prodSampleSV, novaSeqSampleSV, overlaps)


######### Summarise
combinedSomaticsSummary = combinedSomatics %>%
  group_by(sampleId, type, source, scope) %>% 
  count() %>% 
  ungroup() %>%
  mutate(scope = ifelse(scope == 'Shared', 'Shared', "Private")) %>%
  unite(source, source, scope) %>%
  spread(source,n, fill = 0)
  

combinedSVSummary = combinedSV %>% 
  group_by(sampleId, type, source, scope) %>% 
  count() %>% 
  ungroup() %>%
  mutate(scope = ifelse(scope == 'Shared', 'Shared', "Private")) %>%
  unite(source, source, scope) %>%
  spread(source,n, fill = 0)

save(combinedSomaticsSummary, combinedSVSummary, file = "~/hmf/RData/NovaseqValidationSummary.RData")
