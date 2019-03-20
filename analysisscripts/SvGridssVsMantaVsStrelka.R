library(tidyr)
library(dplyr)
library(RMySQL)
library(GenomicRanges)
library(vcfR)

######################## LOAD DATA
load_sv_file <- function(directory, fileName, gridss) {
  chromFactors = c(1:22, "X","Y")
  
  vcf <- read.vcfR(paste0(directory, fileName), verbose = F)
  vcf  = vcfR2tidy(vcf, single_frame = F)$fix 
  if (gridss) {
    vcf$MATEID = vcf$PARID
    vcf$INV3 = F
    vcf$INV5 = F
  }
  
  vcf = vcf %>% filter(is.na(FILTER) | FILTER %in% c(".", "PASS", "")) %>%
    mutate(CHROM = factor(CHROM, chromFactors)) %>%
    filter(!is.na(CHROM)) %>% 
    arrange(CHROM, POS)
  vcf.id = vcf %>% select(ID) %>% mutate(IDRank = row_number())
  vcf.mateid = vcf %>% select(ID = MATEID) %>% mutate(MATEIDRank = row_number())  
  
  vcf = full_join(vcf, vcf.id, by = "ID") %>% 
    left_join(vcf.mateid, by = "ID") %>% 
    mutate(primary = IDRank < MATEIDRank, CHROM = as.character(CHROM)) %>% 
    select(-IDRank, -MATEIDRank) 
  
  vcf$primary = ifelse(vcf$SVTYPE == "BND", vcf$primary, T)
  vcf$orientation = ifelse(grepl('^([A-Z]+)', vcf$ALT), 1, -1)
  
  vcf.start = vcf %>% filter(primary)
  vcf.end = vcf %>% filter(!primary)
  vcf.joined = left_join(vcf.start, vcf.end, by = c("ID" = "MATEID"), suffix = c("_start", '_end')) %>% 
    filter(primary_start) %>% 
    select(-primary_start, -primary_end) %>%
    select(startChromosome = CHROM_start, endChromosome = CHROM_end, startPosition = POS_start, endPosition = POS_end, startOrientation = orientation_start, endOrientation = orientation_end, everything()) 
  
  
  vcf.joined = vcf.joined %>% mutate(
    endChromosome = ifelse(SVTYPE_start == "BND", endChromosome, startChromosome),
    endPosition = ifelse(SVTYPE_start == "BND", endPosition, END_start),
    
    startOrientation = ifelse(SVTYPE_start == "DEL", 1, startOrientation),
    endOrientation = ifelse(SVTYPE_start == "DEL", -1, endOrientation),
    startOrientation = ifelse(SVTYPE_start == "INS", 1, startOrientation),
    endOrientation = ifelse(SVTYPE_start == "INS", -1, endOrientation),
    startOrientation = ifelse(SVTYPE_start == "DUP", -1, startOrientation),
    endOrientation = ifelse(SVTYPE_start == "DUP", 1, endOrientation),
    startOrientation = ifelse(SVTYPE_start == "INV" & INV3_start, 1, startOrientation),
    endOrientation = ifelse(SVTYPE_start == "INV" & INV3_start, 1, endOrientation),
    startOrientation = ifelse(SVTYPE_start == "INV" & INV5_start, -1, startOrientation),
    endOrientation = ifelse(SVTYPE_start == "INV" & INV5_start, -1, endOrientation)
  )
  
  return (vcf.joined)
}

cohort = data.frame(sampleId = c("CPCT02030461T", "CPCT02050327T", "CPCT02130091T", "CPCT02450014T", "CPCT02150016T", "DRUP01330008T", "CPCT02120143T", "CPCT02330102T", "CPCT02160052T", "CPCT02070386T", "CPCT02030516T", "DRUP01010096T", "CPCT02370037T"), stringsAsFactors = F)

query_indels <- function(dbConnect, cohort) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste0("select * from somaticVariant where type = 'INDEL' and sampleId in(", sampleIdString, ") and filter = 'PASS'")
  return (dbGetQuery(dbConnect, query))
}

query_svs <- function(dbConnect, cohort) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste0("select * from structuralVariant where sampleId in (", sampleIdString, ") and filter = 'PASS'")
  cat(query)
  return (dbGetQuery(dbConnect, query))
}

dbProd = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysisWrite")
dbDisconnect(dbProd)
rm(dbProd)

manta = data.frame()

for (sample in cohort$sampleId) {
  cat("Loading", sample, "\n")
  mantaSample = load_sv_file("/Users/jon/hmf/analysis/mantaVgridss/", paste0(sample, ".manta.vcf"), FALSE) %>% mutate(sampleId = sample, scope = "Private") %>% select(sampleId, scope, everything())
  manta = bind_rows(manta, mantaSample)
}
rm(mantaSample)

strelka = read.csv(file = "~/hmf/analysis/mantaVgridss/Indels.csv", stringsAsFactors = F) %>% mutate(length = abs(nchar(ref) - nchar(alt))) %>% filter(length >= 10, length <= 100)
gridss = read.csv(file = "~/hmf/analysis/mantaVgridss/Svs.csv", stringsAsFactors = F) 

save(manta, strelka, gridss, file = "/Users/jon/hmf/analysis/mantaVgridss/rawData.RData")




########################### VALIDATION ALGORITHM
indel_overlaps <- function(svs, indels, maxgap) {
  indelContig =  paste0(indels$chromosome, indels$sampleId)
  svContig =  paste0(svs$startChromosome, svs$sampleId)
  
  indelRanges = GRanges(indelContig, IRanges(indels$position, indels$position + indels$length))
  svRanges = GRanges(svContig, IRanges(svs$startPosition, svs$endPosition))
  ol = as.matrix(findOverlaps(svRanges, indelRanges, type="any", select="all", maxgap))
  
  return (data.frame(ol))
}


sv_overlaps <- function(query, subject, additional = 10) {
  require(tidyr)
  require(dplyr)
  require(GenomicRanges)
  
  attach_position_range <- function(svs) {
    svs$CIPOS_start = ifelse(is.na(svs$CIPOS_start), "-0,0",svs$CIPOS_start)
    svs$CIPOS_end= ifelse(is.na(svs$CIPOS_end), "-0,0",svs$CIPOS_end)
    svsStartPositionCI = data.frame( do.call( rbind, strsplit( svs$CIPOS_start, ',' ) ) ) 
    svsEndPositionCI = data.frame( do.call( rbind, strsplit( svs$CIPOS_start, ',' ) ) ) 
    names(svsStartPositionCI) <- c("startPositionCIBufferLeft", "startPositionCIBufferRight")
    names(svsEndPositionCI) <- c("endPositionCIBufferLeft", "endPositionCIBufferRight")
    svs = bind_cols(svs, svsStartPositionCI)
    svs = bind_cols(svs, svsEndPositionCI)
    svs$startPositionRangeStart = svs$startPosition - abs(as.integer(svs$startPositionCIBufferLeft)) - additional
    svs$startPositionRangeEnd = svs$startPosition + abs(as.integer(svs$startPositionCIBufferRight)) + additional
    svs$endPositionRangeStart = svs$endPosition - abs(as.integer(svs$endPositionCIBufferLeft))
    svs$endPositionRangeEnd = svs$endPosition + abs(as.integer(svs$endPositionCIBufferRight))
    
    svs = svs %>% select(-startPositionCIBufferLeft, -startPositionCIBufferRight, endPositionCIBufferLeft, endPositionCIBufferRight)
    return (svs)
  }
  
  query = attach_position_range(query)
  subject = attach_position_range(subject)
  
  queryStartRange <- GRanges(paste0(query$sampleId, query$startChromosome), IRanges(query$startPositionRangeStart, query$startPositionRangeEnd))
  subjectStartRange <- GRanges(paste0(subject$sampleId, subject$startChromosome), IRanges(subject$startPositionRangeStart, subject$startPositionRangeEnd))
  startOverlaps = data.frame(findOverlaps(queryStartRange, subjectStartRange, type="any", select="all", maxgap = -1))
  
  queryEndRange <- GRanges(paste0(query$sampleId, query$endChromosome), IRanges(query$endPositionRangeStart, query$endPositionRangeEnd))
  subjectEndRange <- GRanges(paste0(subject$sampleId, subject$endChromosome), IRanges(subject$endPositionRangeStart, subject$endPositionRangeEnd))
  endOverlaps = data.frame(findOverlaps(queryEndRange, subjectEndRange, type="any", select="all", maxgap = -1))
  
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

load(file = "/Users/jon/hmf/analysis/mantaVgridss/rawData.RData")
gridss = gridss %>% mutate(
  endChromosome = ifelse(endChromosome == "NULL", NA, endChromosome), 
  endPosition = ifelse(endPosition == "NULL", NA, endPosition), 
  endPosition = as.numeric(endPosition))

strelka$scope <- "Private"
manta$scope <- "Private"
gridss$matchStrelka <- F
gridss$matchManta <- F

########################### STRELKA DELS
strelkaIndex = strelka %>% mutate(originalIndex = row_number(), isEligible = nchar(ref) > nchar(alt)) %>% 
  filter(isEligible) %>% 
  mutate(delIndex = row_number()) %>% select(originalIndex, delIndex)
gridssIndex = gridss %>% mutate(length = endPosition - startPosition - 1, originalIndex = row_number(), isEligible =  type == "DEL" ) %>% 
  filter(isEligible, length >= 10, length <= 100) %>% 
  mutate(delIndex = row_number()) %>% select(originalIndex, delIndex)

strelkaOverlaps = indel_overlaps(gridss[gridssIndex$originalIndex, ], strelka[strelkaIndex$originalIndex, ], 20)
strelkaHits = strelkaIndex[strelkaOverlaps$subjectHits, "originalIndex"]
gridsssHits = gridssIndex[strelkaOverlaps$queryHits, "originalIndex"]

strelka[strelkaHits, "scope"] <- "Shared"
gridss[gridsssHits, "matchStrelka"] <- T

########################### STRELKA INSERTS
strelkaIndex = strelka %>% mutate(originalIndex = row_number(), isEligible = nchar(ref) < nchar(alt)) %>% 
  filter(isEligible) %>% 
  mutate(delIndex = row_number()) %>% select(originalIndex, delIndex)
gridssIndex = gridss %>% mutate(length = endPosition - startPosition + nchar(insertSequence), originalIndex = row_number(), isEligible =  type %in% c("DUP", "INS")) %>% 
  filter(isEligible, length >= 10, length <= 100) %>% 
  mutate(delIndex = row_number()) %>% select(originalIndex, delIndex)

strelkaOverlaps = indel_overlaps(gridss[gridssIndex$originalIndex, ], strelka[strelkaIndex$originalIndex, ], 20)
strelkaHits = strelkaIndex[strelkaOverlaps$subjectHits, "originalIndex"]
gridsssHits = gridssIndex[strelkaOverlaps$queryHits, "originalIndex"]

strelka[strelkaHits, "scope"] <- "Shared"
gridss[gridsssHits, "matchStrelka"] <- T

########################### GRIDDS
gridss = gridss %>% mutate(
  CIPOS_start = paste(startIntervalOffsetStart, startIntervalOffsetEnd, sep = ","),
  CIPOS_end = paste(endIntervalOffsetStart, endIntervalOffsetEnd, sep = ",")) 

gridssIndex = gridss %>% mutate(originalIndex = row_number()) %>% filter(!is.na(endChromosome)) %>% select(originalIndex)

mantaOverlaps = sv_overlaps(gridss[gridssIndex$originalIndex, ], manta, 10)
mantaHits = mantaOverlaps$subjectHits
gridsssHits = gridssIndex[mantaOverlaps$queryHits, "originalIndex"]

manta[mantaHits, "scope"] <- "Shared"
gridss[gridsssHits, "matchManta"] <- T

########################### SINGLES
singles_overlaps <- function(svs, indels, maxgap) {
  indelContig =  paste0(indels$chromosome, indels$sampleId, indels$orientation)
  svContig =  paste0(svs$startChromosome, svs$sampleId, svs$startOrientation)
  
  indelRanges = GRanges(indelContig, IRanges(indels$position, indels$position))
  svRanges = GRanges(svContig, IRanges(svs$startPosition, svs$startPosition))
  ol = as.matrix(findOverlaps(svRanges, indelRanges, type="any", select="all", maxgap))
  return (data.frame(ol))
}


mantaSinglesStart = manta %>% filter(scope == "Private") %>% select(sampleId, chromosome = startChromosome, position = startPosition,  orientation = startOrientation)
mantaSinglesEnd = manta %>% filter(scope == "Private") %>% select(sampleId, chromosome = endChromosome, position = endPosition,  orientation = endOrientation)
mantaSingles = bind_rows(mantaSinglesStart, mantaSinglesEnd)

gridssIndex = gridss %>% mutate(originalIndex = row_number()) %>% filter(type == 'SGL') %>% select(originalIndex)
singlesOverlap = singles_overlaps(gridss[gridssIndex$originalIndex, ], mantaSingles, 20)
gridsssHits = gridssIndex[singlesOverlap$queryHits, "originalIndex"]
gridss[gridsssHits, "matchManta"] <- T

gridss$scope <- "Private"
gridss$scope <- ifelse(gridss$matchManta, "SharedManta", gridss$scope)
gridss$scope <- ifelse(gridss$matchStrelka, "SharedStrelka", gridss$scope)
gridss$scope <- ifelse(gridss$matchManta & gridss$matchStrelka, "SharedBoth", gridss$scope)


summary = gridss %>% group_by(sampleId, type,  scope) %>% count() %>% spread(scope, n)
summary[is.na(summary)] <- 0
View(summary)




