library(tidyr)
library(dplyr)
library(GenomicRanges)
library(MutationalPatterns)
detach("package:purple", unload=TRUE)
library(purple)
library(ggplot2)
library(vcfR)
library(RMySQL)

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


sample = "CPCT02030461T"
mantaSample = load_sv_file("/Users/jon/hmf/analysis/mantaVgridss/", paste0(sample, ".manta.vcf"), FALSE) %>% mutate(sampleId = sample, scope = "Private") %>% select(sampleId, scope, everything())
gridssSample = load_sv_file("/Users/jon/hmf/analysis/mantaVgridss/", paste0(sample, ".purple.sv.ann.vcf.gz"), TRUE) %>% mutate(sampleId = sample, scope = "Private") %>% select(sampleId, scope, everything())  

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
databaseSample = dbGetQuery(dbProd, "select * from structuralVariant where sampleId = 'CPCT02030461T' and filter = 'PASS' and type != 'SGL';")
dbDisconnect(dbProd)
rm(dbProd)

common = inner_join(databaseSample, gridssSample, by = c("startChromosome", "startPosition", "startOrientation", "endChromosome", "endPosition", "endOrientation"))


overlaps = sv_overlaps(mantaSample, gridssSample, 10)
mantaSample[overlaps$queryHits, "scope"] <- "Shared"
gridssSample[overlaps$subjectHits, "scope"] <- "Shared"
View(mantaSample %>% group_by(scope,SVTYPE_start,sampleId) %>% count() %>% spread(scope,n))

directory = "/Users/jon/hmf/analysis/mantaVgridss/"
fileName = "CPCT02030461T.purple.sv.ann.vcf.gz"
gridss = T


cohort = c("CPCT02030461T", "CPCT02050327T", "CPCT02130091T", "CPCT02450014T", "CPCT02150016T", "DRUP01330008T", "CPCT02120143T", "CPCT02330102T", "CPCT02160052T", "CPCT02070386T", "CPCT02030516T", "DRUP01010096T", "CPCT02370037T")
mantaCohort = data.frame()
gridssCohort = data.frame()

for (sample in cohort) {
  cat("Loading", sample, "\n")
  mantaSample = load_sv_file("/Users/jon/hmf/analysis/mantaVgridss/", paste0(sample, ".manta.vcf"), FALSE) %>% mutate(sampleId = sample, scope = "Private") %>% select(sampleId, scope, everything())
  gridssSample = load_sv_file("/Users/jon/hmf/analysis/mantaVgridss/", paste0(sample, ".purple.sv.ann.vcf.gz"), TRUE) %>% mutate(sampleId = sample, scope = "Private") %>% select(sampleId, scope, everything())  
  mantaCohort = bind_rows(mantaCohort, mantaSample)
  gridssCohort = bind_rows(gridssCohort, gridssSample)
}
rm(mantaSample, gridssSample)
save(mantaCohort, gridssCohort, file = "/Users/jon/hmf/analysis/mantaVgridss/mantaVGridss.RData")

unique(mantaCohort$sampleId)
unique(gridssCohort$sampleId)
head(mantaCohort)

######################## FIND SHARED
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

load(file = "/Users/jon/hmf/analysis/mantaVgridss/mantaVGridss.RData")
overlaps = sv_overlaps(mantaCohort, gridssCohort, 10)
mantaCohort[overlaps$queryHits, "scope"] <- "Shared"
gridssCohort[overlaps$subjectHits, "scope"] <- "Shared"
save(mantaCohort, gridssCohort, file = "/Users/jon/hmf/analysis/mantaVgridss/mantaVGridss.RData")

View(mantaCohort %>% group_by(scope,SVTYPE_start,sampleId) %>% count() %>% spread(scope,n))
