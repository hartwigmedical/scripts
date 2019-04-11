library(tidyverse)
library(GenomicRanges)
library(RMySQL)
library(Biostrings)
library(StructuralVariantAnnotation)


query_somatic_structuralVariants = function(dbConnect, table="structuralVariant") {
  query = paste(
    "SELECT id, sampleId, startChromosome, endChromosome, startPosition, endPosition, startOrientation, endOrientation, startHomologySequence, endHomologySequence, insertSequence, startIntervalOffsetStart,startIntervalOffsetEnd, endIntervalOffsetStart,endIntervalOffsetEnd, startLinkedBy,endLinkedBy, startAnchoringSupportDistance, endAnchoringSupportDistance",
    " FROM ", table,
    " WHERE filter = 'PASS'",
    sep = "")
  return (DBI::dbGetQuery(dbConnect, query) %>%
            mutate(id=as.character(id)))
}
sv_gr <- function(dbdf, include.homology=TRUE) {
  grs = GRanges(
    seqnames=dbdf$startChromosome,
    ranges=IRanges(start=dbdf$startPosition + ifelse(include.homology, ifelse(is.na(dbdf$startIntervalOffsetStart), 0, dbdf$startIntervalOffsetStart), 0),
                   end=dbdf$startPosition + ifelse(include.homology, ifelse(is.na(dbdf$startIntervalOffsetEnd), 0, dbdf$startIntervalOffsetEnd), 0)),
    strand=ifelse(dbdf$startOrientation == 1, "+", "-"),
    insertSequence=dbdf$insertSequence,
    partner=ifelse(is.na(dbdf$endChromosome), NA_character_, paste0(dbdf$id, "h")),
    id=dbdf$id,
    sampleid=dbdf$sampleId,
    anchorSupportDistance=dbdf$startAnchoringSupportDistance,
    beid=paste0(dbdf$id, ifelse(is.na(dbdf$endChromosome), "b",  "o")),
    linkedBy=dbdf$startLinkedBy)
  names(grs)=grs$beid
  dbdf = dbdf %>% filter(!is.na(endChromosome))
  rc_insert_sequence = dbdf$insertSequence
  rc_insert_sequence[!str_detect(rc_insert_sequence, "[^ACGTN]")] = as.character(reverseComplement(DNAStringSet(rc_insert_sequence[!str_detect(rc_insert_sequence, "[^ACGTN]")])))
  grh = GRanges(
    seqnames=dbdf$endChromosome,
    ranges=IRanges(start=dbdf$endPosition + ifelse(include.homology, ifelse(is.na(dbdf$endIntervalOffsetStart), 0, dbdf$endIntervalOffsetStart), 0),
                   end=dbdf$endPosition + ifelse(include.homology, ifelse(is.na(dbdf$endIntervalOffsetEnd), 0, dbdf$endIntervalOffsetEnd), 0)),
    strand=ifelse(dbdf$endOrientation == 1, "+", "-"),
    insertSequence=ifelse(dbdf$startOrientation != dbdf$endOrientation, dbdf$insertSequence, rc_insert_sequence),
    partner=paste0(dbdf$id, "o"),
    id=dbdf$id,
    sampleid=dbdf$sampleId,
    anchorSupportDistance=dbdf$endAnchoringSupportDistance,
    beid=paste0(dbdf$id, "h"),
    linkedBy=dbdf$endLinkedBy)
  names(grh)=grh$beid
  return(c(grs, grh))
}
db = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
dbdf = query_somatic_structuralVariants(db)
#save(dbdf, file = "d:/hartwig/anchorsupport.RData")
#dbdf = load(file = "d:/hartwig/anchorsupport.RData")
gr = sv_gr(dbdf)
#pgr = sgr[ifelse(is.na(sgr$partner), names(sgr), sgr$partner)]

asm_linked_breakends = as.data.frame(gr) %>%
  dplyr::select(sampleid, beid, linkedBy) %>%
  dplyr::mutate(linkedBy = str_split(as.character(linkedBy), stringr::fixed(","))) %>%
  tidyr::unnest(linkedBy) %>%
  group_by(sampleid) %>%
  dplyr::filter(str_detect(linkedBy, "^asm"))
asm_links = asm_linked_breakends %>%
  inner_join(asm_linked_breakends, by=c("sampleid", "linkedBy"), suffix = c("1", "2")) %>%
  filter(beid1 != beid2) %>%
  group_by(sampleid, beid1, beid2) %>%
  summarise(linkedBy=paste0(linkedBy, collapse=","))

maxdistance=2000
anchor_df = data.frame(sampleId=unique(gr$sampleid)) %>%
  group_by(sampleId) %>% do({
    sgr = gr[gr$sampleid==.$sampleId]
    hitdf = findOverlaps(sgr, sgr, maxgap=maxdistance, ignore.strand=TRUE) %>%
      as.data.frame() %>%
      filter(
        sgr$id[queryHits] != sgr$id[subjectHits],
        as.logical(strand(sgr)[queryHits] == "-"),
        as.logical(strand(sgr)[subjectHits] == "+"),
        start(sgr)[queryHits] <= start(sgr)[subjectHits]) %>%
      mutate(
        beid1=sgr$beid[queryHits],
        beid2=sgr$beid[subjectHits],
        distance=abs(start(sgr[queryHits])-start(sgr[subjectHits])),
        anchorSupportDistance1=sgr$anchorSupportDistance[queryHits],
        anchorSupportDistance2=sgr$anchorSupportDistance[subjectHits]) %>%
      dplyr::select(-queryHits, -subjectHits)
    hitdf
  })


















