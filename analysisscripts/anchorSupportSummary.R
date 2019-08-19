library(tidyverse)
library(GenomicRanges)
library(RMySQL)
library(Biostrings)
library(StructuralVariantAnnotation)
library(cowplot)


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
#load(file = "d:/hartwig/anchorsupport.RData")
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

maxdistance=1000
anchor_df = findOverlaps(gr, gr, maxgap=maxdistance, ignore.strand=TRUE) %>%
  as.data.frame() %>%
  filter(
    gr$id[queryHits] != gr$id[subjectHits],
    gr$sampleid[queryHits] == gr$sampleid[subjectHits],
    as.logical(strand(gr)[queryHits] == "-"),
    as.logical(strand(gr)[subjectHits] == "+"),
    start(gr)[queryHits] <= start(gr)[subjectHits]) %>%
  mutate(
    beid1=gr$beid[queryHits],
    beid2=gr$beid[subjectHits],
    distance=abs(start(gr[queryHits])-start(gr[subjectHits])),
    anchorSupportDistance1=gr$anchorSupportDistance[queryHits],
    anchorSupportDistance2=gr$anchorSupportDistance[subjectHits]) %>%
  dplyr::select(-queryHits, -subjectHits)


adj_df = anchor_df %>%
  group_by(beid1) %>%
  mutate(isClosest=distance==min(distance)) %>%
  ungroup() %>%
  left_join(asm_links, by=c("beid1"="beid1", "beid2"="beid2")) %>%
  mutate(is_asm_linked=!is.na(linkedBy)) %>%
  mutate(phasing=ifelse(is_asm_linked, "cis", ifelse(pmax(anchorSupportDistance1, anchorSupportDistance2) > distance + 10, "trans", "unphased")))


ggplot(adj_df %>% filter(is_asm_linked)) +
  aes(x=distance, fill=pmax(anchorSupportDistance1, anchorSupportDistance2) > distance + 10) +
  geom_histogram(bins=50) +
  labs(title="Assembly linked variants", fill="supportDistance > distance")

ggplot() +
  aes(x=distance, y=pmax(anchorSupportDistance1, anchorSupportDistance2), colour=is_asm_linked) +
  geom_point(data=adj_df %>% filter(!is_asm_linked) %>% sample_frac(size=0.1), colour="red", size=0.1) +
  geom_point(data=adj_df %>% filter(is_asm_linked), colour="blue", size=0.1) +
  scale_x_continuous(limits=c(0, 600))


ggplot() +
  aes(x=distance, y=pmax(anchorSupportDistance1, anchorSupportDistance2), colour=is_asm_linked) +
  geom_point(data=adj_df %>% filter(!is_asm_linked) %>% sample_frac(size=0.1), colour="red", size=0.1) +
  geom_point(data=adj_df %>% filter(is_asm_linked), colour="blue", size=0.1) +
  scale_x_continuous(limits=c(0, 600))


ggplot(adj_df %>% filter(isClosest | is_asm_linked) %>% filter(distance > 35)) +
  aes(x=distance, fill=phasing) +
  geom_histogram(bin=30) +
  scale_x_continuous(limits = c(50, 800))
# TODO: fix x axis labels

adj_df %>% group_by(beid1) %>%
  summarise(state=max(ifelse(phasing=="cis", 2, ifelse(phasing=="trans", 1, 0)))) %>%
  group_by(state) %>%
  summarise(n=n()) %>%
  mutate(percentage=n/length(gr))
# TODO: fix percentage as we're:
# - only counting the left side
# - not filtering to > 50bp

nearby_summary_df = gr %>% as.data.frame() %>% dplyr::select(beid, sampleid) %>%
  left_join(rbind(
        adj_df %>% filter(distance >= 50 & distance <= 1000 & str_sub(adj_df$beid1, end=8) != str_sub(adj_df$beid2, end=8)) %>% mutate(beid=beid1) %>% dplyr::select(beid, distance),
        adj_df %>% filter(distance >= 50 & distance <= 1000 & str_sub(adj_df$beid1, end=8) != str_sub(adj_df$beid2, end=8)) %>% mutate(beid=beid2) %>% dplyr::select(beid, distance)) %>%
      group_by(beid) %>% summarise(distance=min(distance)),
  by="beid") %>%
group_by(sampleid) %>%
  summarise(n=n(), hasNearby=sum(!is.na(distance)))

nearby_summary_df = nearby_summary_df %>%
  filter(hasNearby > 1) %>%
  group_by(sampleid, n, hasNearby) %>%
  do({data.frame(
      pvalue=prod(.$n - seq(1, .$hasNearby - 1)) * ((1000/3000000000) ** (.$hasNearby - 1)),
      log10pvalue=sum(log10(.$n - seq(1, .$hasNearby - 1))) + (.$hasNearby - 1) * log10(1000/3000000000)
  )})
sum(nearby_summary_df$log10pvalue)

