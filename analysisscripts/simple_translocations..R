library(purple)
detach("package:purple", unload=TRUE)
library(purple)
library(RMySQL)
library(tidyverse)
source("../gridss/libgridss.R")

db = dbConnect(MySQL(), dbname = "gridss_test")
#
# Detects simple insertion-like translocation events
#
maxgap=35 # manta uses 35bp

# svdf = dbGetQuery(db,"SELECT *FROM structuralVariant")

gr = query_structural_variants_as_GRanges(db, data.frame(sampleId=svdf$sampleId))
bpgr = gr[!is.na(gr$partner)]
begr = gr[is.na(gr$partner)]

# Insertion-like events
hits = findOverlaps(begr, begr, maxgap=maxgap, ignore.strand=TRUE) %>%
  as.data.frame() %>%
  filter(begr$sampleId[queryHits] == begr$sampleId[subjectHits]) %>% # intra-sample
  filter(as.logical(strand(begr)[queryHits] != strand(begr)[subjectHits])) %>% # opposite strand
  # matching pairs with the best qual
  mutate(qqual=begr$QUAL[queryHits], squal=begr$QUAL[subjectHits]) %>%
  group_by(queryHits) %>%
  top_n(1, squal) %>%
  group_by(subjectHits) %>%
  top_n(1, qqual) %>%
  ungroup() %>%
  filter(queryHits %in% subjectHits & subjectHits %in% queryHits)

begr$bestSimpleInsertionBreakendPartner = NA_character_
begr$bestSimpleInsertionBreakendPartner[hits$queryHits] = names(begr)[hits$queryHits]

hits = findOverlaps(bpgr, begr, maxgap=maxgap, ignore.strand=TRUE) %>%
  as.data.frame() %>%
  filter(bpgr$sampleId[queryHits] == begr$sampleId[subjectHits]) %>% # intra-sample
  filter(as.logical(strand(bpgr)[queryHits] != strand(begr)[subjectHits])) %>% # opposite strand
  # matching pairs with the best qual
  mutate(qqual=bpgr$QUAL[queryHits], squal=begr$QUAL[subjectHits]) %>%
  group_by(queryHits) %>%
  top_n(1, squal) %>%
  group_by(subjectHits) %>%
  top_n(1, qqual) %>%
  ungroup()

bpgr$bestPartialTranslocationPartner = NA_character_
begr$bestPartialTranslocationPartner = NA_character_
bpgr$bestPartialTranslocationPartner[hits$queryHits] = names(begr)[hits$subjectHits]
begr$bestPartialTranslocationPartner[hits$subjectHits] = names(bpgr)[hits$queryHits]


table(!is.na(begr$bestSimpleInsertionBreakendPartner), !is.na(begr$bestPartialTranslocationPartner))
table(!is.na(bpgr$bestPartialTranslocationPartner))


query_structural_variants_as_GRanges(db, data.frame(sampleId=svdf$sampleId))

