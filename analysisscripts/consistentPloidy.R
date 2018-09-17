library(purple)
library(RMySQL)
library(tidyverse)
library(Biostrings)
library(StructuralVariantAnnotation)
library(testthat)

dbConn = dbConnect(MySQL(), dbname = "hmfpatients", port=3307)
currentSampleId = "CPCT02140004T"
#currentSampleId = "COLO829T"


query_manta_structural_variants_as_GRanges <- function(dbConnect, cohort) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT id,sampleId,startChromosome,endChromosome,startPosition,endPosition,startOrientation,endOrientation,
    startHomologySequence, endHomologySequence, startAF, endAF, ploidy, adjustedStartAF, adjustedEndAF,adjustedStartCopyNumber,adjustedStartCopyNumberChange,adjustedEndCopyNumberChange,
    insertSequence,type,filter,somaticScore,mantaPrecise ",
    "FROM structuralVariant",
    "WHERE sampleId in (",sampleIdString, ")",
    sep = " ")
  dbdf = dbGetQuery(dbConnect, query)
  grs = GRanges(
    seqnames=dbdf$startChromosome,
    ranges=IRanges(start=dbdf$startPosition, width=1),
    strand=ifelse(dbdf$startOrientation == 1, "+", "-"),
    QUAL=NA_real_,
    FILTER=dbdf$filter,
    sampleId=dbdf$sampleId,
    ploidy=dbdf$ploidy,
    insertSequence=dbdf$insertSequence,
    type=dbdf$type,
    af=dbdf$startAF,
    homseq=dbdf$startHomologySequence,
    adjustedaf=dbdf$adjustedStartAF,
    adjustedcn=dbdf$adjustedStartCopyNumber,
    adjustedcn_delta=dbdf$adjustedStartCopyNumberChange,
    partner=ifelse(is.na(dbdf$endChromosome), NA_character_, paste0(dbdf$id, "h")),
    tumourVariantFragmentCount=NA,
    tumourReferenceFragmentCount=NA,
    normalVariantFragmentCount=NA,
    normalReferenceFragmentCount=NA,
    ihomlen=NA,
    somaticScore=dbdf$somaticScore,
    imprecise=dbdf$mantaPrecise =="false",
    event=NA,
    id=dbdf$id,
    vcfid=NA)
  names(grs)=paste0(dbdf$id, ifelse(is.na(dbdf$endChromosome), "b",  "o"))
  dbdf = dbdf %>% filter(!is.na(endChromosome))
  grh = GRanges(
    seqnames=dbdf$endChromosome,
    ranges=IRanges(start=dbdf$endPosition, width=1),
    strand=ifelse(dbdf$endOrientation == 1, "+", "-"),
    QUAL=NA_real_,
    FILTER=dbdf$filter,
    sampleId=dbdf$sampleId,
    ploidy=dbdf$ploidy,
    insertSequence=dbdf$insertSequence,
    type=dbdf$type,
    af=dbdf$endAF,
    homseq=dbdf$endHomologySequence,
    adjustedaf=dbdf$adjustedEndAF,
    adjustedcn=dbdf$adjustedEndCopyNumber,
    adjustedcn_delta=dbdf$adjustedEndCopyNumberChange,
    partner=paste0(dbdf$id, "o"),
    tumourVariantFragmentCount=NA,
    tumourReferenceFragmentCount=NA,
    normalVariantFragmentCount=NA,
    normalReferenceFragmentCount=NA,
    ihomlen=NA,
    somaticScore=dbdf$somaticScore,
    imprecise=dbdf$mantaPrecise =="false",
    event=NA,
    id=dbdf$id,
    vcfid=NA)
  names(grh)=paste0(dbdf$id, "h")
  return(c(grs, grh))
}

query_all_segments = function(dbConnect, sampleId) {
  query = paste(
    "SELECT * ",
    " FROM copyNumber ",
    "WHERE sampleId = '",sampleId, "'",
    sep = "")
  return (DBI::dbGetQuery(dbConnect, query))
}
query_cnv_gr = function(dbConnect, sampleId) {
  df = query_all_segments(dbConnect, sampleId)
  gr = with(df, GRanges(
    seqnames=chromosome,
    ranges=IRanges(
      start=start,
      end=end),
    strand="*",
    sampleId=sampleId,
    segmentStartSupport=segmentStartSupport,
    segmentEndSupport=segmentEndSupport,
    bafCount=bafCount,
    observedBaf=observedBaf,
    actualBaf=actualBaf,
    copyNumber=copyNumber,
    copyNumberMethod=copyNumberMethod
  ))
  names(gr) = df$id
  return(gr)
}
annotate_sv_with_cnv_id = function(cnv_gr, sv_gr, ...) {
  shits = as.data.frame(findOverlaps(query=sv_gr, subject=cnv_gr, type="start", select="all", ignore.strand=TRUE, ...)) %>%
    filter(as.logical(strand(sv_gr[queryHits]) == "-")) %>%
    mutate(distance=abs((start(sv_gr[queryHits]) + end(sv_gr[queryHits])) / 2 - start(cnv_gr[subjectHits]))) %>%
    # match to closest
    group_by(queryHits) %>%
    top_n(1, -distance) %>%
    ungroup()
  ehits = as.data.frame(findOverlaps(query=sv_gr, subject=cnv_gr, type="end", select="all", ignore.strand=TRUE, ...)) %>%
    filter(as.logical(strand(sv_gr[queryHits]) == "+")) %>%
    mutate(distance=abs((start(sv_gr[queryHits]) + end(sv_gr[queryHits])) / 2 - end(cnv_gr[subjectHits]))) %>%
    group_by(queryHits) %>%
    top_n(1, -distance) %>%
    ungroup()
  sv_gr$cnv_id = NA_character_
  sv_gr$cnv_id[shits$queryHits] = names(cnv_gr)[shits$subjectHits]
  sv_gr$cnv_id[ehits$queryHits] = names(cnv_gr)[ehits$subjectHits]
  return(sv_gr$cnv_id)
}
induced_edge_gr = function (cnv_gr, ...) {
  start_cnv_gr = cnv_gr
  start(start_cnv_gr) = start(start_cnv_gr) - 1
  width(start_cnv_gr) = 1
  hits = findOverlaps(cnv_gr, start_cnv_gr)
  induced_gr_left = GRanges(
    seqnames=seqnames(cnv_gr[queryHits(hits)]),
    ranges=IRanges(start=end(cnv_gr[queryHits(hits)]), width = 1),
    strand="+",
    id=paste0("ref", names(cnv_gr[queryHits(hits)])))
  names(induced_gr_left) = paste0("end_", names(cnv_gr[queryHits(hits)]))
  induced_gr_left$cnv_id = names(cnv_gr[queryHits(hits)])
  induced_gr_right = GRanges(
    seqnames=seqnames(cnv_gr[queryHits(hits)]),
    ranges=IRanges(start=end(cnv_gr[queryHits(hits)])+1, width = 1),
    strand="-",
    id=paste0("ref", names(cnv_gr[queryHits(hits)])))
  names(induced_gr_right) = paste0("start_", names(start_cnv_gr[subjectHits(hits)]))
  induced_gr_right$cnv_id = names(start_cnv_gr[subjectHits(hits)])
  induced_gr_left$partner = names(induced_gr_right)
  induced_gr_right$partner = names(induced_gr_left)
  return( c(induced_gr_left, induced_gr_right))
}

annotate_reference_fragment_count = function(induced_sv_gr, sv_gr) {
  as.data.frame(induced_sv_gr) %>% dplyr::select(strand, cnv_id) %>%
    left_join(as.data.frame(sv_gr) %>% dplyr::select(strand, cnv_id, tumourReferenceFragmentCount),
      by = c("strand", "cnv_id")) %>%
    distinct(strand, cnv_id, .keep_all=TRUE) %>%
    pull(tumourReferenceFragmentCount)
}


cnv_gr = query_cnv_gr(dbConn, currentSampleId)
sv_gr = query_manta_structural_variants_as_GRanges(dbConn, data.frame(sampleId=currentSampleId))
sv_gr$cnv_id = annotate_sv_with_cnv_id(cnv_gr, sv_gr, maxgap=1000)
if (any(is.na(sv_gr$cnv_id))) {
  stop("Missing CNV end point for SV")
}
induced_sv_gr = induced_edge_gr(cnv_gr)
induced_sv_gr$fragment_count = annotate_reference_fragment_count(induced_sv_gr, sv_gr)


puritydf = dbGetQuery(dbConn, paste0(
  "SELECT purity, normFactor, ploidy ",
  "FROM purity ",
  "WHERE sampleId ='",currentSampleId, "'"))

cndf = data.frame(
  seg_id=names(cnv_gr),
  # reverse the purple purity adjustment
  depth=(2 * (1-puritydf$purity) + puritydf$purity * cnv_gr$copyNumber),
  length=end(cnv_gr)-start(cnv_gr),
  baf=cnv_gr$observedBaf,
  baf_count=cnv_gr$bafCount,
  start=cnv_gr$segmentStartSupport,
  end=cnv_gr$segmentEndSupport,
  stringsAsFactors = FALSE)

svdf = data.frame(
  sv_id = as.character(sv_gr$id),
  be_id = names(sv_gr),
  be_id_partner = sv_gr$partner,
  orientation=as.character(strand(sv_gr)),
  seg_id = sv_gr$cnv_id,
  fragment_count = sv_gr$tumourReferenceFragmentCount,
  is_ref = FALSE,
  stringsAsFactors = FALSE) %>%
  bind_rows(data.frame(
    sv_id = induced_sv_gr$id,
    be_id = names(induced_sv_gr),
    be_id_partner = induced_sv_gr$partner,
    orientation=as.character(strand(induced_sv_gr)),
    seg_id = induced_sv_gr$cnv_id,
    fragment_count = induced_sv_gr$fragment_count,
    is_ref = TRUE,
    stringsAsFactors = FALSE)
  )
source("D:/dev/flowcnsv/libGurobi.R")
library(gurobi)
library(ggplot2)
results = list()
for (c1purity in seq(0, 1, 0.05)) {
  if (is.null(results[[paste0("purity", c1purity)]])) {
    model_filename = paste0(currentSampleId, ".model.lp")
    writeLPmodel_for_fixed_purity(model_filename, cndf, svdf, subclonal_purity=c(c1purity))
    model = gurobi_read(model_filename)
    results[[paste0("purity", c1purity)]] = run_model(model_filename, params=list(TimeLimit=600))
  }
}
resultdf = data.frame(
  name=names(results),
  purity=as.numeric(str_replace(names(results), "purity", "")),
  objval=unlist(lapply(results, function(x) ifelse(is.null(x$objval), NA_real_, x$objval))),
  cn_c1=unlist(lapply(results, function(x) ifelse(is.null(x$x), NA_real_, x$x[which(model$varnames == "avg_ploidy_c1")])))
)
ggplot(resultdf) + aes(x=purity, y=objval) + geom_point()
resultdf %>% group_by(purity) %>% filter(!is.na(objval)) %>% do({
  fit = results[[.$name]]
  cndf %>% mutate(
    c1_ploidy = fit$x[match(paste0("ploidy_seg", seg_id, "_c1_major"), model$varnames)]  +
      fit$x[match(paste0("ploidy_seg", seg_id, "_c1_minor"), model$varnames)],
    cn_delta = fit$x[match(paste0("cn_delta_seg_", seg_id), model$varnames)],
    start=cumsum(as.numeric(length)) - length,
    end=cumsum(as.numeric(length)))
}) %>%
  mutate(total_cn = 2 * (1 - purity) + purity * c1_ploidy) %>%
  mutate(fit_depth = total_cn * sum(as.numeric(depth * length)) / sum(as.numeric(total_cn * length))) %>%
  mutate(depth_delta = fit_depth - depth) %>%
  mutate(purple_copy_number=cnv_gr$copyNumber) %>%
ggplot() +
  aes(x=start) +
  geom_step(aes(y=purple_copy_number), colour="grey") +
  geom_step(aes(y=c1_ploidy)) +
  #geom_step(aes(y=cn_delta), colour="blue") +
  #geom_step(aes(y=fit_depth), colour="red") +
  #geom_step(aes(y=depth), colour="green") +
  coord_cartesian(ylim=c(0, 5)) +
  facet_wrap(~ purity, scales="free")
