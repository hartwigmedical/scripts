library(tidyverse)

CN_ROUNDING= 0.2
CN_DIFF_MARGIN = 0.25
CN_CHANGE_MIN = 0.8
DB_MAX_LENGTH = 1000
MIN_LOH_CN = 0.5

query_all_copy_numer = function(dbConnect) {
  query = paste(
    "SELECT * ",
    " FROM copyNumber ",
    sep = "")
  return (DBI::dbGetQuery(dbConnect, query))
}
query_somatic_structuralVariants = function(dbConnect) {
  query = paste(
    "SELECT * ",
    " FROM structuralVariant ",
    " WHERE filter = 'PASS'",
    sep = "")
  return (DBI::dbGetQuery(dbConnect, query))
}
to_cn_gr = function(df) {
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
to_sv_gr <- function(svdf, include.homology=TRUE) {
  dbdf = svdf
  grs = GRanges(
    seqnames=dbdf$startChromosome,
    ranges=IRanges(start=dbdf$startPosition + ifelse(include.homology, ifelse(is.na(dbdf$startIntervalOffsetStart), 0, dbdf$startIntervalOffsetStart), 0),
                   end=dbdf$startPosition + ifelse(include.homology, ifelse(is.na(dbdf$startIntervalOffsetEnd), 0, dbdf$startIntervalOffsetEnd), 0)),
    strand=ifelse(dbdf$startOrientation == 1, "+", "-"),
    QUAL=dbdf$qualScore,
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
    tumourVariantFragmentCount=dbdf$startTumourVariantFragmentCount,
    tumourReferenceFragmentCount=dbdf$startTumourReferenceFragmentCount,
    normalVariantFragmentCount=dbdf$startNormalVariantFragmentCount,
    normalReferenceFragmentCount=dbdf$startNormalReferenceFragmentCount,
    ihomlen=dbdf$inexactHomologyOffsetEnd-dbdf$inexactHomologyOffsetStart,
    somaticScore=dbdf$somaticScore,
    imprecise=dbdf$imprecise != 0,
    event=dbdf$event,
    id=dbdf$id,
    vcfid=dbdf$vcfId,
    beid=paste0(dbdf$id, ifelse(is.na(dbdf$endChromosome), "b",  "o")))
  names(grs)=grs$beid
  dbdf = dbdf %>% filter(!is.na(endChromosome))
  grh = GRanges(
    seqnames=dbdf$endChromosome,
    ranges=IRanges(start=dbdf$endPosition + ifelse(include.homology, ifelse(is.na(dbdf$endIntervalOffsetStart), 0, dbdf$endIntervalOffsetStart), 0),
                   end=dbdf$endPosition + ifelse(include.homology, ifelse(is.na(dbdf$endIntervalOffsetEnd), 0, dbdf$endIntervalOffsetEnd), 0)),
    strand=ifelse(dbdf$endOrientation == 1, "+", "-"),
    QUAL=dbdf$qualScore,
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
    tumourVariantFragmentCount=dbdf$endTumourVariantFragmentCount,
    tumourReferenceFragmentCount=dbdf$endTumourReferenceFragmentCount,
    normalVariantFragmentCount=dbdf$endNormalVariantFragmentCount,
    normalReferenceFragmentCount=dbdf$endNormalReferenceFragmentCount,
    ihomlen=dbdf$inexactHomologyOffsetEnd-dbdf$inexactHomologyOffsetStart,
    somaticScore=dbdf$somaticScore,
    imprecise=dbdf$imprecise != 0,
    event=dbdf$event,
    id=dbdf$id,
    vcfid=dbdf$vcfId,
    beid=paste0(dbdf$id, "h"))
  names(grh)=grh$beid
  return(c(grs, grh))
}
annotate_sv_with_cnv_id = function(cnv_gr, sv_gr, ...) {
  shits = as.data.frame(findOverlaps(query=sv_gr, subject=cnv_gr, type="start", select="all", ignore.strand=TRUE, ...)) %>%
    filter(sv_gr$sampleId[queryHits] == cnv_gr$sampleId[subjectHits] &
             as.logical(strand(sv_gr[queryHits]) == "-")) %>%
    mutate(
      distance=abs((start(sv_gr[queryHits]) + end(sv_gr[queryHits])) / 2 - start(cnv_gr[subjectHits])),
      qual=sv_gr$QUAL[queryHits]) %>%
    # match the closest SV
    group_by(queryHits) %>%
    arrange(distance, -qual) %>%
    filter(row_number() == 1) %>%
    # only match the best hit
    #group_by(subjectHits) %>%
    #top_n(1, qual) %>%
    #filter(row_number() == 1) %>%
    ungroup()
  ehits = as.data.frame(findOverlaps(query=sv_gr, subject=cnv_gr, type="end", select="all", ignore.strand=TRUE, ...)) %>%
    filter(sv_gr$sampleId[queryHits] == cnv_gr$sampleId[subjectHits] &
             as.logical(strand(sv_gr[queryHits]) == "+")) %>%
    mutate(
      distance=abs((start(sv_gr[queryHits]) + end(sv_gr[queryHits])) / 2 - end(cnv_gr[subjectHits])),
      qual=sv_gr$QUAL[queryHits]) %>%
    group_by(queryHits) %>%
    arrange(distance, -qual) %>%
    filter(row_number() == 1) %>%
    #group_by(subjectHits) %>%
    #top_n(1, qual) %>%
    #filter(row_number() == 1) %>%
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

gridss_sv_links = function(svdf) {
  linkdf = svdf %>% dplyr::select(sampleId, id, linkedBy) %>%
    dplyr::filter(linkedBy != "" & linkedBy != "." ) %>%
    # remove suffix from transitive calls
    dplyr::mutate(linkedBy = str_replace_all(linkedBy, "o", "")) %>%
    dplyr::mutate(linkedBy = str_replace_all(linkedBy, "h", "")) %>%
    dplyr::mutate(linkedBy = str_split(as.character(linkedBy), stringr::fixed(","))) %>%
    tidyr::unnest(linkedBy)
  linkedsvs = linkdf %>% inner_join(linkdf, by=c("sampleId"="sampleId", "linkedBy"="linkedBy"), suffix=c("1", "2")) %>%
    filter(id1 != id2) %>%
    mutate(
      id1=as.character(id1),
      id2=as.character(id2))
  #pairwise = linkedsvs %>% group_by(sampleId, id1, id2) %>%
  #  summarise(linkedBy=paste(linkedBy, collapse=","))
  return(linkedsvs)
}

cluster_from_links = function(svdf, links) {
  id_cluster = seq_len(nrow(svdf))
  names(id_cluster) = svdf$id
  last_id_cluster = -1
  while (any(last_id_cluster != id_cluster)) {
    last_id_cluster = id_cluster
    id_cluster[links$id1] = pmin(id_cluster[links$id1], id_cluster[links$id2])
  }
  return(id_cluster)
}


find_cn_loh_links = function(cndf) {
  loh_bounds_df = cndf %>% group_by(sampleId, chromosome) %>%
    arrange(sampleId, chromosome, start) %>%
    # Q: why no min bafCount on the LOH segment in svanalyser?
    filter(bafCount > 0) %>%
    mutate(
      minorCN = (1 - actualBaf) * copyNumber,
      is_loh = minorCN < MIN_LOH_CN,
      is_loh_start_flank = !is_loh & lead(is_loh) & !is.na(lead(is_loh)),
      is_loh_end_flank = !is_loh & lag(is_loh) & !is.na(lag(is_loh))) %>%
    # just the bounding non-LOH segments - these should have the SV breakpoints
    filter(is_loh_start_flank | is_loh_end_flank) %>%
    # Remove start/end indicator if we don't pair up
    # this happens when a LOH continues to the  start/end of the chromosome
    mutate(
      is_loh_start_flank=is_loh_start_flank & lead(is_loh_end_flank) & !is.na(lead(is_loh_end_flank)),
      is_loh_end_flank=is_loh_end_flank & lag(is_loh_start_flank) & !is.na(lag(is_loh_start_flank)))
  loh_cn_pair_df = data.frame(
      cnid_start_flank = loh_bounds_df %>% filter(is_loh_start_flank) %>% pull(id) %>% as.character(),
      cnid_end_flank = loh_bounds_df %>% filter(is_loh_end_flank) %>% pull(id)  %>% as.character(),
      stringsAsFactors=TRUE) %>%
    mutate(linked_by=paste0("loh", row_number()))
  return(loh_cn_pair_df)
}
find_sv_loh_links = function(cndf, cngr, svgr, maxgap=1000, ...) {
  loh_cn_pair_df = find_cn_loh_links(cndf) %>%
    mutate(
      beid_start_flank = find_closest_sv_to_segment(cnid_start_flank, cngr, svgr, "end", maxgap=maxgap, ...),
      beid_end_flank = find_closest_sv_to_segment(cnid_end_flank, cngr, svgr, "start", maxgap=maxgap, ...)) %>%
    filter(is.na(beid_start_flank) | is.na(beid_end_flank) |
             # ifelse is just a placeholder so we don't index by NA
             svgr[ifelse(is.na(beid_start_flank), names(svgr)[1], beid_start_flank)]$partner !=
             names(svgr[ifelse(is.na(beid_end_flank), names(svgr)[1], beid_end_flank)]))
  return(loh_cn_pair_df)
}
find_closest_sv_to_segment = function(cnid, cngr, svgr, position, ...) {
  result = rep(NA, length(cnid))
  cngr = cngr[as.character(cnid)]
  expected_sv_strand = ifelse(position=="start", "-", "+")
  hits = findOverlaps(query=cngr, subject=svgr, type=position, select="all", ignore.strand=TRUE, ...) %>%
    as.data.frame() %>%
    filter(as.logical(strand(svgr)[subjectHits] == expected_sv_strand)) %>%
    mutate(
      distance=abs((start(svgr[subjectHits]) + end(svgr[subjectHits])) / 2 - ifelse(position=="start", start(cngr[queryHits]), end(cngr[queryHits]))),
      QUAL=svgr$QUAL[subjectHits]) %>%
    group_by(queryHits) %>%
    arrange(distance, -QUAL) %>%
    filter(row_number() == 1)
  result[hits$queryHits] = svgr$beid[hits$subjectHits]
  return(result)
}
loh_sv_links = function(cndf, svdf) {
  cndf = cndf %>%
    group_by(sampleId, chromosome) %>%
    arrange(start) %>%
    mutate(
      is_loh_start_flank = lead(is_loh_start),
      start_loh_id = lead(loh_id),
      is_loh_end_flank = lag(is_loh_end),
      end_loh_id = lag(loh_id))

  start_sv = cndf %>%
    filter(is_loh_start_flank) %>%
    inner_join(svdf %>%
      filter(strand=="+") %>%
      dplyr::select(sampleId, id, cnid, strand),
      by=c("id"="cnid", "sampleId"="sampleId"),
      suffix=c("", ".sv")) %>%
    mutate(bound="start_sv_id")

    bind_rows(cndf %>% dplyr::select(sampleId, id, loh_id, is_loh_end) %>%
      filter(is_loh_end) %>%
      inner_join(svdf %>% dplyr::select(sampleId, id, cnid, strand) %>% filter(strand=="+"), by=c("id"="cnid", "sampleId"="sampleId"), suffix=c("", ".sv")) %>%
      mutate(bound="end_sv_id")) %>%
    dplyr::select(sampleId, loh_id, sv_id=id.sv, bound) %>%
    group_by(sampleId, loh_id, bound) %>%
    spread(bound, sv_id)
}



