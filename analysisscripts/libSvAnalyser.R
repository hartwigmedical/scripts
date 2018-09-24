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
  return (DBI::dbGetQuery(dbConnect, query) %>%
    mutate(id=as.character(id)))
}
query_somatic_structuralVariants = function(dbConnect) {
  query = paste(
    "SELECT * ",
    " FROM structuralVariant ",
    " WHERE filter = 'PASS'",
    sep = "")
  return (DBI::dbGetQuery(dbConnect, query) %>%
    mutate(id=as.character(id)))
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
    imprecise=dbdf$imprecise != 0,
    event=dbdf$event,
    id=dbdf$id,
    vcfid=dbdf$vcfId,
    beid=paste0(dbdf$id, ifelse(is.na(dbdf$endChromosome), "b",  "o")),
    linkedBy=dbdf$startLinkedBy)
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
    insertSequence=ifelse(dbdf$startOrientation != dbdf$endOrientation, dbdf$insertSequence, as.character(reverseComplement(DNAStringSet(dbdf$insertSequence)))),
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
    imprecise=dbdf$imprecise != 0,
    event=dbdf$event,
    id=dbdf$id,
    vcfid=dbdf$vcfId,
    beid=paste0(dbdf$id, "h"),
    linkedBy=dbdf$endLinkedBy)
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

gridss_sv_links = function(svdf, svgr) {
  bind_rows(
    # gridss links
    svgr %>% as.data.frame() %>%
      dplyr::select(sampleId, id, beid, linkedBy) %>%
      filter(!is.na(linkedBy) & linkedBy != "." & linkedBy != "") %>%
      dplyr::mutate(linkedBy = str_split(as.character(linkedBy), stringr::fixed(","))) %>%
      tidyr::unnest(linkedBy) %>%
      group_by(sampleId, linkedBy) %>%
      filter(n() == 2) %>%
      arrange(beid) %>%
      summarise(
        id1=paste0(ifelse(row_number() == 1, id, ""), collapse=""),
        id2=paste0(ifelse(row_number() == 2, id, ""), collapse=""),
        beid1=paste0(ifelse(row_number() == 1, beid, ""), collapse=""),
        beid2=paste0(ifelse(row_number() == 2, beid, ""), collapse="")),
    # breakpoint partner links
    svdf %>%
      filter(!is.na(endChromosome)) %>%
      mutate(
        linkedBy=paste0("partner", id),
        id1=id,
        id2=id,
        beid1=paste0(id, "o"),
        beid2=paste0(id, "h")) %>%
      dplyr::select(sampleId, linkedBy, id1, id2, beid1, beid2)) %>%
    ungroup()
}

cluster_links = function(svdf, linked_breakends) {
  linkdf = linked_breakends %>% ungroup() %>% dplyr::select(id1, id2)
  cluster_id = 1:nrow(svdf)
  names(cluster_id) = svdf$id
  last_cluster_id = -1
  while (any(cluster_id != last_cluster_id)) {
    last_cluster_id = cluster_id
    # update cluster_id
    linkdf = linkdf %>%
      mutate(
        cid1 = cluster_id[id1],
        cid2 = cluster_id[id2],
        cid=pmin(cluster_id[id1], cluster_id[id2]))
    update_ids = bind_rows(
        linkdf %>% filter(cid1 > cid) %>% dplyr::select(id=id1, cid),
        linkdf %>% filter(cid2 > cid) %>% dplyr::select(id=id2, cid)) %>%
      group_by(id) %>%
      mutate(cid=min(cid))
    cluster_id[update_ids$id] = update_ids$cid
  }
  return(cluster_id)
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

hg19_centromeres = function() {
  GRanges(
    seqnames=c(1:22, "X", "Y"),
    ranges=IRanges(start=c(121535434, 92326171, 90504854, 49660117, 46405641, 58830166, 58054331, 43838887, 47367679, 39254935,
                           51644205, 34856694, 16000000, 16000000, 17000000, 35335801, 22263006, 15460898, 24681782, 26369569, 11288129, 13000000, 58632012, 10104553),
                   end=c(124535434, 95326171, 93504854, 52660117, 49405641, 61830166, 61054331, 46838887, 50367679, 42254935,
                         54644205, 37856694, 19000000, 19000000, 20000000, 38335801, 25263006, 18460898, 27681782, 29369569, 14288129, 16000000, 61632012, 13104553)))
}
hg19_primary_seqinfo = function() {
  require(R.cache)
  seqinfo = addMemoization(SeqinfoForUCSCGenome)("hg19")
  seqlevelsStyle(seqinfo) = "NCBI"
  seqinfo = seqinfo[c(1:22, "X", "Y")]
  return(seqinfo)
}
hg19_arms = function() {
  centomeres = hg19_centromeres()
  parm = GRanges(seqnames=seqnames(centomeres), IRanges(start=1, end=start(centomeres) - 1))
  qarm = GRanges(seqnames=seqnames(centomeres), IRanges(start=end(centomeres) + 1, end=seqlengths(hg19_primary_seqinfo())))
  centomeres$arm = paste0(seqnames(centomeres), "C")
  parm$arm = paste0(seqnames(parm), "P") # short arm, and lower position
  qarm$arm = paste0(seqnames(qarm), "Q")
  return(c(parm, qarm, centomeres))
}
on_hg19_arm = function(gr) {
  hg19_arms()$arm[findOverlaps(gr, hg19_arms(), select="first", ignore.strand=TRUE)]
}
line_elements = function()

cluster_consistency = function(svgr) {
  svgr$arm = on_hg19_arm(svgr)
  svgr %>% as.data.frame() %>%
    mutate(towardsCentromere = (str_detect(arm, "P") & strand == "+") | (str_detect(arm, "Q") & strand == "-")) %>%
    group_by(sampleId, cluster, arm) %>%
    # overall cluster info
    mutate(cluster_calls=length(unique(id)),
           cluster_breakend_call=sum(str_detect(beid, "b"))) %>%
    group_by(sampleId, cluster, arm, cluster_calls, cluster_breakend_call) %>%
    # arm-level consistency
    summarise(
      toward_centromere_count=sum(towardsCentromere),
      toward_telomere_count=n() - toward_centromere_count,
      toward_centromere_ploidy=sum(towardsCentromere * ploidy),
      toward_telomere_ploidy=sum(ploidy) - toward_centromere_ploidy)
}
load_line_elements = function(file) {
  linedf = readr::read_csv(file)
  linegr = with(linedf, GRanges(
    seqnames=Chromosome,
    ranges=IRanges(start=PosStart, end=PosEnd),
    type=Type))
}
load_fragile_sites = function(file) {
  df = readr::read_csv(file)
  with(df, GRanges(
    seqnames=Chromosome,
    ranges=IRanges(start=PosStart, end=PosEnd),
    type=Type,
    gene_name=GeneName,
    cfs_name=CFSName))
}
svids_of_overlap = function(gr, annotation_gr, maxgap=1, ignore.strand=TRUE) {
  unique(gr$id[overlapsAny(gr, annotation_gr, maxgap=maxgap, ignore.strand=ignore.strand)])
}

find_simple_events = function(cndf, svdf, svgr) {

}
add_prev_next_cnid = function(cndf) {
  cndf %>%
    group_by(sampleId, chromosome) %>%
    arrange(start) %>%
    mutate(
      next_cnid = lead(id),
      prev_cnid = lag(id)) %>%
    ungroup()
}
#    ------ CN1 -> CN3
#   /      \
# CN1  CN2  CN3
find_simple_deletions = function(cndf, svdf, svgr) {
  cndf = add_prev_next_cnid(cndf) %>% as.data.frame()
  svdf = svdf %>% as.data.frame()
  row.names(cndf) = cndf$id
  row.names(svdf) = svdf$id
  bpgr = with(mcols(svgr), svgr[!is.na(partner)])
  dels = bpgr[strand(bpgr) == "+" & strand(partner(bpgr)) == "-" &
      cndf[bpgr$cnid,]$next_cnid == cndf[partner(bpgr)$cnid,]$prev_cnid]
  data.frame(
      simple_event_type="DEL",
      svid=dels$id,
      left_flank_cnid=dels$cnid,
      cnid=cndf[dels$cnid,]$next_cnid,
      right_flank_cnid=partner(bpgr)[dels$partner]$cnid,
      stringsAsFactors=FALSE) %>%
    mutate(
      left_flank_ploidy=cndf[left_flank_cnid,]$copyNumber,
      right_flank_ploidy=cndf[right_flank_cnid,]$copyNumber,
      ploidy=cndf[cnid,]$copyNumber,
      svploidy=svdf[svid,]$ploidy) %>%
    mutate(
      flanking_ploidy_delta=left_flank_ploidy-right_flank_ploidy,
      ploidy_inconsistency_delta=(left_flank_ploidy + right_flank_ploidy) / 2 - svploidy - ploidy)
}
#    ------ CN2 -> CN2
#    \    /
# CN1  CN2  CN3
find_simple_duplications = function(cndf, svdf, svgr) {
  cndf = add_prev_next_cnid(cndf) %>% as.data.frame()
  svdf = svdf %>% as.data.frame()
  row.names(cndf) = cndf$id
  row.names(svdf) = svdf$id
  bpgr = with(mcols(svgr), svgr[!is.na(partner)])
  dups = bpgr[strand(bpgr) == "-" & strand(partner(bpgr)) == "+" & bpgr$cnid == partner(bpgr)$cnid]
  data.frame(
    simple_event_type="DUP",
    svid=dups$id,
    cnid=dups$cnid,
    left_flank_cnid=cndf[dups$cnid,]$prev_cnid,
    right_flank_cnid=cndf[dups$cnid,]$next_cnid,
    stringsAsFactors=FALSE) %>%
    mutate(
      left_flank_ploidy=cndf[left_flank_cnid,]$copyNumber,
      right_flank_ploidy=cndf[right_flank_cnid,]$copyNumber,
      ploidy=cndf[cnid,]$copyNumber,
      svploidy=svdf[svid,]$ploidy) %>%
    mutate(
      flanking_ploidy_delta=left_flank_ploidy-right_flank_ploidy,
      ploidy_inconsistency_delta=(left_flank_ploidy + right_flank_ploidy) / 2 + svploidy - ploidy)
}
#    a    b     c   d
#    -----------------
#    |    |     |    |
# CN1  CN2  CN3  CN4  CN5   (CN2 and CN4 can be 0bp in size)
find_simple_inversions = function(cndf, svdf, svgr, max.breakend.gap=35) {
  cndf = add_prev_next_cnid(cndf) %>% as.data.frame()
  svdf = svdf %>% as.data.frame()
  row.names(cndf) = cndf$id
  row.names(svdf) = svdf$id
  bpgr = svgr[!is.na(svgr$partner)]
  bpgr = bpgr[seqnames(bpgr) == seqnames(partner(bpgr)) & strand(bpgr) == strand(partner(bpgr))]
  findBreakpointOverlaps(bpgr, bpgr, maxgap=max.breakend.gap, ignore.strand=TRUE, sizemargin=NULL, restrictMarginToSizeMultiple=NULL) %>%
    filter(
      bpgr[queryHits]$sampleId == bpgr[subjectHits]$sampleId &
      as.logical(strand(bpgr[queryHits]) != strand(bpgr[subjectHits])) &
      end(bpgr[queryHits]) < start(partner(bpgr)[queryHits]) &
        end(bpgr[subjectHits]) < start(partner(bpgr)[subjectHits]) &
      start(bpgr[queryHits]) < start(bpgr[subjectHits])) %>%
    mutate(
      beida=names(bpgr[queryHits]),
      beidb=names(bpgr[subjectHits]),
      beidc=ifelse(start(partner(bpgr)[queryHits]) <= start(partner(bpgr)[subjectHits]), bpgr$partner[queryHits], bpgr$partner[subjectHits]),
      beidd=ifelse(start(partner(bpgr)[queryHits]) <= start(partner(bpgr)[subjectHits]), bpgr$partner[subjectHits], bpgr$partner[queryHits])) %>%
    mutate(
      cnid_left_a=ifelse(strand(bpgr[beida]) == "+", bpgr[beida]$cnid, cndf[bpgr[beida]$cnid,]$prev_cnid),
      cnid_right_a=ifelse(strand(bpgr[beida]) == "-", bpgr[beida]$cnid, cndf[bpgr[beida]$cnid,]$next_cnid),

      cnid_left_b=ifelse(strand(bpgr[beidb]) == "+", bpgr[beidb]$cnid, cndf[bpgr[beidb]$cnid,]$prev_cnid),
      cnid_right_b=ifelse(strand(bpgr[beidb]) == "-", bpgr[beidb]$cnid, cndf[bpgr[beidb]$cnid,]$next_cnid),

      cnid_left_c=ifelse(strand(bpgr[beidc]) == "+", bpgr[beidc]$cnid, cndf[bpgr[beidc]$cnid,]$prev_cnid),
      cnid_right_c=ifelse(strand(bpgr[beidc]) == "-", bpgr[beidc]$cnid, cndf[bpgr[beidc]$cnid,]$next_cnid),

      cnid_left_d=ifelse(strand(bpgr[beidd]) == "+", bpgr[beidd]$cnid, cndf[bpgr[beidd]$cnid,]$prev_cnid),
      cnid_right_d=ifelse(strand(bpgr[beidd]) == "-", bpgr[beidd]$cnid, cndf[bpgr[beidd]$cnid,]$next_cnid)) %>%
    mutate(
      # simple inversions have nothing happening within the inversion
      # ie: CN3 = right(b) == left(c)
      simple_event_type="INV",
      is_simple_inversion=cnid_right_b==cnid_left_c,
      left_flank_cnid=cnid_left_a,
      left_overlap_cnid=ifelse(cnid_right_a == cnid_left_b, cnid_right_a, NA_character_),
      cnid=ifelse(cnid_right_b==cnid_left_c, cnid_right_b, NA_character_),
      right_overlap_cnid=ifelse(cnid_right_c == cnid_left_d, cnid_right_c, NA_character_),
      right_flank_cnid=cnid_right_d) %>%
    mutate(
      left_flank_ploidy=cndf[left_flank_cnid,]$copyNumber,
      left_overlap_ploidy=cndf[left_overlap_cnid,]$copyNumber,
      ploidy=cndf[cnid,]$copyNumber,
      right_overlap_ploidy=cndf[right_overlap_cnid,]$copyNumber,
      right_flank_ploidy=cndf[right_flank_cnid,]$copyNumber,
      svploidy_a=bpgr[beida]$ploidy,
      svploidy_b=bpgr[beidb]$ploidy) %>%
    mutate(
      sv_delta = svploidy_a - svploidy_b,
      flanking_ploidy_delta=left_flank_ploidy-right_flank_ploidy,
      ploidy_left_flank_delta=left_flank_ploidy - ploidy,
      ploidy_right_flank_delta=left_flank_ploidy - ploidy,
      ploidy_sv_delta = (svploidy_a + svploidy_b) / 2 - ploidy)
}
chain_events = function(links) {
  #chains = data.frame(
   # beid=unique(links$beid
    #index=1
  stop("TODO")
}



