detach("package:purple", unload=TRUE)
library(purple)
library(RMySQL)
library(tidyverse)
library(Biostrings)
library(StructuralVariantAnnotation)
library(testthat)

dbConn = dbConnect(MySQL(), dbname = "gridss_test")
currentSampleId = "COLO829T"

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
    strand="+")
  names(induced_gr_left) = paste0("end_", names(cnv_gr[queryHits(hits)]))
  induced_gr_left$cnv_id = names(cnv_gr[queryHits(hits)])
  induced_gr_right = GRanges(
    seqnames=seqnames(cnv_gr[queryHits(hits)]),
    ranges=IRanges(start=end(cnv_gr[queryHits(hits)])+1, width = 1),
    strand="-")
  names(induced_gr_right) = paste0("start_", names(start_cnv_gr[subjectHits(hits)]))
  induced_gr_right$cnv_id = names(start_cnv_gr[subjectHits(hits)])
  induced_gr_left$partner = names(induced_gr_right)
  induced_gr_right$partner = names(induced_gr_left)
  return( c(induced_gr_left, induced_gr_right))
}

#' 1-based ordinal of the given allele specific ploidy for the given subclone
cnv_x_ordinal = function(cnv_gr, id, subclone_ordinal, total_subclones) {
  name_ordinal = match(id, names(cnv_gr))
  return(subclone_ordinal + (name_ordinal - 1) * total_subclones)
}
test_that("cnv_x_ordinal", {
  test_cnv = c("a", "b", "c")
  names(test_cnv) = test_cnv
  expect_that(cnv_x_ordinal(test_cnv, "a", 1, 3), equals(1))
  expect_that(cnv_x_ordinal(test_cnv, "a", 2, 3), equals(2))
  expect_that(cnv_x_ordinal(test_cnv, "a", 3, 3), equals(3))
  expect_that(cnv_x_ordinal(test_cnv, "b", 1, 3), equals(4))
  expect_that(cnv_x_ordinal(test_cnv, "b", 2, 3), equals(5))
  expect_that(cnv_x_ordinal(test_cnv, "b", 3, 3), equals(6))
  expect_that(cnv_x_ordinal(test_cnv, "c", 1, 3), equals(7))
  expect_that(cnv_x_ordinal(test_cnv, "c", 2, 3), equals(8))
})
#' 1-based ordinal of the given allele specific ploidy for the given subclone
sv_x_ordinal = function(cnv_gr, sv_gr, id, subclone_ordinal, total_subclones) {
  snv_start_offset = total_subclones * length(cnv_gr)
  name_ordinal = match(id, names(sv_gr))
  return(snv_start_offset - 1 + subclone_ordinal + (name_ordinal - 1) * total_subclones)
}

test_that("x_ordinal", {
  test_cnv = c("a", "b", "c")
  names(test_cnv) = test_cnv
  expect_that(cnv_x_ordinal(test_cnv, "a", 1, 3), equals(1))
  expect_that(cnv_x_ordinal(test_cnv, "a", 2, 3), equals(2))
  expect_that(cnv_x_ordinal(test_cnv, "a", 3, 3), equals(3))
  expect_that(cnv_x_ordinal(test_cnv, "b", 1, 3), equals(4))
  expect_that(cnv_x_ordinal(test_cnv, "b", 2, 3), equals(5))
  expect_that(cnv_x_ordinal(test_cnv, "b", 3, 3), equals(6))
  expect_that(cnv_x_ordinal(test_cnv, "c", 1, 3), equals(7))
  expect_that(cnv_x_ordinal(test_cnv, "c", 2, 3), equals(8))
  test_sv = c("aa", "bb", "cc", "dd")
  names(test_sv) = test_sv
  expect_that(sv_x_ordinal(test_cnv, test_sv, "aa", 1, 3), equals(9))
  expect_that(sv_x_ordinal(test_cnv, test_sv, "aa", 2, 3), equals(10))
  expect_that(sv_x_ordinal(test_cnv, test_sv, "aa", 3, 3), equals(11))
  expect_that(sv_x_ordinal(test_cnv, test_sv, "bb", 1, 3), equals(12))
  expect_that(sv_x_ordinal(test_cnv, test_sv, "bb", 2, 3), equals(13))
})

annotate_reference_fragment_count = function(induced_sv_gr, sv_gr) {
  as.data.frame(induced_sv_gr) %>% dplyr::select(strand, cnv_id) %>%
    left_join(as.data.frame(sv_gr) %>% dplyr::select(strand, cnv_id, tumourReferenceFragmentCount),
      by = c("strand", "cnv_id")) %>%
    distinct(strand, cnv_id, .keep_all=TRUE) %>%
    pull(tumourReferenceFragmentCount)
}


cnv_gr = query_cnv_gr(dbConn, currentSampleId)
sv_gr = query_structural_variants_as_GRanges(dbConn, data.frame(sampleId=currentSampleId))
sv_gr$cnv_id = annotate_sv_with_cnv_id(cnv_gr, sv_gr, maxgap=1000)
if (any(is.na(sv_gr$cnv_id))) {
  stop("Missing CNV end point for SV")
}
induced_sv_gr = induced_edge_gr(cnv_gr)
induced_sv_gr$fragment_count = annotate_reference_fragment_count(induced_sv_gr, sv_gr)


cndf = data.frame(
  seg_id=names(cnv_gr),
  # TODO: do we need to reverse the purity adjustment?
  depth=cnv_gr$copyNumber,
  length=end(cnv_gr)-start(cnv_gr),
  baf=cnv_gr$observedBaf,
  baf_count=cnv_gr$bafCount,
  start=cnv_gr$segmentStartSupport,
  end=cnv_gr$segmentEndSupport,
  stringsAsFactors = FALSE)

svdf = data.frame(
  sv_id = names(sv_gr),
  orientation=as.character(strand(sv_gr)),
  seg_id = sv_gr$cnv_id,
  partner = sv_gr$partner,
  fragment_count = sv_gr$tumourReferenceFragmentCount,
  stringsAsFactors = FALSE) %>%
  bind_rows(data.frame(
    sv_id = names(induced_sv_gr),
    orientation=as.character(strand(induced_sv_gr)),
    seg_id = induced_sv_gr$cnv_id,
    partner = induced_sv_gr$partner,
    fragment_count = induced_sv_gr$fragment_count,
    stringsAsFactors = FALSE)
  )

# Generate LP file
writeLPmodel = function(file, cndf, subclones=2) {
  model = list(
    clones = c("normal", paste0("c", seq(subclones))),
    alleles = c("J", "n"))
  model$ascn = expand.grid(allele=paste0("_", model$alleles), subclone = paste0("_", model$clones)) %>%
      mutate(ordinal = seq(nrow(.)))
  model$q_ascn = df_cross_product(model$ascn, model$ascn, suffix=c("", "2"))

  # BAF must be minor allele frequency
  cndf = cndf %>% mutate(baf = pmin(baf, 1 - baf))

  # add null placeholder
  cndf_with_null = data.frame(
      seg_id="NULL",
      depth=NA,
      length=NA,
      baf=NA,
      baf_count=NA,
      start="NULL",
      end="NULL",
      stringsAsFactors = FALSE) %>%
    bind_rows(cndf)
  svdf_with_null = bind_rows(svdf,
    cndf %>% dplyr::select(seg_id) %>% mutate(
      partner=paste0("null_to_start", seg_id),
      sv_id=paste0("start_to_null", seg_id),
      orientation="-",
      fragment_count=NA),
    cndf %>% dplyr::select(seg_id) %>% mutate(
      partner=paste0("start_to_null", seg_id),
      sv_id=paste0("null_to_start", seg_id),
      seg_id="NULL",
      orientation="*",
      fragment_count=NA),
    cndf %>% dplyr::select(seg_id) %>% mutate(
      partner=paste0("null_to_end", seg_id),
      sv_id=paste0("end_to_null", seg_id),
      orientation="+",
      fragment_count=NA),
    cndf %>% dplyr::select(seg_id) %>% mutate(
      partner=paste0("end_to_null", seg_id),
      sv_id=paste0("null_to_end", seg_id),
      seg_id="NULL",
      orientation="*",
      fragment_count=NA))

  # need normalisation multipliers between read depth and coverage
  # and fragment counts and coverage
  # that is, read counts to effective coverage
  # for each segment:
  #   normalisation * sum(purity * ploidy) = observed read depth
  # and we have a circular definition since the normalisation factor
  # actually depends on the solution estimate.
  # need to do an initial estimate: assume ploidy = 2 then re-estimate

  # KISS!
  "Minimise"
  # segment read depth error
  cndf %>% mutate(
      var = paste0("rd_delta_seg", seg_id),
      weight = sqrt(length)) %>%
    filter(!is.na(weight)) %>%
    mutate(
      eqn = paste(weight, var)
    ) %>%
    pull(eqn) %>%
    paste0(collapse = " + ")
  # segment BAF error
  # edge fragment count error
  "Subject To"
  # sum ploidy = 1
  paste0("total_ploidy: ", model$ascn %>%
    mutate(var=paste0("ploidy", subclone, allele)) %>%
    pull(var) %>%
    paste0(collapse = " + ") %>%
    paste(" - 1 = 0", collapse = ""))
  # total cn = sum purity * ploidy
  cndf %>% mutate(
      eqn = model$ascn %>%
        mutate(
          var1=paste0("purity", subclone),
          var2=paste0("ploidy_seg", "[SEG_ID]", subclone, allele)) %>%
        mutate(
          eqn=paste0(var1, " * ", var2)
        ) %>%
        pull(eqn) %>%
        paste0(collapse = " + ") %>%
        paste0(" - cn_[SEG_ID] = 0", collapse = "") %>%
        str_replace_all(stringr::fixed("[SEG_ID]"), seg_id)) %>%
    mutate(eqn=paste0("copy_number_seg", seg_id, ": ", eqn)) %>%
    pull(eqn)
  # sum edge ASCN in = segment ASCN (if start != telemore/centromere)
  cndf %>% filter(!(start %in% c("TELOMERE", "NULL"))) %>%
    left_join(svdf_with_null %>% filter(orientation=="-"), by=c("seg_id"), suffix=c("", ".sv")) %>%
    df_cross_product(model$ascn) %>%
    mutate(eqn = paste("ploidy_sv", sv_id, subclone, allele, sep="")) %>%
    group_by(seg_id, subclone, allele) %>%
    summarise(eqn = paste0(eqn, collapse=" + ")) %>%
    ungroup() %>%
    mutate(eqn = paste0(eqn, " - ", "ploidy_seg", seg_id, subclone, allele, " = 0"))
  # sum edge ASCN out = segment ASCN (if start != telemore/centromere)
  cndf %>% filter(!(start %in% c("TELOMERE", "NULL"))) %>%
    left_join(svdf_with_null %>% filter(orientation=="+"), by=c("seg_id"), suffix=c("", ".sv")) %>%
    df_cross_product(model$ascn) %>%
    mutate(eqn = paste("ploidy_sv", sv_id, subclone, allele, sep="")) %>%
    group_by(seg_id, subclone, allele) %>%
    summarise(eqn = paste0(eqn, collapse=" + ")) %>%
    ungroup() %>%
    mutate(eqn = paste0(eqn, " - ", "ploidy_seg", seg_id, subclone, allele, " = 0"))
  #TODO: SV_ID must be common for both breakends (need a second id or can we get away with 1 and no partner?)

  # NULL ASCNs are a special case since they don't have a direction
  # we could constrain the NULL copy number based on how confident we are that we have all/most SVs

  # TODO: since we can get centromeric swaps, we can instead constrain
  cndf %>% mutate(
    eqn =
  )

  # for each segment

  #   sum edge ASCN out = ASCN (if end != telemore/centromere)

  "Bounds"
  # ploidies <= 1
  # Variables:
  # total cn depth
  # total edge depth
  # add null edges
  # segment ASCN
  # edge ASCN
  # Variables
  # ploidy

  "Integers"
  # segment total CN
  # segment ASCN
  # edge ASCN

  # General Constraints

  # relative weightings
  cndf = cndf %>%
    mutate(
      read_depth_weight = sqrt(length),
      baf_weight = baf_count)

  "Minimize"
   cndf %>% mutate(
    objname = paste0("OBJ_rd_", name),
    weight = 1)

   # ploidies sum to 1
   paste(paste(paste0(c("normal", model$clones), sep="_ploidy"), collapse=" + "), "= 1")
}
df_cross_product = function(df1, df2, ...) {
  full_join(df1 %>% mutate(df_cross_product_placeholder=1), df2 %>% mutate(df_cross_product_placeholder=1), by="df_cross_product_placeholder", ...) %>%
    dplyr::select(-df_cross_product_placeholder)
}
segment_read_depth_objective_function = function(model, cndf, ascn) {
  # (depth - depth_hat) ^ 2
  # = depth ^ 2 - 2 * depth_hat * depth + (constant we can ignore)
  # depth = sum over i in subclones ploidy(i) * (depth(i, major) + depth(i, minor))
  # + normal_ploidy) * 2
  quartic_terms = cndf %>% top_n(1) %>%
    mutate(weight = read_depth_weight) %>%
    df_cross_product(model$q_ascn) %>%
    mutate(
      p1=paste0("ploidy", subclone),
      p2=paste0("ploidy", subclone2),
      cn1=paste0(name, subclone, allele),
      cn2=paste0(name, subclone2, allele2)) %>%
    dplyr::select(weight, p1, p2, cn1, cn2)


  quadradic_terms = cndf %>% mutate(weight = - 2 * read_depth_weight * depth) %>%
    df_cross_product(model$ascn) %>%
    mutate(
      ploidy=paste0("ploidy", subclone),
      cn=paste0(name, subclone, allele)) %>%
    dplyr::select(weight, ploidy, allele)
}

fileConn = file(paste0(currentSampleId, ".model.lp"))
writeLines(writeLPmodel(cndf, svdf), fileConn)
close(fileConn)

