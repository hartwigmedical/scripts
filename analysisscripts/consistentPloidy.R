detach("package:purple", unload=TRUE)
library(purple)
library(RMySQL)
library(tidyverse)
library(Biostrings)
library(StructuralVariantAnnotation)
library(testthat)

dbConn = dbConnect(MySQL(), dbname = "gridss_test")
sampleId = "COLO829T"

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

#' Determines copy number segments that should be the same
#'
#' Simple deletion:
equivalent_cn_segments = function(cnv_gr, sv_gr) {
  cnv_gr
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


cnv_gr = query_cnv_gr(dbConn, sampleId)
sv_gr = query_structural_variants_as_GRanges(dbConn, data.frame(sampleId=sampleId))
sv_gr$cnv_id = annotate_sv_with_cnv_id(cnv_gr, sv_gr, maxgap=1000)
if (any(is.na(sv_gr$cnv_id))) {
  stop("Missing CNV end point for SV")
}
induced_sv_gr = induced_edge_gr(cnv_gr)





