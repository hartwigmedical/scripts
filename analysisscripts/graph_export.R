# Outputs reduced
library(RMySQL)
library(tidyverse)
library(Biostrings)
library(StructuralVariantAnnotation)
library(readr)
library(stringr)
tmpwd=getwd()
setwd("../gridss/")
source("libgridss.R")
setwd(tmpwd)
remove(tmpwd)
source("libSvAnalyser.R")

graph_dir = "../../graphs/"

samples = str_replace(list.files(path=graph_dir, pattern=".*_cn_reduced.purple.cnv"), stringr::fixed("_cn_reduced.purple.cnv"), "")

export_to_visjs = function(sample, graph_dir, cn_file_suffix="_cn_reduced.purple.cnv", sv_file_suffix="_sv_remaining.csv.sv") {
  cndf = read_tsv(paste0(graph_dir, "/", sample, cn_file_suffix), na = c("NA", "null"))
  svdf = read_tsv(paste0(graph_dir, "/", sample, sv_file_suffix), na = c("NA", "null"),
    col_names = c("id",	"startChromosome",	"endChromosome",	"startPosition",	"endPosition",	"startOrientation",	"endOrientation",	"startHomologySequence",	"endHomologySequence",	"startAF",	"endAF",	"ploidy",	"adjustedStartAF",	"adjustedEndAF",	"adjustedStartCopyNumber",	"adjustedEndCopyNumber",	"adjustedStartCopyNumberChange",	"adjustedEndCopyNumberChange",	"insertSequence",	"type",	"filter",	"imprecise",	"qualScore",	"event",	"startTumourVariantFragmentCount",	"startTumourReferenceFragmentCount",	"startNormalVariantFragmentCount",	"startNormalReferenceFragmentCount",	"endTumourVariantFragmentCount",	"endTumourReferenceFragmentCount",	"endNormalVariantFragmentCount",	"endNormalReferenceFragmentCount",	"startIntervalOffsetStart",	"startIntervalOffsetEnd",	"endIntervalOffsetStart",	"endIntervalOffsetEnd",	"inexactHomologyOffsetStart",	"inexactHomologyOffsetEnd",	"startLinkedBy",	"endLinkedBy"),
    col_types=cols(
      id = col_character(),
      startChromosome = col_character(),
      endChromosome  = col_character(),
      startPosition = col_integer(),
      endPosition =  col_integer(),
      startOrientation = col_integer(),
      endOrientation = col_integer(),
      startHomologySequence = col_character(),
      endHomologySequence = col_character(),
      startAF = col_double(),
      endAF = col_double(),
      ploidy = col_double(),
      adjustedStartAF = col_double(),
      adjustedEndAF = col_double(),
      adjustedStartCopyNumber = col_double(),
      adjustedEndCopyNumber = col_double(),
      adjustedStartCopyNumberChange = col_double(),
      adjustedEndCopyNumberChange = col_double(),
      insertSequence = col_character(),
      type = col_character(),
      filter = col_character(),
      imprecise = col_logical(),
      qualScore = col_double(),
      event = col_character(),
      startTumourVariantFragmentCount = col_integer(),
      startTumourReferenceFragmentCount = col_integer(),
      startNormalVariantFragmentCount = col_integer(),
      startNormalReferenceFragmentCount = col_integer(),
      endTumourVariantFragmentCount = col_integer(),
      endTumourReferenceFragmentCount = col_integer(),
      endNormalVariantFragmentCount = col_integer(),
      endNormalReferenceFragmentCount = col_integer(),
      startIntervalOffsetStart = col_integer(),
      startIntervalOffsetEnd = col_integer(),
      endIntervalOffsetStart = col_integer(),
      endIntervalOffsetEnd = col_integer(),
      inexactHomologyOffsetStart = col_integer(),
      inexactHomologyOffsetEnd = col_integer(),
      startLinkedBy = col_character(),
      endLinkedBy = col_character()
    ),
    comment = "#")
  names(cndf) = str_replace(names(cndf), "#", "")
  cndf = cndf %>% mutate(
    id = paste0("Unknown", row_number()),
    sampleId = sample,
    observedBaf = observedBAF,
    actualBaf = actualBAF,
    copyNumberMethod = method)
  svdf = svdf %>% mutate(
    sampleId = sample,
    vcfId = id)
  cngr = to_cn_gr(cndf)
  svgr = to_sv_gr(svdf, include.homology=FALSE)
  svgr$cnid = annotate_sv_with_cnv_id(cngr, svgr, maxgap=1000)
  if (any(is.na(svgr$cnid))) {
    bad_svid = svgr$id[is.na(svgr$cnid)]
    warning(paste("Missing matching segment boundary for", length(bad_svid), "SVs. Removing from analysis. TODO: confirm these are all ALT contigs, missing sex or missing short arm events"))
    svdf = svdf %>% filter(!(id %in% bad_svid))
    svgr = svgr[!(svgr$id %in% bad_svid)]
  }
  tmpwd=getwd()
  setwd(graph_dir)
  export_to_visNetwork(cndf, svdf, svgr, sample, file=paste0("breakpointgraph.", sample, "simplified.html"))
  setwd(tmpwd)
}

for (sample in samples) {
  export_to_visjs(sample, graph_dir)
}
