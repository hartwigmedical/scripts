#!/usr/bin/env Rscript
#
# Incorporates the given VCF files into the Panel Of Normals
#
library(tidyverse)
library(stringr)
library(rtracklayer)
library(R.cache)
source("libgridss.R")

min_pon_normal_qual = 150

usage = "Usage: Rscript create_gridss_pon.R <pon directory> <input VCFs>"
args = commandArgs(TRUE)
if (str_detect(args[1], "create_gridss_pon")) {
  args = args[-1]
}
if (length(args) < 2) {
  write(usage, stderr())
  q(save="no", status=1)
}
pon_dir = args[1]
args = args[-1]
vcf_list = args
if (!dir.exists(pon_dir)) {
  write("PON directory not found", stderr())
  write(usage, stderr())
  q(save="no", status=1)
}
if (!any(file.exists(vcf_list))) {
  write("VCF input file not found", stderr())
  write(usage, stderr())
  q(save="no", status=1)
}
#pon_dir = "C:/hartwig/pon"
#vcf_list = c(
#"C:/hartwig/down/COLO829R_COLO829T.gridss.vcf",
#"C:/hartwig/down/CPCT02100013R_CPCT02100013T.gridss.vcf",
#"C:/hartwig/down/CPCT02100013R_CPCT02100013TII.gridss.vcf")

setCacheRootPath(paste0(pon_dir, "/Rcache"))
options("R.cache.compress"=TRUE)
load_germline_pon_calls = addMemoization(function(vcf_file, sampleId) {
  write(paste("Start load", vcf_file), stderr())
  full_vcf = readVcf(vcf_file, "hg19")
  bpgr = breakpointRanges(full_vcf, unpartneredBreakends=FALSE)
  begr = breakpointRanges(full_vcf, unpartneredBreakends=TRUE)

  bpgr = bpgr[geno(full_vcf[bpgr$vcfId])$QUAL[,1] > min_pon_normal_qual | geno(full_vcf[bpgr$partner])$QUAL[,1] > min_pon_normal_qual]
  begr = begr[geno(full_vcf[begr$vcfId])$BQ[,1] > min_pon_normal_qual * gridss.single_breakend_multiplier]

  minimal_bpgr = bpgr
  mcols(minimal_bpgr) = NULL
  minimal_bpgr$vcf = rep(sampleId, length(bpgr))
  minimal_bpgr$IMPRECISE = info(full_vcf[names(minimal_bpgr)])$IMPRECISE
  names(minimal_bpgr) = paste(minimal_bpgr$vcf, names(minimal_bpgr), sep="_")
  minimal_bpgr$partner = paste(minimal_bpgr$vcf, bpgr$partner, sep="_")

  minimal_begr = begr
  mcols(minimal_begr) = NULL
  minimal_begr$vcf = rep(sampleId, length(begr))
  minimal_begr$IMPRECISE = info(full_vcf[names(minimal_begr)])$IMPRECISE
  names(minimal_begr) = NULL
  write(paste("End load", vcf_file), stderr())
  return(list(bp=minimal_bpgr, be=minimal_begr))
})
full_bp = list()
full_be = list()
for (vcf_file in vcf_list) {
  sampleId = str_replace(basename(vcf_file), ".gridss.vcf", "")
  calls = load_germline_pon_calls(vcf_file, sampleId)
  full_bp[[sampleId]] = calls$bp
  full_be[[sampleId]] = calls$be
}

bpdf = bind_rows(lapply(full_bp, function(x) {
  data.frame(
    seqnames = seqnames(x),
    start = start(x),
    end = end(x),
    strand = strand(x),
    vcf = x$vcf,
    name = names(x),
    partner = x$partner,
    IMPRECISE = x$IMPRECISE
  )})) %>%
  # preferentially call the precise call with the largest homology
  # name included purely to ensure stable sort order
  arrange(IMPRECISE, desc(end - start), name)
bpgr = GRanges(
  seqnames=bpdf$seqnames,
  ranges=IRanges(
    start=bpdf$start,
    end=bpdf$end),
  strand=bpdf$strand,
  vcf=bpdf$vcf,
  partner=bpdf$partner,
  IMPRECISE = bpdf$IMPRECISE)
names(bpgr) = bpdf$name
# PON requires two or more hits
hits = findBreakpointOverlaps(bpgr, bpgr) %>%
  filter(queryHits < subjectHits) %>%
  filter(bpdf$vcf[queryHits] != bpdf$vcf[subjectHits]) %>%
  group_by(queryHits) %>%
  summarise(n=n()) %>%
  filter(n >= 2)
ponbp = bpgr[hits$queryHits]
ponbp$hits = hits$n
bedpe = data.frame(
  chrom1=seqnames(ponbp),
  start1=start(ponbp) - 1,
  end1=end(ponbp),
  chrom2=seqnames(partner(bpgr)[hits$queryHits]),
  start2=start(partner(bpgr)[hits$queryHits]) - 1,
  end2=end(partner(bpgr)[hits$queryHits]),
  name=".",
  score=".",
  strand1=strand(ponbp),
  strand2=strand(partner(bpgr)[hits$queryHits]),
  IMPRECISE=ponbp$IMPRECISE,
  hits=ponbp$hits
) %>% filter(as.numeric(chrom1) < as.numeric(chrom2) | (chrom1 == chrom2 & start1 <= start2)) %>%
  arrange(chrom1, start1, chrom2, start2)
write.table(bedpe, paste(pon_dir, "gridss_pon_breakpoint.bedpe", sep="/"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)


bedf = bind_rows(lapply(full_be, function(x) {
  data.frame(
    seqnames = seqnames(x),
    start = start(x),
    end = end(x),
    strand = strand(x),
    vcf = x$vcf,
    IMPRECISE = x$IMPRECISE
  )})) %>%
  arrange(vcf, seqnames, start, desc(end - start))
begr = GRanges(
  seqnames=bedf$seqnames,
  ranges=IRanges(
    start=bedf$start,
    end=bedf$end),
  strand=bedf$strand,
  vcf=bedf$vcf,
  IMPRECISE = bedf$IMPRECISE)
behits = findOverlaps(begr, begr) %>%
  as.data.frame() %>%
  filter(queryHits < subjectHits) %>%
  filter(begr$vcf[queryHits] != begr$vcf[subjectHits]) %>%
  group_by(queryHits) %>%
  summarise(n=n()) %>%
  filter(n >= 2)
ponbe = begr[behits$queryHits]

export(ponbe, con=paste(pon_dir, "gridss_pon_single_breakend.bed", sep="/"), format="bed")
