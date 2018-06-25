#!/usr/bin/env Rscript
#
# Incorporates the given VCF files into the Panel Of Normals
#
library(tidyverse)
library(readr)
library(stringr)
usage = "Usage: Rscript add_to_pon.R <pon directory> <input VCFs>"
args = commandArgs(TRUE)
if (str_detect(args[1], "add_to_pon")) {
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
source("libgridss.R")
pon_dir = "C:/hartwig/pon"
vcf_list = c(
"C:/hartwig/down/COLO829R_COLO829T.gridss.vcf",
"C:/hartwig/down/CPCT02100013R_CPCT02100013T.gridss.vcf",
"C:/hartwig/down/CPCT02100013R_CPCT02100013TII.gridss.vcf")

full_bp = list()
full_be = list()
for (vcf_file in vcf_list) {
  sampleId = str_replace(basename(vcf_file), ".gridss.vcf", "")
  full_vcf = readVcf(vcf_file, "hg19")
  bpgr = breakpointRanges(full_vcf, unpartneredBreakends=FALSE)
  begr = breakpointRanges(full_vcf, unpartneredBreakends=TRUE)

  bpgr = bpgr[geno(full_vcf[bpgr$vcfId])$QUAL[,1] > gridss.min_qual / 2 | geno(full_vcf[bpgr$partner])$QUAL[,1] > gridss.min_qual / 2]
  begr = begr[geno(full_vcf[bpgr$vcfId])$QUAL[,1] > gridss.min_qual * gridss.single_breakend_multiplier / 2]

  minimal_bpgr = bpgr
  mcols(minimal_bpgr) = NULL
  minimal_bpgr$vcf = sampleId
  minimal_bpgr$IMPRECISE = info(full_vcf[names(minimal_bpgr)])$IMPRECISE
  names(minimal_bpgr) = paste(minimal_bpgr$vcf, names(minimal_bpgr), sep="_")
  minimal_bpgr$partner = paste(minimal_bpgr$vcf, bpgr$partner, sep="_")

  minimal_begr = begr
  mcols(minimal_begr) = NULL
  minimal_begr$vcf = sampleId
  minimal_begr$IMPRECISE = info(full_vcf[names(minimal_begr)])$IMPRECISE
  names(minimal_begr) = NULL

  full_bp[[sampleId]] = minimal_bpgr
  full_be[[sampleId]] = minimal_begr
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
  # name included purely forto ensure stable sort order
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

bedpe <- data.frame(
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
  IMPRECISE=ponbp$IMPRECISE
)
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
  arrange(IMPRECISE, desc(end - start), vcf, seqnames, start)
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
