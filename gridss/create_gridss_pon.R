#!/usr/bin/env Rscript
#
# Incorporates the given VCF files into the Panel Of Normals
#
library(argparser, quietly=TRUE)
argp = arg_parser("Filters a raw GRIDSS VCF into somatic call subsets.")
argp = add_argument(argp, "--pondir", default=NA, help="Directory to write PON to.")
argp = add_argument(argp, "--scriptdir", default=ifelse(sys.nframe() == 0, "./", dirname(sys.frame(1)$ofile)), help="Path to libgridss.R script")
argp = add_argument(argp, "--cache-only", flag=TRUE, help="Only generate cache objects for input files. Useful for parallel processing of inputs.")
argp = add_argument(argp, "--normalordinal", type="integer", default=1, help="Ordinal of normal sample in the VCF")
argp = add_argument(argp, "--input", nargs=Inf, help="Input VCFs normal")
# argv = parse_args(argp, argv=c("--pondir", "~/pon/", "--scriptdir", "/data/common/repos/scripts/gridss/", "--input", "a.vcf", "b.vcf"))
argv = parse_args(argp)

if (!dir.exists(argv$pondir)) {
  write("PON directory not found", stderr())
  write(usage, stderr())
  q(save="no", status=1)
}
vcf_list = argv$input
if (!any(file.exists(vcf_list))) {
  write("VCF input file not found", stderr())
  q(save="no", status=1)
}
libgridssfile = paste0(argv$scriptdir, "/", "libgridss.R")
if (file.exists(libgridssfile)) {
  tmpwd = getwd()
  setwd(argv$scriptdir)
  library(tidyverse, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(rtracklayer, quietly=TRUE)
  library(R.cache, quietly=TRUE)
  source("libgridss.R")
  setwd(tmpwd)
} else {
  msg = paste("Could not find libgridss.R in", argv$scriptdir, " - please specify a --scriptdir path to a directory containing the required scripts")
  write(msg, stderr())
  print(argp)
  stop(msg)
}

setCacheRootPath(paste0(argv$pondir, "/Rcache"))
options("R.cache.compress"=TRUE)
load_germline_pon_calls = function(vcf_file, sampleId) {
  vcf_file = normalizePath(vcf_file)
  key=list(vcf_file=vcf_file)
  cached = loadCache(key=key)
  if (!is.null(cached)) {
    write(paste("Using cached data for ", vcf_file), stderr())
    return(cached)
  }
  write(paste("Start load", vcf_file), stderr())
  full_vcf = readVcf(vcf_file, "hg19")
  bpgr = breakpointRanges(full_vcf, unpartneredBreakends=FALSE)
  begr = breakpointRanges(full_vcf, unpartneredBreakends=TRUE)

  bpgr = bpgr[geno(full_vcf[bpgr$vcfId])$QUAL[,argv$normalordinal] > gridss.pon.min_normal_qual | geno(full_vcf[bpgr$partner])$QUAL[,argv$normalordinal] > gridss.pon.min_normal_qual]
  begr = begr[geno(full_vcf[begr$vcfId])$BQ[,argv$normalordinal] > gridss.pon.min_normal_qual * gridss.single_breakend_multiplier]

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
  result = list(bp=minimal_bpgr, be=minimal_begr)
  saveCache(result, key=key)
  write(paste("End load", vcf_file), stderr())
  return(result)
}
full_bp = list()
full_be = list()
for (vcf_file in vcf_list) {
  sampleId = str_replace(str_replace(basename(vcf_file), ".gridss.vcf.gz", ""), ".gridss.vcf", "")
  calls = load_germline_pon_calls(vcf_file, sampleId)
  full_bp[[sampleId]] = calls$bp
  full_be[[sampleId]] = calls$be
}
if (argv$cache_only) {
  write("Caching complete ", stderr())
  q(save="no", status=0)
}

bpdf = bind_rows(lapply(full_bp, function(x) {
  data.frame(
    seqnames = as.character(seqnames(x)),
    start = start(x),
    end = end(x),
    strand = as.character(strand(x)),
    name = names(x),
    score=1,
    partner = x$partner,
    IMPRECISE = x$IMPRECISE,
    stringsAsFactors=FALSE
  )})) %>%
  # preferentially call the precise call with the largest homology
  # name included purely to ensure stable sort order
  arrange(IMPRECISE, desc(end - start), name)
bpgr = as(bpdf, "GRanges")
names(bpgr) = bpdf$name
hits = findBreakpointOverlaps(bpgr, bpgr) %>%
  group_by(queryHits) %>%
  summarise(n=n())
bpdf$score[hits$queryHits] = hits$n
bedpe = bpdf %>% mutate(
    chrom1=seqnames,
    start1=start,
    end1=end,
    chrom2=as.character(GenomeInfoDb::seqnames(partner(bpgr))),
    start2=BiocGenerics::start(partner(bpgr)),
    end2=BiocGenerics::end(partner(bpgr)),
    score=score,
    strand1=strand,
    strand2=as.character(BiocGenerics::strand(partner(bpgr))),
    IMPRECISE=IMPRECISE) %>%
  mutate(
    # switch from 1-based [ ] to BED 0-based [ ) indexing
    start1=start1 - 1,
    start2=start2 - 1) %>%
  group_by(chrom1, start1, end1, strand1, IMPRECISE, chrom2, start2, end2, strand2) %>%
  summarise(score=max(score)) %>%
  filter(score >= gridss.pon.min_samples) %>%
  filter(chrom1 < chrom2 | (chrom1 == chrom2 & start1 <= start2)) %>%
  arrange(chrom1, start1, chrom2, start2) %>%
  mutate(name=".") %>%
  dplyr::select(chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2, IMPRECISE)
withr::with_options(c(scipen = 10), write.table(bedpe, paste(argv$pondir, "gridss_pon_breakpoint.bedpe", sep="/"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE))

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
  strand=bedf$strand)
bedf$score = countOverlaps(begr, begr)
bedf = bedf %>%
  group_by(seqnames, start, end, strand, IMPRECISE) %>%
  summarise(score=max(score)) %>%
  filter(score >= gridss.pon.min_samples)
ponbe = GRanges(
    seqnames=bedf$seqnames,
    ranges=IRanges(
      start=bedf$start,
      end=bedf$end),
    strand=bedf$strand,
    score=bedf$score,
    IMPRECISE=bedf$IMPRECISE)
export(ponbe, con=paste(argv$pondir, "gridss_pon_single_breakend.bed", sep="/"), format="bed")
