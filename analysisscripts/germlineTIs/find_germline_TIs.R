#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)

library(argparser)
argp = arg_parser("Identifies germline templated insertions")
argp = add_argument(argp, "--input", help="Input GRIDSS VCF")
argp = add_argument(argp, "--output", help="Candidate templated insertions")
argp = add_argument(argp, "--normal_ordinal", type="integer", default=1, help="Ordinal of normal sample in the VCF")
argp = add_argument(argp, "--normal_qual", type="integer", default=350, help="Minimum QUAL score of the normal sample")
argp = add_argument(argp, "--insert_site_max_gap", type="integer", default=35, help="Maximum distance between insert site breakends")
argp = add_argument(argp, "--template_insert_max_length", type="integer", default=10000, help="Maximum size of templated insertion")
# argv = parse_args(argp, c("--input", "D:/hartwig/down/pre.vcf", "--output", "D:/hartwig/down/tis.vcf"))
argv = parse_args(argp)

library(StructuralVariantAnnotation)
library(VariantAnnotation)
library(tidyverse)

vcf = readVcf(argv$input)
vcf = vcf[geno(vcf)$QUAL[,argv$normal_ordinal] >= argv$normal_qual]
gr = breakpointRanges(vcf)
gr = gr[gr$partner %in% names(gr)]
# filter out small events
gr = gr[!(as.logical(seqnames(gr) == seqnames(partner(gr))) & abs(start(gr) - start(partner(gr))) <= max(argv$insert_site_max_gap, argv$template_insert_max_length))]
# insert site matches
pgr = partner(gr)
hits = findOverlaps(gr, gr, maxgap=argv$insert_site_max_gap, ignore.strand=TRUE) %>%
  as.data.frame() %>%
  filter(as.logical(
    strand(gr)[subjectHits] != strand(gr)[queryHits] &
    # partner breakends make a TI
    seqnames(pgr)[subjectHits] == seqnames(pgr)[queryHits] &
    strand(pgr)[subjectHits] == "-" & strand(pgr)[queryHits] == "+") &
    start(pgr)[subjectHits] < start(pgr)[queryHits]) %>%
  mutate(
    insdellen=ifelse(as.logical(strand(gr)[subjectHits] == "-"), 1, -1) * start(gr)[subjectHits] - start(gr)[queryHits],
    tilen = abs(start(pgr)[subjectHits] - start(pgr)[queryHits])) %>%
  filter(abs(insdellen) <= argv$insert_site_max_gap & tilen <= argv$template_insert_max_length) %>%
  # ignore local events
  filter(as.logical(seqnames(gr)[subjectHits] != seqnames(pgr)[subjectHits]) | abs(start(gr)[subjectHits] - start(pgr)[subjectHits]) > argv$template_insert_max_length)

out_vcf = vcf[unique(c(names(gr)[c(hits$queryHits, hits$subjectHits)], gr$partner[c(hits$queryHits, hits$subjectHits)]))]
writeVcf(out_vcf, argv$output)

