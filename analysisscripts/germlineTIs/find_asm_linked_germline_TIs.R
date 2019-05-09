#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)

library(argparser)
argp = arg_parser("Identifies germline templated insertions")
argp = add_argument(argp, "--input", nargs=Inf, help="Input GRIDSS VCFs")
argp = add_argument(argp, "--output_vcf", nargs=Inf, help="Candidate templated insertions")
#argp = add_argument(argp, "--normal_ordinal", type="integer", default=1, help="Ordinal of normal sample in the VCF")
argp = add_argument(argp, "--qual", type="integer", default=350, help="Minimum QUAL score of variants")
argp = add_argument(argp, "--scriptdir", default=ifelse(sys.nframe() == 0, "./", dirname(sys.frame(1)$ofile)), help="Path to libgridss.R script")
#argv = parse_args(argp, c("--scriptdir", "D:/hartwig/scripts/gridss/", "--input", "D:/hartwig/down/CPCT02180037R_CPCT02180037T.gridss.vcf.gz", "--output_vcf", "D:/hartwig/down/test_ti.vcf"))
argv = parse_args(argp)

libgridssfile = paste0(argv$scriptdir, "/", "libgridss.R")
if (file.exists(libgridssfile)) {
  tmpwd = getwd()
  setwd(argv$scriptdir)
  source("libgridss.R")
  setwd(tmpwd)
} else {
  msg = paste("Could not find libgridss.R in", argv$scriptdir, " - please specify a --scriptdir path to a directory containing the required scripts")
  write(msg, stderr())
  print(argp)
  stop(msg)
}

library(StructuralVariantAnnotation)
library(VariantAnnotation)
library(tidyverse)
find_asm_linked_tis = function(vcf, bpgr, sample_name) {
  linkeddf = linked_assemblies(vcf)
  linked_pairs = linkeddf %>%
    arrange(vcfId) %>%
    group_by(linked_by) %>%
    summarise(
      vcfId1=dplyr::first(vcfId),
      vcfId2=dplyr::last(vcfId)) %>%
    group_by(vcfId1, vcfId2) %>%
    summarise(linked_by=paste0(linked_by, collapse=","))
  linked_bp = linked_pairs %>%
    ungroup() %>%
    filter(vcfId1 %in% names(bpgr) & vcfId2 %in% names(bpgr)) %>%
    mutate(
      left_vcfId = ifelse(strand(gr[vcfId1]) == "+", vcfId2, vcfId1),
      right_vcfId = ifelse(strand(gr[vcfId1]) == "+", vcfId1, vcfId2))
  geno_left = geno(vcf[linked_bp$left_vcfId])
  geno_right = geno(vcf[linked_bp$right_vcfId])
  gr_left = gr[linked_bp$left_vcfId]
  gr_right = gr[linked_bp$right_vcfId]
  linked_bp = linked_bp %>% mutate(
      tilen=start(gr_right) - start(gr_left),
      left_hom=end(gr_left) - start(gr_left),
      right_hom=end(gr_right) - start(gr_right),
      normal_qual = geno_left$QUAL[,1] + geno_right$QUAL[,1],
      somatic_qual = geno_left$QUAL[,2] + geno_right$QUAL[,2],
      status=ifelse(normal_qual > 350, "germline", ifelse(somatic_qual < 350, "low_qual", "somatic")),
      ti_desc=paste0(
        seqnames(gr[gr_left$partner]), ":", start(gr[gr_left$partner]),
        "<->",
        seqnames(gr_left), ":", start(gr_left), "-", start(gr_right),
        "<->",
        seqnames(gr[gr_right$partner]), ":", start(gr[gr_right$partner]))
    )
  return(linked_bp %>% mutate(sample_name=sample_name))
}
all_df = NULL
for (i in seq_along(argv$input)) {
  input_filename = argv$input[i]
  output_filename = argv$output_vcf[i]
  vcf = readVcf(input_filename)
  vcf = vcf[VariantAnnotation::fixed(vcf)$QUAL >= argv$qual]
  gr = breakpointRanges(vcf)
  gr = gr[gr$partner %in% names(gr)]
  linkeddf = find_asm_linked_tis(vcf, gr, input_filename)
  if (!is.na(output_filename)) {
    writeVcf(filename=output_filename, vcf[names(vcf) %in% c(
      linkeddf$left_vcfId,
      linkeddf$right_vcfId,
      gr[linkeddf$left_vcfId]$partner,
      gr[linkeddf$right_vcfId]$partner)])
  }
  all_df = bind_rows(all_df, linkeddf)
}
