#!/usr/bin/env Rscript
#
# Filters a raw GRIDSS VCF to high quality somatic calls
#
library(stringr)
usage <- "Usage: Rscript gridss_somatic_filter.R <input VCF> <output VCF>"
args <- commandArgs(TRUE)
if (length(args) == 3 & str_detect(args[1], "gridss_somatic_filter")) {
  args = args[-1]
}
if (length(args) != 2) {
  write(usage, stderr())
  q(save="no", status=1)
}
if (!file.exists(args[1])) {
  write(paste(args[1], "not found"), stderr())
  q(save="no", status=1)
}
input_vcf <- args[1]
output_vcf <- args[2]
source("libgridss.R")

# Filter to somatic calls
vcf <- readVcf(input_vcf, "")
gr <- breakpointRanges(vcf)
filtered <- gridss_filter(gr, vcf)

vcf <- vcf[names(gr)[!filtered]]
gr <- breakpointRanges(vcf)
vcf <- vcf[names(gr)] # sanity reordering in case of asymetrical filtering
gr$af <- gridss_af(gr, vcf, 2)
info(vcf)$BPI_AF <- paste(gr$af, partner(gr)$af, sep=",")
VariantAnnotation::fixed(vcf)$FILTER <- "PASS"
writeVcf(vcf, output_vcf)
