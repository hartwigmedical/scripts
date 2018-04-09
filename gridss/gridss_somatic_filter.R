#!/usr/bin/env Rscript
#
# Filters a raw GRIDSS VCF to high quality somatic calls
#
usage <- "Usage: gridss_somatic_filter <input VCF> <output VCF>"
args <- commandArgs(TRUE)
if (length(args) != 3) {
  write(usage, stderr())
  q(save="no", status=1)
}
if (!file.exists(args[2])) {
  write(paste(args[2], "not found"), stderr())
  q(save="no", status=1)
}
input_vcf <- args[2]
output_vcf <- args[3]
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
