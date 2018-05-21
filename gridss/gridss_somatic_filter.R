#!/usr/bin/env Rscript
#
# Filters a raw GRIDSS VCF to high quality somatic calls
#
library(stringr)
usage = "Usage: Rscript gridss_somatic_filter.R <input VCF> <output VCF>"
args = commandArgs(TRUE)
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
input_vcf = args[1]
output_vcf = args[2]
source("libgridss.R")

# Filter to somatic calls
vcf = readVcf(input_vcf, "")
bpgr = breakpointRanges(vcf, unpartneredBreakends=FALSE)
begr = breakpointRanges(vcf, unpartneredBreakends=TRUE)
filtered = gridss_filter(bpgr, vcf)

# filter out 'shadow' calls of strong multi-mapping calls
# bwa overestimates the MAPQ of some multimapping reads
# and makes FP calls to alternate locations.
combinedgr <- bpgr
combinedgr$partner <- NULL
combinedgr <- c(combinedgr, begr)
# require exact overlaps
bestOverlap = findOverlaps(bpgr, combinedgr) %>%
  as.data.frame() %>%
  # bpgr offsets are the same as combinedgr ofsets
  filter(queryHits != subjectHits) %>%
  mutate(beQUAL = combinedgr$QUAL[subjectHits]) %>%
  group_by(queryHits) %>%
  summarise(beQUAL = max(beQUAL))
bpgr$overlapQUAL = 0
bpgr$overlapQUAL[bestOverlap$queryHits] = bestOverlap$beQUAL
i <- info(vcf[bpgr$vcfId])
better_call_filter = !is_short_event(bpgr) &
  bpgr$overlapQUAL > 3 * bpgr$QUAL &
  i$BANSRQ + i$BANRPQ > i$ASQ

filtered = filtered | better_call_filter

vcf = vcf[names(bpgr)[as.logical(!filtered)]]
bpgr = breakpointRanges(vcf)
vcf = vcf[names(bpgr)] # sanity reordering in case of asymetrical filtering
bpgr$af = gridss_af(bpgr, vcf, 2)
info(vcf)$BPI_AF = paste(bpgr$af, partner(bpgr)$af, sep=",")
VariantAnnotation::fixed(vcf)$FILTER = "PASS"
writeVcf(vcf, output_vcf)
