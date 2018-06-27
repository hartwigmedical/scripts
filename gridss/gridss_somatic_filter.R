#!/usr/bin/env Rscript
#
# Filters a raw GRIDSS VCF to high quality somatic calls
#
library(tidyverse)
library(readr)
library(stringr)
usage = "Usage: Rscript gridss_somatic_filter.R <pon directory> <input VCF> <output VCF>"
args = commandArgs(TRUE)
if (length(args) == 3 & str_detect(args[1], "gridss_somatic_filter")) {
  args = args[-1]
}
if (length(args) != 3) {
  write(usage, stderr())
  q(save="no", status=1)
}
if (!file.exists(args[1])) {
  write(paste(args[1], "not found"), stderr())
  q(save="no", status=1)
}
if (!file.exists(args[2])) {
  write(paste(args[2], "not found"), stderr())
  q(save="no", status=1)
}
input_vcf = args[1]
pon_dir = args[2]
output_vcf = args[3]
source("libgridss.R")

# Filter to somatic calls
full_vcf = readVcf(input_vcf, "hg19")
# hard filter variants that are over the final calling threshold in the normal - they're definitely not somatic
full_vcf = full_vcf[geno(full_vcf)$QUAL[,1] < gridss.min_qual]
bpgr = breakpointRanges(full_vcf, unpartneredBreakends=FALSE)
begr = breakpointRanges(full_vcf, unpartneredBreakends=TRUE)
bpfiltered = gridss_breakpoint_filter(bpgr, full_vcf, pon_dir=pon_dir)
befiltered = gridss_breakend_filter(begr, full_vcf, pon_dir=pon_dir)
bpfiltered = .addFilter(bpfiltered, "shadow", is_shadow_breakpoint(bpgr, begr, full_vcf))

#bpfiltered = .addFilter(bpfiltered, "LOW_Qual", bpgr$QUAL < gridss.min_qual)
#som_llr = gridss_breakpoint_somatic_llr(full_vcf, bpgr=bpgr, contamination_rate=gridss.allowable_normal_contamination)

# - filter to only decent length assemblies?
#begr$calls_1k_window = countOverlaps(begr, rowRanges(full_vcf), ignore.strand=TRUE, maxgap=1000)

bpfiltered = bpfiltered != ""
befiltered = befiltered != ""

bp_vcf = full_vcf[names(bpgr)[!bpfiltered]]
bpgr = breakpointRanges(bp_vcf) # fix any asymetrical filtering
begr = begr[!befiltered]
vcf = full_vcf[names(full_vcf) %in% c(names(bpgr), names(begr))]
bpgr$af = gridss_somatic_bp_af(bpgr, vcf)
bpgr$af_str = paste(bpgr$af, partner(bpgr)$af, sep=",")
begr$af = gridss_somatic_be_af(begr, vcf)
begr$af_str = as.character(begr$af)
info(vcf)$BPI_AF = rep("", length(vcf))
info(vcf[names(bpgr)])$BPI_AF = bpgr$af_str
info(vcf[names(begr)])$BPI_AF = begr$af_str
VariantAnnotation::fixed(vcf)$FILTER = "PASS"

# Assembly-based event linking
asm_linked_df = linked_assemblies(vcf)
# transitive calling reduction
transitive_df = transitive_calls(vcf, bpgr)
link_df = asm_linked_df %>%
  inner_join(data.frame(vcfid=names(vcf), event=info(vcf)$EVENT), by="event") %>%
  dplyr::select(vcfid, linked_by) %>%
  bind_rows(transitive_df) %>%
  group_by(vcfid) %>%
  summarise(linked_by=paste0(linked_by, collapse=","))

#####
# EVENT COMPLETION
#

# TODO: complete chains
# TODO: simplify simple insertions

#####
# Final filtering
#
info(vcf)$LINKED_BY = ""
info(vcf[link_df$vcfid])$LINKED_BY = link_df$linked_by
vcf = vcf[!(names(vcf) %in% transitive_df$linked_by)] # remove transitive calls
vcf = vcf[rowRanges(vcf)$QUAL >= gridss.min_qual | info(vcf)$LINKED_BY != ""]
vcf = vcf[rowRanges(vcf)$QUAL >= gridss.min_qual * gridss.single_breakend_multiplier | info(vcf)$LINKED_BY != ""]

bpgr = breakpointRanges(vcf, unpartneredBreakends=FALSE)
begr = breakpointRanges(vcf, unpartneredBreakends=TRUE)

vcf = vcf[names(vcf) %in% c(names(bpgr), names(begr))]
# ggplot(as.data.frame(begr)) + aes(x=QUAL, fill=info(vcf[names(begr)])$LINKED_BY == "") + geom_histogram(bins=100) + scale_x_log10()
# TODO: filter breakends within 100bp of microsatellites; we don't want purple to segment at every one

vcf = align_breakpoints(vcf)
writeVcf(vcf, output_vcf)




