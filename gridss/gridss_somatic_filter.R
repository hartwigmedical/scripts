#!/usr/bin/env Rscript
#
# Filters a raw GRIDSS VCF to high quality somatic calls
#
library(tidyverse)
library(readr)
library(stringr)
usage = "Usage: Rscript gridss_somatic_filter.R <pon directory> <input VCF> <output QUAL filter VCF> <output full somatic VCF> <filter reason CSV>"
args = commandArgs(TRUE)
if (str_detect(args[1], "gridss_somatic_filter")) {
  args = args[-1]
}
if (length(args) != 4 & length(args) != 5) {
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
pon_dir = args[1]
input_vcf = args[2]
output_vcf = args[3]
output_full_vcf = args[4]
filter_reason_csv = args[5]
source("libgridss.R")

# Filter to somatic calls
write(paste0("Reading ", input_vcf), stderr())
full_vcf = readVcf(input_vcf, "hg19")
# hard filter variants that are over the final calling threshold in the normal - they're definitely not somatic
full_vcf = full_vcf[geno(full_vcf)$QUAL[,1] < gridss.min_qual]
write(paste0("Parsing ", input_vcf), stderr())
bpgr = breakpointRanges(full_vcf, unpartneredBreakends=FALSE)
begr = breakpointRanges(full_vcf, unpartneredBreakends=TRUE)
write(paste0("Filtering pass 1 ", input_vcf), stderr())
bpfiltered = gridss_breakpoint_filter(bpgr, full_vcf)
befiltered = gridss_breakend_filter(begr, full_vcf)
# shadow breakpoint removed due to initial mapq20 filter reducing FP rate
# bpfiltered = .addFilter(bpfiltered, "shadow", is_shadow_breakpoint(bpgr, begr, full_vcf))

#bpfiltered = .addFilter(bpfiltered, "LOW_Qual", bpgr$QUAL < gridss.min_qual)
#som_llr = gridss_breakpoint_somatic_llr(full_vcf, bpgr=bpgr, contamination_rate=gridss.allowable_normal_contamination)

# - filter to only decent length assemblies?
#begr$calls_1k_window = countOverlaps(begr, rowRanges(full_vcf), ignore.strand=TRUE, maxgap=1000)
if (!is.null(filter_reason_csv)) {
  write(paste0("Writing ", filter_reason_csv), stderr())
  readr::write_csv(
    x=data.frame(vcfId=c(names(bpgr), names(begr)), filters=c(bpfiltered, befiltered)),
    path=filter_reason_csv)
}
bpfiltered = bpfiltered != ""
befiltered = befiltered != ""
bp_vcf = full_vcf[names(bpgr)[!bpfiltered]]
bpgr = breakpointRanges(bp_vcf) # fix any asymetrical filtering
begr = begr[!befiltered]
vcf = full_vcf[names(full_vcf) %in% c(names(bpgr), names(begr))]

write(paste0("Calculating transitive links", input_vcf), stderr())
# transitive calling
transitive_df = transitive_calls(vcf, bpgr, report="max2") %>%
  # only make transitive calls were we actually know the path
  filter(!has_multiple_paths) %>%
  mutate(type="transitive")
# now we filter imprecise variants
vcf = vcf[is.na(info(vcf)$IMPRECISE) | !info(vcf)$IMPRECISE]
vcf = vcf[info(vcf)$ASSR + info(vcf)$SR + info(vcf)$IC > 0] # not unanchored
bpgr = breakpointRanges(vcf, unpartneredBreakends=FALSE)
begr = breakpointRanges(vcf, unpartneredBreakends=TRUE)
vcf = vcf[names(vcf) %in% c(names(bpgr), names(begr))]
bpgr = breakpointRanges(vcf, unpartneredBreakends=FALSE)

write(paste0("Calculating VAF ", input_vcf), stderr())
bpgr$af = gridss_somatic_bp_af(bpgr, vcf)
bpgr$af_str = paste(bpgr$af, partner(bpgr)$af, sep=",")
begr$af = gridss_somatic_be_af(begr, vcf)
begr$af_str = as.character(begr$af)
info(vcf)$BPI_AF = rep("", length(vcf))
info(vcf[names(bpgr)])$BPI_AF = bpgr$af_str
info(vcf[names(begr)])$BPI_AF = begr$af_str
VariantAnnotation::fixed(vcf)$FILTER = ifelse(names(vcf) %in% c(
  names(bpgr)[gridss_overlaps_breakpoint_pon(bpgr, pon_dir=pon_dir)],
  names(begr)[gridss_overlaps_breakend_pon(begr, pon_dir=pon_dir)]), "PON", "PASS")

write(paste0("Calculating assembly links ", input_vcf), stderr())
# Assembly-based event linking
asm_linked_df = linked_assemblies(vcf) %>%
  mutate(type="asm")

link_df = bind_rows(asm_linked_df, transitive_df) %>%
  mutate(linking_group=str_replace(linked_by, "/.*$", "")) %>%
  mutate(pass=passes_final_filters(vcf[vcfId])) %>%
  group_by(linking_group) %>%
  mutate(pass=any(pass)) %>%
  ungroup() %>%
  filter(pass)

write(paste0("Calculating bebe insertion links ", input_vcf), stderr())
# Insertion linkage
bebeins_link_df = linked_by_breakend_breakend_insertion_classification(begr) %>%
  group_by(linked_by) %>%
  mutate(pass=passes_final_filters(vcf[vcfId])) %>%
  mutate(pass=any(pass)) %>%
  ungroup() %>%
  filter(pass) %>%
  mutate(type="bebeins")
write(paste0("Calculating bebp insertion links ", input_vcf), stderr())
bebpins_link_df = linked_by_breakpoint_breakend_insertion_classification(bpgr, begr) %>%
  group_by(linked_by) %>%
  mutate(pass=passes_final_filters(vcf[vcfId])) %>%
  mutate(pass=any(pass)) %>%
  ungroup() %>%
  filter(pass) %>%
  mutate(type="bebpins")
# Inversion linkage
write(paste0("Calculating simple inversions ", input_vcf), stderr())
inv_link_df = linked_by_simple_inversion_classification(bpgr) %>%
  group_by(linked_by) %>%
  mutate(pass=passes_final_filters(vcf[vcfId])) %>%
  mutate(pass=any(pass)) %>%
  ungroup() %>%
  filter(pass) %>%
  mutate(type="inv")
# Deletion bridge linkage
# TODO: do we want to do this?
# I'm suspicious of the model used in ChainFinder PMC3690918
# Notably: I'm suspicous that "repair with major DNA loss" is actually a thing
# given the catastrophic nature of chromo*, a more reasonable explaination is
# an additional DSB with the subsequent loss of that DNA fragment.
# Given the focal nature of chromoplexy, ChainFinder works because it just
# finds the focal events, not because the model is correct.
# TODO: show this by modelling additional focal DSBs
write(paste0("Calculating dsb links ", input_vcf), stderr())
dsb_link_df = linked_by_dsb(bpgr) %>%
  group_by(linked_by) %>%
  mutate(pass=passes_final_filters(vcf[vcfId])) %>%
  mutate(pass=any(pass)) %>%
  ungroup() %>%
  filter(pass) %>%
  mutate(type="dsb")


write(paste0("Removing duplicated/conflicting links ", input_vcf), stderr())
# linking priorities:
# - asm independent of other linkages
# - transitive independent of other linkages
# - ins, inv, dsb linkages
event_link_df = bind_rows(
  bebeins_link_df,
  bebpins_link_df,
  inv_link_df,
  dsb_link_df) %>%
  dplyr::select(vcfId, linked_by) %>%
  mutate(QUAL=rowRanges(vcf)[vcfId]$QUAL)
# Only keep the best QUAL event linkage
event_link_df = event_link_df %>%
  group_by(linked_by) %>%
  mutate(linkQUAL = pmin(QUAL)) %>%
  group_by(vcfId) %>%
  filter(QUAL == linkQUAL) %>%
  group_by(linked_by)
# Don't event link to PON filtered variants
event_link_df = event_link_df %>%
  filter(VariantAnnotation::fixed(vcf[vcfId])$FILTER != "PON")
# Fix up pairing
event_link_df = event_link_df %>%
  filter(n() == 2) %>%
  ungroup()

link_summary_df = bind_rows(link_df, event_link_df) %>%
  group_by(vcfId) %>%
  summarise(linked_by=paste0(linked_by, collapse=","))

# include both breakends of any linked breakpoints
# as linkage can be breakend specific (e.g. assembly, bpbeins)
linked_vcfIds = c(link_summary_df$vcfId,
  link_summary_df %>%
    filter(vcfId %in% names(bpgr)) %>%
    mutate(partner_vcfId=bpgr[vcfId]$partner) %>%
    pull(partner_vcfId))

#####
# Final QUAL filtering
#
write(paste0("Writing ", output_full_vcf), stderr())
info(vcf)$LOCAL_LINKED_BY = ""
info(vcf)$REMOTE_LINKED_BY = ""
info(vcf[link_summary_df$vcfId])$LOCAL_LINKED_BY = link_summary_df$linked_by
info(vcf[!is.na(info(vcf)$PARID)])$REMOTE_LINKED_BY = info(vcf[info(vcf[!is.na(info(vcf)$PARID)])$PARID])$LOCAL_LINKED_BY

vcf = vcf[!(names(vcf) %in% c(transitive_df$transitive_start, transitive_df$transitive_end))] # remove transitive calls

writeVcf(align_breakpoints(vcf), output_full_vcf)


write(paste0("Writing ", output_vcf), stderr())
vcf = vcf[passes_final_filters(vcf) | names(vcf) %in% linked_vcfIds]
bpgr = breakpointRanges(vcf, unpartneredBreakends=FALSE)
begr = breakpointRanges(vcf, unpartneredBreakends=TRUE)

vcf = vcf[names(vcf) %in% c(names(bpgr), names(begr))]
# ggplot(as.data.frame(begr)) + aes(x=QUAL, fill=info(vcf[names(begr)])$LINKED_BY == "") + geom_histogram(bins=100) + scale_x_log10()
# TODO: filter breakends within 100bp of microsatellites; we don't want purple to segment at every one

vcf = align_breakpoints(vcf)
writeVcf(vcf, output_vcf)





