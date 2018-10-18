#!/usr/bin/env Rscript
#
# Filters a raw GRIDSS VCF to high quality somatic calls
#
library(stringr)
usage = "Usage: Rscript gridss_somatic_filter.R <pon directory> <input VCF> <output QUAL filter VCF> <output full somatic VCF>"
args = commandArgs(TRUE)
if (str_detect(args[1], "gridss_somatic_filter")) {
  args = args[-1]
}
if (length(args) != 4) {
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
if (pon_dir == "") {
  pon_dir = NULL
}
library(tidyverse)
library(readr)
source("libgridss.R")

# Filter to somatic calls
write(paste0("Reading ", input_vcf), stderr())
full_vcf = readVcf(input_vcf, "hg19")
# hard filter unpaired breakpoints (caused by inconsistent scoring across the two breakends)
full_vcf = full_vcf[is.na(info(full_vcf)$PARID) | info(full_vcf)$PARID %in% names(full_vcf)]
full_vcf = align_breakpoints(full_vcf)
write(paste0("Parsing SVs in ", input_vcf), stderr())
full_bpgr = breakpointRanges(full_vcf, unpartneredBreakends=FALSE)
full_begr = breakpointRanges(full_vcf, unpartneredBreakends=TRUE)
write(paste0("Calculating VAF ", input_vcf), stderr())
full_bpgr$af = gridss_somatic_bp_af(full_bpgr, full_vcf)
full_bpgr$af_str = paste(full_bpgr$af, partner(full_bpgr)$af, sep=",")
full_begr$af = gridss_somatic_be_af(full_begr, full_vcf)
full_begr$af_str = as.character(full_begr$af)
info(full_vcf)$BPI_AF = ""
info(full_vcf[names(full_bpgr)])$BPI_AF = full_bpgr$af_str
info(full_vcf[names(full_begr)])$BPI_AF = full_begr$af_str

write(paste0("Filtering pass 1 ", input_vcf), stderr())
bpfiltered = gridss_breakpoint_filter(full_bpgr, full_vcf, pon_dir=pon_dir)
befiltered = gridss_breakend_filter(full_begr, full_vcf, pon_dir=pon_dir)
# shadow breakpoint removed due to initial mapq20 filter reducing FP rate
# bpfiltered = .addFilter(bpfiltered, "shadow", is_shadow_breakpoint(bpgr, begr, full_vcf))

#bpfiltered = .addFilter(bpfiltered, "LOW_Qual", bpgr$QUAL < gridss.min_qual)
#som_llr = gridss_breakpoint_somatic_llr(full_vcf, bpgr=bpgr, contamination_rate=gridss.allowable_normal_contamination)

# - filter to only decent length assemblies?
#begr$calls_1k_window = countOverlaps(begr, rowRanges(full_vcf), ignore.strand=TRUE, maxgap=1000)

filters = rep("", length(full_vcf))
names(filters) = names(full_vcf)
filters[names(full_bpgr)] = bpfiltered
filters[names(full_begr)] = befiltered

vcf = full_vcf[passes_soft_filters(filters)]
vcf = vcf[is.na(info(vcf)$PARID) | info(vcf)$PARID %in% names(vcf)]
bpgr = full_bpgr[names(full_bpgr) %in% names(vcf)]
begr = full_begr[names(full_begr) %in% names(vcf)]

write(paste0("Calculating transitive links", input_vcf), stderr())
# transitive calling
transitive_df = transitive_calls(vcf, bpgr, report="max2") %>%
  # only make transitive calls were we actually know the path
  filter(!has_multiple_paths) %>%
  mutate(type="transitive")
# now we filter imprecise variants
is_imprecise = !(is.na(info(vcf)$IMPRECISE) | !info(vcf)$IMPRECISE) |
  !((!is.na(info(vcf)$PARID) & info(vcf)$ASSR + info(vcf)$SR + info(vcf)$IC > 0) |
      (is.na(info(vcf)$PARID) & info(vcf)$BASSR + info(vcf)$BSC > 0))
filters[names(vcf)[is_imprecise]] = paste0(filters[names(vcf)[is_imprecise]], ";imprecise")
filters[transitive_df$linked_by] = paste0(filters[transitive_df$linked_by], ";transitive")

vcf = full_vcf[passes_soft_filters(filters)]
vcf = vcf[is.na(info(vcf)$PARID) | info(vcf)$PARID %in% names(vcf)]
bpgr = full_bpgr[names(full_bpgr) %in% names(vcf)]
begr = full_begr[names(full_begr) %in% names(vcf)]

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
  filter(!str_detect(filters[vcfId], "PON"))
# Fix up pairing
event_link_df = event_link_df %>%
  filter(n() == 2) %>%
  ungroup()

# include both breakends of any linked breakpoints
# as linkage can be breakend specific (e.g. assembly, bpbeins)
link_rescue = bind_rows(link_df, event_link_df) %>% pull(vcfId) %>% unique()
link_rescue = c(link_rescue, bpgr[link_rescue[link_rescue %in% names(bpgr)]]$partner)

# Note that we don't rescue equivalent events
begr$partner = NA
eqv_link_df = linked_by_equivalent_variants(full_vcf, as(rbind(as.data.frame(bpgr), as.data.frame(begr)), "GRanges")) %>%
  filter(passes_final_filters(vcf[vcfId]) | vcfId %in% link_rescue) %>%
  group_by(linked_by) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  mutate(type="eqv")

link_summary_df = bind_rows(link_df, event_link_df, eqv_link_df) %>%
  group_by(vcfId) %>%
  summarise(linked_by=paste0(linked_by, collapse=","))

# Add linking information
info(full_vcf)$LOCAL_LINKED_BY = ""
info(full_vcf)$REMOTE_LINKED_BY = ""
info(full_vcf[link_summary_df$vcfId])$LOCAL_LINKED_BY = link_summary_df$linked_by
info(full_vcf[!is.na(info(full_vcf)$PARID)])$REMOTE_LINKED_BY = info(full_vcf[info(full_vcf[!is.na(info(full_vcf)$PARID)])$PARID])$LOCAL_LINKED_BY

# final qual filtering
fails_qual_without_rescue = !passes_final_filters(full_vcf, include.existing.filters=FALSE) & !(names(full_vcf) %in% link_rescue)
filters[names(full_vcf)[fails_qual_without_rescue]] = paste0(filters[names(full_vcf)[fails_qual_without_rescue]], ";qual")

################
# Write outputs
VariantAnnotation::fixed(full_vcf)$FILTER = ifelse(str_remove(filters, "^;") == "", "PASS", str_remove(filters, "^;"))

write(paste0("Writing ", output_vcf), stderr())
vcf = full_vcf[passes_soft_filters(filters)]
vcf = vcf[is.na(info(vcf)$PARID) | info(vcf)$PARID %in% names(vcf)]
writeVcf(vcf, output_vcf, index=TRUE)

write(paste0("Writing ", output_full_vcf), stderr())
vcf = full_vcf[passes_very_hard_filters(filters)]
vcf = vcf[is.na(info(vcf)$PARID) | info(vcf)$PARID %in% names(vcf)]
writeVcf(vcf, output_full_vcf, index=TRUE)






