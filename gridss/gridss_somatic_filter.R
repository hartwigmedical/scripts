#!/usr/bin/env Rscript
#
# Filters a raw GRIDSS VCF to high quality somatic calls
#
library(tidyverse)
library(readr)
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
full_vcf = readVcf(input_vcf, "hg19")
# work-around for https://github.com/Bioconductor/VariantAnnotation/issues/8
read_alt_directly_from_vcf = function() {
  alt = read_tsv(input_vcf, comment="#", col_names=c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", seq_len(ncol(geno(full_vcf)[[1]]))), cols_only(ALT=col_character()))$ALT
  return(CharacterList(lapply(as.character(alt), function(x) x)))
}
VariantAnnotation::fixed(full_vcf)$ALT = read_alt_directly_from_vcf()
bpgr = breakpointRanges(full_vcf, unpartneredBreakends=FALSE)
begr = breakpointRanges(full_vcf, unpartneredBreakends=TRUE)
bpfiltered = gridss_breakpoint_filter(bpgr, full_vcf)
befiltered = gridss_breakend_filter(begr, full_vcf)

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
i <- info(full_vcf[bpgr$vcfId])
better_call_filter = !is_short_event(bpgr) &
  bpgr$overlapQUAL > 3 * bpgr$QUAL &
  i$BANSRQ + i$BANRPQ > i$ASQ

bpfiltered = bpfiltered | better_call_filter

# - filter to only decent length assemblies?
begr$calls_1k_window = countOverlaps(begr, rowRanges(full_vcf), ignore.strand=TRUE, maxgap=1000)

bp_vcf = full_vcf[names(bpgr)[as.logical(!bpfiltered)]]
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
asm_linked_df = data.frame(
  event=info(vcf)$EVENT,
  beid=sapply(info(vcf)$BEID, function(x) paste0(c("", unlist(x)), collapse=";"))) %>%
  separate_rows(beid, sep=";") %>%
  filter(!is.na(beid) & nchar(beid) > 0) %>%
  group_by(event, beid) %>%
  distinct() %>%
  group_by(beid) %>%
  filter(n() > 1) %>%
  mutate(linked_by=beid) %>%
  dplyr::select(event, linked_by=beid)

# transitive calling reduction
transitive_df = transitive_breakpoints(bpgr, min_segment_length=20, report="all")
transitive_df = transitive_df %>%
  filter(info(vcf[transitive_df$transitive])$IMPRECISE) %>%
  mutate(full_path=bp_path)
if (nrow(transitive_df) != 0) {
  transitive_df = transitive_df %>%
    separate_rows(bp_path, sep=";") %>%
    mutate(imprecise=info(vcf[bp_path])$IMPRECISE) %>%
    group_by(transitive, full_path) %>%
    summarise(
      imprecise=sum(imprecise),
      hops=n(),
      min_length=min(min_length),
      max_length=max(max_length)) %>%
    # for now we just link imprecise transitive events
    # via precise calls
    filter(imprecise == 0) %>%
    group_by(transitive) %>%
    top_n(1, min_length) %>%
    ungroup() %>%
    dplyr::select(linked_by = transitive, vcfid = full_path) %>%
    separate_rows(vcfid, sep=";")
} else {
  # work-around for https://github.com/tidyverse/tidyr/issues/470
  transitive_df = data.frame(linked_by="placeholder", vcfid="placeholder", stringsAsFactors=FALSE) %>% filter(FALSE)
}

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
vcf = vcf[rowRanges(vcf)$QUAL >= gridss.min_breakpoint_qual | info(vcf)$LINKED_BY != ""]
vcf = vcf[rowRanges(vcf)$QUAL >= gridss.min_breakend_qual | info(vcf)$LINKED_BY != ""]

bpgr = breakpointRanges(vcf, unpartneredBreakends=FALSE)
begr = breakpointRanges(vcf, unpartneredBreakends=TRUE)

vcf = vcf[names(vcf) %in% c(names(bpgr), names(begr))]
# ggplot(as.data.frame(begr)) + aes(x=QUAL, fill=info(vcf[names(begr)])$LINKED_BY == "") + geom_histogram(bins=100) + scale_x_log10()
# TODO: filter breakends within 100bp of microsatellites; we don't want purple to segment at every one

vcf = align_breakpoints(vcf)
writeVcf(vcf, output_vcf)




