#!/usr/bin/env Rscript
#
# Filters a raw GRIDSS VCF to high quality somatic calls
#
library(tidyverse)
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
library(data.table)
fixed(full_vcf)$ALT = CharacterList(lapply(fread(file=input_vcf, sep="\t", sep2=NULL, header=FALSE, stringsAsFactors=FALSE, select=5, skip=120)$V5, function(x) x))
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

vcf = full_vcf[names(bpgr)[as.logical(!bpfiltered)]]
bpgr = breakpointRanges(vcf)
vcf = vcf[names(bpgr)] # sanity reordering in case of asymetrical filtering
bpgr$af = gridss_af(bpgr, vcf, 2)
info(vcf)$BPI_AF = paste(bpgr$af, partner(bpgr)$af, sep=",")
VariantAnnotation::fixed(vcf)$FILTER = "PASS"
writeVcf(vcf, output_vcf)


# Assembly-based event linking
asm_linked_df <- data.frame(
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

# TODO:
# transitive calling reduction


# event completion

# Simple insertion events







