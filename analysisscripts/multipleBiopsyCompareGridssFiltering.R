library(purple)
detach("package:purple", unload=TRUE)
library(purple)
library(StructuralVariantAnnotation)
library(tidyverse)
library(ggplot2)
source("libgridss.R")

vcf1 = readVcf("../../down/manta_compare/CPCT02010567R_CPCT02010567T.gridss.vcf", "hg19")
vcf2 = readVcf("../../down/manta_compare/CPCT02010567R_CPCT02010567TII.gridss.vcf", "hg19")


## Breakpoint comparison

bpgr1 = breakpointRanges(vcf1, unpartneredBreakends=FALSE)
bpgr2 = breakpointRanges(vcf2, unpartneredBreakends=FALSE)
bpgr1$filters_applied = gridss_breakpoint_filter(bpgr1, vcf1)
bpgr1$filters_applied[bpgr1$filters_applied == ""] = "PASS"
bpgr2$filters_applied = gridss_breakpoint_filter(bpgr2, vcf2)
bpgr2$filters_applied[bpgr2$filters_applied == ""] = "PASS"
hits = findBreakpointOverlaps(bpgr1, bpgr2, maxgap=20)
bpgr1$bestHit = NA_character_
bpgr2$bestHit = NA_character_
hits = hits %>% arrange(bpgr2$QUAL[hits$subjectHits])
bpgr1$bestHit[hits$queryHits] = names(bpgr2)[hits$subjectHits]
hits = hits %>% arrange(bpgr1$QUAL[hits$queryHits])
bpgr2$bestHit[hits$subjectHits] = names(bpgr1)[hits$queryHits]

bpgr1$matched_filter = "No Match"
bpgr1$matched_filter[!is.na(bpgr1$bestHit)] = bpgr2[bpgr1$bestHit[!is.na(bpgr1$bestHit)]]$filters_applied


table(bpgr1$filters_applied, bpgr1$matched_filter) %>%
  as.data.frame() %>%
  filter(Freq != 0) %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  filter(Var1 == "PASS" | Var2 == "PASS") %>%
  View()

pass_normal = bpgr1[bpgr1$filters_applied == ";normalSupport" & bpgr1$matched_filter == "PASS"]
normal_pass = bpgr1[bpgr1$filters_applied == "PASS" & bpgr1$matched_filter == ";normalSupport"]
gr = c(pass_normal, normal_pass)
data.frame(
  n1VF = .genosum(geno(vcf1[gr$vcfId])$VF,1),
  n1REF = .genosum(geno(vcf1[gr$vcfId])$REF,1),
  t1VF = .genosum(geno(vcf1[gr$vcfId])$VF,2),
  t1REF = .genosum(geno(vcf1[gr$vcfId])$REF,2),
  t1svlen = gr$svLen,

  n2VF = .genosum(geno(vcf2[gr$bestHit])$VF,1),
  n2REF = .genosum(geno(vcf2[gr$bestHit])$REF,1),
  t2VF = .genosum(geno(vcf2[gr$bestHit])$VF,2),
  t2REF = .genosum(geno(vcf2[gr$bestHit])$REF,2),
  t2svlen = bpgr2[gr$bestHit]$svLen
) %>% View()

ggplot(bpgr1[bpgr1$filters_applied == "PASS" & bpgr1$matched_filter == "No Match"] %>% as.data.frame(), aes(x=QUAL)) + geom_histogram(bins=100)
