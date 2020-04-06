library(StructuralVariantAnnotation)
library(cowplot)
library(tidyverse)
library(rtracklayer)
theme_set(theme_bw())

#' Replaces the NA values in a with corresponding values in b
#' @param a,b objects to be tested or coerced.
#' @return The altered object.
'%na%' <- function(a, b) {
  if (is.null(a) || length(a) == 0) return(b)
  if (is.null(b) || length(b) == 0) return(a)
  return(ifelse(is.na(a), b, a))
}

datadir="D:/hartwig/svtoolkit/germline/"


giab_tier1_regions = import(paste0(datadir, "HG002_SVs_Tier1_v0.6.bed"))
giab_tier1_vcf = readVcf(paste0(datadir, "HG002_SVs_Tier1_v0.6.vcf.gz"))
giab_tier1_gr = breakpointRanges(giab_tier1_vcf)
giab_tier1_gr$caller = "GIAB Tier 1 Truth"
giab_tier1_pass_gr = giab_tier1_gr[giab_tier1_gr$FILTER == "PASS"]


cgrs = list()
for (caller_file in list.files(datadir, pattern="*.vcf$")) {
  caller_name = str_replace(caller_file, ".vcf", "")
  write(paste("Processing", caller_name), stderr())
  caller_vcf = readVcf(paste0(datadir, caller_file))
  VariantAnnotation::fixed(caller_vcf)$QUAL = score_for_caller(caller_name, caller_vcf)
  caller_bpgr = breakpointRanges(caller_vcf)
  caller_begr = breakendRanges(caller_vcf)
  caller_gr = c(caller_bpgr, caller_begr)
  caller_gr$caller = caller_name
  cgrs[[caller_name]] = caller_gr
}

score_for_caller = function(caller_name, vcf) {
  qual = switch(
    caller_name,
    "crest"=(info(vcf)$right_softclipped_read_count %na% 0) + (info(vcf)$left_softclipped_read_count %na% 0),
    "pindel_0.2.5b6"=geno(vcf)$AD,
    "delly_0.7.6"=(info(vcf)$PE %na% 0) + (info(vcf)$SR %na% 0),
    #"breakdancer_1.4.5"
    #"gridss_1.6.1",
    #"gridss_2.8.1",
    #"hydra_master20160129"=,
    #"manta_1.1.1",
    rowRanges(vcf)$QUAL)
  qual[is.na(qual)] = 0
  return(qual)
}

calc_roc_pass_all = function(truth_gr, caller_gr, ...) {
  all_calls=calc_roc(truth_gr, caller_gr, ...)
  all_calls$roc$subset = "All calls"
  all_calls$gr$subset = "All calls"
  pass_calls=calc_roc(truth_gr, caller_gr[is.na(caller_gr$FILTER) | caller_gr$FILTER %in% c("PASS", "", ".")], ...)
  pass_calls$roc$subset = "PASS only"
  pass_calls$gr$subset = "PASS only"
  return(list(roc=bind_rows(all_calls$roc, pass_calls$roc), gr=c(all_calls$gr, pass_calls$gr)))
}
calc_roc = function(truth_gr, caller_gr, filter_to_region_gr=NULL, bpmaxgap=100, bemaxgap=5, minsize=50, additional_filter=function(gr) { gr } ) {
  caller_name=unique(caller_gr$caller)
  write(paste("Processing", caller_name), stderr())
  if (!is.null(filter_to_region_gr)) {
    truth_gr = truth_gr[overlapsAny(truth_gr, filter_to_region_gr, ignore.strand=TRUE)]
    caller_gr = caller_gr[overlapsAny(caller_gr, filter_to_region_gr, ignore.strand=TRUE)]
  }
  truth_gr = additional_filter(truth_gr)
  caller_gr = additional_filter(caller_gr)
  caller_gr = caller_gr[is.na(caller_gr$partner) | caller_gr$partner %in% names(caller_gr)]
  truth_gr = truth_gr[is.na(truth_gr$partner) | truth_gr$partner %in% names(truth_gr)]
  caller_bpgr = caller_gr[!is.na(caller_gr$partner)]
  caller_begr = caller_gr[is.na(caller_gr$partner)]
  
  truth_gr = truth_gr[abs(truth_gr$svLen) >= minsize]
  caller_bpgr = caller_bpgr[!is.na(caller_bpgr$svLen) & abs(caller_bpgr$svLen) >= minsize]
  
  bphits = findBreakpointOverlaps(caller_bpgr, truth_gr, maxgap=bpmaxgap) %>%
    as.data.frame() %>%
    mutate(QUAL=caller_bpgr$QUAL[queryHits]) %>%
    group_by(subjectHits) %>%
    arrange(desc(QUAL)) %>%
    mutate(isBestHit=row_number() == 1) %>%
    ungroup()
  bestbphits = bphits %>% filter(isBestHit)
  
  behits = findOverlaps(caller_begr, truth_gr, maxgap=bemaxgap) %>%
    as.data.frame() %>%
    mutate(QUAL=caller_begr$QUAL[queryHits]) %>%
    group_by(subjectHits) %>%
    arrange(desc(QUAL)) %>%
    mutate(isBestHit=row_number() == 1) %>%
    ungroup()
  bestbehits = behits %>% filter(isBestHit)

  caller_bpgr$tp = FALSE
  caller_bpgr$tpdup = FALSE
  caller_begr$tp = rep(TRUE, length(caller_begr))
  caller_begr$tpdup = rep(TRUE, length(caller_begr))
  truth_gr$tp = FALSE
  truth_gr$QUAL = -1
  truth_gr$betp = FALSE
  truth_gr$beQUAL = -1

  # Match to breakpoint hits
  truth_gr$tp[bphits$subjectHits] = TRUE
  truth_gr$QUAL[bestbphits$subjectHits] = bestbphits$QUAL
  caller_bpgr$tpdup[bphits$queryHits] = TRUE
  caller_bpgr$tp[bestbphits$queryHits] = TRUE
  caller_bpgr$tpdup = caller_bpgr$tpdup & !caller_bpgr$tp
  
  caller_begr$tpdup[behits$queryHits] = TRUE
  caller_begr$tp[bestbehits$queryHits] = TRUE
  caller_begr$tpdup = caller_begr$tpdup & !caller_begr$tp
  truth_gr$tp[behits$subjectHits] = TRUE
  truth_gr$QUAL[bestbehits$subjectHits] = pmax(truth_gr$QUAL[bestbehits$subjectHits], bestbehits$QUAL)
  
  truth_gr$fn = !truth_gr$tp
  roc_gr = c(truth_gr, caller_bpgr[!caller_bpgr$tp], caller_begr[!caller_begr$tp])
  roc_gr$caller = caller_name
  
  roc = roc_gr %>% as.data.frame() %>%
    replace_na(list(fn=FALSE, tpdup=FALSE)) %>%
    mutate(QUAL=round(QUAL)) %>%
    dplyr::select(QUAL, tp, tpdup, fn, caller) %>%
    group_by(caller, QUAL) %>%
    summarise(tp=sum(tp), tpdup=sum(tpdup), fn=sum(fn), ncalls=n()) %>%
    group_by(caller) %>%
    arrange(desc(QUAL)) %>%
    mutate(cumtp=cumsum(tp), cumtpdup=cumsum(tpdup), fn=sum(fn),cumncalls=cumsum(ncalls)) %>%
    mutate(
      recall=cumtp/(max(cumtp) + fn),
      precision=cumtp/cumncalls) %>%
    ungroup() %>%
    filter(QUAL >= 0)
  return(list(roc=roc, gr=roc_gr))
}

cgrs[["HG002_SVs_Tier1_v0.6"]] = NULL
# Add GRIDSS breakpoint only
bponly=cgrs[["gridss_2.8.1"]]
bponly$caller=paste0(bponly$caller, " breakpoints only")
bponly = bponly[!is.na(bponly$partner)]
cgrs[[unique(bponly$caller)]] = bponly
roclist = lapply(cgrs, function(gr) {calc_roc_pass_all(giab_tier1_pass_gr, gr, giab_tier1_regions)})
combined_rocgr = unlist(GRangesList(lapply(roclist, function(x) x$gr)))
rocdf = bind_rows(lapply(roclist, function(x) x$roc))

ggplot(rocdf) +
  aes(x=recall, y=precision, colour=caller, shape=subset) +
  geom_point() +
  scale_color_brewer(palette="Dark2") +
  labs(title="HG002 60x\nGIAB Tier 1 truth set")

ggplot(as.data.frame(combined_rocgr)) +
  aes(x=abs(svLen), fill=tp) + 
  geom_histogram() +
  facet_grid(caller ~. , scales="free") +
  scale_x_log10()

ggplot(as.data.frame(combined_rocgr) %>% filter(tp | (fn %na% FALSE), str_detect(caller, "gridss") | str_detect(caller, "manta"))) +
  aes(x=abs(svLen), fill=ifelse(tp, "TP", ifelse(fn, "FN", "FP"))) + 
  geom_histogram(bins=100) +
  scale_x_log10() +
  facet_grid(caller ~. , scales="free")




xgr = combined_rocgr[!is.na(combined_rocgr)]










