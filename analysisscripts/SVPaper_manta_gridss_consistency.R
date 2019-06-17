source("libSVPaper.R")
library(StructuralVariantAnnotation)

data_dir=paste0(dropbox, "/gridss_manta_consistency")

bpgr = NULL
begr = NULL

ann_caller_sample = function(filename, caller, sample, befct) {
  gr = befct(readVcf(filename))
  if (length(gr) > 0) {
    gr$caller = caller
    gr$sample = sample
  }
  return(gr)
}
for (manta_file in list.files(path=data_dir, pattern="*.manta.vcf.gz", full.names=TRUE)) {
  sample=str_extract(manta_file, "[^/]+(.manta.vcf.gz)$")
  sample=str_extract(sample, "^[^.]+")
  write(paste("Processing", sample), stderr())
  gridss_file = str_replace(manta_file, ".manta.vcf.gz", "T.gridss.somatic.vcf.gz")
  purple_file = str_replace(manta_file, ".manta.vcf.gz", "T.purple.sv.ann.vcf.gz")
  bpgr = c(bpgr,
    ann_caller_sample(manta_file, "manta", sample, breakpointRanges),
    ann_caller_sample(gridss_file, "gridss", sample, breakpointRanges),
    ann_caller_sample(purple_file, "purple", sample, breakpointRanges))
  begr = c(begr,
    ann_caller_sample(gridss_file, "gridss", sample, breakendRanges),
    ann_caller_sample(purple_file, "purple", sample, breakendRanges))
}
hits_by_caller = function(bpgr) {
  bpgr$hits = 0
  for (caller in unique(bpgr$caller)) {
    callergr = bpgr[bpgr$caller == caller]
    for (sample in unique(callergr$sample)) {
      hits = bpgr$hits[bpgr$caller == caller & bpgr$sample == sample]
      for (otherSample in unique(callergr$sample)) {
        if (sample != otherSample) {
          hits = hits + (countBreakpointOverlaps(
              bpgr[bpgr$caller == caller & bpgr$sample == sample],
              bpgr[bpgr$caller == caller & bpgr$sample == otherSample]) > 0)
        }
      }
      bpgr$hits[bpgr$caller == caller & bpgr$sample == sample] = hits
    }
  }
  return(bpgr)
}
passingbpgr = bpgr[bpgr$FILTER == "PASS"]
bpna12878 = hits_by_caller(passingbpgr[passingbpgr$sample %in% c("GIABvsSELFv001", "GIABvsSELFv002", "GIABvsSELFv003")])
bpcolo829 = hits_by_caller(passingbpgr[passingbpgr$sample %in% c("COLO829v001", "COLO829v002", "COLO829v003")])

summary_na12878 = bpna12878 %>% data.frame() %>%
    group_by(caller, hits) %>%
    summarise(n=n()) %>%
    mutate(
      called_in_samples=hits + 1,
      events=n / called_in_samples)
summary_colo829 =  bpcolo829 %>% data.frame() %>%
  group_by(caller, hits) %>%
  summarise(n=n()) %>%
  mutate(
    called_in_samples=hits + 1,
    events=n / called_in_samples)

ggplot_sample_consistency = function(df) {
  return(
    ggplot(df %>% data.frame() %>%
        filter(caller %in% c("manta", "purple")) %>%
        mutate(
          caller=factor(ifelse(caller=="purple", "gridss", caller), labels=c("gridss", "manta")),
          sample=paste0("Run ", str_extract(sample, ".$"))) %>%
        group_by(caller, sample, hits) %>%
        summarise(n=n())) +
      aes(fill=factor(hits, levels=c(0,1,2), labels=c("1 sample", "2 samples", "3 samples")), y=n, x=sample) +
      geom_bar(stat="identity") +
      facet_grid(~ caller) +
      scale_fill_brewer(palette="Dark2") +
      scale_y_continuous(expand=expand_scale(mult=c(0, 0.02))) +
      theme(
        strip.text.x=element_text(margin=margin(3,0,3,0, "mm")),
        strip.background=element_rect(fill="white")) +
    labs(x="Sequencing run", y="Calls", fill="Called in")
  )
}

ggplot_sample_consistency(bpna12878) +
  labs(title="NA12878 vs self", fill="Called in ", y="Calls")
figsave("manta_vs_gridss_na12878", width=8, height=5)

ggplot_sample_consistency(bpcolo829) +
  labs(title="COLO829 Consistency", fill="Called in ", y="Calls")
figsave("manta_vs_gridss_colo829", width=8, height=5)

