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
  bpgr = c(
    ann_caller_sample(manta_file, "manta", sample, breakpointRanges),
    ann_caller_sample(gridss_file, "gridss", sample, breakpointRanges),
    ann_caller_sample(purple_file, "purple", sample, breakpointRanges),
    bpgr)
  begr = c(
    ann_caller_sample(gridss_file, "gridss", sample, breakendRanges),
    ann_caller_sample(purple_file, "purple", sample, breakendRanges),
    begr)
}
hits_by_caller = function(bpgr, countFunction) {
  bpgr$hits = rep(0, length(bpgr))
  for (caller in unique(bpgr$caller)) {
    callergr = bpgr[bpgr$caller == caller]
    if (length(callergr) > 0) {
      for (sample in unique(callergr$sample)) {
        isCurrent = bpgr$caller == caller & bpgr$sample == sample
        if (any(isCurrent)) {
          hits = bpgr$hits[isCurrent]
          for (otherSample in unique(callergr$sample)) {
            if (sample != otherSample) {
              othergr = bpgr[bpgr$caller == caller & bpgr$sample == otherSample]
              hits = hits + (countFunction(bpgr[isCurrent], othergr) > 0)
            }
          }
          bpgr$hits[isCurrent] = hits
        }
      }
    }
  }
  return(bpgr)
}
passingbpgr = bpgr[bpgr$FILTER == "PASS"]
bpna12878 = hits_by_caller(passingbpgr[passingbpgr$sample %in% c("GIABvsSELFv001", "GIABvsSELFv002", "GIABvsSELFv003")], countBreakpointOverlaps)
bpcolo829 = hits_by_caller(passingbpgr[passingbpgr$sample %in% c("COLO829v001", "COLO829v002", "COLO829v003")], countBreakpointOverlaps)
passingbegr = begr[begr$FILTER == "PASS"]
passingbegr$partner = NA_character_
bena12878 = hits_by_caller(passingbegr[passingbegr$sample %in% c("GIABvsSELFv001", "GIABvsSELFv002", "GIABvsSELFv003")], countOverlaps)
becolo829 = hits_by_caller(passingbegr[passingbegr$sample %in% c("COLO829v001", "COLO829v002", "COLO829v003")], countOverlaps)
grna12878 = c(bpna12878, bena12878)
grcolo829 = c(bpcolo829, becolo829)

summary_na12878 = grna12878 %>%
    group_by(caller, hits) %>%
    summarise(n=n()) %>%
    mutate(
      called_in_samples=hits + 1,
      events=n / called_in_samples)
summary_colo829 = grcolo829 %>%
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

ggplot_sample_consistency(grna12878) +
  labs(title="NA12878 vs self", fill="Called in ", y="Calls")
figsave("manta_vs_gridss_na12878", width=8, height=5)

ggplot_sample_consistency(grcolo829) +
  labs(title="COLO829 Consistency", fill="Called in ", y="Calls")
figsave("manta_vs_gridss_colo829", width=8, height=5)

ggplot(passingbpgr %>%
    as.data.frame(row.names=NULL) %>%
    mutate(precise=end - start - passingbpgr$HOMLEN == 0) %>%
    group_by(caller, sample) %>%
    filter(
      sample %in% c("COLO829v001", "COLO829v002", "COLO829v003") &
      caller %in% c("manta", "purple")) %>%
    ungroup() %>%
    mutate(caller=ifelse(caller=="purple", "gridss", caller)) %>%
    group_by(caller) %>%
    summarise(prec=sum(precise)/n())) +
  aes(x=caller, y=prec) +
  geom_bar(stat="identity") +
  scale_y_continuous(labels = percent_format()) +
  labs(title="COLO829", x="", y="Base pair precise")
figsave("manta_vs_gridss_call_percent_precise", width=3, height=3)

#
# Probe validation results
#

load(paste0(data_dir, "/probeResult.RData"))


rbind(
    probeResult %>% filter(scope!="SharedBoth"),
    probeResult %>% filter(scope=="SharedBoth") %>% mutate(scope="SharedManta"),
    probeResult %>% filter(scope=="SharedBoth") %>% mutate(scope="SharedStrelka")) %>%
  filter(probeQuality >= 20) %>%
  filter(!is.na(source)) %>%
  group_by(source, callset, scope, supported) %>%
  summarise(n=n()) %>%
  mutate(category=factor(paste(callset, scope),
    levels=c("Manta Private", "Gridss SharedManta", "Gridss Private", "Gridss SharedStrelka", "Strelka Private"),
    labels=c("manta", "gridss+manta", "gridss", "gridss+strelka", "strelka (32bp+)"))) %>%
ggplot() +
  aes(x=source, fill=supported, y=n) +
  geom_bar(stat="identity") +
  facet_grid(category ~ ., scales="free") +
  theme_bw() +
  scale_fill_manual(values=c("#d95f02", "#1b9e77")) +
  labs(fill="Validated", x="", y="breakpoint calls") +
  theme(axis.text.x = element_text(angle = 90))
figsave("probe_results_raw", width=5, height=6)


