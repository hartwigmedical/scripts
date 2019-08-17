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
  if (length(gr) > 0) {
    names(gr) = paste0(sample, caller, names(gr))
    gr$partner = paste0(sample, caller, gr$partner)
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
passingbegr = begr[begr$FILTER == "PASS"]
passingbegr$partner = NA_character_

bpna12878 = hits_by_caller(passingbpgr[passingbpgr$sample %in% c("GIABvsSELFv001", "GIABvsSELFv002", "GIABvsSELFv003")], countBreakpointOverlaps)
bena12878 = hits_by_caller(passingbegr[passingbegr$sample %in% c("GIABvsSELFv001", "GIABvsSELFv002", "GIABvsSELFv003")], countOverlaps)
grna12878 = c(bpna12878, bena12878)
grna12878$betp = FALSE
grna12878$bptp = FALSE

becolo829 = hits_by_caller(passingbegr[passingbegr$sample %in% c("COLO829v001", "COLO829v002", "COLO829v003")], countOverlaps)
bpcolo829 = hits_by_caller(passingbpgr[passingbpgr$sample %in% c("COLO829v001", "COLO829v002", "COLO829v003")], countBreakpointOverlaps)
colo829truth = readVcf(paste0(data_dir, "/../colo829/COLO829.somatic.overlapped.vcf"))
bpcolo829t = breakpointRanges(colo829truth, inferMissingBreakends=TRUE)
becolo829t = breakendRanges(colo829truth)
bpcolo829t$bpid = ifelse(names(bpcolo829t) < bpcolo829t$partner, names(bpcolo829t), bpcolo829t$partner)
becolo829$betp = overlapsAny(becolo829, c(bpcolo829t, becolo829t), maxgap=4)
bpcolo829$betp = overlapsAny(bpcolo829, c(bpcolo829t, becolo829t), maxgap=4)
becolo829$bptp = FALSE
bpcolo829$bptp = countBreakpointOverlaps(bpcolo829, bpcolo829t, maxgap=4) > 0
becolo829$bpid = bpcolo829t$bpid[findOverlaps(becolo829, bpcolo829t, select="first")]
bpcolo829$bpid = bpcolo829t$bpid[findOverlaps(bpcolo829, bpcolo829t, select="first")]
grcolo829 = c(bpcolo829, becolo829)

summary_na12878 = grna12878 %>%
  as.data.frame() %>%
  group_by(caller, hits) %>%
  summarise(n=n()) %>%
  mutate(
    called_in_samples=hits + 1,
    events=n / called_in_samples)
summary_colo829 = grcolo829 %>%
  as.data.frame() %>%
  group_by(caller, hits) %>%
  summarise(n=n()) %>%
  mutate(
    called_in_samples=hits + 1,
    events=n / called_in_samples)
ggplot_fp_df = function(gr) {
  gr %>% data.frame() %>%
    filter(caller %in% c("manta", "purple")) %>%
    filter(!betp & !bptp) %>%
    filter(!is.na(partner)) %>%
    mutate(
      caller=factor(ifelse(caller=="purple", "gridss", caller), labels=c("gridss", "manta")),
      sample=paste0("Run ", str_extract(sample, ".$"))) %>%
    group_by(caller, sample, hits) %>%
    summarise(n=n() / 2)
}
ggplot_fp_consistency = function(df) {
  return(
    ggplot(df) +
      aes(fill=factor(hits, levels=c(0,1,2), labels=c("1 sample", "2 samples", "3 samples")), y=n, x=sample) +
      geom_bar(stat="identity") +
      facet_grid(~ caller) +
      scale_fill_brewer(palette="Set2") +
      scale_y_continuous(expand=expand_scale(mult=c(0, 0.02))) +
      theme(
        strip.text.x=element_text(margin=margin(3,0,3,0, "mm")),
        strip.background=element_rect(fill="white")) +
    labs(x="Sequencing run", y="False Positives", fill="Called in")
  )
}
plot_data = ggplot_fp_df(grna12878)
ggplot_fp_consistency(plot_data) +
  labs(title="NA12878 vs self")
figsave("manta_vs_gridss_fp_na12878", width=8, height=5)

plot_data = ggplot_fp_df(grcolo829)
ggplot_fp_consistency(plot_data) +
  labs(title="COLO829 Consistency")
figsave("manta_vs_gridss_fp_colo829", width=8, height=5)

grcolo829 %>% as.data.frame() %>%
  group_by(sample, caller, bpid) %>%
  summarise(hits=max(hits)) %>%
  group_by(sample, caller, hits) %>%
  # score a breakend that matches a breakpoint as a match
  summarise(n=length(unique(bpid))) %>%
  ungroup() %>%
  filter(caller %in% c("manta", "purple")) %>%
  mutate(
    caller=factor(ifelse(caller=="purple", "gridss", caller), labels=c("gridss", "manta")),
    sample=paste0("Run ", str_extract(sample, ".$"))) %>%
  # hack in a zero row so the legend has three categories
  rbind(data.frame(sample="Run 1", caller="gridss", hits=0, n=0)) %>%
ggplot() +
  aes(fill=factor(hits, levels=c(0,1,2), labels=c("1 sample", "2 samples", "3 samples")), y=n, x=sample) +
  geom_bar(stat="identity") +
  facet_grid(~ caller) +
  scale_fill_brewer(palette="Set2") +
  scale_y_continuous(
    limits=c(0, length(bpcolo829t) / 2),
    expand=expand_scale(mult=c(0, 0.02)),
    sec.axis=sec_axis(
      trans=function(x) x / length(bpcolo829t) * 2 * 100,
      breaks=1:5*20,
      label=paste0(1:5*20, "%")))+
  theme(
    strip.text.x=element_text(margin=margin(3,0,3,0, "mm")),
    strip.background=element_rect(fill="white")) +
  labs(x="Sequencing run", y="True Positives", fill="Called in")
figsave("manta_vs_gridss_tp_colo829", width=8, height=5)

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
figsave("probe_results_by_sample", width=5, height=6)


# GRIDSS vs manta
rbind(
  probeResult %>% filter(scope!="SharedBoth"),
  probeResult %>% filter(scope=="SharedBoth") %>% mutate(scope="SharedManta"),
  probeResult %>% filter(scope=="SharedStrelka") %>% mutate(callset="Gridss", scope="Private")) %>%
  filter(probeQuality >= 20) %>%
  filter(!is.na(source)) %>%
  filter(source!="CPCT02450014T") %>%
  mutate(category=factor(paste(callset, scope),
                         levels=c("Gridss Private", "Gridss SharedManta", "Manta Private"),
                         labels=c("gridss", "gridss+manta", "manta"))) %>%
  mutate(under50bp=!is.na(length) & abs(length) <= 50) %>%
  filter(!under50bp & !is.na(category)) %>%
  group_by(category, supported) %>%
  summarise(n=n()) %>%
ggplot() +
  aes(x=category, fill=supported, y=n) +
  geom_bar(stat="identity") +
  theme_bw() +
  scale_fill_manual(values=c("#d95f02", "#1b9e77")) +
  labs(fill="Probe\nValidated", x="", y="Variant calls (50bp+)") +
  theme(axis.text.x = element_text(angle = 90))
figsave("probe_results_vs_manta_over_50bp", width=5, height=4)
rbind(
  probeResult %>% filter(scope!="SharedBoth"),
  probeResult %>% filter(scope=="SharedBoth") %>% mutate(scope="SharedStrelka"),
  probeResult %>% filter(scope=="SharedManta") %>% mutate(callset="Gridss", scope="Private")) %>%
  filter(probeQuality >= 20) %>%
  filter(!is.na(source)) %>%
  filter(source!="CPCT02450014T") %>%
  mutate(category=factor(paste(callset, scope),
                         levels=c("Gridss Private", "Gridss SharedStrelka", "Strelka Private"),
                         labels=c("gridss", "gridss+strelka", "strelka"))) %>%
  mutate(under50bp=!is.na(length) & abs(length) <= 50) %>%
  filter(under50bp & !is.na(category)) %>%
  group_by(category, supported) %>%
  summarise(n=n()) %>%
ggplot() +
  aes(x=category, fill=supported, y=n) +
  geom_bar(stat="identity") +
  theme_bw() +
  scale_fill_manual(values=c("#d95f02", "#1b9e77")) +
  labs(fill="Probe\nValidated", x="", y="Variant calls (32-50bp)") +
  theme(axis.text.x = element_text(angle = 90))
figsave("probe_results_vs_stelka_under_50bp", width=5, height=4)

sumdf = probeResult %>%
  filter(probeQuality >= 20) %>%
  filter(source!="CPCT02450014T") %>%
  group_by(callset, scope) %>%
  summarise(calls=n(), validated=sum(supported)) %>%
  ungroup()
rbind(
  sumdf %>%
    filter(str_detect(paste0(callset,scope), "Manta")) %>%
    summarise(calls=sum(calls), validated=sum(validated)) %>%
    mutate(caller="Manta"),
  sumdf %>%
    filter(str_detect(paste0(callset,scope), "Gridss")) %>%
    summarise(calls=sum(calls), validated=sum(validated)) %>%
    mutate(caller="gridss"),
  sumdf %>%
    filter(str_detect(paste0(callset,scope), "Strelka")) %>%
    summarise(calls=sum(calls), validated=sum(validated)) %>%
    mutate(caller="Strelka")) %>%
  mutate(
    prec=validated/calls,
    fdr=1-validated/calls)

shortsumdf = probeResult %>%
  filter(probeQuality >= 20) %>%
  filter(source!="CPCT02450014T") %>%
  mutate(under50bp=!is.na(length) & abs(length) <= 50) %>%
  filter(under50bp) %>%
  group_by(callset, scope) %>%
  summarise(calls=n(), validated=sum(supported)) %>%
  ungroup()
rbind(
  shortsumdf %>%
    filter(str_detect(paste0(callset,scope), "Manta")) %>%
    summarise(calls=sum(calls), validated=sum(validated)) %>%
    mutate(caller="Manta"),
  shortsumdf %>%
    filter(str_detect(paste0(callset,scope), "Gridss")) %>%
    summarise(calls=sum(calls), validated=sum(validated)) %>%
    mutate(caller="gridss"),
  shortsumdf %>%
    filter(str_detect(paste0(callset,scope), "Strelka")) %>%
    summarise(calls=sum(calls), validated=sum(validated)) %>%
    mutate(caller="Strelka")) %>%
  mutate(
    prec=validated/calls,
    fdr=1-validated/calls)
