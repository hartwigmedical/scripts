library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(StructuralVariantAnnotation)
options(stringsAsFactors = FALSE)
theme_set(theme_bw())

#### PCAWG ####
data_dir="D:/dev/gpl_benchmarking/pcawgcnsv"
pcawg_evaluate_cn_transitions = function(sampleId, ...) {
  write(paste("Processing ", sampleId), stderr())
  cnfile=paste0(data_dir, "/", sampleId, ".consensus.20170119.somatic.cna.txt")
  svfile=paste0(data_dir, "/", sampleId, ".pcawg_consensus_1.6.161116.somatic.sv.bedpe")
  if (!file.exists(cnfile)) {
    write(paste("Missing ", cnfile), stderr())
    return(data.frame())
  }
  if (!file.exists(svfile)) {
    write(paste("Missing ", svfile), stderr())
    return(data.frame())
  }
  sv_bedpe = read_delim(svfile, delim="\t", col_names=TRUE, col_types="cnncnncncccc")
  cndf = read_delim(cnfile, delim="\t", col_names=TRUE, col_types="cnnnnnn", na="NA")
  if (nrow(sv_bedpe) == 0 || nrow(cndf) ==0) {
    return(data.frame())
  }
  svgr = with(sv_bedpe, GRanges(
    seqnames=c(chrom1, chrom2),
    ranges=IRanges(start=c(start1 + 1, start2 + 1), end=c(end1, end2)),
    strand=c(strand1, strand2),
    sampleId=sampleId,
    sourceId=c(paste0(sampleId, sv_id, "_o"), paste0(sampleId, sv_id, "_h")),
    partner=c(paste0(sampleId, sv_id, "_h"), paste0(sampleId, sv_id, "_o"))
  ))
  names(svgr) = svgr$sourceId
  cngr = with(cndf, GRanges(
    seqnames=chromosome,
    ranges=IRanges(start=start, end=end),
    sampleId=sampleId,
    cn=total_cn,
    cn_major=major_cn,
    cn_minor=minor_cn,
    star=star))
  cn_consistency = evaluate_cn_transitions(cngr, svgr, ...)
  cn_consistency$sampleId = sampleId
  as.data.frame(cn_consistency)
}

sampleIds = unique(str_extract(list.files(path=data_dir), "^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"))
pcawg_cn_transitions_list = lapply(sampleIds, pcawg_evaluate_cn_transitions, distance="cn_transition")
pcawg_sv_transitions_list = lapply(sampleIds, pcawg_evaluate_cn_transitions, distance="sv")
# batched version so we don't lose everything when a file is bad
#pcawg_cn_transitions_list = list()
#batch_size = 20
#while (length(sampleIds) - length(pcawg_cn_transitions_list) > 0) {
#  batch_sampleIds = sampleIds[!(sampleIds %in% names(pcawg_cn_transitions_list))]
#  batch_sampleIds = batch_sampleIds[1:min(length(batch_sampleIds), batch_size)]
#  batch_pcawg_evaluate_cn_transitions = lapply(batch_sampleIds, pcawg_evaluate_cn_transitions)
#  names(batch_pcawg_evaluate_cn_transitions) = batch_sampleIds
#  pcawg_cn_transitions_list = c(pcawg_cn_transitions_list, batch_pcawg_evaluate_cn_transitions)
#}
pcawg_cn_transitions = bind_rows(pcawg_cn_transitions_list[unlist(sapply(pcawg_cn_transitions_list, nrow)) > 0])
pcawg_sv_transitions = bind_rows(pcawg_sv_transitions_list[unlist(sapply(pcawg_sv_transitions_list, nrow)) > 0])

#### Hartwig ####
library(tidyverse)
library(readr)
lnx_svs = read_csv("D:/hartwig/sv/paper_hpc/LNX_SVS.csv")
lnx_cns = read_tsv("D:/hartwig/sv/paper_hpc/LNX_VIS_COPY_NUMBER.tsv")
hartwig_svgr = lnx_svs %>%
  filter(!(ChrEnd == "0" & is.na(InsertSeq))) %>% # strip purple placeholder calls
  lnx_to_gr()
hartwig_cngr = with(lnx_cns, GRanges(
  seqnames=Chromosome,
  ranges=IRanges(start=Start, end=End),
  cn=CopyNumber,
  cn_major=CopyNumber * BAF,
  cn_minor=CopyNumber * (1 - BAF),
  sampleId=SampleId))
hartwig_evaluate_cn_transitions = function(sampleId, ...) {
  write(paste("Processing ", sampleId), stderr())
  svgr = hartwig_svgr[hartwig_svgr$SampleId == sampleId,]
  cngr = hartwig_cngr[hartwig_cngr$sampleId == sampleId,]
  cn_consistency = evaluate_cn_transitions(cngr, svgr, ...)
  cn_consistency$sampleId = sampleId
  as.data.frame(cn_consistency)
}
hartwig_cn_transitions_list = lapply(unique(lnx_svs$SampleId), hartwig_evaluate_cn_transitions)
hartwig_cn_transitions = bind_rows(hartwig_cn_transitions_list[unlist(sapply(hartwig_cn_transitions_list, nrow)) > 0])
hartwig_sv_transitions_list = lapply(unique(lnx_svs$SampleId), hartwig_evaluate_cn_transitions, distance="sv")
hartwig_sv_transitions = bind_rows(hartwig_sv_transitions_list[unlist(sapply(hartwig_sv_transitions_list, nrow)) > 0])


#### Analysis ####

cn_transitions = bind_rows(
    pcawg_cn_transitions %>% mutate(cohort="PCAWG"),
    hartwig_cn_transitions %>% mutate(cohort="Hartwig")) %>%
  mutate(
    inGap = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_gaps, maxgap=100000),
    inCentromere = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_centromeres, maxgap=100000),
    nearGapOrCentromere=inGap | inCentromere)

sv_transitions = bind_rows(
  pcawg_sv_transitions %>% mutate(cohort="PCAWG"),
  hartwig_sv_transitions %>% mutate(cohort="Hartwig")) %>%
  mutate(
    inGap = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_gaps, maxgap=100000),
    inCentromere = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_centromeres, maxgap=100000),
    nearGapOrCentromere=inGap | inCentromere,
    distance=ifelse(is.na(distance), 100000, distance))

ggplot(sv_transitions) +
  aes(x=distance, fill=cohort) +
  geom_histogram(position=position_dodge()) +
  scale_x_log10() +
  scale_y_log10()

sv_transitions %>%
  filter(!nearGapOrCentromere) %>%
  mutate(match_cn=distance <= 2 & !is.na(distance)) %>%
  group_by(cohort) %>%
  summarise(
    percent=1-sum(match_cn) / n(),
    match_cn=sum(match_cn),
    n=n())


missed_summary =
  cn_transitions %>%
  filter(!is.na(cn_left) & !is.na(cn_right)) %>%
  mutate(missed_sv = sv_matches == 0) %>%
  group_by(cohort, nearGapOrCentromere) %>%
  summarise(portion_missed=sum(missed_sv) / n(), transitions=n())

missed_clonal_summary =
  cn_transitions %>%
  filter(!is.na(cn_left) & !is.na(cn_right) & abs(cn_left - cn_right) > 0.5) %>%
  mutate(missed_sv = sv_matches == 0) %>%
  group_by(cohort, nearGapOrCentromere) %>%
  summarise(portion_missed=sum(missed_sv) / n(), transitions=n())

cn_transitions %>%
  group_by(cohort) %>%
  filter(!is.na(distance)) %>%
  summarise(percent_matching_position=sum(distance <= 1) / n())

ggplot(cn_transitions %>% filter(can_evaluate_cn_error)) +
  aes(x=pmin(cn_error, 8)) +
  facet_wrap(~ cohort, scales="free") +
  geom_histogram() +
  labs("Magnitude of copy number transition inconsistency")

cn_transitions %>%
  filter(can_evaluate_cn_error) %>%
  filter(!is.na(cn_error)) %>%
  mutate(rounded_cn_error=pmin(round(cn_error), 8)) %>%
  group_by(cohort, rounded_cn_error) %>%
  summarise(count=n()) %>%
  group_by(cohort) %>%
  mutate(percent=1-count/sum(count))

ggplot(cn_transitions %>%
    filter(can_evaluate_cn_error) %>%
    mutate(rounded_cn_error=pmin(round(cn_error), 8)) %>%
    group_by(cohort, rounded_cn_error) %>%
    summarise(count=n()) %>%
    group_by(cohort) %>%
    mutate(percent=1-count/sum(count))) +
  aes(x=rounded_cn_error, y=percent, fill=cohort) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(
    title="Copy number consistency",
    x="Magnitude of copy number transition inconsistency",
    y="Portion of copy number transitions\nwith a single supporting structural variant")

ggplot(cn_transitions %>% filter(!is.na(distance))) +
  aes(x=distance) +
  facet_wrap(~ cohort, scales="free") +
  geom_histogram() +
  scale_x_log10()

cn_transitions %>%
  filter(!is.na(cn_left) & !is.na(cn_right)) %>%
  group_by(cohort, sv_matches, nearGapOrCentromere) %>%
  summarise(n=n()) %>%
  mutate(portion=n/sum(n))

cn_transitions %>%
  filter(!is.na(cn_left) & !is.na(cn_right) & abs(cn_delta) > 0) %>%
  group_by(sv_matches, cohort) %>%
  summarise(n=n()) %>%
  group_by(cohort) %>%
  mutate(portion=n/sum(n))

missed_by_sample = cn_transitions %>%
  filter(!is.na(cn_left) & !is.na(cn_right)) %>%
  mutate(missed_sv = sv_matches == 0) %>%
  group_by(sampleId) %>%
  summarise(portion_missed=sum(missed_sv) / n(), transitions=n()) %>%
  arrange(desc(portion_missed))

ggplot(missed_by_sample) +
  aes(x=portion_missed) +
  geom_histogram() +
  labs(title="Portion of CN transitions with no corresponding SV")
ggplot(missed_by_sample) +
  aes(x=portion_missed, y=transitions) +
  geom_point() +
  scale_y_log10() +
  labs(title="Portion of CN transitions with no corresponding SV")

missed_by_sample_excluding_cndelta0 = pcawg_cn_transitions %>%
  filter(!is.na(cn_left) & !is.na(cn_right) & abs(cn_delta) > 0) %>%
  mutate(missed_sv = sv_matches == 0) %>%
  group_by(sampleId) %>%
  summarise(portion_missed=sum(missed_sv) / n(), transitions=n()) %>%
  arrange(desc(portion_missed))
ggplot(missed_by_sample_excluding_cndelta0) +
  aes(x=portion_missed, y=transitions) +
  geom_point() +
  scale_y_log10() +
  labs(title="Portion of CN transitions with no corresponding SV")

ggplot(pcawg_cn_transitions %>% filter(sv_matches)) +
  aes(x=distance, fill="") +
  geom_histogram() +
  scale_y_log10() +
  scale_x_log10()


