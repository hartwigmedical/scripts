library(GenomicRanges)
library(tidyverse)
library(stringr)
library(grid)
library(gridExtra)
library(cowplot)
library(RMySQL)
source("libSvAnalyser.R")

# Paper common
dropbox = "~/../Dropbox (HMF Australia)/HMF Australia team folder/Structural Variant Analysis"
figdir=paste0(dropbox, "/figures/")
load(paste0(dropbox, "/paper/cohort.RData"))
db = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
#sample_cancer_type_df = query_cancer_type_by_sample(db)
sva_svs = read_csv('../../sv/SVA_SVS.csv')
sva_links = read_csv('../../sv/SVA_LINKS.csv')
sva_clusters = read_csv('../../sv/SVA_CLUSTERS.csv')

cancerTypeOrder = sort(unique(cohort$cancerType))
cancerTypeOrdinal = seq_along(cancerTypeOrder)
names(cancerTypeOrdinal) = cancerTypeOrder

insdf = sva_svs %>%
  filter(SampleId %in% (cohort %>% filter(hpc) %>% pull(sampleId))) %>%
  filter(!is.na(InsertSeq)) %>%
  mutate(
    beSeq = ifelse(OrientStart == 1, InsertSeq, as.character(reverseComplement(DNAStringSet(InsertSeq)))),
    uid=paste0(SampleId, "_", Id))
insdf %>%
  mutate(fq=paste0(">", uid, "\n", beSeq)) %>%
  pull(fq) %>%
  writeLines(con="ins.fa")
# RepeatMasker: RMBlast, default speed, report query aln, rep lower case, skip bacterial, no contamination check, matrix: RM choice
rminsdf = bind_rows(lapply(list.files(path = "../../sv/sgl/", pattern=".*ins.*[.]out$", full.names=TRUE), function(x)
  import.repeatmasker.insertseq(x)))
rminssumdf = repeatmasker.insertseq.summarise(rminsdf)
inslen = str_length(insdf$beSeq)
names(inslen) = insdf$uid
insperbaserm = rminsdf %>%
  group_by(query) %>%
  arrange(qend - qstart) %>%
  do({
    baseRepeats = rep("No repeat", inslen[.$query[1]])
    for (i in seq(nrow(.))) {
      baseRepeats[seq(.$qstart[i], .$qend[i])] = .$repeatClass[i]
    }
    data.frame(
      offset = seq_along(baseRepeats),
      repeatClass = baseRepeats)
  })
norepeathit = insdf %>%
  filter(!(uid %in% rminsdf$query)) %>%
  group_by(uid) %>%
  do({
    data.frame(
      offset=seq(str_length(.$InsertSeq)),
      repeatClass="No repeat")
  })
insperbase = bind_rows(insperbaserm %>% dplyr::select(uid=query, offset, repeatClass), norepeathit) %>%
  mutate(length_bin = ceiling(inslen[insperbase$uid]/binsize)*binsize)

################
# SGL analysis
insperbase_summary = insperbase %>%
  semi_join(insdf %>% filter(Type=="SGL"), by="uid") %>%
  group_by(offset, repeatClass) %>%
  summarise(n=n())
ggplot(insperbase_summary) +
  aes(x=offset, y=n, fill=getRepeatAnn(repeatClass)) +
  geom_bar(stat="identity") +
  scale_fill_brewer(palette="Dark2") +
  labs(y="Breakend count", x="Distance from breakend", fill="RepeatMasker\nAnnotation")
ggsave(paste0(figdir, "sgl_rm_per_base_overall.pdf"), width=8, height=5)

binsize=100
insperbase_binned = insperbase %>%
  semi_join(insdf %>% filter(Type=="SGL"), by="uid") %>%
  mutate(repeatAnn=getRepeatAnn(repeatClass)) %>%
  group_by(offset, repeatAnn, length_bin) %>%
  summarise(n=n())
ggplot(insperbase_binned) +
  aes(x=offset, y=n, fill=repeatAnn) +
  geom_bar(stat="identity") +
  facet_wrap(~ length_bin, scales="free_y") +
  scale_fill_brewer(palette="Dark2") +
  labs(y="Breakend count", x="Distance from breakend", fill="RepeatMasker\nAnnotation")
ggsave(paste0(figdir, "sgl_rm_per_base_binned.pdf"), width=8, height=5)

sgllongrepeat = left_join(
    insdf %>%
      mutate(inslen=str_length(InsertSeq)) %>%
      filter(Type=="SGL") %>%
      dplyr::select(uid, inslen),
    rminssumdf %>%
      mutate(longest_reapeatAnn=getRepeatAnn(longest_repeatClass)) %>%
      dplyr::select(uid=query, longest_reapeatAnn),
    by="uid") %>%
  replace_na(list(longest_reapeatAnn="No repeat"))
ggplot(sgllongrepeat) +
  aes(x=inslen, fill=longest_reapeatAnn) +
  geom_histogram(bins=50) +
  scale_fill_brewer(palette="Dark2") +
  labs(y="Single breakend count", x="Assembled breakend sequence length", fill="RepeatMasker annotation of\nlongest repeat in\nbreakend sequence")
ggsave(paste0(figdir, "sgl_longest_repeat.pdf"), width=8, height=5)
# viral integration?


################
# LINE analysis
linedf = sva_svs %>% filter(ResolvedType=="LINE") %>%
  # filter to highest purity sample per patient
  filter(SampleId %in% (cohort %>% filter(hpc) %>% pull(sampleId)))
# Issue: we shouldn't classify lone polyA SGL as LINE without
# any other supporting evidence
# export all LINE SGL events then run RepeatMasker to see
# how widespread this issue is
line_gr =SVA_SVS_to_gr(linedf)

# site level analysis
max_distance = 5000
cluster_gr = GenomicRanges::reduce(flank(line_gr[line_gr$isLine], max_distance, both=TRUE, ignore.strand=TRUE), ignore.strand=TRUE)
cluster_gr$hits = countOverlaps(cluster_gr, line_gr, ignore.strand=TRUE)

# Insertion sites

per_sample_df = linedf %>%
  group_by(SampleId) %>%
  summarise(
    sites=length(unique(ClusterId)),
    breaks=n()) %>%
  rbind(data.frame(SampleId=cohort$sampleId, sites=0, breaks=0)) %>%
  group_by(SampleId) %>%
  summarise(sites=sum(sites), breaks=sum(breaks)) %>%
  left_join(cohort, by=c("SampleId"="sampleId"))

ggplot(per_sample_df) +
  aes(x=cancerType, y=log10(breaks+1)) +
  geom_violin(scale="width")
  #geom_violin(aes(y=log10(breaks+1)), scale="width", colour="red", alpha=0.5)
  #geom_jitter()

ggplot(per_sample_df) +
  aes(x=sites, y=breaks/sites, colour=cancerType) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  labs("LINE insertions")

# Clean point-based CDF
per_sample_df %>%
  arrange(breaks) %>%
  group_by(cancerType) %>%
  mutate(x_offset = cancerTypeOrdinal[cancerType] + row_number() / n()) %>%
ggplot() +
  geom_segment(aes(
      x=x_offset + 0.1,
      xend=x_offset + 0.9,
      y=median,
      yend=median),
    data=per_sample_df %>%
      group_by(cancerType) %>%
      summarise(median=median(breaks+1)) %>%
      mutate(x_offset = cancerTypeOrdinal[cancerType]),
    size=1,
    alpha=0.4) +
  geom_point(aes(x=x_offset, y=breaks + 1), size=0.1) +
  scale_y_log10(breaks=c(2,11,101,1001), labels=c(1,10,100,1000)) +
  scale_x_continuous(expand=c(0.001,0), breaks = cancerTypeOrdinal, labels=cancerTypeOrder, minor_breaks=cancerTypeOrdinal) +
  labs(x="", y="LINE insertions") +
  theme(
    axis.text.x = element_text(angle = 90, hjust=1, vjust=1.6),
    panel.grid.major.x = element_line(colour=alpha("black", 0.05)))
ggsave(paste0(figdir, "/line_cancer_type_cdf.pdf"), width=8, height=3.2)


