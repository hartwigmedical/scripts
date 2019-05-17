library(GenomicRanges)
library(tidyverse)
library(stringr)
library(grid)
library(gridExtra)
library(cowplot)
library(RMySQL)
library(igraph)
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



#
# LINE analysis
linedf = sva_svs %>% filter(ResolvedType=="LINE") %>%
  # filter to highest purity sample per patient
  filter(SampleId %in% (cohort %>% filter(hpc) %>% pull(sampleId)))
# Issue: we shouldn't classify lone polyA SGL as LINE without
# any other supporting evidence
# export all LINE SGL events then run RepeatMasker to see
# how widespread this issue is
linedf %>%
  filter(PosEnd==-1) %>%
  mutate(fq=paste0(">", SampleId, "_", Id, "\n", InsertSeq)) %>%
  pull(fq) %>%
  writeLines(con="line_single_breakends.fa")
line_gr = SVA_SVS_to_gr(linedf)

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

