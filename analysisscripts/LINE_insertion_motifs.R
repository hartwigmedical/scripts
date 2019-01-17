library(tidyverse)
library(cowplot)
library(GenomicRanges)
library(StructuralVariantAnnotation)

cluster_raw_df = read_csv('~/../Dropbox (HMF Australia)/HMF Australia team folder/Structural Variant Analysis/CLUSTER.csv')
full_gr = c(with(cluster_raw_df, GRanges(
  seqnames=ChrStart,
  ranges=IRanges(start=PosStart, width=1),
  strand=ifelse(OrientStart == -1, "-", "+"),
  Id=Id,
  beid=paste0(Id, ifelse(ChrEnd != 0, "o", "b")),
  SampleId=SampleId,
  isLine=LEStart != "None",
  partner=ifelse(ChrEnd != 0, paste0(Id, "h"), NA)
)), with(line_raw_df %>% filter(ChrEnd != 0), GRanges(
  seqnames=ChrEnd,
  ranges=IRanges(start=PosEnd, width=1),
  strand=ifelse(OrientEnd == -1, "-", "+"),
  Id=Id,
  beid=paste0(Id, "h"),
  SampleId=SampleId,
  isLine=LEEnd != "None",
  partner=paste0(Id, "o")
)))
line_raw_df = cluster_raw_df %>% filter(LEStart != "None" | LEEnd != "None")
line_gr = c(with(line_raw_df, GRanges(
    seqnames=ChrStart,
    ranges=IRanges(start=PosStart, width=1),
    strand=ifelse(OrientStart == -1, "-", "+"),
    Id=Id,
    beid=paste0(Id, ifelse(ChrEnd != 0, "o", "b")),
    SampleId=SampleId,
    isLine=LEStart != "None",
    partner=ifelse(ChrEnd != 0, paste0(Id, "h"), NA)
  )), with(line_raw_df %>% filter(ChrEnd != 0), GRanges(
    seqnames=ChrEnd,
    ranges=IRanges(start=PosEnd, width=1),
    strand=ifelse(OrientEnd == -1, "-", "+"),
    Id=Id,
    beid=paste0(Id, "h"),
    SampleId=SampleId,
    isLine=LEEnd != "None",
    partner=paste0(Id, "o")
  )))
names(line_gr) = line_gr$beid
bp_line_gr = line_gr[!is.na(line_gr$partner)]

intra_line = line_raw_df %>%
  filter(LEStart != "None" & LEEnd != "None") %>%
  filter(ChrStart==ChrEnd) %>%
  mutate(distance_between_bp=abs(PosStart - PosEnd)) %>%
  filter(distance_between_bp < 10000)
table(intra_line$Type)
ggplot(intra_line) +
  aes(x=distance_between_bp, fill=Type) +
  geom_histogram() +
  scale_x_log10() +
  labs("Distance between breakends on intra-LINE inversions")
# TODO: do these intra-LINE elements have corresponding insertion sites that inidicate a -- or ++ insertion site?

bpbp_line_df = findOverlaps(bp_line_gr, bp_line_gr, maxgap=100, ignore.strand=TRUE) %>%
  as.data.frame() %>%
  filter(as.logical(
    bp_line_gr$SampleId[queryHits] == bp_line_gr$SampleId[subjectHits] &
    strand(bp_line_gr)[queryHits] == "+" & strand(bp_line_gr)[subjectHits] == "-" &
    !bp_line_gr$isLine[queryHits] & !bp_line_gr$isLine[subjectHits] &
    seqnames(partner(bp_line_gr)[queryHits]) == seqnames(partner(bp_line_gr)[subjectHits]) &
    abs(start(partner(bp_line_gr)[queryHits]) - start(partner(bp_line_gr)[subjectHits])) < 50000)) %>%
  mutate(deleted_bases=start(bp_line_gr)[subjectHits] - start(bp_line_gr)[queryHits] - 1) %>%
  mutate(SampleId=bp_line_gr$SampleId[queryHits])

ggplot(bpbp_line_df) +
  aes(x=deleted_bases) +
  geom_histogram(bins=60) +
  scale_x_continuous(limits=c(-30, 30))

ggplot(bpbp_line_df %>%
  group_by(SampleId) %>%
  summarise(
    mean=mean(deleted_bases),
    n=n(),
    lt5=sum(deleted_bases < -5),
    gt5=sum(deleted_bases >= -5),
    lt0=sum(deleted_bases < 0),
    gt0=sum(deleted_bases >= 0))) +
  aes(x=n, y=gt0) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="LINE counts by sample")

ggplot(bpbp_line_df %>% filter(SampleId %in% (bpbp_line_df %>% group_by(SampleId) %>% summarise(n=n()) %>% filter(n >= 100) %>% pull(SampleId)))) +
  aes(x=deleted_bases, fill=SampleId) +
  geom_histogram(bins=60) +
  scale_x_continuous(limits=c(-30, 30)) +
  labs(title="LINE counts by sample (min 20 events)")







