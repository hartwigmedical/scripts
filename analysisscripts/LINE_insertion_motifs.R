library(tidyverse)
library(cowplot)
library(GenomicRanges)
library(StructuralVariantAnnotation)

#cluster_raw_df = read_csv('~/../Dropbox (HMF Australia)/HMF Australia team folder/Structural Variant Analysis/SVA_CLUSTERS.csv')
sv_raw_df = read_csv('~/../Dropbox (HMF Australia)/HMF Australia team folder/Structural Variant Analysis/SVA_SVS.csv')
row.names(sv_raw_df) = sv_raw_df$Id
full_gr = c(with(sv_raw_df, GRanges(
  seqnames=ChrStart,
  ranges=IRanges(start=PosStart, width=1),
  strand=ifelse(OrientStart == -1, "-", "+"),
  Id=Id,
  beid=paste0(Id, ifelse(ChrEnd != 0, "o", "b")),
  SampleId=SampleId,
  isLine=LEStart != "None",
  partner=ifelse(ChrEnd != 0, paste0(Id, "h"), NA),
  Homology=Homology,
  ihomlen=InexactHOEnd-InexactHOStart,
  insSeq=InsertSeq,
  qual=QualScore,
  cn=Ploidy,
  refContext=RefContextStart
)), with(sv_raw_df %>% filter(ChrEnd != 0), GRanges(
  seqnames=ChrEnd,
  ranges=IRanges(start=PosEnd, width=1),
  strand=ifelse(OrientEnd == -1, "-", "+"),
  Id=Id,
  beid=paste0(Id, "h"),
  SampleId=SampleId,
  isLine=LEEnd != "None",
  partner=paste0(Id, "o"),
  Homology=Homology,
  ihomlen=InexactHOEnd-InexactHOStart,
  insSeq=InsertSeq,
  qual=QualScore,
  cn=Ploidy,
  refContext=RefContextEnd
)))
names(full_gr) = full_gr$beid
line_raw_df = sv_raw_df %>% filter(LEStart != "None" | LEEnd != "None")
line_gr = full_gr[full_gr$Id %in% line_raw_df$Id]
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

bpbp_line_df = findOverlaps(bp_line_gr, bp_line_gr, maxgap=50, ignore.strand=TRUE) %>%
  as.data.frame() %>%
  filter(as.logical(
    bp_line_gr$SampleId[queryHits] == bp_line_gr$SampleId[subjectHits] &
    strand(bp_line_gr)[queryHits] == "+" & strand(bp_line_gr)[subjectHits] == "-" &
    !bp_line_gr$isLine[queryHits] & !bp_line_gr$isLine[subjectHits] &
    seqnames(partner(bp_line_gr)[queryHits]) == seqnames(partner(bp_line_gr)[subjectHits]) &
    abs(start(partner(bp_line_gr)[queryHits]) - start(partner(bp_line_gr)[subjectHits])) < 50000)) %>%
  mutate(deleted_bases=start(bp_line_gr)[subjectHits] - start(bp_line_gr)[queryHits] - 1) %>%
  mutate(SampleId=bp_line_gr$SampleId[queryHits]) %>%
  mutate(line_chr=as.character(seqnames(partner(bp_line_gr)[queryHits])),
    line_pos1 = start(partner(bp_line_gr)[queryHits]),
    line_pos2 = start(partner(bp_line_gr)[subjectHits])) %>%
  mutate(line_length=abs(line_pos1-line_pos2))

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

ggplot(bpbp_line_df %>% filter(SampleId %in% (bpbp_line_df %>% group_by(SampleId) %>% summarise(n=n()) %>% filter(n >= 200) %>% pull(SampleId)))) +
  aes(x=deleted_bases) +
  geom_histogram(bins=60) +
  facet_wrap(~SampleId) +
  scale_x_continuous(limits=c(-20, 11)) +
  labs(title="LINE counts by sample (min 100 events)")

density = countOverlaps(line_gr, line_gr, maxgap=1000)
max_gr = line_gr[density == max(density)][1]
ggplot(bpbp_line_df %>% filter(line_chr == as.character(seqnames(max_gr)) &
    ((abs(line_pos1 - start(max_gr)) <= 2000) | abs(line_pos2 - start(max_gr)) <= 2000)) %>%
      mutate(deleted_bases_jittered=deleted_bases-0.5+row_number()/n())) +
  aes(x=pmin(line_pos1, line_pos2), xend=pmax(line_pos1, line_pos2), y=deleted_bases_jittered, yend=deleted_bases_jittered) +
  geom_segment() +
  coord_cartesian(xlim=c(start(max_gr) - 1000, end(max_gr) + 1000), ylim=c(-20, 11)) +
  labs(title="Snapshot of most inserted LINE", y="deleted bases", x="Genomic position")

ggplot(bpbp_line_df) +
  aes(x=deleted_bases, y=abs(line_pos1-line_pos2)) +
  geom_jitter(width = 0.5, height = 0.5) +
  coord_cartesian(xlim=c(-20, 11), ylim=c(0,2000)) +
  geom_density2d() +
  geom_marginal(aes(group=deleted_bases, color="red")) +
  labs(title="vs length of LINE insertion")

ggplot(bpbp_line_df %>% filter(deleted_bases >= -20 & deleted_bases < 12)) +
  aes(x=line_length, fill=as.factor(deleted_bases)) +
  geom_histogram(bins=25) +
  scale_x_continuous(limits=c(0,1000)) +
  labs(title="Length of LINE insertions")
ggplot(bpbp_line_df %>% filter(deleted_bases >= -20 & deleted_bases < 12)) +
  aes(x=line_length) +
  facet_wrap(~ deleted_bases, scales="free") +
  geom_histogram(bins=25) +
  scale_x_continuous(limits=c(0,1000)) +
  labs(title="Length of LINE insertion compared to insertion site overlap")
ggplot(bpbp_line_df %>% filter(deleted_bases >= -20 & deleted_bases < 12) %>%
         filter(line_length <= 1000) %>%
         mutate(bin=floor(line_length / (1000 / 20)) * (1000/20)) %>%
         group_by(bin, deleted_bases) %>%
         summarise(n=n()) %>%
         group_by(bin) %>%
         mutate(portion=n/sum(n))) +
  aes(x=bin, y= portion, fill=as.factor(deleted_bases)) +
  geom_bar(stat="identity") +
  labs(title="vs LINE insertion length")

cbind_hitdf = function(hitdf, subject_gr, query_gr, suffix=c("1", "2")) {
  todf = function(gr, suffix="") {
    names(gr) = NULL
    x = as.data.frame(gr)
    names(x) = paste0(names(x), suffix)
    return(x)
  }
  return(bind_cols(hitdf, todf(subject_gr[hitdf$subjectHits], suffix[1]), todf(query_gr[hitdf$queryHits], suffix[2])))
}
bpbp_df = cbind_hitdf(bpbp_line_df, bp_line_gr, bp_line_gr) %>%
  replace_na(list(insSeq1="", Homology1="", insSeq2="", Homology2="")) %>%
  mutate(
    ins_length = str_length(insSeq1) + str_length(insSeq2),
    homlen=str_length(Homology1) + str_length(Homology2),
    hasHom=Homology1 != "" | Homology2 != "",
    hasNoIns = ins_length == 0,
    hasOneIns = pmax(str_length(insSeq1), str_length(insSeq2)) == ins_length,
    hasTwoIns = !hasNoIns & !hasOneIns,
    isPolyA = str_detect(refContext1, "AAAA") | str_detect(refContext1, "TTTT") | str_detect(refContext2, "AAAA") | str_detect(refContext2, "TTTT"),
    ihomlen=ihomlen1 + ihomlen2,
    cn=(cn1 + cn2) / 2,
    ins1As = str_count(insSeq1, stringr::fixed("A")),
    ins2As = str_count(insSeq2, stringr::fixed("A")),
    ins1ATs = str_count(insSeq1, stringr::fixed("A")) + str_count(insSeq1, stringr::fixed("T")),
    ins2ATs = str_count(insSeq2, stringr::fixed("A")) + str_count(insSeq2, stringr::fixed("T")))
clusterdf = bpbp_df %>%
  dplyr::select(
    deleted_bases,
    line_length,
    ins_length,
    hasHom,
    hasNoIns,
    hasOneIns,
    hasTwoIns,
    isPolyA,
    ihomlen,
    homlen,
    cn)
clustermat=scale(as.matrix(clusterdf %>% dplyr::select(-deleted_bases)))
require(Rtsne)
tsnedf = Rtsne(clustermat)
ggplot(data.frame(
    x=tsnedf$Y[,1],
    y=tsnedf$Y[,2],
    deleted_bases=clusterdf$deleted_bases)) +
  aes(x=x, y=y, col=cut(deleted_bases, c(-20, -16, -11, -7, 0, 5))) +
  geom_point()

require(rpart)
fit = rpart(isStrandInvasion ~ line_length + ins_length + hasHom + hasNoIns + hasOneIns + hasTwoIns + isPolyA + ihomlen + cn + homlen,
            clusterdf %>% mutate(isStrandInvasion=ifelse(deleted_bases < -7, "StrandInvasion", "Clean")),
            method="class")
rsq.rpart(fit)

ggplot(clusterdf) +
  aes(x=deleted_bases, fill=ifelse(hasNoIns, "Both Clean", ifelse(hasOneIns, "One Clean", "Both breaks have inserted sequence"))) +
  geom_histogram(bins=45) +
  scale_x_continuous(limits=c(-25, 20))


bpbp_df %>% filter(hasTwoIns) %>% dplyr::select(insSeq1, refContext1, insSeq2, refContext2) %>% View()

ggplot(bpbp_df) +
  aes(x=str_length(insSeq1), y=str_length(insSeq2), color=(ins1As + ins2As)/ins_length > 0.8) +
  geom_point()
# non-polyA insert sequence?
ggplot(bpbp_df) +
  aes(x=deleted_bases, y=(ins1ATs + ins2ATs)/ins_length, color=hasTwoIns, size=ins_length) +
  geom_jitter() +
  scale_x_continuous(limits=c(-25, 20))




