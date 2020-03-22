library(tidyverse)
library(stringi)
library(stringr)
library(grid)
library(gridExtra)
library(cowplot)
library(GenomicRanges)
library(StructuralVariantAnnotation)
options(stringsAsFactors=FALSE)

basedir="~/../Dropbox (HMF Australia)/HMF Australia team folder/Structural Variant Analysis/"

svdf = read_csv(file=paste0(basedir, "/RNA/SVA_SVS.csv"))

bpdf = svdf %>% filter(PosEnd != -1)
# TODO: include single breakends in this GRanges as they should be included when determining the closest event
sva_svs_to_bpgr <- function(svdf) {
  grs = GRanges(
    seqnames=svdf$ChrStart,
    ranges=IRanges(svdf$PosStart, width=1),
    strand=ifelse(svdf$OrientStart == 1, "+", "-"),
    partner=paste0(svdf$Id, "h"),
    beid=paste0(svdf$Id, "o"),
    Id=svdf$Id,
    SampleId=svdf$SampleId)
  names(grs)=grs$beid
  gre = GRanges(
    seqnames=svdf$ChrEnd,
    ranges=IRanges(svdf$PosEnd, width=1),
    strand=ifelse(svdf$OrientStart == 1, "+", "-"),
    partner=paste0(svdf$Id, "o"),
    beid=paste0(svdf$Id, "h"),
    Id=svdf$Id,
    SampleId=svdf$SampleId)
  names(gre)=gre$beid
  return(c(grs, gre))
}
bpgr = sva_svs_to_bpgr(bpdf)
hits = findOverlaps(bpgr, bpgr, maxgap=10000) %>%
  as.data.frame() %>%
  filter(queryHits != subjectHits) %>%
  filter(bpgr[queryHits]$SampleId == bpgr[subjectHits]$SampleId) %>%
  mutate(distance=abs(start(bpgr[queryHits]) - start(bpgr[subjectHits]))) %>%
  mutate(Id=bpgr[queryHits]$Id) %>%
  group_by(Id) %>%
  summarise(distance=min(distance))

bpdf = bpdf %>% left_join(hits %>% dplyr::select(Id, distance), by="Id")


errdf = bpdf %>%
  mutate(distanceBin=cut(distance, c(-1000, 1, 2, 10, 50, 100, 500, 1000, 1000000000))) %>%
  mutate(cnError=abs(abs(CNChgStart)-abs(CNChgEnd))) %>%
  mutate(svError=abs(abs(CNChgStart+CNChgEnd) / 2 - Ploidy))

ggplot(errdf) +
  aes(x=cnError, color=distanceBin) +
  geom_density() +
  scale_x_continuous(expand = c(0,0), limits=c(0, 2)) +
  scale_y_continuous(expand = c(0,0))

