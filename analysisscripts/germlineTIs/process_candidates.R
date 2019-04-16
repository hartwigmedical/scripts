library(StructuralVariantAnnotation)
library(VariantAnnotation)
library(rtracklayer)
library(tidyverse)
library(stringr)
library(cowplot)

rmfile = "D:/hartwig/hg19.fa.out"

candidate_to_df = function(file, insert_site_max_gap=35, template_insert_max_length=10000) {
  vcf = readVcf(file)
  gr = breakpointRanges(vcf, nominalPosition=TRUE)
  gr = gr[gr$partner %in% names(gr)]
  gr$repeatClass = grrm$repeatClass[ findOverlaps(gr, grrm, select="first", ignore.strand=TRUE)]
  gr$repeatType = grrm$repeatType[ findOverlaps(gr, grrm, select="first", ignore.strand=TRUE)]
  pgr = partner(gr)
  # matches logic in find_germline_TIs.Rs
  hits = findOverlaps(gr, gr, maxgap=insert_site_max_gap, ignore.strand=TRUE) %>%
    as.data.frame() %>%
    filter(as.logical(
      strand(gr)[subjectHits] != strand(gr)[queryHits] &
        # partner breakends make a TI
        seqnames(pgr)[subjectHits] == seqnames(pgr)[queryHits] &
        strand(pgr)[subjectHits] == "-" & strand(pgr)[queryHits] == "+") &
        start(pgr)[subjectHits] < start(pgr)[queryHits]) %>%
    mutate(
      insdellen=ifelse(as.logical(strand(gr)[subjectHits] == "-"), 1, -1) * start(gr)[subjectHits] - start(gr)[queryHits],
      tilen = abs(start(pgr)[subjectHits] - start(pgr)[queryHits])) %>%
    filter(abs(insdellen) <= insert_site_max_gap & tilen <= template_insert_max_length) %>%
    # ignore local events
    filter(as.logical(seqnames(gr)[subjectHits] != seqnames(pgr)[subjectHits]) | abs(start(gr)[subjectHits] - start(pgr)[subjectHits]) > template_insert_max_length)

  # repeatmasker annotate
  inslocgr = GRanges(
    seqnames=seqnames(gr)[hits$subjectHits],
    ranges=IRanges(
      start=pmin(start(gr)[hits$subjectHits], start(gr)[hits$queryHits]),
      end=pmax(start(gr)[hits$subjectHits], start(gr)[hits$queryHits])))
  tigr = GRanges(
    seqnames=seqnames(pgr)[hits$subjectHits],
    ranges=IRanges(
      start=pmin(start(pgr)[hits$subjectHits], start(pgr)[hits$queryHits]),
      end=pmax(start(pgr)[hits$subjectHits], start(pgr)[hits$queryHits])))
  # TODO annotate with portion of insertion
  # Get logic from BEALN annotation code

  df = hits %>% mutate(
      sampleId=str_match(file, "([^/\\\\]*).gridss.vcf.gz.ti.vcf")[,2],
      vcfid1=names(gr)[subjectHits],
      vcfid2=names(gr)[queryHits],
      hom1=as.integer(info(vcf[names(gr)[subjectHits]])$HOMLEN),
      hom2=as.integer(info(vcf[names(gr)[queryHits]])$HOMLEN),
      rcins1=gr$repeatClass[subjectHits],
      rcins2=gr$repeatClass[queryHits],
      rcti1=pgr$repeatClass[subjectHits],
      rcti2=pgr$repeatClass[queryHits],
      rtins1=gr$repeatType[subjectHits],
      rtins2=gr$repeatType[queryHits],
      rtti1=pgr$repeatType[subjectHits],
      rtti2=pgr$repeatType[queryHits]) %>%
    replace_na(list(
      "hom1"=0, "hom2"=0,
      "rcins1"="", "rcins2"="", "rcti1"="", "rcti2"="",
      "rtins1"="", "rtins2"="", "rtti1"="", "rtti2"=""))
  return(df)
}
candidate_to_gr = function(file, insert_site_max_gap=35, template_insert_max_length=10000) {
  vcf = readVcf(file)
  gr = breakpointRanges(vcf, nominalPosition=TRUE)
  gr$sampleId=str_match(file, "([^/\\\\]*).gridss.vcf.gz.ti.vcf")[,2]
  return(gr)
}


# from http://github.com/PapenfussLab/sv_benchmark
import.repeatmasker.fa.out <- function(repeatmasker.fa.out) {
  rmdt <- read_table2(repeatmasker.fa.out, col_names=FALSE, skip=3)
  grrm <- GRanges(
    seqnames=rmdt$X5,
    ranges=IRanges(start=rmdt$X6 + 1, end=rmdt$X7),
    strand=ifelse(rmdt$X9=="C", "-", "+"),
    repeatType=rmdt$X10,
    repeatClass=rmdt$X11)
  return(grrm)
}
# wget http://www.repeatmasker.org/genomes/hg19/RepeatMasker-rm330-db20120124/hg19.fa.out.gz
cache_filename = paste0(rmfile, ".grrm.rds")
if (file.exists(cache_filename)) {
  grrm = readRDS(cache_filename)
} else {
  grrm = import.repeatmasker.fa.out(rmfile)
  saveRDS(grrm, file=cache_filename)
}
seqlevelsStyle(grrm) = "NCBI"

fulldf = bind_rows(lapply(list.files(path="../analysisscripts/germlineTIs/candidates/", pattern="*.gridss.vcf.gz.ti.vcf", full.names=TRUE), candidate_to_df))
fullgr = lapply(list.files(path="../analysisscripts/germlineTIs/candidates/", pattern="*.gridss.vcf.gz.ti.vcf", full.names=TRUE), candidate_to_gr)
fullgr = unlist(GRangesList(fullgr))
names(fullgr) = paste0(fullgr$sampleId, names(fullgr))
fullgr$partner = paste0(fullgr$sampleId, fullgr$partner)

# can't use findBreakpointOverlaps() because there's too many pairwise hits
pfullgr = partner(fullgr)
merged = as.data.frame(fullgr) %>%
  dplyr::mutate(
    seqnames2=as.character(seqnames(pfullgr)),
    start2=start(pfullgr),
    strand2=as.character(strand(pfullgr))) %>%
  group_by(seqnames, seqnames2, start, start2, strand, strand2) %>%
  summarise(
    vcfIds=paste(vcfId, collapse=";"),
    sampleIds=paste(sampleId, collapse=";"),
    uids=paste(paste0(fullgr$sampleId, fullgr$vcfId), collapse=";"),
    n=n(),
    id=paste0("merged:", seq_along(.)))
merged %>% separate_rows(uid, ";") %>%
  inner_join(fulldf, by=c("uid"=


fulldf = fulldf %>%
  mutate(
    repins1=str_extract(rcins1, "^[^/]*"),
    repins2=str_extract(rcins2, "^[^/]*"),
    repti1=str_extract(rcti1, "^[^/]*"),
    repti2=str_extract(rcti2, "^[^/]*"),
    repins=ifelse(repins1 == "", repins2, repins1),
    repti=ifelse(repti1 == "", repti2, repti1))

# TODO:
# 0) Generated merged VCF by adding a new annotation
# 1) add TI top level repeat class and %coverage
# 2) dedup variants common across multiple samples


ggplot(fulldf) +
  aes(x=insdellen) +
  geom_histogram(bins=70) +
  scale_x_continuous(limits=c(-35, 35))
  facet_wrap(. ~ repti)

plot_scatter_main = ggplot(fulldf) +
  aes(x=insdellen, y=tilen, colour=repti) +
  geom_point() +
  scale_y_log10()
plot_scatter_xhist = axis_canvas(plot_scatter_main, axis="x") +
  geom_histogram(data=fulldf, bins=70) +
  aes(x=insdellen, fill=repti)
plot_scatter_ybox = axis_canvas(plot_scatter_main, axis = "y", coord_flip=TRUE) +
  geom_density(data=fulldf) +
  aes(x=tilen, fill=repti) +
  scale_x_log10() +
  coord_flip()
plot_scatter = plot_scatter_main %>%
  insert_xaxis_grob(plot_scatter_xhist, position = "bottom") %>%
  insert_yaxis_grob(plot_scatter_ybox, position = "right")
ggdraw(plot_scatter)


ggplot(fulldf) +
  aes(x=tilen) +
  geom_histogram(bins=70) +
  scale_x_log10() +
  facet_wrap(. ~ repti)

