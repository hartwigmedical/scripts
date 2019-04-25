library(StructuralVariantAnnotation)
library(VariantAnnotation)
library(rtracklayer)
library(tidyverse)
library(readr)
library(stringr)
library(cowplot)

rmfile = "D:/hartwig/hg19.fa.out"

vcf_to_breakpoints = function(file, sampleId=str_match(file, "([^/\\\\]*).gridss.vcf.gz.ti.vcf")[,2]) {
  vcf = readVcf(file)
  gr = breakpointRanges(vcf, nominalPosition=TRUE)
  gr = gr[gr$partner %in% names(gr)]
  gr$sampleId = str_match(file, "([^/\\\\]*).gridss.vcf.gz.ti.vcf")[,2]
  gr$asmlinks=info(vcf[names(gr)])$BEID
  return(gr)
}

simple_ti_candidates = function(gr, insert_site_max_gap=35, template_insert_max_length=10000) {
  gr$repeatClass = grrm$repeatClass[ findOverlaps(gr, grrm, select="first", ignore.strand=TRUE)]
  gr$repeatType = grrm$repeatType[ findOverlaps(gr, grrm, select="first", ignore.strand=TRUE)]
  pgr = partner(gr)
  # matches logic in find_germline_TIs.Rs
  hits = findOverlaps(gr, gr, maxgap=insert_site_max_gap, ignore.strand=TRUE) %>%
    as.data.frame() %>%
    filter(as.logical(
      strand(gr)[subjectHits] != strand(gr)[queryHits] &
      gr$sampleId[subjectHits] == gr$sampleId[queryHits] &
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
      asmlinked=elementNROWS(intersect(gr[queryHits]$asmlinks, gr[subjectHits]$asmlinks)) > 0,
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

germline_raw_grl = lapply(list.files(path="../analysisscripts/germlineTIs/candidates/", pattern="*.gridss.vcf.gz.ti.vcf", full.names=TRUE), vcf_to_breakpoints)
fulldf = bind_rows(lapply(germline_raw_grl, simple_ti_candidates))
fullgr = unlist(GRangesList(germline_raw_grl))

source("../analysisscripts/libSvAnalyser.R")
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
svTable = query_somatic_structuralVariants(dbProd)
somatic_raw_gr = to_sv_gr(svTable, FALSE)
somatic_raw_gr = somatic_raw_gr[!is.na(somatic_raw_gr$partner)]
somatic_raw_gr$asmlinks = setdiff(CharacterList(paste(somatic_raw_gr$linkedBy, partner(somatic_raw_gr)$linkedBy, sep=",")), ".")


somatic_list_gr = lapply(unique(somatic_raw_gr$sampleId))

fullgr = unlist(GRangesList(fullgr))
names(fullgr) = paste0(fullgr$sampleId, names(fullgr))
fullgr$partner = paste0(fullgr$sampleId, fullgr$partner)
fulldf$uid1 = paste0(fulldf$sampleId, fulldf$vcfid1)
fulldf$uid2 = paste0(fulldf$sampleId, fulldf$vcfid2)

# can't use findBreakpointOverlaps() because there's too many pairwise hits
pfullgr = partner(fullgr)
merged = as.data.frame(fullgr) %>%
  dplyr::mutate(
    seqnames2=as.character(seqnames(pfullgr)),
    start2=start(pfullgr),
    strand2=as.character(strand(pfullgr))) %>%
  group_by(seqnames, seqnames2, start, start2, strand, strand2) %>%
  mutate(uid=paste0(sampleId, vcfId)) %>%
  summarise(
    vcfIds=paste(vcfId, collapse=";"),
    sampleIds=paste(sampleId, collapse=";"),
    uids=paste(uid, collapse=";"),
    n=n()) %>%
  ungroup() %>%
  mutate(id=paste0("merged", row_number()))
mlookup = dplyr::select(merged, id, uids) %>%
  separate_rows(uids, sep=";")
fulldf = fulldf %>%
  left_join(mlookup, by=c("uid1"="uids")) %>%
  left_join(mlookup, by=c("uid2"="uids"), suffix=c("1", "2"))

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

dedupdf = fulldf %>%
  group_by(id1, id2, repins1, repins2, repti1, repti2, repins, repti) %>%
  summarise(
    insdellen=mean(insdellen),
    tilen=mean(tilen),
    hom1=mean(hom1),
    hom2=mean(hom2),
    count=n()) %>%
  ungroup()

ggplot(dedupdf) +
  aes(x=insdellen, fill=factor(pmin(count, 10), levels=1:10, labels=c(1:9, "10+"))) +
  geom_histogram(bins=70) +
  scale_x_continuous(limits=c(-35, 35)) +
  facet_wrap(. ~ repti) +
  labs(title="Germline templated insertions", fill="Samples", x="Overlap at insertion site(bp)")

ggplot(dedupdf) +
  aes(x=count) +
  geom_histogram() +
  scale_x_log10()

dedupdf_3 = dedupdf %>% filter(count > 3)

plot_scatter_main = ggplot(dedupdf_3) +
  aes(x=insdellen, y=tilen, colour=repti1, size=log10(count)) +
  geom_jitter() +
  scale_y_log10()
plot_scatter_xhist = axis_canvas(plot_scatter_main, axis="x") +
  geom_histogram(data=dedupdf_3, bins=70) +
  aes(x=insdellen, fill=repti1)
plot_scatter_ybox = axis_canvas(plot_scatter_main, axis = "y", coord_flip=TRUE) +
  geom_density(data=dedupdf_3, alpha=0.5) +
  aes(x=tilen, fill=repti1) +
  scale_x_log10() +
  coord_flip()
plot_scatter = plot_scatter_main %>%
  insert_xaxis_grob(plot_scatter_xhist, position = "bottom") %>%
  insert_yaxis_grob(plot_scatter_ybox, position = "right")
ggdraw(plot_scatter)


ggplot(dedupdf) +
  aes(x=tilen) +
  geom_histogram(bins=70) +
  scale_x_log10() +
  facet_wrap(. ~ repti)

