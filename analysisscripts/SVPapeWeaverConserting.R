library(IRanges)
library(tidyverse)
library(reshape2)
library(stringi)
library(stringr)
library(grid)
library(gridExtra)
library(cowplot)
library(GenomicRanges)
library(StructuralVariantAnnotation)
options(stringsAsFactors=FALSE)

basedir="D:/hartwig/colo829/"


# CHR BEGIN END ALLELE_1_CN ALLELE_2_CN
weaver_cn = with(read_tsv(
  file=paste0(basedir, "weaver/REGION_CN_PHASE"),
  col_names=c("chr", "start", "end", "cn1", "cn2"),
  col_type= "ciinn"),
  GRanges(seqnames=chr, ranges=IRanges(start=start, end=end), cn1=cn1, cn2=cn2))
# CHR_1 POS_1 ORI_1 ALLELE_ CHR_2 POS_2 ORI_2 ALLELE_ CN germline/somatic_post_aneuploidy/somatic_pre_aneuploidy
weaver_bp = with(read_tsv(
  file=paste0(basedir, "weaver/SV_CN_PHASE"),
  col_names=c("chr1", "start1", "ori1", "allele1", "chr2", "start2", "ori2", "allele2", "cn1", "cn2", "cncn", "type"),
  col_type= "ciciciciiicc"), {
    bp_name = paste0("bp", seq_along(chr1))
    gro = GRanges(seqnames="chr1", ranges=IRanges(start=start1, width=1), strand=ori1, allele=allele1, cn=cn1, type=type, partner=paste0(bp_name, "h"))
    grh = GRanges(seqnames="chr2", ranges=IRanges(start=start2, width=1), strand=ori2, allele=allele2, cn=cn2, type=type, partner=paste0(bp_name, "o"))
    names(gro) = paste0(bp_name, "o")
    names(grh) = paste0(bp_name, "h")
    return(c(gro, grh))
  })
seqlevelsStyle(weaver_bp) = "NCBI"

gridss_bp_gr = breakpointRanges(readVcf(paste0(basedir, "purple/COLO829T.purple.sv.ann.vcf.gz")))
gridss_be_gr = breakendRanges(readVcf(paste0(basedir, "purple/COLO829T.purple.sv.ann.vcf.gz")))
gridss_gr = c(gridss_bp_gr, gridss_be_gr)
purple_cn = with(read_tsv(paste0(basedir, "purple/COLO829T.purple.cnv")) %>% rename("chromosome"="#chromosome"),
  GRanges(seqnames=chromosome, ranges=IRanges(start=start, end=end),
          cn=copyNumber,
          bafCount=bafCount,
          observedBAF=observedBAF,
          segmentStartSupport=segmentStartSupport,
          segmentEndSupport=segmentEndSupport,
          method=method,
          depthWindowCount=depthWindowCount,
          gcContent=gcContent,
          minStart=minStart,
          maxStart=maxStart))
purple_germline_cn = with(read_tsv(paste0(basedir, "purple/COLO829T.purple.germline.cnv")) %>% rename("chromosome"="#chromosome"),
  GRanges(seqnames=chromosome, ranges=IRanges(start=start, end=end),
    cn=copyNumber,
    bafCount=bafCount,
    observedBAF=observedBAF,
    segmentStartSupport=segmentStartSupport,
    segmentEndSupport=segmentEndSupport,
    method=method,
    depthWindowCount=depthWindowCount,
    gcContent=gcContent,
    minStart=minStart,
    maxStart=maxStart))

ggplot(bind_rows(
  weaver_cn %>% as.data.frame() %>%
    mutate(seqnames=as.character(seqnames), caller="weaver", total_cn=cn1+cn2) %>%
    dplyr::select(seqnames, start, end, total_cn, caller),
  purple_cn %>% as.data.frame() %>%
      mutate(seqnames=as.character(seqnames), caller="purple", total_cn=cn) %>%
    dplyr::select(seqnames, start, end, total_cn, caller),
  purple_germline_cn %>% as.data.frame() %>%
      mutate(seqnames=as.character(seqnames), caller="purple_germline", total_cn=cn) %>%
    dplyr::select(seqnames, start, end, total_cn, caller))) +
  aes(x=start, xend=end, y=total_cn, yend=total_cn, colour=caller) +
  geom_segment() +
  facet_wrap(~ seqnames, scales="free") +
  coord_cartesian(ylim=c(0, 16))


sv_distance_to_segmentation = function(cngr, svgr) {
  #TODO work out why distanceToNearest() function is wrong
  distanceToNearest(
    svgr,
    c(resize(cngr, width=1, fix="start"), resize(cngr, width=1, fix="end")),
    ignore.strand=TRUE)
}
ggplot(bind_rows(
  sv_distance_to_segmentation(weaver_cn, weaver_bp) %>% as.data.frame() %>% mutate(caller="weaver"),
  sv_distance_to_segmentation(purple_cn, gridss_gr) %>% as.data.frame() %>% mutate(caller="purple-gridss"))) +
  aes(x=distance + 1, fill=caller) +
  scale_x_log10() +
  geom_histogram()









