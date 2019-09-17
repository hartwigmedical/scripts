library(tidyverse)
library(stringi)
library(stringr)
library(grid)
library(gridExtra)
library(cowplot)
library(scales)
library(GenomicRanges)
library(StructuralVariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
options(stringsAsFactors=FALSE)
source("libSVPaper.R")

basedir="~/../Dropbox (HMF Australia)/HMF Australia team folder/Structural Variant Analysis/colo829/"

# UCSC table export
hg19_gaps = with(read_tsv(paste0(basedir,"hg19_gap")),
                 GRanges(seqnames=str_replace(chrom, "chr", ""), ranges=IRanges(start=chromStart, end=chromEnd), type=type))
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
hg19_cytobands =  with(read_tsv(
  file=paste0(basedir, "cytoband.txt"),
  col_names=c("chr", "start", "end", "band", "type"),
  col_type= "ciicc"),
  GRanges(seqnames=chr, ranges=IRanges(start=start+1, end=end), band=band, type=type))
seqlevelsStyle(hg19_cytobands) = "NCBI"
hg19_centromeres = hg19_cytobands[hg19_cytobands$type == "acen"]

.lp_lookup_df = hg19_cytobands %>% as.data.frame() %>%
  group_by(seqnames) %>%
  summarise(len=max(end)) %>%
  ungroup() %>%
  mutate(chrint=as.integer(ifelse(seqnames == "X", 23, ifelse(seqnames == "Y", 24, as.character(seqnames))))) %>%
  arrange(chrint) %>%
  mutate(offset=cumsum(as.numeric(len))-len) %>%
  mutate(chr=as.character(seqnames))
.lp_lookup = .lp_lookup_df$offset
names(.lp_lookup) = .lp_lookup_df$chr
as_linear_pos = function(chr, x) { .lp_lookup[as.character(chr)] + x }

######
# Weaver load
# CHR BEGIN END ALLELE_1_CN ALLELE_2_CN
weaver_cn = with(read_tsv(
  file=paste0(basedir, "weaver/REGION_CN_PHASE"),
  col_names=c("chr", "start", "end", "cn1", "cn2"),
  col_type= "ciinn"),
  GRanges(seqnames=chr, ranges=IRanges(start=start, end=end), cn1=cn1, cn2=cn2, caller="weaver"))
# CHR_1 POS_1 ORI_1 ALLELE_ CHR_2 POS_2 ORI_2 ALLELE_ CN germline/somatic_post_aneuploidy/somatic_pre_aneuploidy
weaver_bp = with(read_tsv(
  file=paste0(basedir, "weaver/SV_CN_PHASE"),
  col_names=c("chr1", "start1", "ori1", "allele1", "chr2", "start2", "ori2", "allele2", "cn1", "cn2", "cncn", "type"),
  col_type= "ciciciciiicc"), {
    bp_name = paste0("bp", seq_along(chr1))
    gro = GRanges(seqnames=chr1, ranges=IRanges(start=start1, width=1), strand=ori1, allele=allele1, cn=cn1, type=type, partner=paste0(bp_name, "h"), caller="weaver")
    grh = GRanges(seqnames=chr2, ranges=IRanges(start=start2, width=1), strand=ori2, allele=allele2, cn=cn2, type=type, partner=paste0(bp_name, "o"), caller="weaver")
    names(gro) = paste0(bp_name, "o")
    names(grh) = paste0(bp_name, "h")
    return(c(gro, grh))
  })

######
# GRIDSS/purple load
gridss_vcf = readVcf(paste0(basedir, "purple/COLO829T.purple.sv.ann.vcf"))
gridss_bp_gr = breakpointRanges(gridss_vcf, nominalPosition=TRUE)
gridss_be_gr = breakendRanges(gridss_vcf, nominalPosition=TRUE)
gridss_gr = c(gridss_bp_gr, gridss_be_gr)
gridss_gr = gridss_gr[!str_detect(names(gridss_gr), "purple")] # strip out placeholder purple breakends
gridss_gr$vaf =
gridss_gr$caller="purple"
purple_cn = with(read_tsv(paste0(basedir, "purple/COLO829T.purple.cnv")) %>% rename("#chromosome"="chromosome"),
  GRanges(seqnames=chromosome, ranges=IRanges(start=start, end=end),
          cn=copyNumber,
          bafCount=bafCount,
          observedBAF=observedBAF,
          actualBAF=actualBAF,
          segmentStartSupport=segmentStartSupport,
          segmentEndSupport=segmentEndSupport,
          method=method,
          depthWindowCount=depthWindowCount,
          gcContent=gcContent,
          minStart=minStart,
          maxStart=maxStart,
          caller="purple"))
purple_germline_cn = with(read_tsv(paste0(basedir, "purple/COLO829T.purple.germline.cnv")) %>% rename("#chromosome"="chromosome"),
  GRanges(seqnames=chromosome, ranges=IRanges(start=start, end=end),
    cn=copyNumber,
    bafCount=bafCount,
    observedBAF=observedBAF,
    actualBAF=actualBAF,
    segmentStartSupport=segmentStartSupport,
    segmentEndSupport=segmentEndSupport,
    method=method,
    depthWindowCount=depthWindowCount,
    gcContent=gcContent,
    minStart=minStart,
    maxStart=maxStart,
    caller="purple"))

######
# Conserting/CREST load
conserting_qual_merge = read_tsv(paste0(basedir, "conserting/colo829_CONSERTING_Mapability_100.txt.QualityMerge"))
conserting_cna_calls = read_tsv(paste0(basedir, "conserting/colo829_CONSERTING_Mapability_100.txt.CNAcalls"))
conserting_conflict = read_tsv(paste0(basedir, "conserting/colo829_CONSERTING_Mapability_100_potential_conflict_segment.txt"))
conserting_crest = read_tsv(paste0(basedir, "conserting/colo829_CREST-map_final_report.txt.QualityMerge"),
  col_names=c("chr1", "pos1", "ori1", "score1", "chr2", "pos2", "ori2", "score2"),
  col_type= "cicicici")
conserting_chr_to_xy = function(chr) { ifelse(chr == 23, "X", ifelse(chr == 24, "Y", chr)) }
conserting_cna = with(conserting_cna_calls, GRanges(
  seqnames=conserting_chr_to_xy(chrom),
  ranges=IRanges(start=loc.start, end=loc.end),
  num.mark=num.mark,
  Log2Ratio=Log2Ratio,
  caller="conserting"))
conserting_cn = with(conserting_qual_merge, GRanges(
  seqnames=conserting_chr_to_xy(chrom),
  ranges=IRanges(start=loc.start, end=loc.end),
  num.mark=num.mark,
  length.ratio=length.ratio,
  seg.mean=seg.mean,
  GMean=GMean,
  DMean=DMean,
  LogRatio=LogRatio,
  QualityScore=QualityScore,
  SV_Matched=SV_Matched,
  caller="conserting"))
crest_bp = with(conserting_crest, {
    bp_name = paste0("bp", seq_along(chr1))
    gro = GRanges(seqnames=conserting_chr_to_xy(chr1), ranges=IRanges(start=pos1, width=1), strand=ori1, score=score1, partner=paste0(bp_name, "h"), caller="conserting")
    grh = GRanges(seqnames=conserting_chr_to_xy(chr2), ranges=IRanges(start=pos2, width=1), strand=ifelse(ori2 == "+", "-", "+"), score=score2, partner=paste0(bp_name, "o"), caller="conserting")
    names(gro) = paste0(bp_name, "o")
    names(grh) = paste0(bp_name, "h")
    return(c(gro, grh))
  })
######
# Ascat
ascat_cn = with(read_tsv(paste0(basedir, "ascat/COLO829T.segments.txt")),
  GRanges(seqnames=chr, ranges=IRanges(start=startpos, end=endpos), nMajor=nMajor, nMinor=nMinor, caller="ascat"))


ascat_cn$cn = ascat_cn$nMajor + ascat_cn$nMinor
ascat_cn$cn_major = ascat_cn$nMajor
ascat_cn$cn_minor = ascat_cn$nMinor
weaver_cn$cn = weaver_cn$cn1 + weaver_cn$cn2
weaver_cn$cn_major = pmax(weaver_cn$cn1, weaver_cn$cn2)
weaver_cn$cn_minor = pmin(weaver_cn$cn1, weaver_cn$cn2)
purple_cn$cn_major = pmax(purple_cn$cn * purple_cn$actualBAF, purple_cn$cn * (1 - purple_cn$actualBAF))
purple_cn$cn_minor = pmin(purple_cn$cn * purple_cn$actualBAF, purple_cn$cn * (1 - purple_cn$actualBAF))
# use the purple average ploidy so result scales match
conserting_dg_normalisation =
  (sum(width(conserting_cn) * conserting_cn$GMean)/sum(width(conserting_cn))) /
  (sum(width(conserting_cn) * conserting_cn$DMean)/sum(width(conserting_cn)))
purple_cn_average = sum(width(purple_cn) * purple_cn$cn)/sum(width(purple_cn))
ascat_cn_average = sum(width(ascat_cn) * ascat_cn$cn)/sum(width(ascat_cn))
conserting_cn$cn = ascat_cn_average * conserting_dg_normalisation *
  (conserting_cn$DMean/conserting_cn$GMean)
cn = c(ascat_cn, purple_cn, conserting_cn, weaver_cn)
# Consistency with SV calls:
# segments
# segments supported by SV
# distance to SV
evaluate_cn_transitions = function (cngr, svgr, margin=100000, distance=c("cn_transition", "sv")) {
  distance <- match.arg(distance)
  cn_transitions = with(cngr %>% as.data.frame(), IRanges::reduce(c(
    GRanges(seqnames=seqnames, ranges=IRanges(start=start, width=1)),
    GRanges(seqnames=seqnames, ranges=IRanges(start=end + 1, width=1)))))
  cn_transitions$caller = unique(cngr$caller)
  cn_transitions$cn_left = NA
  cn_transitions$cn_right = NA
  hits = findOverlaps(cn_transitions, with(cngr %>% as.data.frame(), GRanges(seqnames=seqnames, ranges=IRanges(start=start, width=1))))
  cn_transitions$cn_right[queryHits(hits)] = cngr[subjectHits(hits)]$cn
  hits = findOverlaps(cn_transitions, with(cngr %>% as.data.frame(), GRanges(seqnames=seqnames, ranges=IRanges(start=end + 1, width=1))))
  cn_transitions$cn_left[queryHits(hits)] = cngr[subjectHits(hits)]$cn
  cn_transitions$cn_delta = cn_transitions$cn_right - cn_transitions$cn_left

  hits = findOverlaps(svgr, cn_transitions, maxgap=margin, ignore.strand=TRUE) %>% as.data.frame() %>%
    mutate(distance=pmax(1, abs(start(svgr)[queryHits] - start(cn_transitions)[subjectHits])))
  best_sv_hit = hits %>% group_by(queryHits) %>% filter(distance==min(distance)) %>% ungroup()
  best_cn_hit = hits %>% group_by(subjectHits) %>% filter(distance==min(distance)) %>% ungroup()
  cn_transitions$distance = NA
  cn_transitions[best_cn_hit$subjectHits]$distance = best_cn_hit$distance
  svgr$distance = NA
  svgr[best_sv_hit$queryHits]$distance = best_sv_hit$distance
  exact_bp_hit = best_sv_hit %>%
    mutate(queryHits.p=match(svgr[queryHits]$partner, names(svgr))) %>%
    filter(!is.na(queryHits.p)) %>%
    inner_join(best_sv_hit, by=c("queryHits.p"="queryHits"), suffix=c("", ".p")) %>%
    filter(distance == 1 & distance.p == 1) %>%
    mutate(
      ori_local = as.character(strand(svgr))[queryHits],
      ori_remote = as.character(strand(svgr))[queryHits.p],
      sv_cn = cn_transitions$cn_delta[subjectHits.p] * ifelse(ori_remote == "+", -1, 1),
      cn_left = cn_transitions$cn_left[subjectHits] + ifelse(ori_local == "+", 0, sv_cn),
      cn_right = cn_transitions$cn_right[subjectHits] + ifelse(ori_local == "+", sv_cn, 0),
      cn_error=abs(cn_left - cn_right))
  sv_hit_count = best_sv_hit %>% group_by(subjectHits) %>% summarise(n=n())
  cn_transitions$sv_matches = 0
  cn_transitions$sv_matches[sv_hit_count$subjectHits] = sv_hit_count$n
  cn_transitions$sv_partner_matches = 0
  cn_transitions$sv_partner_matches[exact_bp_hit$subjectHits] = cn_transitions$sv_matches[exact_bp_hit$subjectHits.p]
  cn_transitions$can_evaluate_cn_error = cn_transitions$sv_matches == 1 & cn_transitions$sv_partner_matches == 1
  cn_transitions$cn_error = NA
  cn_transitions$cn_error[exact_bp_hit$subjectHits] = exact_bp_hit$cn_error
  cn_transitions$cn_error[!cn_transitions$can_evaluate_cn_error] = NA

  if (distance == "cn_transition") {
    return(cn_transitions)
  } else {
    return(svgr)
  }
}
cn_transistions = c(
  evaluate_cn_transitions(ascat_cn, gridss_gr),
  evaluate_cn_transitions(purple_cn, gridss_gr),
  evaluate_cn_transitions(conserting_cn, crest_bp),
  evaluate_cn_transitions(weaver_cn, weaver_bp)
)

cn_transistions$inGap = overlapsAny(cn_transistions, hg19_gaps, maxgap=100000)
cn_transistions$inCentromere = overlapsAny(cn_transistions, hg19_centromeres, maxgap=100000)
cn_transistions$isFirstOrLast = is.na(lead(as.character(seqnames(cn_transistions)))) |
  is.na(lag(as.character(seqnames(cn_transistions)))) |
  seqnames(cn_transistions) != lead(as.character(seqnames(cn_transistions))) |
  seqnames(cn_transistions) != lag(as.character(seqnames(cn_transistions)))


########
# Ploidy estimates
sum(width(ascat_cn) * ascat_cn$cn)/sum(width(ascat_cn))
sum(width(purple_cn) * purple_cn$cn)/sum(width(purple_cn))
sum(width(weaver_cn) * weaver_cn$cn)/sum(width(weaver_cn))
sum(width(conserting_cn) * conserting_cn$cn)/sum(width(conserting_cn))

########
# Standardised CN export
ascat_cn$score = ascat_cn$cn
purple_cn$score = purple_cn$cn
conserting_cn$score = conserting_cn$cn
weaver_cn$score = weaver_cn$cn
export(ascat_cn, paste0(basedir, "/out/ascat_cn.bed"))
export(purple_cn, paste0(basedir, "/out/purple_cn.bed"))
export(conserting_cn, paste0(basedir, "/out/conserting_cn.bed"))
export(weaver_cn, paste0(basedir, "/out/weaver_cn.bed"))
library('openxlsx')
suppfigurewb <- createWorkbook()
addWorksheet(suppfigurewb, "cn")
writeData(suppfigurewb, "cn", c(ascat_cn, purple_cn, weaver_cn, conserting_cn) %>% as.data.frame() %>% dplyr::select(1:11))
addWorksheet(suppfigurewb, "sv")
crest_bp$beid = names(crest_bp)
weaver_bp$beid = names(weaver_bp)
gridss_gr$beid = names(gridss_gr)
writeData(suppfigurewb, "sv", c(gridss_gr, crest_bp, weaver_bp) %>% as.data.frame(row.names=NULL) %>% dplyr::select(-c("vaf","score","allele","cn","type")))
addWorksheet(suppfigurewb, "cntransitions")
writeData(suppfigurewb, "cntransitions", cn_transistions %>% as.data.frame())
saveWorkbook(suppfigurewb, file = paste0(basedir, "../figures/supptable_purple_conserting_weaver_cn_comparison.xlsx"), overwrite = TRUE)

purple_first = function(x) {factor(x, levels=c("purple", "conserting", "weaver", "ascat")) }
cn_transistions$caller = purple_first(cn_transistions$caller)
########
# Plots
cn_size_distribution = c(ascat_cn, purple_cn, conserting_cn, weaver_cn) %>%
  as.data.frame() %>%
  dplyr::select(seqnames, start, end, caller) %>%
  mutate(length=end-start)
ggplot(cn_size_distribution) +
  aes(x=length) +
  geom_histogram() +
  scale_x_log10() +
  facet_grid(caller ~ ., scales="free_y") +
  labs(title="Copy number segment size distribution")

cn_transition_plot = cn_transistions %>%
  as.data.frame() %>%
  replace_na(list(distance=100000)) %>%
  filter(!(inGap | isFirstOrLast | inCentromere))
ggplot(cn_transition_plot) +
  aes(x=distance) +
  geom_histogram(bins=10) +
  scale_x_log10(breaks=c(1,10,100,1000, 10000, 100000), labels=c("0", "10", "100", "1000", "10000", "No matching SV")) +
  facet_grid(caller ~ ., scales="free_y") +
  labs(y="CN transitions", x="Distance to nearest SV")

fig3_sv_at_cn_transition = cn_transition_plot %>%
  mutate(distance_bin = cut(
    distance,
    labels=c("Less than 10 bp", "10-99 bp", "100-999 bp", "1000-9999 bp", "10000-99999 bp", "No matching SV"),
    breaks=c(-100, 1,10,100,1000, 10000, 100000))) %>%
  group_by(caller) %>%
  mutate(caller_n=n()) %>%
  group_by(caller, distance_bin) %>%
  summarise(portion=n()/max(caller_n)) %>%
  ungroup() %>%
ggplot() +
  aes(x=caller, y=portion, fill=distance_bin) +
  geom_col() +
  scale_fill_manual(values=c("#ABD9E9", "#B4B3C0", "#BD8C97","#C5666E", "#CE3F45", "#D7191C"), name="Distance to SV", drop=FALSE) +
  scale_y_continuous(labels=percent) +
  labs(
    x="",
    y="Percent of copy number transitions") +
  theme(axis.text.x = element_text(angle = 90))
fig3_sv_at_cn_transition
ggsave(filename=paste0(basedir, "/sv_at_cn_transition.pdf"), width=6, height=4)


cn_transistions %>%
  as.data.frame() %>%
  filter(!is.na(cn_error)) %>%
  mutate(cn_error=pmin(2, cn_error)) %>%
ggplot() +
  aes(x=cn_error) +
  geom_histogram(bins=20) +
  facet_grid(caller ~ ., scales="free_y") +
  scale_x_continuous(breaks=seq(0, 2, 0.25), labels=c(head(seq(0, 2, 0.25), -1), "2+")) +
  labs(
    title="Copy number consistency of breakpoint-connected segments",
    x="Magnitude of copy number inconsistency",
    y="Count of copy number segmentation transitions")

fig3_cn_consistency = cn_transistions %>%
  as.data.frame() %>%
  filter(!is.na(cn_error)) %>%
  mutate(cn_error=pmin(2, cn_error)) %>%
  mutate(cn_error_bin=cut(
    cn_error,
    labels=c("0.0 - 0.05","0.05 - 0.1", "0.1 - 0.2", "0.2 - 0.5", "1.0 - 2.0", "2.0+"),
    breaks=c(-100, 0.05, 0.1, 0.2, 0.5, 1, 2))) %>%
  group_by(caller) %>%
  mutate(caller_n=n()) %>%
  group_by(caller, cn_error_bin) %>%
  summarise(portion=n()/max(caller_n)) %>%
ggplot() +
  aes(x=caller, y=portion, fill=cn_error_bin) +
  geom_col() +
  scale_fill_manual(values=c("#ABD9E9", "#B4B3C0", "#BD8C97","#C5666E", "#CE3F45", "#D7191C"), name="Magnitude of CN\ninconsistency") +
  scale_y_continuous(labels=percent) +
  labs(
    x="",
    y="Portion of breakpoints") +
  theme(axis.text.x = element_text(angle = 90))
fig3_cn_consistency
ggsave(filename=paste0(basedir, "/cn_consistency.pdf"), width=5, height=4)

fig3_sv_at_cn_transition_percentage = cn_transistions %>%
  as.data.frame() %>%
  filter(!is.na(cn_error)) %>%
  mutate(cn_percentage_error=cn_error/pmax(cn_left, cn_right)) %>%
  mutate(cn_percentage_error_bin=cut(
    cn_percentage_error,
    labels=c("0-1","1-5", "5-10", "10-20", "20-50", "50+"),
    breaks=c(-1, 0.01, 0.05, 0.1, 0.2, 0.5, 1000))) %>%
  group_by(caller) %>%
  mutate(caller_n=n()) %>%
  group_by(caller, cn_percentage_error_bin) %>%
  summarise(portion=n()/max(caller_n)) %>%
ggplot() +
  aes(x=caller, y=portion, fill=cn_percentage_error_bin) +
  geom_col() +
  scale_fill_manual(values=c("#ABD9E9", "#B4B3C0", "#BD8C97","#C5666E", "#CE3F45", "#D7191C"), name="Magnitude of CN\ninconsistency (%)") +
  scale_y_continuous(labels=percent) +
  scale_x_discrete(drop=FALSE) +
  labs(
    title="",
    y="Portion of breakpoints") +
  theme(axis.text.x = element_text(angle = 90))
fig3_sv_at_cn_transition_percentage
ggsave(filename=paste0(basedir, "/cn_consistency_percentage.pdf"), width=5, height=4)

fig3_comparison = plot_grid(
  fig3_sv_at_cn_transition + theme(legend.position="none"),
  get_legend(fig3_sv_at_cn_transition),
  fig3_sv_at_cn_transition_percentage + theme(legend.position="none"),
  get_legend(fig3_sv_at_cn_transition_percentage),
  ncol=2, nrow=2, rel_widths=c(3.2,2))
figsave("../colo829/purple_conserting_weaver_comparisons", fig3_comparison, width=4, height=8)


cn_transistions %>%
  as.data.frame() %>%
  filter(!is.na(cn_error)) %>%
  mutate(cn_error=pmin(2, cn_error)) %>%
ggplot() +
  aes(x=cn_error, y=(cn_left+cn_right)/2) +
  geom_point() +
  #geom_jitter(width=0.02) +
  facet_grid(caller ~ ., scales="free_y") +
  scale_x_continuous(breaks=seq(0, 2, 0.5), labels=c(head(seq(0, 2, 0.5), -1), "2+")) +
  labs(title="Copy number consistency of breakpoint-connected segments", x="Magnitude of copy number inconsistency", y="Average copy number")

with(cn %>% as.data.frame() %>% mutate(chr=as.character(seqnames)), bind_rows(
    data.frame(caller=caller, chr=chr, pos=start, cn=cn),
    data.frame(caller=caller, chr=chr, pos=end, cn=cn))) %>%
  mutate(linearpos=as_linear_pos(chr, pos)) %>%
  arrange(caller, linearpos) %>%
ggplot() +
  geom_segment(aes(x=x, xend=xend, y=cn, yend=cn, colour=caller), data=cn %>%
    as.data.frame() %>%
    mutate(
      x=as_linear_pos(seqnames, start),
      xend=as_linear_pos(seqnames, end))) +
  geom_segment(aes(x=x, xend=xend, y=cn, yend=cn_next, colour=caller), data=cn %>%
                 as.data.frame() %>%
                 mutate(
                   x=as_linear_pos(seqnames, end) + 1,
                   xend=lead(as_linear_pos(seqnames, start)),
                   cn_next=lead(cn)) %>%
                 filter(x == xend & seqnames==lead(seqnames))) +
  geom_vline(aes(xintercept=pos), data=data.frame(pos=as_linear_pos(c(1:22, "X", "Y"), 0))) +
  geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=10), fill="gray", data=data.frame(
    start=as_linear_pos(seqnames(hg19_centromeres), start(hg19_centromeres)),
    end=as_linear_pos(seqnames(hg19_centromeres), end(hg19_centromeres)))) +
  coord_cartesian(ylim=c(0, 8)) +
  scale_x_continuous(breaks=with(.lp_lookup_df, offset+len/2), labels=c(1:22, "X", "Y")) +
  theme(axis.ticks.x=element_blank()) +
  facet_grid(caller ~ . ) +
  labs(x="", y="Copy Number")







