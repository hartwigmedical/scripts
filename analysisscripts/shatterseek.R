library(devtools)
#install_github("parklab/ShatterSeek")
options(stringsAsFactors = FALSE)

acceptable_chr=c("X", 1:22)
# TODO: https://github.com/parklab/ShatterSeek/blob/master/tutorial.pdf
# docs do not match example
sssvdf = with(
  lnx_svs %>% filter(PosEnd >= 1),
  data.frame(
    sampleId=SampleId,
    chrom1=ChrStart,
    Pos1=PosStart,
    chrom2=ChrEnd,
    Pos2=PosEnd,
    strand1=ifelse(OrientStart > 0, "+", "-"),
    strand2=ifelse(OrientEnd > 0, "+", "-"))) %>%
  mutate(SVtype=
           ifelse(chrom1 != chrom2, "TRA",
           ifelse(strand1 == "+" & strand2 == "-", "DEL",
           ifelse(strand1 == "-" & strand2 == "+", "DUP",
           ifelse(strand1 == "+" & strand2 == "-", "h2hINV", "t2tINV"))))) %>%
  filter(chrom1 %in% acceptable_chr & chrom2 %in% acceptable_chr)

sscndf = with(lnx_cns, data.frame(
  sampleId=SampleId,
  chrom=Chromosome,
  start=Start,
  end=End,
  total_cn=round(CopyNumber))

# Issues
# - Does not use ASCN
# - Does not handle subclonality
# - Does not handle multiple distinct instances of chromothripsis on the same chromosome
# - Add all SVs in chromothripsis interval to
# - Does it break up chromothripsis+BFB events?
# -


