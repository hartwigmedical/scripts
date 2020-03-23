library(devtools)
#install_github("parklab/ShatterSeek")
library(ShatterSeek)
library(tidyverse)
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

sscndf = with(lnx_cns,
  data.frame(
    sampleId=SampleId,
    chrom=Chromosome,
    start=Start,
    end=End,
    total_cn=round(CopyNumber))) %>%
  filter(chrom %in% acceptable_chr)

# Issues
# - Does not use ASCN
# - Does not handle subclonality
# - Does not handle multiple distinct instances of chromothripsis on the same chromosome
# - Add all SVs in chromothripsis interval to
# - Does it break up chromothripsis+BFB events?
# -
#sid = "CPCT02020542T"
library(doParallel)
registerDoParallel(cores=16)
foreach(sid=unique(sssvdf$sampleId)) %dopar% {
#for (sid in unique(sssvdf$sampleId)) {
  require(tidyverse)
  require(ShatterSeek)
  file_name = paste0("D:/hartwig/svtoolkit/shatterseek/shatterseek_", sid, ".RData")
  if (!file.exists(file_name)) {
    SV_data <- with(sssvdf %>% filter(sampleId == sid), SVs(
      chrom1=chrom1,
      pos1=Pos1,
      chrom2=chrom2,
      pos2=Pos2,
      SVtype=SVtype,
      strand1=strand1,
      strand2=strand2))
    CN_data =  with(sscndf %>% filter(sampleId == sid), CNVsegs(
      chrom=chrom,
      start=start,
      end=end,
      total_cn=total_cn))

    chromothripsis = shatterseek(SV.sample=SV_data, seg.sample=CN_data)
    saveRDS(chromothripsis, file_name)
  }
}
