
# awk ' { if ($8 > 2) print ; } ' < full.gridss_pon_breakpoint.bedpe > gridss_pon_breakpoint.bedpe
# awk ' { if ($5 > 2) print ; } ' < full.gridss_pon_single_breakend.bed > gridss_pon_single_breakend.bed
setwd("../gridss")
source("libgridss.R")
pon_dir="D:/hartwig/pon"
ponbpgr = read_gridss_breakpoint_pon(paste(pon_dir, "gridss_pon_breakpoint.bedpe", sep="/"))
ponbegr = import(paste(pon_dir, "gridss_pon_single_breakend.bed", sep="/"))

ponbpgr = ponbpgr[ponbpgr$score > 2]

dbload <- function(dbConnect) {
  dbdf = DBI::dbGetQuery(dbConnect, "SELECT id, startChromosome, startPosition, startOrientation, endChromosome, endPosition, endOrientation, qualScore, filter FROM structuralVariant")
  grs = GRanges(
    seqnames=dbdf$startChromosome,
    ranges=IRanges(start=dbdf$startPosition,
                   end=dbdf$startPosition),
    strand=ifelse(dbdf$startOrientation == 1, "+", "-"),
    QUAL=dbdf$qualScore,
    partner=ifelse(is.na(dbdf$endChromosome), NA_character_, paste0(dbdf$id, "h")),
    id=dbdf$id,
    beid=paste0(dbdf$id, ifelse(is.na(dbdf$endChromosome), "b",  "o")))
  names(grs)=grs$beid
  dbdf = dbdf %>% filter(!is.na(endChromosome))
  grh = GRanges(
    seqnames=dbdf$endChromosome,
    ranges=IRanges(start=dbdf$endPosition,
                   end=dbdf$endPosition),
    strand=ifelse(dbdf$endOrientation == 1, "+", "-"),
    QUAL=dbdf$qualScore,
    partner=paste0(dbdf$id, "o"),
    FILTER=dbdf$filter,
    id=dbdf$id,
    beid=paste0(dbdf$id, "h"))
  names(grh)=grh$beid
  return(c(grs, grh))
}
library(RMySQL)
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
dbgr = dbload(dbProd)
dbbpgr = dbgr[!is.na(dbgr$partner)]
dbepgr = dbgr[is.na(dbgr$partner)]


query = dbbpgr
subject = ponbpgr
hits <- bind_rows(
  as.data.frame(findOverlaps(query, subject, maxgap=4, type="any", select="all", ignore.strand=FALSE), row.names=NULL),
  as.data.frame(findOverlaps(partner(query), partner(subject), maxgap=4, type="any", select="all", ignore.strand=FALSE), row.names=NULL))
hits = hits %>% arrange(queryHits, subjectHits) %>%
  filter(!is.na(lead(queryHits)) & !is.na(lead(subjectHits)) & lead(queryHits) == queryHits & lead(subjectHits) == subjectHits) %>%
  mutate(count=as.numeric(ponbpgr$score[subjectHits]))
length(unique(hits$subjectHits))
median(as.numeric(ponbpgr$score)[unique(hits$subjectHits)])
median(as.numeric(ponbpgr$score))
mean(as.numeric(ponbpgr$score)[unique(hits$subjectHits)])
mean(as.numeric(ponbpgr$score))

hits = hits %>%
  group_by(queryHits) %>%
  summarise(hits=sum(count))
dbbpgr$ponhits = 0
dbbpgr$ponhits[hits$queryHits] = hits$hits

ggplot(as.data.frame(dbbpgr)) +
  aes(x=ponhits) +
  geom_histogram(bins=100) +
  scale_y_log10() +
  theme_bw() +
  labs("prod variants by PON hits\n
       3 variant PON inclusion")

ggplot(as.data.frame(dbbpgr)) +
  aes(x=ponhits) +
  geom_histogram(bins=100) +
  scale_x_continuous(limits=c(1,100)) +
  theme_bw() +
  labs("prod variants with PON hits\n
       3 variant PON inclusion")


ggplot(as.data.frame(dbbpgr) %>% filter(ponhits > 0)) +
  aes(x=ponhits, y=QUAL) +
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10()


median(dbbpgr$QUAL[dbbpgr$ponhits == 0])
median(dbbpgr$QUAL[dbbpgr$ponhits > 0])

median

