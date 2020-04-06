library(tidyverse)
library(readr)
library(BSgenome.Hsapiens.UCSC.hg19, quietly=TRUE)
lnx_svs = read_csv("D:/hartwig/sv/paper_hpc/LNX_SVS.csv")
lnx_to_gr <- function(lnx_svs) {
  lnx_svs = lnx_svs %>% replace_na(list(InsertSeq=""))
  grs = GRanges(
    seqnames=lnx_svs$ChrStart,
    ranges=IRanges(start=lnx_svs$PosStart, width=1),
    strand=ifelse(lnx_svs$OrientStart == 1, "+", "-"),
    InsertSeq=lnx_svs$InsertSeq,
    partner=ifelse(lnx_svs$ChrEnd == 0, NA_character_, paste0(lnx_svs$Id, "h")),
    Id=lnx_svs$Id,
    SampleId=lnx_svs$SampleId,
    beid=paste0(lnx_svs$Id, ifelse(is.na(lnx_svs$ChrEnd), "b",  "o")))
  names(grs)=grs$beid
  lnx_svs = lnx_svs %>% filter(ChrEnd != 0)
  rc_insert_sequence = lnx_svs$InsertSeq
  rc_insert_sequence[!str_detect(rc_insert_sequence, "[^ACGTN]")] = as.character(reverseComplement(DNAStringSet(rc_insert_sequence[!str_detect(rc_insert_sequence, "[^ACGTN]")])))
  grh = GRanges(
    seqnames=lnx_svs$ChrEnd,
    ranges=IRanges(start=lnx_svs$PosEnd, width=1),
    strand=ifelse(lnx_svs$OrientEnd == 1, "+", "-"),
    insertSequence=ifelse(lnx_svs$OrientStart != lnx_svs$OrientEnd, lnx_svs$InsertSeq, rc_insert_sequence),
    partner=paste0(lnx_svs$Id, "o"),
    Id=lnx_svs$Id,
    SampleId=lnx_svs$SampleId,
    beid=paste0(lnx_svs$Id, "h"))
  names(grh)=grh$beid
  return(c(grs, grh))
}
svgr = lnx_to_gr(lnx_svs)
seqlevelsStyle(svgr) = "UCSC"
refseq = getSeq(BSgenome.Hsapiens.UCSC.hg19, names=flank(svgr, 50, both=TRUE), as.character=TRUE)
revSeq = as.character(reverseComplement(DNAStringSet(refseq)))
fwdseq = ifelse(as.logical(strand(svgr) == "+"), refseq, revSeq)
dnass = DNAStringSet(fwdseq)
names(dnass) = names(svgr)
writeXStringSet(dnass, "W:/projects/hartwig/motif/lnx_svs_flank50bp_noinsseq.fa")
