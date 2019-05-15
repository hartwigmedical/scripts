library(GenomicRanges)
library(tidyverse)
library(stringr)
library(grid)
library(gridExtra)
library(cowplot)

figdir="~/../Dropbox (HMF Australia)/HMF Australia team folder/Structural Variant Analysis/figures/"
sampleCancerTypes= (read.csv(paste(figdir,'../sample_cancer_types.csv',sep=''), header = T, stringsAsFactors = F))

# from http://github.com/PapenfussLab/sv_benchmark
import.repeatmasker.fa.out <- function(repeatmasker.fa.out) {
  rmdt <- read_table2(repeatmasker.fa.out, col_names=FALSE, skip=3)
  grrm <- GRanges(
    seqnames=rmdt$X5,
    ranges=IRanges(start=rmdt$X6 + 1, end=rmdt$X7),
    repeatType=rmdt$X10,
    repeatClass=rmdt$X11)
  #grrm$repeatClass <- str_replace(str_replace(grrm$repeatClass, "[?]", ""), "/.*", "")
  return(grrm)
}
# wget http://www.repeatmasker.org/genomes/hg19/RepeatMasker-rm330-db20120124/hg19.fa.out.gz
grrm = import.repeatmasker.fa.out("D:/hartwig/hg19.fa.out")
seqlevelsStyle(grrm) = "NCBI"

#svData = read_csv('C:/Users/Daniel/Dropbox (HMF Australia)/HMF Australia team folder/Structural Variant Analysis/SVA_SVS.csv')
svData = read_csv('../../sv/SVA_SVS.csv')

# All (!) non-template inserted sequence 16 or longer have an alignment to either the reference genome, or our viral db
svData %>% filter(str_length(InsertSeq) >= 16) %>% mutate(hasAln = InsSeqAlignments != "") %>% pull(hasAln) %>% table()

insseqgr = with(svData %>% filter(!is.na(InsSeqAlignments)) %>%
    dplyr::select(Id, InsSeqAlignments) %>%
    mutate(InsSeqAlignments=str_split(InsSeqAlignments, ";")) %>%
    unnest(InsSeqAlignments) %>%
    # many mapping locations can result in  truncation of the insertion alignment field
    # just drop the last one if it doesn't parse
    filter(str_detect(InsSeqAlignments, "[^:]+:[0-9]+[|][-+][|]([0-9]+[MIDNSHPX=])*[|]")) %>%
    separate(InsSeqAlignments, sep="[:|]", into=c("chr", "start", "orientation", "cigar", "maqp")) %>%
    mutate(
      start=as.integer(start),
      end=start+GenomicAlignments::cigarWidthAlongReferenceSpace(cigar)),
  GRanges(seqnames=chr, ranges=IRanges(start=start, end=end), Id=Id))
hits = findOverlaps(insseqgr, grrm, select="all")
hits = hits %>% as.data.frame() %>%
  mutate(
    overlap=pmin(end(insseqgr[queryHits]), end(grrm[subjectHits])) - pmax(start(insseqgr[queryHits]), start(grrm[subjectHits])),
    repeatOverlapPercentage=overlap/(end(insseqgr[queryHits])-start(insseqgr[queryHits])),
    repeatType=grrm[subjectHits]$repeatType,
    repeatClass=grrm[subjectHits]$repeatClass,
    Id=insseqgr[queryHits]$Id)
insrmdf = hits %>%
  group_by(Id) %>%
  top_n(1, repeatOverlapPercentage) %>%
  distinct(Id, .keep_all=TRUE) %>%
  ungroup() %>%
  dplyr::select(Id, repeatType, repeatClass, repeatOverlapPercentage)
svDataWithIns = svData %>% left_join(insrmdf, by="Id")
ggplot(insrmdf) +
  aes(x=repeatOverlapPercentage) +
  geom_histogram(bins=100) +
  scale_y_log10()

svDataWithIns = svDataWithIns %>%
  mutate(ins_class=ifelse(is.na(InsSeqAlignments),
                          "No alignment",
                          ifelse(!str_detect(InsSeqAlignments, "(^([0-9]?[0-9])|X|Y|(MT))[:]"), "Viral Insertion",
                                 ifelse(!is.na(repeatOverlapPercentage) & repeatOverlapPercentage > 0,
                                        "Repeat",
                                        "Not repeat"))))
ggplot(svDataWithIns %>% filter(!is.na(InsertSeq))) +
  aes(x=str_length(InsertSeq), fill=ins_class) +
  geom_histogram(bins=10) +
  scale_x_log10() +
  scale_fill_brewer(palette="Dark2") +
  facet_grid(ifelse(ChrEnd == 0, "Breakend", "Breakpoint") ~ .) +
  labs(
    x="Assembled sequence length",
    y="count",
    fill="Sequence\nAnnotation"
    )
ggsave(filename=paste0(figdir, "insert_sequence_annotation_histogram.pdf"), width=4, height=6)
svDataWithIns %>%
  filter(ins_class=="Viral Insertion") %>%
  pull(SampleId) %>%
  unique() %>%
  length()


# Viral insertions
virushostdb.tsv = read_tsv("../virus/virushostdb.tsv")
names(virushostdb.tsv) = str_replace_all(names(virushostdb.tsv), " ", "")
virushostlookup = virushostdb.tsv %>%
  dplyr::select(virustaxid, virusname, refseqid) %>%
  mutate(refseqid=str_replace_all(refseqid, " ", "")) %>%
  separate_rows(refseqid, sep = ",") %>%
  filter(!duplicated(refseqid))
get_best_viral_hit = function(x) {
  df = filter(x, !duplicated(virusname))
  if (nrow(df) == 1) return (df)
  df2 = df %>%
    mutate(virusname=Biobase::lcPrefix(virusname)) %>%
    filter(!duplicated(virusname))
  if (nrow(df2) == 1 & all(str_length(df2$virusname) > 10)) return (df2)
  df = df %>% filter(!(virusname %in% c("Alphapapillomavirus 7", "HBV recombinant C/G", "Hepatitis B virus subtype adr")))
  if (nrow(df) == 1) return (df)
  df2 = top_n(df, 1, match_length)
  if (nrow(df2) == 1) return(df2)
  browser()
  return(df)
}
viral_insertions = svDataWithIns %>%
  filter(ins_class=="Viral Insertion") %>%
  dplyr::select(SampleId, ChrStart, PosStart, InsSeqAlignments) %>%
  separate_rows(InsSeqAlignments, sep=";") %>%
  mutate(bin=paste(ChrStart, ceiling(PosStart / 100000))) %>%
  mutate(viral_contig=str_match(InsSeqAlignments, "^([^:]+):[0-9]+")[,2]) %>%
  mutate(match_length=str_match(InsSeqAlignments, "([0-9]+)M")[,2]) %>%
  inner_join(virushostlookup, by=c("viral_contig"="refseqid")) %>%
  group_by(SampleId, bin) %>%
  do(get_best_viral_hit(.))

viral_translation=c(
  "Human papillomavirus type 16"="HPV",
  "Human papillomavirus type 18"="HPV",
  "Alphapapillomavirus 7"="HPV",
  "Human papillomavirus type 45"="HPV",
  "Merkel cell polyomavirus"="Merkel cell",
  "HBV genotype D"="HepB",
  "Human betaherpesvirus 5"="Herpes",
  "Human alphaherpesvirus 1"="Herpes",
  "Adeno-associated virus - 2"="AAV",
  "Adeno-associated virus 3B"="AAV",
  "Adeno-associated virus - 6"="AAV",
  "HBV genotype C"="HepB",
  "HBV recombinant C/G"="HepB",
  "Hepatitis B virus subtype adr"="HepB",
  "HBV genotype B"="HepB",
  "Hepatitis B virus"="HepB")

viral_insertions = svDataWithIns %>%
  filter(ins_class=="Viral Insertion") %>%
  dplyr::select(SampleId, ChrStart, PosStart, InsSeqAlignments) %>%
  separate_rows(InsSeqAlignments, sep=";") %>%
  mutate(bin=paste(ChrStart, ceiling(PosStart / 100000))) %>%
  mutate(viral_contig=str_match(InsSeqAlignments, "^([^:]+):[0-9]+")[,2]) %>%
  mutate(match_length=str_match(InsSeqAlignments, "([0-9]+)M")[,2]) %>%
  inner_join(virushostlookup, by=c("viral_contig"="refseqid")) %>%
  group_by(SampleId, bin) %>%
  do(get_best_viral_hit(.))

viral_insertions_simple_name = svDataWithIns %>%
  filter(ins_class=="Viral Insertion") %>%
  dplyr::select(SampleId, ChrStart, PosStart, InsSeqAlignments) %>%
  separate_rows(InsSeqAlignments, sep=";") %>%
  mutate(bin=paste(ChrStart, ceiling(PosStart / 100000))) %>%
  mutate(viral_contig=str_match(InsSeqAlignments, "^([^:]+):[0-9]+")[,2]) %>%
  mutate(match_length=str_match(InsSeqAlignments, "([0-9]+)M")[,2]) %>%
  inner_join(virushostlookup, by=c("viral_contig"="refseqid")) %>%
  mutate(virusname=viral_translation[virusname]) %>%
  group_by(SampleId, bin) %>%
  summarise(virusname=unique(virusname)) %>%
  dplyr::select(-bin)

viral_summary = viral_insertions_simple_name %>%
  left_join(sampleCancerTypes) %>%
  replace_na(list("CancerType"="Unknown")) %>%
  mutate(CancerType=str_to_lower(CancerType)) %>%
  group_by(CancerType, SampleId, virusname) %>%
  summarise(sites=n()) %>%
  group_by(CancerType, virusname) %>%
  summarise(samples=n(), sites=sum(sites))






ins16df = svDataWithIns %>%
  mutate(
    PosStart=as.integer(PosStart),
    PosEnd=as.integer(PosEnd)) %>%
  filter(str_length(InsertSeq) >= 16) %>%
  mutate(
    firstInsAlignment_chr=str_extract(InsSeqAlignments, "^[^:]+"),
    firstInsAlignment_start=as.integer(str_match(InsSeqAlignments, "^[^:]+:([0-9]+)[|][-+][|]([0-9]+[MIDNSHPX=])*[|]")[,2]),
    firstInsAlignment_cigar=str_match(InsSeqAlignments, "^[^:]+:[0-9]+[|][-+][|]([0-9]+[MIDNSHPX=])*[|]")[,2],
    firstInsAlignment_orientation=str_match(InsSeqAlignments, "^[^:]+:[0-9]+[|]([-+])[|]([0-9]+[MIDNSHPX=])*[|]")[,2],
    firstInsAlignment_end=firstInsAlignment_start + GenomicAlignments::cigarWidthAlongReferenceSpace(firstInsAlignment_cigar),
    insAlnDistance1=ifelse(firstInsAlignment_chr == ChrStart, pmin(abs(PosStart - firstInsAlignment_start), abs(PosStart - firstInsAlignment_end)), NA_integer_),
    insAlnDistance2=ifelse(firstInsAlignment_chr == ChrEnd, pmin(abs(PosEnd - firstInsAlignment_start), abs(PosEnd - firstInsAlignment_end)), NA_integer_),
    insAlnDistance=ifelse(is.na(insAlnDistance1), insAlnDistance2, ifelse(is.na(insAlnDistance2), insAlnDistance1, pmin(insAlnDistance1, insAlnDistance1))))

ggplot(ins16df) +
  aes(x=insAlnDistance + 1, fill=Type, color=is.na(repeatClass)) +
  geom_histogram() +
  scale_x_log10()


ins16df %>%
  filter(Type %in% c("DEL", "DUP", "INS")) %>%
  mutate(
    size=str_length(InsertSeq)+abs(as.integer(PosEnd)-as.integer(PosStart)),
    sizeBin=pmin(2**10, 2**floor(log2(size+1))),
    basesInEvent=ifelse(firstInsAlignment_start > PosEnd | firstInsAlignment_end < PosStart, 0, pmin(firstInsAlignment_end, PosEnd) - pmax(firstInsAlignment_start, PosStart)),
    portionOfInsertionInEvent=basesInEvent / str_length(InsertSeq)) %>%
  filter(!is.na(firstInsAlignment_orientation)) %>%
ggplot() +
  aes(x=insAlnDistance + 1, fill=as.factor(cut(portionOfInsertionInEvent, seq(0, 1, 0.25)))) +
  geom_histogram() +
  scale_x_log10() +
  facet_grid(sizeBin ~ Type + firstInsAlignment_orientation, scales="free")

fqNoAln100bp = svDataWithIns %>% filter(!is.na(InsertSeq) & is.na(InsSeqAlignments) & str_length(InsertSeq) > 100) %>% mutate(fastq=paste0(">", Id, "\n", InsertSeq)) %>% summarise(fq=paste0(fastq, collapse="\n")) %>% pull(fq)
write(fqNoAln100bp, file="D:/hartwig/fqNoAln100bp.fastq")
