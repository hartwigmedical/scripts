library(tidyr)
library(dplyr)
library(RMySQL)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

load(file = "/Users/jon/hmf/analysis/mantaVgridss/scopedData.RData")

#gridss %>% filter(type != 'SGL') %>% count()
#gridss %>% filter(type != 'SGL') %>% group_by(type) %>% count()

privateManta = manta %>% filter(scope == 'Private') %>% select(sampleId, startChromosome, endChromosome, startPosition, endPosition, insertSequence = SVINSSEQ_start, startOrientation, endOrientation, type = SVTYPE_start, scope)
privateManta[is.na(privateManta)] <- ""
rm(manta, strelka)

probeBeforeBreakend <- function(svs) {
  probes = svs %>%
    mutate(
      orientation = as.numeric(orientation),
      probeStart = ifelse(orientation == 1, position - 140, position + 21), 
      probeEnd = ifelse(orientation == 1, position - 21, position + 140),
      probe = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, paste0("chr", chromosome), probeStart, probeEnd)),
      length = nchar(probe)
      ) 
    return (probes)
}

probeOverlaps <- function(svs) {
  overlaps = svs %>% 
    mutate(lengthInsertSequence = nchar(insertSequence)) %>%
    filter(lengthInsertSequence < 120) %>%
    mutate(
      firstBreakendStart = ifelse(startOrientation == 1, startPosition - 59 + ceiling(lengthInsertSequence / 2), startPosition), 
      firstBreakendEnd = ifelse(startOrientation == 1, startPosition, startPosition + 59 -ceiling(lengthInsertSequence / 2)),
      secondBreakendStart = ifelse(endOrientation == 1, endPosition - 59 + floor(lengthInsertSequence / 2), endPosition), 
      secondBreakendEnd = ifelse(endOrientation == 1, endPosition, endPosition + 59 - floor(lengthInsertSequence / 2)),
      length = firstBreakendEnd - firstBreakendStart + secondBreakendEnd - secondBreakendStart + 2 + lengthInsertSequence
    ) %>% select(sampleId, lengthInsertSequence,
                 firstBreakendChromosome = startChromosome, firstBreakendStart, firstBreakendEnd, firstBreakendOrientation = startOrientation,
                 insertSequence, 
                 secondBreakendChromosome = endChromosome, secondBreakendStart, secondBreakendEnd, secondBreakendOrientation = endOrientation)
  
  overlaps = overlaps %>%
    mutate(
      anchorStart = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, paste0("chr", firstBreakendChromosome), firstBreakendStart, firstBreakendEnd)),
      anchorEnd = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, paste0("chr", secondBreakendChromosome), secondBreakendStart, secondBreakendEnd)),
      anchorEnd = ifelse(firstBreakendOrientation == secondBreakendOrientation, reverse(chartr("ATGC", "TACG", anchorEnd)), anchorEnd),
      probe = ifelse(firstBreakendOrientation == 1, paste0(anchorStart, insertSequence, anchorEnd), paste0(anchorEnd, insertSequence, anchorStart)),
      length = nchar(probe)
    )
  
  return (overlaps)
}


gridssStartProbes = probeBeforeBreakend(gridss %>% select(sampleId, chromosome = startChromosome, position = startPosition, orientation = startOrientation, scope))
gridssEndProbes = probeBeforeBreakend(gridss%>% filter(!is.na(endChromosome)) %>% select(sampleId, chromosome = endChromosome, position = endPosition, orientation = endOrientation, scope))
gridssStartAndEndProbes = bind_rows(gridssStartProbes, gridssEndProbes)
rm(gridssStartProbes, gridssEndProbes)

mantaStartProbes = probeBeforeBreakend(privateManta %>% select(sampleId, chromosome = startChromosome, position = startPosition, orientation = startOrientation))
mantaEndProbes = probeBeforeBreakend(privateManta %>% select(sampleId, chromosome = endChromosome, position = endPosition, orientation = endOrientation, scope))
mantaStartAndEndProbes = bind_rows(mantaStartProbes, mantaEndProbes)
rm(mantaStartProbes, mantaEndProbes)


gridssOverlaps = probeOverlaps(gridss %>% filter(type != 'SGL') %>% select(sampleId, startChromosome, endChromosome, startPosition, endPosition, insertSequence, startOrientation, endOrientation))
mantaOverlaps = probeOverlaps(privateManta %>% select(sampleId, startChromosome, endChromosome, startPosition, endPosition, insertSequence, startOrientation, endOrientation))

gridssSingles = gridss %>% filter(type %in% c('SGL')) %>%
  mutate(
    lengthInsertSequence = nchar(insertSequence), 
    insertSequenceStart = substr(insertSequence, 1, 60),
    insertSequenceEnd = substr(insertSequence, lengthInsertSequence - 59, lengthInsertSequence),
    truncatedInsertSequence = ifelse(startOrientation == 1, insertSequenceStart, insertSequenceEnd),
    truncatedInsertSequenceLength = nchar(truncatedInsertSequence),
    probeStart = ifelse(startOrientation == 1, startPosition - 119 + truncatedInsertSequenceLength, startPosition), 
    probeEnd = ifelse(startOrientation == 1, startPosition, startPosition + 119 - truncatedInsertSequenceLength), 
    length = probeEnd - probeStart + 1 + truncatedInsertSequenceLength
    ) %>% 
  select(sampleId, chromosome = startChromosome, probeStart, probeEnd, truncatedInsertSequence, orientation = startOrientation)

gridssSingles = gridssSingles %>%
  mutate(
    probe = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, paste0("chr", chromosome), probeStart, probeEnd)),
    probe = ifelse(orientation == 1, paste0(probe, truncatedInsertSequence), paste0(truncatedInsertSequence, probe)),
    length = nchar(probe)
  )

probes = data.frame()
probes = bind_rows(probes, gridssStartAndEndProbes %>% select(sampleId, probe) %>% mutate(source = "GRIDSS"))
probes = bind_rows(probes, gridssOverlaps %>% select(sampleId, probe)%>% mutate(source = "GRIDSS"))
probes = bind_rows(probes, gridssSingles %>% select(sampleId, probe)%>% mutate(source = "GRIDSS"))
probes = bind_rows(probes, mantaStartAndEndProbes %>% select(sampleId, probe)%>% mutate(source = "MANTA"))
probes = bind_rows(probes, mantaOverlaps %>% select(sampleId, probe)%>% mutate(source = "MANTA"))
probes = probes %>% distinct(probe)

#PASS 
verification = gridssOverlaps %>% 
  filter(
    firstBreakendChromosome == secondBreakendChromosome, 
    sampleId == 'CPCT02030461T', 
    firstBreakendChromosome == 1, 
    firstBreakendOrientation != secondBreakendOrientation, 
    firstBreakendOrientation == 1)

#PASS 
verification2 = gridssOverlaps %>% 
  filter(
    firstBreakendChromosome != secondBreakendChromosome, 
    sampleId == 'CPCT02030461T', 
    firstBreakendChromosome == 1, 
    firstBreakendOrientation != secondBreakendOrientation, 
    firstBreakendOrientation == 1)

#PASS 
verification3 = gridssOverlaps %>% 
  filter(
    firstBreakendChromosome == secondBreakendChromosome, 
    sampleId == 'CPCT02030461T', 
    firstBreakendChromosome == 1, 
    firstBreakendOrientation != secondBreakendOrientation, 
    firstBreakendOrientation == -1)

#PASS 
verification4 = gridssOverlaps %>% 
  filter(
    firstBreakendChromosome != secondBreakendChromosome, 
    sampleId == 'CPCT02030461T', 
    firstBreakendChromosome == 1, 
    firstBreakendOrientation != secondBreakendOrientation, 
    firstBreakendOrientation == -1)

#PASS 
verification5 = gridssOverlaps %>% 
  filter(
    firstBreakendChromosome == secondBreakendChromosome, 
    sampleId == 'CPCT02030461T', 
    firstBreakendChromosome == 1, 
    firstBreakendOrientation == secondBreakendOrientation, 
    firstBreakendOrientation == 1)

#PASS 
verification6 = gridssOverlaps %>% 
  filter(
    firstBreakendChromosome == secondBreakendChromosome, 
    sampleId == 'CPCT02030461T', 
    firstBreakendChromosome == 1, 
    firstBreakendOrientation == secondBreakendOrientation, 
    firstBreakendOrientation == -1)

overlapCommmands = verification4 %>%
  mutate(
    sam= paste0("/data/common/tools/samtools_v1.2/samtools view CPCT02030461R_CPCT02030461T.assembly.bam.sv.bam ", firstBreakendChromosome, ":" ,firstBreakendStart, "-", firstBreakendEnd, "> verification.sam"),
    grep = paste0("grep ", probe, " verification.sam")
  ) %>%
  select(sam, grep)

#PASS 
verification7 = gridssStartAndEndProbes %>% 
  filter(
    sampleId == 'CPCT02030461T', 
    chromosome == 1, 
    orientation == -1)

#PASS 
verification8 = gridssStartAndEndProbes %>% 
  filter(
    sampleId == 'CPCT02030461T', 
    chromosome == 1, 
    orientation == 1)

startEndCommmands = verification9 %>%
  mutate(
    sam= paste0("/data/common/tools/samtools_v1.2/samtools view CPCT02030461R_CPCT02030461T.assembly.bam.sv.bam ", chromosome, ":" ,probeStart, "-", probeEnd, "> verification.sam"),
    grep = paste0("grep ", probe, " verification.sam")
  ) %>%
  select(sam, grep)


#PASS
verification9 = gridssSingles %>% 
  filter(
    sampleId == 'CPCT02050327T', orientation == 1)

#PASS
verification10 = gridssSingles %>% 
  filter(
    sampleId == 'CPCT02050327T', orientation == -1)


SGLCommmands = verification10 %>%
    mutate(
      sam= paste0("/data/common/tools/samtools_v1.2/samtools view CPCT02050327R_CPCT02050327T.assembly.bam.sv.bam ", chromosome, ":" ,probeStart, "-", probeEnd, "> verification.sam"),
      grep = paste0("grep ", probe, " verification.sam")
    ) %>%
    select(sam, grep)

  