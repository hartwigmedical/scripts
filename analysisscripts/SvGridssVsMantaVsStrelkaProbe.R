library(tidyr)
library(dplyr)
library(RMySQL)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(seqinr)

load(file = "/Users/jon/hmf/analysis/mantaVgridss/scopedData.RData")

#gridss %>% filter(type != 'SGL') %>% count()
#gridss %>% filter(type != 'SGL') %>% group_by(type) %>% count()

privateStrelka = strelka %>% filter(scope == 'Private') %>% 
  mutate(
    startChromosome = chromosome, 
    endChromosome = chromosome, 
    endPosition = position + nchar(ref), 
    startOrientation = 1,
    endOrientation = -1,
    insertSequence = substring(alt, 2), 
    insertSequenceLength = nchar(insertSequence)) %>%
  select(sampleId, startChromosome, endChromosome, startPosition = position, endPosition, startOrientation, endOrientation, ref, alt, insertSequence, insertSequenceLength, scope)

privateManta = manta %>% filter(scope == 'Private') %>% 
  select(sampleId, startChromosome, endChromosome, startPosition, endPosition, insertSequence = SVINSSEQ_start, startOrientation, endOrientation, type = SVTYPE_start, scope,IMPRECISE_start, IMPRECISE_end, CIPOS_start, CIPOS_end)
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

probeOverlaps <- function(svs, contigLength = 120) {
  positionOffset = contigLength / 2 - 1
  
  if (!"startHomologySequence" %in% names(svs)) {
    svs[ , c("ploidy","startHomologySequence","endHomologySequence", "qualScore", "startLinkedBy", "endLinkedBy", "recoveryFilter")] <- NA
  }
  
  if (!"IMPRECISE_start" %in% names(svs)) {
    svs[ , c("IMPRECISE_start","IMPRECISE_end","CIPOS_start", "CIPOS_end")] <- NA
  }
  
  overlaps = svs %>% 
    mutate(lengthInsertSequence = nchar(insertSequence)) %>%
    filter(lengthInsertSequence < 120) %>%
    mutate(
      firstBreakendStart = ifelse(startOrientation == 1, startPosition - positionOffset + ceiling(lengthInsertSequence / 2), startPosition), 
      firstBreakendEnd = ifelse(startOrientation == 1, startPosition, startPosition + positionOffset -ceiling(lengthInsertSequence / 2)),
      secondBreakendStart = ifelse(endOrientation == 1, endPosition - positionOffset + floor(lengthInsertSequence / 2), endPosition), 
      secondBreakendEnd = ifelse(endOrientation == 1, endPosition, endPosition + positionOffset - floor(lengthInsertSequence / 2)),
      length = firstBreakendEnd - firstBreakendStart + secondBreakendEnd - secondBreakendStart + 2 + lengthInsertSequence
    ) %>% select(sampleId, id, event, lengthInsertSequence,
                 firstBreakendChromosome = startChromosome, firstBreakendStart, firstBreakendEnd, firstBreakendOrientation = startOrientation,
                 insertSequence, 
                 secondBreakendChromosome = endChromosome, secondBreakendStart, secondBreakendEnd, secondBreakendOrientation = endOrientation,
                 startPosition, endPosition, scope, type,
                 ploidy, startHomologySequence, endHomologySequence, qualScore, startLinkedBy, endLinkedBy, insertSequence, recoveryFilter,
                 IMPRECISE_start, IMPRECISE_end, CIPOS_start, CIPOS_end) %>%
    mutate(firstBreakendOrientation = as.numeric(firstBreakendOrientation), secondBreakendOrientation = as.numeric(secondBreakendOrientation))
  
  overlaps = overlaps %>%
    mutate(
      anchorStart = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, paste0("chr", firstBreakendChromosome), firstBreakendStart, firstBreakendEnd)),
      anchorEnd = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, paste0("chr", secondBreakendChromosome), secondBreakendStart, secondBreakendEnd)),
      anchorEnd = ifelse(firstBreakendOrientation == secondBreakendOrientation, reverse(chartr("ATGC", "TACG", anchorEnd)), anchorEnd),
      contigPosition = ifelse(firstBreakendOrientation == 1, nchar(anchorStart), nchar(anchorEnd)),
      probe = ifelse(firstBreakendOrientation == 1, paste0(anchorStart, insertSequence, anchorEnd), paste0(anchorEnd, insertSequence, anchorStart)),
      length = nchar(probe)
    )
  
  return (overlaps)
}

probe_single <- function(svs, maxProbeLength = 120) {
  maxSequenceSize = maxProbeLength / 2
  
  result = svs %>% filter(type %in% c('SGL')) %>%
    mutate(
      lengthInsertSequence = nchar(insertSequence), 
      insertSequenceStart = substr(insertSequence, 1, maxSequenceSize),
      insertSequenceEnd = substr(insertSequence, lengthInsertSequence - maxSequenceSize + 1, lengthInsertSequence),
      truncatedInsertSequence = ifelse(startOrientation == 1, insertSequenceStart, insertSequenceEnd),
      truncatedInsertSequenceLength = nchar(truncatedInsertSequence),
      probeStart = ifelse(startOrientation == 1, startPosition - maxProbeLength + 1 + truncatedInsertSequenceLength, startPosition), 
      probeEnd = ifelse(startOrientation == 1, startPosition, startPosition + maxProbeLength - 1 - truncatedInsertSequenceLength), 
      length = probeEnd - probeStart + 1 + truncatedInsertSequenceLength
    ) %>% 
    select(sampleId, chromosome = startChromosome, probeStart, probeEnd, truncatedInsertSequence, orientation = startOrientation, lengthInsertSequence, event, id, startPosition, scope, type, 
           ploidy, startHomologySequence, endHomologySequence, qualScore, startLinkedBy, endLinkedBy, insertSequence, recoveryFilter)
  
  result = result %>%
    mutate(
      probe = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, paste0("chr", chromosome), probeStart, probeEnd)),
      contigPosition = ifelse(orientation == 1, nchar(probe), nchar(truncatedInsertSequence)),
      probe = ifelse(orientation == 1, paste0(probe, truncatedInsertSequence), paste0(truncatedInsertSequence, probe)),
      length = nchar(probe)
    )
  
  return (result)
}

gridssSingles = probe_single(gridss %>% filter(type == 'SGL'), maxProbeLength = 2000) %>% mutate(contig = paste0(event, "_", id))
gridssOverlaps = probeOverlaps(gridss %>% filter(type != 'SGL'), contigLength = 2000) %>% mutate(contig = paste0(event, "_", id)) 
mantaOverlaps = probeOverlaps(privateManta %>% mutate(id = row_number(), event = "MANTA"), contigLength = 2000) %>% mutate(contig = paste0(event, "_", id))
strelkaOverlaps = probeOverlaps(privateStrelka %>% mutate(id = row_number(), event = "STRELKA", type = ifelse(nchar(ref) > nchar(alt), "DEL", "INS")), contigLength = 2000) %>% mutate(contig = paste0(event, "_", id))
allOverlaps = bind_rows(gridssOverlaps, mantaOverlaps) %>% bind_rows(strelkaOverlaps)
save(allOverlaps, file = "~/hmf/analysis/probes/allOverlaps.RData")

##### CREATE SAGE HOTSPOT FILE
allOverlaps = allOverlaps %>% mutate(uniqueId = contig)
#%>% filter(sampleId == 'CPCT02370037T')

hotspotLeft = allOverlaps %>% select(sampleId, chromosome = firstBreakendChromosome, position = startPosition, uniqueId, scope)
hotspotRight = allOverlaps %>% select(sampleId, chromosome = secondBreakendChromosome, position = endPosition, uniqueId, scope)
hotspotOverlap = allOverlaps %>% select(sampleId, chromosome = contig, position = contigPosition, uniqueId, scope)

gridssSingles = gridssSingles %>% mutate(endPosition = startPosition +  orientation, uniqueId = contig)
save(gridssSingles, file = "~/hmf/analysis/probes/gridssSingles.RData")

sglLeft = gridssSingles %>% select(sampleId, chromosome, position = startPosition, uniqueId, scope)
sglRight = gridssSingles %>% select(sampleId, chromosome, position = endPosition, uniqueId, scope)
sglOverlap = gridssSingles %>% select(sampleId, chromosome = contig, position = contigPosition, uniqueId, scope)
sglHotspots = bind_rows(sglLeft, sglRight) %>% bind_rows(sglOverlap) 

hotspots = bind_rows(hotspotLeft, hotspotRight) %>% 
  bind_rows(hotspotOverlap) %>% 
  bind_rows(sglHotspots) %>% 
  mutate(ref = 'N', alt = 'A')

save(hotspots, file = "~/hmf/analysis/probes/hotspots.RData")

uniqueHotspots = hotspots %>% select(chromosome, position, ref, alt) %>% distinct(chromosome, position, ref, alt)
write.table(uniqueHotspots, file = "~/hmf/analysis/probes/probeHotspots.tsv", row.names = F, col.names = F, sep = "\t", quote = F)

##### CREATE CONTIGS FOR SYNTHETIC REF GENOME
gridssSglContigs =  gridssSingles %>% select(contig, probe)
gridssContigs =  gridssOverlaps %>% select(contig, probe)
mantaContigs = mantaOverlaps %>% select(contig, probe)
strelkaContigs =  strelkaOverlaps %>% select(contig, probe)

contigs = bind_rows(gridssContigs, mantaContigs) %>% bind_rows(strelkaContigs) %>% bind_rows(gridssSglContigs) %>% mutate(len = nchar(probe))
write.fasta(sequences = as.list(contigs$probe), names = (contigs$contig), file.out = "~/hmf/analysis/svPaper/probes.fasta", nbchar = 60)

##### CREATE 120 BASE PROBES
gridssStartProbes = probeBeforeBreakend(gridss %>% select(sampleId, chromosome = startChromosome, position = startPosition, orientation = startOrientation, scope))
gridssEndProbes = probeBeforeBreakend(gridss%>% filter(!is.na(endChromosome)) %>% select(sampleId, chromosome = endChromosome, position = endPosition, orientation = endOrientation, scope))
gridssStartAndEndProbes = bind_rows(gridssStartProbes, gridssEndProbes)
rm(gridssStartProbes, gridssEndProbes)

mantaStartProbes = probeBeforeBreakend(privateManta %>% select(sampleId, chromosome = startChromosome, position = startPosition, orientation = startOrientation))
mantaEndProbes = probeBeforeBreakend(privateManta %>% select(sampleId, chromosome = endChromosome, position = endPosition, orientation = endOrientation, scope))
mantaStartAndEndProbes = bind_rows(mantaStartProbes, mantaEndProbes)
rm(mantaStartProbes, mantaEndProbes)

strelkaStartProbes = probeBeforeBreakend(privateStrelka %>% select(sampleId, chromosome = startChromosome, position = startPosition, orientation = startOrientation))
strelkaEndProbes = probeBeforeBreakend(privateStrelka %>% select(sampleId, chromosome = endChromosome, position = endPosition, orientation = endOrientation, scope))
strelkaStartAndEndProbes = bind_rows(strelkaStartProbes, strelkaEndProbes)
rm(strelkaStartProbes, strelkaEndProbes)

gridssOverlaps = probeOverlaps(gridss %>% filter(type != 'SGL') %>% select(sampleId, startChromosome, endChromosome, startPosition, endPosition, insertSequence, startOrientation, endOrientation))
mantaOverlaps = probeOverlaps(privateManta %>% select(sampleId, startChromosome, endChromosome, startPosition, endPosition, insertSequence, startOrientation, endOrientation))
strelkaOverlaps = probeOverlaps(privateStrelka %>% select(sampleId, startChromosome, endChromosome, startPosition, endPosition, insertSequence, startOrientation, endOrientation))

gridssSingles = probe_single(gridss, 120)

probes = data.frame()
probes = bind_rows(probes, gridssStartAndEndProbes %>% select(sampleId, probe) %>% mutate(source = "GRIDSS"))
probes = bind_rows(probes, gridssOverlaps %>% select(sampleId, probe)%>% mutate(source = "GRIDSS"))
probes = bind_rows(probes, gridssSingles %>% select(sampleId, probe)%>% mutate(source = "GRIDSS"))
probes = bind_rows(probes, mantaStartAndEndProbes %>% select(sampleId, probe)%>% mutate(source = "MANTA"))
probes = bind_rows(probes, mantaOverlaps %>% select(sampleId, probe)%>% mutate(source = "MANTA"))
probes = bind_rows(probes, strelkaStartAndEndProbes %>% select(sampleId, probe)%>% mutate(source = "STRELKA"))
probes = bind_rows(probes, strelkaOverlaps %>% select(sampleId, probe)%>% mutate(source = "STRELKA"))
probes = probes %>% distinct(probe) %>% filter(!grepl(',', probe), !grepl('N', probe))


write.csv(probes, file = "/Users/jon/hmf/analysis/mantaVgridss/probes.csv", row.names = F)

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

strelkaVerification1 = strelkaOverlaps %>% filter(sampleId == 'CPCT02030461T')

overlapCommmands = strelkaVerification1 %>%
  mutate(
    sam= paste0("/data/common/tools/samtools_v1.2/samtools view CPCT02030461R_CPCT02030461T.assembly.bam.sv.bam ", firstBreakendChromosome, ":" ,firstBreakendStart, "-", firstBreakendEnd, "> verification.sam"),
    grep = paste0("grep ", probe, " verification.sam")
  ) %>%
  select(sam, grep)

strelkaVerification2 = strelkaStartAndEndProbes %>% filter(sampleId == 'CPCT02030461T')
startEndCommmands = strelkaVerification2 %>%
  mutate(
    sam= paste0("/data/common/tools/samtools_v1.2/samtools view CPCT02030461R_CPCT02030461T.assembly.bam.sv.bam ", chromosome, ":" ,probeStart, "-", probeEnd, "> verification.sam"),
    grep = paste0("grep ", probe, " verification.sam")
  ) %>%
  select(sam, grep)

gridss %>% select(sampleId) %>% distinct()
  