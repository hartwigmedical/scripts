library(VariantAnnotation)
library(tidyr)
library(dplyr)

vcf_data_frame<- function(sample) {
  
  file = paste0("~/hmf/analysis/probes/", sample, ".vcf")
  vcf = readVcf(file)
  
  vcf.info = info(vcf)
  vcf.geno = geno(vcf)
  
  vcf.df = data.frame(
    chromosome = seqnames(vcf), 
    position = start(vcf), 
    reads = vcf.geno$DP[,1])
  
  return (vcf.df)
}

probe_summary<- function(sample, hotspots) {

  df = vcf_data_frame(sample) %>% mutate(sampleId = sample)
  
  combined = left_join(df, hotspots %>% select(chromosome, position, uniqueId, source), by = c("chromosome", "position")) %>% 
    group_by(uniqueId) %>%
    mutate(
      last = position == max(position), 
      standard = chromosome %in% c(1:22, "X","Y"),
      type = ifelse(last, "Right", "Left"),
      type = ifelse(standard, type, "SV"),
      reads = as.numeric(reads),
      type = factor(type, levels = c("Left", "Right", "SV"))) %>%
    select(source, sampleId, uniqueId, reads, type, chromosome, position) %>%
    ungroup() %>%
    distinct()

  
  # NOTES 
  #   By definition, cannot get 0 reads in an alt contig... it will just be absent.
  #   If there is an entry here at all, it is okay to set refDepth to 0
  
  result = combined %>% 
    group_by(uniqueId) %>%  
    spread(type, reads, fill = 0) %>%
    mutate(
      refDepth = pmax(Left, Right, na.rm = T),
      altDepth = SV, 
      sourceRefDepth = ifelse(sampleId == source, refDepth, 0),
      sourceAltDepth = ifelse(sampleId == source, altDepth, 0),
      otherRefDepth = ifelse(sampleId != source, refDepth, 0),
      otherAltDepth = ifelse(sampleId != source, altDepth, 0),
      otherAltCount = ifelse(sampleId != source & altDepth > 0, 1, 0)
    ) %>%
    select(source, uniqueId, sourceRefDepth, sourceAltDepth, otherRefDepth, otherAltDepth, otherAltCount)
  
  #result = combined %>% 
  #  group_by(uniqueId) %>%  
  #  spread(type, reads) %>%
  #  mutate(
  #    refDepth = pmax(Left, Right, na.rm = T),
  #    refDepth = ifelse(is.na(refDepth), 0, refDepth),
  #    support = paste0(SV, "|", refDepth),
  #    support = ifelse(is.na(SV), "None", support)
  #  ) %>% 
  #  select(source, uniqueId, support)

  #names(result) <- c("source", "uniqueId", sample)

  return(result)
}

load(file = "~/hmf/analysis/probes/gridssSingles.RData")
load(file = "~/hmf/analysis/probes/allOverlaps.RData")
svData = allOverlaps %>% select(uniqueId = contig, startChromosome = firstBreakendChromosome, startPosition, startOrientation = firstBreakendOrientation, endChromosome = secondBreakendChromosome, endPosition, endOrientation = secondBreakendOrientation, type, scope, ploidy, startHomologySequence, endHomologySequence, qualScore, startLinkedBy, endLinkedBy, insertSequence, recoveryFilter, IMPRECISE_start, IMPRECISE_end, CIPOS_start, CIPOS_end)
svSGLData = gridssSingles %>% select(uniqueId, startChromosome = chromosome, startPosition, startOrientation = orientation, type, scope, ploidy, startHomologySequence, endHomologySequence, qualScore, startLinkedBy, endLinkedBy, insertSequence, recoveryFilter)
svData = bind_rows(svData, svSGLData)


load(file = "~/hmf/analysis/probes/hotspots.RData")
hotspots = hotspots %>%
  mutate(
    source = sampleId,
    caller = ifelse(grepl("gridss", uniqueId), "Gridss", "Manta"),
    caller = ifelse(grepl("STRELKA", uniqueId), "Strelka", caller)
    ) %>% 
  select(-sampleId)

#sample = "CPCT02120143T"
samples = c("CPCT02030461T","CPCT02120143T","CPCT02330102T","CPCT02450014T","CPCT02030516T","CPCT02130091T","CPCT02370037T","DRUP01010096T","CPCT02050327T","CPCT02150016T","DRUP01330008T","CPCT02070386T","CPCT02160052T")

result = data.frame()
for (sample in samples) {
  cat ("Processing", sample, "\n")
  sampleResult = probe_summary(sample, hotspots)
  result = bind_rows(result, sampleResult)
}

resultSummary = result %>% group_by(source, uniqueId) %>% 
  summarise(
    sourceRefDepth = sum(sourceRefDepth), 
    sourceAltDepth = sum(sourceAltDepth), 
    otherRefDepth = sum(otherRefDepth),
    otherAltDepthMax = max(otherAltDepth),
    otherAltDepth = sum(otherAltDepth),
    otherAltCount = sum(otherAltCount)
    ) 
  
probeResultCondensed = full_join(resultSummary, svData, by = "uniqueId") %>%
  mutate(p=ppois(sourceAltDepth,otherAltDepthMax,FALSE,FALSE),length=endPosition-startPosition,supported=p<0.001&otherAltDepth<40&sourceAltDepth>1,callset=substr(uniqueId,1,5))

save(probeResultCondensed, file = "~/hmf/analysis/probes/probeResultCondensed.RData")


load( file = "~/hmf/analysis/probes/probeResultCondensed.RData")
probeQuality = read.table(file = "/Users/jon/hmf/analysis/probes/ProbeQuality.tsv", sep = "\t", header = F) %>% select(probe = V2, probeQuality = V3) %>% filter(probeQuality > 0) #%>%  filter(!grepl("gridss", probe), !grepl("MANTA", probe) )
probeResult = probeResultCondensed %>% left_join(probeQuality, by = c("uniqueId" = "probe")) %>% mutate(probeQuality = ifelse(is.na(probeQuality), 0, probeQuality))
probeResult = probeResult %>%
  mutate(
    callset = ifelse(callset == "grids", "Gridss", callset),
    callset = ifelse(callset == "STREL", "Strelka", callset),
    callset = ifelse(callset == "MANTA", "Manta", callset),
    simpleScope = ifelse(scope == "Private", "Private", "Shared"),
    supported = ifelse(is.na(supported), F, supported)
  )

save(probeResult, file = "~/hmf/analysis/probes/probeResult.RData")
load(file = "~/hmf/analysis/probes/probeResult.RData")

shortDup = probeResult %>% filter(source == 'CPCT02450014T', length < 1000, type == 'DUP') %>% 
  mutate(sequence = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, paste0("chr", startChromosome), startPosition, endPosition)))

resultSummary = probeResult %>% group_by(callset, simpleScope, supported, type) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(supported = ifelse(supported, "Found", "NotFound")) %>%
  unite(united, callset, simpleScope, supported) %>% spread(united, n)





