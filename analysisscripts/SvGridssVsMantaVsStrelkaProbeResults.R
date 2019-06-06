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
  
  combined = left_join(df, hotspots, by = c("chromosome", "position")) %>% 
    group_by(uniqueId) %>%
    mutate(
      last = position == max(position), 
      standard = chromosome %in% c(1:22, "X","Y"),
      type = ifelse(last, "Right", "Left"),
      type = ifelse(standard, type, "SV"),
      reads = as.numeric(reads),
      type = factor(type, levels = c("Left", "Right", "SV"))) %>%
    select(source, sampleId, uniqueId, reads, type) %>%
    ungroup() %>%
    distinct()

  result = combined %>% 
    group_by(uniqueId) %>%  
    spread(type, reads) %>%
    mutate(
      refDepth = pmax(Left, Right, na.rm = T),
      refDepth = ifelse(is.na(refDepth), 0, refDepth),
      support = paste0(SV, "|", refDepth),
      support = ifelse(is.na(SV), "None", support)
    ) %>% 
    select(source, uniqueId, support)

  names(result) <- c("source", "uniqueId", sample)

  return(result)
}

load(file = "~/hmf/analysis/probes/gridssSingles.RData")
load(file = "~/hmf/analysis/probes/allOverlaps.RData")
svData = allOverlaps %>% select(uniqueId = contig, startChromosome = firstBreakendChromosome, startPosition, startOrientation = firstBreakendOrientation, endChromosome = secondBreakendChromosome, endPosition, endOrientation = secondBreakendOrientation, type, scope)
svSGLData = gridssSingles %>% select(uniqueId, startChromosome = chromosome, startPosition, startOrientation = orientation, type, scope)
svData = bind_rows(svData, svSGLData)


load(file = "~/hmf/analysis/probes/hotspots.RData")
hotspots = hotspots %>%
  mutate(
    source = sampleId,
    caller = ifelse(grepl("gridss", uniqueId), "Gridss", "Manta"),
    caller = ifelse(grepl("STRELKA", uniqueId), "Strelka", caller)
    ) %>% 
  select(-sampleId)

samples = c("CPCT02030461T","CPCT02120143T","CPCT02330102T","CPCT02450014T","CPCT02030516T","CPCT02130091T","CPCT02370037T","DRUP01010096T","CPCT02050327T","CPCT02150016T","DRUP01330008T","CPCT02070386T","CPCT02160052T")

result = data.frame(source = as.character(), uniqueId = as.character())
for (sample in samples) {
  cat ("Processing", sample, "\n")
  sampleResult = probe_summary(sample, hotspots)
  result = full_join(result, sampleResult, by = c("source", "uniqueId"))
}

probeResultCondensed = full_join(result, svData, by = "uniqueId")
save(probeResultCondensed, file = "~/hmf/analysis/probes/probeResultCondensed.RData")
