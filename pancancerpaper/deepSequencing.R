#library(BiocManager)
#install("VariantAnnotation")
#install("dplyr")
#install("ggplot2")

library(GenomicRanges)
library(VariantAnnotation)
library(dplyr)
library(cowplot)
theme_set(theme_bw())

vcf_data_frame<- function(vcf) {
  vcf.rowRanges = rowRanges(vcf)
  vcf.info = info(vcf)
  vcf.alt = CharacterList(alt(vcf))
  vcf.alt[elementNROWS(vcf.alt) > 1L ] <- lapply(vcf.alt[elementNROWS(vcf.alt) > 1L ], paste0, collapse=",")
  
  vcf.df = data.frame(
    chr = seqnames(vcf), 
    pos = start(vcf), 
    ref = ref(vcf), 
    alt = as.character(vcf.alt),  
    filter = as.character(vcf.rowRanges$FILTER)
  )
  
  return (vcf.df)
}

sv_data_frame<- function(vcf) {
  vcf.rowRanges = rowRanges(vcf)
  vcf.info = info(vcf)
  vcf.alt = CharacterList(alt(vcf))
  vcf.alt[elementNROWS(vcf.alt) > 1L ] <- lapply(vcf.alt[elementNROWS(vcf.alt) > 1L ], paste0, collapse=",")
  
  vcf.df = data.frame(
    id = names(vcf.rowRanges),
    chr = seqnames(vcf), 
    pos = start(vcf), 
    ref = ref(vcf), 
    alt = as.character(vcf.alt),  
    filter = as.character(vcf.rowRanges$FILTER),
    pairId = vcf.info$PARID
  )
  
  return (vcf.df)
}


somatic_dataframe <- function(sampleId, category, vcf) {
  result  = vcf_data_frame(vcf) 
  result$type <- "INDEL"
  result$sampleId <- sampleId
  result$category <- category
  
  result = result %>% mutate(
    alt = as.character(alt), 
    type = ifelse(nchar(ref) == nchar(alt) & nchar(ref) == 1, "SNV", type),
    type = ifelse(nchar(ref) == nchar(alt) & nchar(ref) != 1, "MNV", type))
  
  return (result %>% filter(filter == 'PASS'))
}

sv_dataframe <- function(sampleId, category, vcf) {
  result  = sv_data_frame(vcf) 
  result$type <- "SV"
  result$sampleId <- sampleId
  result$category <- category

  vcf.id = result %>% dplyr::select(id) %>% mutate(IDRank = row_number())
  vcf.mateid = result %>% dplyr::select(id = pairId) %>% mutate(MATEIDRank = row_number())  
  
  result = left_join(result, vcf.id, by = "id") %>% 
    left_join(vcf.mateid, by = "id") %>% 
    mutate(primary = IDRank < MATEIDRank) %>% 
    filter(primary | is.na(primary), filter == 'PASS') %>%
    dplyr::select(-IDRank, -MATEIDRank) 
  
  
  return (result)
}


purity_dataframe <- function(sample, category, file) {
  df = read.table(file = file, sep = "\t", header = T, comment.char = "?")
  result = data.frame(sampleId = sample, category = category, type = "Purity", value = df$X.Purity)
  result = bind_rows(result, data.frame(sampleId = sample, category = category, type = "Ploidy", value = df$Ploidy))
}

sample_data <- function(sample) {
  
  combined = somatic_dataframe(sample, "Combined", readVcf(file = paste0("~/hmf/analysis/deep/", sample, "T/", sample, "C.purple.somatic.vcf.gz"), 'hg19'))
  rep1 = somatic_dataframe(sample, "Replicate1", readVcf(file = paste0("~/hmf/analysis/deep/", sample, "T/", sample, "A.purple.somatic.vcf.gz"), 'hg19'))
  rep2 = somatic_dataframe(sample, "Replicate2", readVcf(file = paste0("~/hmf/analysis/deep/", sample, "T/", sample, "B.purple.somatic.vcf.gz"), 'hg19'))
  somatics = bind_rows(bind_rows(combined, rep1), rep2) %>% group_by(sampleId, category, type) %>% dplyr::summarise(value = n())
  
  combined_sv = sv_dataframe(sample, "Combined", readVcf(file = paste0("~/hmf/analysis/deep/", sample, "T/", sample, "C.purple.sv.vcf.gz"), 'hg19'))
  rep1_sv = sv_dataframe(sample, "Replicate1", readVcf(file = paste0("~/hmf/analysis/deep/", sample, "T/", sample, "A.purple.sv.vcf.gz"), 'hg19'))
  rep2_sv = sv_dataframe(sample, "Replicate2", readVcf(file = paste0("~/hmf/analysis/deep/", sample, "T/", sample, "B.purple.sv.vcf.gz"), 'hg19'))
  svs = bind_rows(bind_rows(combined_sv, rep1_sv), rep2_sv) %>% group_by(sampleId, category, type) %>% dplyr::summarise(value = n())

  combined_purity = purity_dataframe(sample, "Combined", paste0("~/hmf/analysis/deep/", sample, "T/", sample, "C.purple.purity"))
  rep1_purity = purity_dataframe(sample, "Replicate1", paste0("~/hmf/analysis/deep/", sample, "T/", sample, "A.purple.purity"))
  rep2_purity = purity_dataframe(sample, "Replicate2", paste0("~/hmf/analysis/deep/", sample, "T/", sample, "B.purple.purity"))
  purity = bind_rows(bind_rows(combined_purity, rep1_purity), rep2_purity) 
  
  return (bind_rows(bind_rows(somatics, svs), purity))
}

sample1 = sample_data("CPCT02050339")
sample2 = sample_data("CPCT02090044")
deepSequencingSummary = bind_rows(sample1, sample2)
save(deepSequencingSummary, file = "~/hmf/RData/Reference/deepSequencingSummary.RData")








