library(VariantAnnotation)
library(tidyr)
library(dplyr)


vcf_data_frame <- function(sample, vcf) {
  vcf.info = info(vcf)
  vcf.fixed = fixed(vcf)
  vcf.alt = CharacterList(alt(vcf))
  vcf.alt[elementNROWS(vcf.alt) > 1L ] <- lapply(vcf.alt[elementNROWS(vcf.alt) > 1L ], paste0, collapse=",")
  
  vcf.df = data.frame(
      sample = sample,
      chromosome = seqnames(vcf),
      pos = start(vcf),
      ref = as.character(ref(vcf)), 
      alt = as.character(vcf.alt),
      filter = vcf.fixed$FILTER,
      repeatCount = vcf.info$REP_C,
      biallelic = vcf.info$BIALLELIC, stringsAsFactors = F) %>%
    mutate(repeatCount = ifelse(is.na(repeatCount), 0, repeatCount))
    
    return (vcf.df)
}

vcf = readVcf("/Users/jon/hmf/analysis/COLO829T/purple/COLO829T.purple.somatic.vcf.gz")
vcf.df = vcf_data_frame("COLO829T", vcf)


