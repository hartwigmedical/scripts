library(RMySQL)

library(IRanges)
detach("package:purple", unload=TRUE); 
library(purple);
library(Biostrings)
library(GenomicRanges)
library(MASS)
library(seqinr)
library(dplyr)
detach("package:dndscv", unload=TRUE); 
library(dndscv)

annotate_somatics <- function(annotmuts, somatics) {
  result = left_join(annotmuts, somatics, by = c("sampleID","chr","pos","ref", "mut")) %>% 
    select(sampleId = sampleID, chromosome = chr, position = pos, ref = ref, alt = mut, gene, type, impact, worstCodingEffect, canonicalCodingEffect, hotspot, biallelic)
  
  result = result %>%
    mutate(impact = ifelse(impact == "Essential_Splice", "Splice", impact)) %>%
    mutate(impact = ifelse(impact == "no-SNV", paste0(substr(canonicalCodingEffect, 1, 1), tolower(substring(canonicalCodingEffect, 2))), impact)) %>%
    mutate(impact = ifelse(impact == "None", paste0(substr(worstCodingEffect, 1, 1), tolower(substring(worstCodingEffect, 2))), impact)) %>%
    mutate(impact = ifelse(impact %in% c("Stop_loss","Nonsense_or_frameshift"), "Nonsense", impact)) %>%
    mutate(impact = ifelse(type == "INDEL" & impact == "Missense", "Inframe", impact)) %>%
    mutate(impact = ifelse(type == "INDEL" & impact == "Nonsense", "Frameshift", impact)) %>%
    filter(!is.na(impact), impact != "None")
 
    return (result) 
}

data("cancergenes_cgc81", package="dndscv")
data("covariates_hg19", package="dndscv")
refdb="~/hmf/RData/HmfRefCDS.RData"
load(refdb)
newgenes = sapply(RefCDS, function(x) x$gene_name) # New list of genes
oldgenes = sapply(strsplit(newgenes,split=":"), function(x) x[1])
cv = covs[oldgenes,]
rownames(cv) = newgenes
kc = newgenes[oldgenes %in% known_cancergenes]

load(file = "~/hmf/RData/highestPurityExonicSomatics.RData")

somaticEnrichment = highestPurityExonicSomatics %>% 
  select(sampleID = sampleId, chr = chromosome, pos = position, ref = ref, mut = alt, type, worstCodingEffect, canonicalCodingEffect, hotspot, biallelic)

######## FILTERED
somatics = highestPurityExonicSomatics %>%  
  filter(type == "INDEL" | type == "SNP") %>%
  select(sampleId, chr = chromosome, pos = position, ref = ref, alt = alt)

output = dndscv(somatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)

highestPurityExonicSomaticsFilteredImpact = annotate_somatics(output$annotmuts, somaticEnrichment)
save(highestPurityExonicSomaticsFilteredImpact, file = "~/hmf/RData/highestPurityExonicSomaticsFilteredImpact.RData")

########### UNFILTERED
somatics = highestPurityExonicSomatics %>%  
  select(sampleId, chr = chromosome, pos = position, ref = ref, alt = alt)

output = dndscv(somatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE, max_muts_per_gene_per_sample = 30000000, max_coding_muts_per_sample = 30000000)

highestPurityExonicSomaticsImpact = annotate_somatics(output$annotmuts, somaticEnrichment)
save(highestPurityExonicSomaticsImpact, file = "~/hmf/RData/highestPurityExonicSomaticsImpact.RData")
