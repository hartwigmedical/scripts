library(dplyr)
library(tidyr)
library(GenomicRanges)

cgi = read.table(file = "~/hmf/resources/cgi_biomarkers_per_variant.tsv", stringsAsFactors = F, header = T, sep = "\t") %>%
  select(Biomarker) %>% 
  separate(Biomarker, c("gene","alteration"), sep = " ") %>%
  filter(alteration %in% c("amplification", "deletion")) %>%
  distinct(gene, alteration)

onco = read.table(file = "~/hmf/resources/oncoKb.tsv", stringsAsFactors = F, header = T, sep = "\t") %>% 
  select(gene = Gene, alteration = Alteration) %>%
  mutate(alteration = tolower(alteration)) %>%
  filter(alteration %in% c("amplification", "deletion")) %>%
  distinct(gene, alteration)

civic = read.table(file = "~/hmf/resources/civic_variants.tsv", stringsAsFactors = T, header = T, sep = "\t", quote = "") %>%
  select(gene, alteration = variant) %>%
  mutate(alteration = tolower(alteration)) %>%
  filter(alteration %in% c("amplification", "deletion")) %>%
  distinct(gene, alteration)

knownAmpsDels = rbind(rbind(cgi, onco), civic) %>% 
  select(gene_name = gene, alteration)  %>%
  distinct(gene_name = gene, alteration) %>% 
  mutate(value = T) %>% spread(alteration, value) 

knownAmpsDels[is.na(knownAmpsDels)] <- F
rm(cgi, onco, civic)

load("~/hmf/RData/fragileGenes.RData")
load("~/hmf/RData/genePanel.RData")
genePanel = merge(genePanel, fragileGenes, by = "gene_name", all=T)
genePanel = merge(genePanel, knownAmpsDels, by = "gene_name", all=T)
rm(fragileGenes)


save(genePanel, file = "~/hmf/RData/ampsDelsGenePanel.RData")
