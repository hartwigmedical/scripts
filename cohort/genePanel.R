library(RMySQL)
library(dplyr)
library(tidyr)
library(GenomicRanges)


########################### Cosmic Genes
cosmicCurated = read.csv("~/hmf/resources/CosmicCurated.csv", stringsAsFactors = F)
cosmicCurated$cosmicCurated <- TRUE
cosmicCurated = cosmicCurated[, c("Genes","cosmicCurated")]
colnames(cosmicCurated) <- c("gene_name", "cosmicCurated")

cosmicCensus = read.csv("~/hmf/resources/CosmicCensus.csv", stringsAsFactors = F)
colnames(cosmicCensus) <- c("gene_name", "cosmicOncogene", "cosmicTsg")
cosmicCensus = cosmicCensus %>% filter(cosmicOncogene | cosmicTsg)

cosmicGenes = merge(cosmicCurated, cosmicCensus, by = "gene_name", all = T)
rm(cosmicCurated, cosmicCensus)

########################### KNOWN AMPS AND DELS
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

manual = data.frame(gene = "ZNF703", alteration = "amplification")
knownAmpsDels = bind_rows(cgi, onco) %>%
    bind_rows(civic) %>%
    bind_rows(manual) %>%
    select(gene_name = gene, alteration) %>%
    distinct(gene_name, alteration) %>%
    mutate(value = T) %>% spread(alteration, value)
rm(cgi, onco, civic, manual)


########################### Gene Panel
load(file="~/hmf/RData/reference/PcawgRefCDSCv.RData")
load(file="~/hmf/analysis/genePanel/HmfRefCDSCv.RData")
sig = 0.01
hmfSignificant =  HmfRefCDSCv %>% filter(qglobal_cv < sig) %>% distinct(gene_name)
hmfSignificant$hmf <- TRUE
martincorenaSignificant =  PcawgRefCDSCv %>% filter(qglobal < sig) %>% distinct(gene_name)
martincorenaSignificant$martincorena <- TRUE

genePanel = merge(martincorenaSignificant, hmfSignificant, by = "gene_name", all = T)
genePanel = merge(genePanel, cosmicGenes, by = "gene_name", all = T)
genePanel = merge(genePanel, knownAmpsDels, by = "gene_name", all=T)
genePanel = genePanel %>% filter(!gene_name %in% c("POM121L12","TRIM49B","LPCAT2"))
save(genePanel, file="~/hmf/analysis/genePanel/genePanel.RData")

#genePanel = merge(genePanel, fragileGenes, by = "gene_name", all=T)
jon = PcawgRefCDSCv %>% filter(gene_name == "RACGAP1")

