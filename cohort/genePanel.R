library(RMySQL)
library(dplyr)
library(tidyr)
library(GenomicRanges)

########################### Cosmic Genes
cosmicCurated = read.csv("~/hmf/resources/CosmicCurated.csv", stringsAsFactors = F)
cosmicCurated$cosmicCurated <- TRUE
cosmicCurated = cosmicCurated[, c("Genes","cosmicCurated")]
colnames(cosmicCurated) <- c("gene", "cosmicCurated")

cosmicCensus = read.csv("~/hmf/resources/CosmicCensus.csv", stringsAsFactors = F)
colnames(cosmicCensus) <- c("gene", "cosmicOncogene", "cosmicTsg")
cosmicCensus = cosmicCensus %>% filter(cosmicOncogene | cosmicTsg)

cosmicGenes = merge(cosmicCurated, cosmicCensus, by = "gene", all = T)
rm(cosmicCurated, cosmicCensus)

########################### Actionable Genes
actionableGenes = read.table("~/hmf/analysis/actionable/actionablePanel.tsv", header = T, stringsAsFactors = F) 
names(actionableGenes)
colnames(actionableGenes) <- c("gene", "actionableAmplification", "actionableDeletion", "actionableFusion", "actionableVariant", "actionableDrup", "actionableDrupCategory", "actionableResponse", "actionableResponseSource", "actionableResistance", "actionableResistanceSource")
actionableGenes$actionableAmplification <- ifelse(actionableGenes$actionableAmplification  == "true", T, NA)
actionableGenes$actionableDeletion <- ifelse(actionableGenes$actionableDeletion  == "true", T, NA)
actionableGenes$actionableFusion <- ifelse(actionableGenes$actionableFusion  == "true", T, NA)
actionableGenes$actionableVariant <- ifelse(actionableGenes$actionableVariant  == "true", T, NA)
actionableGenes$actionableDrup <- ifelse(actionableGenes$actionableDrup  == "true", T, NA)


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
colnames(knownAmpsDels) <- c("gene", "knownAmplification", "knownDeletion")


########################### Gene Panel
load(file="~/hmf/RData/reference/PcawgRefCDSCv.RData")
load(file="~/hmf/analysis/cohort/processed/HmfRefCDSCv.RData")
sig = 0.01
hmfSignificant =  HmfRefCDSCv %>% filter(qglobal_cv < sig) %>% distinct(gene_name) %>% select(gene = gene_name)
hmfSignificant$hmfDnds <- TRUE
martincorenaSignificant =  PcawgRefCDSCv %>% filter(qglobal < sig) %>% distinct(gene_name) %>% select(gene = gene_name)
martincorenaSignificant$martincorenaDnds <- TRUE

genePanel = merge(martincorenaSignificant, hmfSignificant, by = "gene", all = T)
genePanel = merge(genePanel, cosmicGenes, by = "gene", all = T)
genePanel = merge(genePanel, actionableGenes, by = "gene", all = T)
genePanel = merge(genePanel, knownAmpsDels, by = "gene", all=T)

genePanel = genePanel %>% filter(!gene %in% c("POM121L12","TRIM49B","LPCAT2"))
save(genePanel, file="~/hmf/analysis/cohort/processed/genePanel.RData")
str(genePanel)


########################### PART 2 WITH TARGET AMPS AND DELS
load(file="~/hmf/analysis/cohort/processed/genePanel.RData")
load(file = "~/hmf/analysis/cohort/processed/geneCopyNumberDeleteTargets.RData")
load(file = "~/hmf/analysis/cohort/processed/geneCopyNumberAmplificationTargets.RData")

completeGenePanel = merge(genePanel, geneCopyNumberAmplificationTargets %>% ungroup() %>% mutate(gene = target, hmfAmplification = T) %>% select(gene, hmfAmplification), by = "gene", all = T)
completeGenePanel = merge(completeGenePanel, geneCopyNumberDeleteTargets %>% ungroup() %>% mutate(gene = target, hmfDeletion = T) %>% select(gene, centromere, telomere, hmfDeletion) %>% distinct(), by = "gene", all = T)
save(completeGenePanel, file="~/hmf/analysis/cohort/processed/completeGenePanel.RData")

