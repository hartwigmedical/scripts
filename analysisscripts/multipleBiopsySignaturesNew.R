#!/usr/bin/Rscript
library(RMySQL)
library(IRanges)
library(MutationalPatterns)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)

## Script Arguments
sample1 = 'XXXT'
sample2 = 'XXXTII'
samples = c(sample1, sample2)

## Define Constants
singleBlue = "#6baed6"
singleRed = "#d94701"
scopeColours = setNames(c(singleRed, singleBlue), c("Private", "Shared"))
cosmicSignatureColours =
c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
"#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a",
"#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
"#dea185","#a0729d","#8a392f")
redGradient = c("#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15")
greenGradient = c("#edf8e9", "#bae4b3", "#74c476", "#31a354", "#006d2c")
blackGradient = c("#f7f7f7", "#cccccc", "#969696", "#636363", "#252525")
indelSignatureColours = c(rev(redGradient), greenGradient)

cosmicSignatures = read.table("https://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt", sep = "\t", header = T)
cosmicSignatures = cosmicSignatures[order(cosmicSignatures[, 1]),]
cosmicSignatures = as.matrix(cosmicSignatures[, 4:33])

## Signature Functions
snp_signature <- function(somaticVariants) {
  standard_mutation <- function(types) {
    types = gsub("G>T", "C>A", types)
    types = gsub("G>C", "C>G", types)
    types = gsub("G>A", "C>T", types)
    types = gsub("A>T", "T>A", types)
    types = gsub("A>G", "T>C", types)
    types = gsub("A>C", "T>G", types)
    return(types)
  }

  standard_context <- function(raw_type, standard_type, context) {
    x = which(raw_type != standard_type)
    context[x] = reverse(chartr("ATGC", "TACG", context[x]))
    return(context)
  }

  create_empty_mutational_signature <- function() {
    cRow = paste0(c(rep("A", 4), rep("C", 4), rep("G", 4), rep("T", 4)), "C", c("A", "C", "G", "T"))
    tRow = paste0(c(rep("A", 4), rep("C", 4), rep("G", 4), rep("T", 4)), "T", c("A", "C", "G", "T"))
    types = c(rep("C>A", 16), rep("C>G", 16), rep("C>T", 16), rep("T>A", 16), rep("T>C", 16), rep("T>G", 16))
    contexts = c(rep(cRow, 3), rep(tRow, 3))
    result = data.frame(type = types, context = contexts, stringsAsFactors = F)
    return (result)
  }

  snps = somaticVariants %>%
    filter(type == 'SNP', !grepl('N', trinucleotideContext), nchar(ref) == 1, nchar(alt) == 1) %>%
    mutate(
    raw_types = paste(ref, alt, sep = ">"),
    type = standard_mutation(raw_types),
    context = standard_context(raw_types, type, trinucleotideContext)) %>%
    select(sampleId, type, context, scope)

  snpSignature = snps %>%
  #top_n(50, type) %>%
    group_by(type, context, scope) %>% count() %>%
    spread(scope, n, fill = 0) %>%
    full_join(create_empty_mutational_signature(), by = c("type", "context")) %>%
    arrange(type, context)


  for (sample in unique(somaticVariants$sampleId)) {
    if (!sample %in% names(snpSignature)) {
      snpSignature[, sample] <- 0
    }
  }

  if (!"Shared" %in% names(snpSignature)) {
    snpSignature$Shared <- 0
  }

  snpSignature[is.na(snpSignature)] <- 0
  snpSignature = snpSignature %>% ungroup() %>% select(type, context, Shared, everything())

  return (snpSignature)
}

indel_signature <- function(somaticVariants) {

  levelFactor = c(-5,-4,-3,-2,-1,1,2,3,4,5)

  indels = somaticVariants %>% filter(type == 'INDEL') %>%
    mutate(
    length = pmax(-5, pmin(5, nchar(alt) - nchar(ref))),
    length = factor(length, levelFactor, ordered = T)
    )

  empty_indel_signature = data.frame(length = factor(levelFactor, levelFactor, ordered = T))

  indelSignature = indels %>%
    group_by(length, scope) %>% count() %>%
    spread(scope, n, fill = 0) %>%
    full_join(empty_indel_signature, by = "length") %>%
    arrange(length)

  indelSignature[is.na(indelSignature)] <- 0

  for (sample in unique(somaticVariants$sampleId)) {
    if (!sample %in% names(indelSignature)) {
      indelSignature[, sample] <- 0
    }
  }

  if (!"Shared" %in% names(indelSignature)) {
    indelSignature$Shared <- 0
  }

  return (indelSignature)
}


## Plot Functions
theme_set(theme_bw() + theme(
axis.text=element_text(size=6),
axis.title=element_text(size=6),
plot.title = element_text(size = 8),
legend.key.size = unit(0.5, "cm"),
legend.text = element_text(size=5)))


plot_ploidy_histogram <- function(df, sample, maxPloidy = 3.5) {
  df = df %>% filter(sampleId == sample) %>% mutate(scope = ifelse(scope == "Shared", "Shared", "Private"))
  ggplot(df , aes(ploidy)) +
    geom_histogram(binwidth = 0.05, aes(fill = scope), position = "identity", alpha = 0.8) +
    xlim(0,maxPloidy) +
    ylab("Count") + xlab("Ploidy") +
    ggtitle(sample) +
    scale_fill_manual(values = scopeColours) +
    theme(legend.position="bottom", legend.title = element_blank())
}

plot_snp_signature <- function(snpSignature) {
  contribution = fit_to_signatures(as.matrix(snpSignature[, -c(1, 2)]), cosmicSignatures)$contribution
  p1 <- plot_contribution(contribution, cosmicSignatures, mode = "absolute") +
    theme(axis.text.x = element_text(size=6))+
    scale_fill_manual(values= cosmicSignatureColours, labels = paste0("Sig ", c(1:30)))+labs(fill="")+
    ggtitle("Mutational Signatures") +
    labs(x = "", y = "Absolute contribution") +
    theme(legend.key.size = unit(0.5, "cm"),  legend.text = element_text(size=5), legend.box.background = element_blank(),
    axis.text=element_text(size=6), axis.title=element_text(size=6), plot.title = element_text(size = 8))
  return (p1)
}

plot_indel_signature <- function(indelSignature) {
  plot_absolute_contribution(indelSignature)+
    ggtitle("Indel Signatures")+
    scale_fill_manual(name = "", values = indelSignatureColours)
}

plot_absolute_contribution <- function(contribution) {
  m_contribution = contribution %>% gather(variable, value, -1) %>% mutate(variable = factor(variable, levels = names(contribution)[-1]))
  colnames(m_contribution) = c("Signature", "Sample", "Contribution")
  plot = ggplot(m_contribution, aes(x = factor(Sample),  y = Contribution, fill = factor(Signature), order = Sample)) +
    geom_bar(stat = "identity", colour = "black") +
    xlim(rev(levels(factor(m_contribution$Sample)))) +
    labs(x = "", y = "Absolute contribution") +
    scale_fill_discrete(name = "") +
    theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) +
    theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
}


####################################### EXECUTION

## Collect data from MySQL
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis", host = "127.0.0.1")
sampleIdString = paste("'", samples, "'", collapse = ",", sep = "")
somaticVariants = dbGetQuery(dbProd, paste0("SELECT * FROM somaticVariant WHERE filter = 'PASS' AND sampleId in (",sampleIdString,")"))
structuralVariants = dbGetQuery(dbProd, paste0("SELECT * FROM structuralVariant WHERE filter = 'PASS' AND sampleId in (",sampleIdString,")"))
dbDisconnect(dbProd)
rm(dbProd, sampleIdString)

##### Add Scope
somaticVariants = somaticVariants %>% group_by(chromosome, position, ref, alt) %>% mutate(scope = ifelse(n() == 1, sampleId, "Shared"), scope = factor(scope, rev(c("Shared", samples)), ordered = T)) %>% ungroup()
structuralVariants = structuralVariants %>% group_by(startChromosome, endChromosome, startPosition, endPosition, startOrientation, endOrientation) %>% mutate(scope = ifelse(n() == 1, "Private", "Shared")) %>% ungroup()


#### Somatic Signature
snpSignature = snp_signature(somaticVariants)
indelSignature = indel_signature(somaticVariants)

p1 = plot_ploidy_histogram(somaticVariants, sample1) + xlab("Somatic Variant Ploidy")
p2 = plot_ploidy_histogram(structuralVariants, sample1) + xlab("Structural Variant Ploidy")
p3 = plot_indel_signature(indelSignature)

p4 = plot_ploidy_histogram(somaticVariants, sample2) + xlab("Somatic Variant Ploidy") + theme(legend.position = "none")
p5 = plot_ploidy_histogram(structuralVariants, sample2) + xlab("Structural Variant Ploidy") + theme(legend.position = "none")
p6 = plot_snp_signature(snpSignature)

pCombined = cowplot::plot_grid(p1, p2, p3, p4, p5, p6)
pCombined




mutationalSignature = snpSignature


ID = commandArgs(trailingOnly = TRUE)

## Get cohort
multipleBiopsyCohort = purple::query_multiple_biopsy_cohort(dbProd, data.frame(sampleId = character(), genesDeleted = integer()))

##### USE FOLLOWING LINE TO FILTER A PARTICULAR PATIENT ID
multipleBiopsyCohort = multipleBiopsyCohort %>% filter(patientId == ID)

## Get somatic and structural variants
multipleBiopsySomaticVariants = purple::query_somatic_variants(dbProd, multipleBiopsyCohort)
multipleBiopsySomaticVariants$ploidy = multipleBiopsySomaticVariants$copyNumber * multipleBiopsySomaticVariants$adjustedVaf
multipleBiopsyStructuralVariants = purple::query_structural_variants(dbProd, multipleBiopsyCohort)

### Clean up
dbDisconnect(dbProd)

### Mutational Signatures
multipleBiopsyMutationalSignature = list()
multipleBiopsyIndelSignature = list()
multipleBiopsySVSignature = list()

multipleBiopsyPatientsId = unique(multipleBiopsyCohort$patientId)
for (patientId in multipleBiopsyPatientsId) {
  cat("Processing", patientId, "\n")

  # Extract patient data
  patientSampleIds = unique(multipleBiopsyCohort[multipleBiopsyCohort$patientId == patientId, c("sampleId")])
  patientSomaticVariants = scope(multipleBiopsySomaticVariants[multipleBiopsySomaticVariants$sample %in% patientSampleIds,])
  patientStructuralVariants = scope(multipleBiopsyStructuralVariants[multipleBiopsyStructuralVariants$sample %in% patientSampleIds,])

  # Signatures
  patientMutationalSignature = purple::mutational_signature_by_scope(patientSomaticVariants)
  patientIndelSignature = purple::indel_signature_by_scope(patientSomaticVariants)
  patientSVSignature = purple::sv_signature_by_scope(patientStructuralVariants)

  # Save Signatures
  multipleBiopsyMutationalSignature[[patientId]] <- patientMutationalSignature
  multipleBiopsyIndelSignature[[patientId]] <- patientIndelSignature
  multipleBiopsySVSignature[[patientId]] <- patientSVSignature

  # Plots
  somaticPloidyPlots = somatic_ploidy_plots(patientSomaticVariants)
  structuralPloidyPlots = structural_ploidy_plots(patientStructuralVariants)
  p7 <- plot_indel_signature(patientIndelSignature)
  p8 <- plot_sv_signature(patientSVSignature)
  p9 <- plot_cosmic_signature(patientMutationalSignature)

  pdf(file = paste(patientId, "MultipleBiopsies.pdf", sep = "_"), height = 14, width = 20)
  multiplot(somaticPloidyPlots[[1]], somaticPloidyPlots[[2]], somaticPloidyPlots[[3]],
  structuralPloidyPlots[[1]], structuralPloidyPlots[[2]], structuralPloidyPlots[[3]],
  p7, p8, p9,
  cols = 3)
  dev.off()
}

#save(multipleBiopsyCohort, file="~/hmf/RData/multipleBiopsyCohort.RData")
#save(multipleBiopsyStructuralVariants, file="~/hmf/RData/multipleBiopsyStructuralVariants.RData")
#save(multipleBiopsySomaticVariants, file="~/hmf/RData/multipleBiopsySomaticVariants.RData")
#save(multipleBiopsyMutationalSignature, multipleBiopsyIndelSignature, multipleBiopsySVSignature, file="~/hmf/multipleBiopsySignatures.RData")


#### LOAD INPUTS
#load(file="~/hmf/multipleBiopysyCohort.RData")
#load(file="~/hmf/multipleBiopsySomaticVariants.RData")
#load(file="~/hmf/multipleBiopsyStructuralVariants.RData")

