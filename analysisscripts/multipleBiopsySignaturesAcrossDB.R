#!/usr/bin/Rscript
library(purple)
library(RMySQL)
library(data.table)
library(IRanges)
library(grid)
library(MutationalPatterns)
library(ggplot2)
library(dplyr)

## To generate a CPCT-DRUP comparison run "multipleBiopsySignatures.R <CPCT patient ID>"
ID = commandArgs(trailingOnly = TRUE)
ID = "CPCT02180045T"

multipleBiopsyCohort = data.frame(sampleId = ID, patientId = ID)

##### PROD DATA
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis", host = "127.0.0.1")
multipleBiopsySomaticVariantsProd = purple::query_somatic_variants(dbProd, multipleBiopsyCohort) %>% mutate(source = "Prod")
multipleBiopsyStructuralVariantsProd = purple::query_structural_variants(dbProd, multipleBiopsyCohort)  %>% mutate(source = "Prod")
dbDisconnect(dbProd)
rm(dbProd)

##### PILOT DATA
dbPilot = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis", host = "127.0.0.1")
multipleBiopsySomaticVariantsPilot = purple::query_somatic_variants(dbPilot, multipleBiopsyCohort) %>% mutate(source = "Pilot")
multipleBiopsyStructuralVariantsPilot = purple::query_structural_variants(dbPilot, multipleBiopsyCohort) %>% mutate(source = "Pilot")
dbDisconnect(dbPilot)

##### COMBINE
multipleBiopsySomaticVariants = bind_rows(multipleBiopsySomaticVariantsProd, multipleBiopsySomaticVariantsPilot) %>% mutate(sampleId = paste0(sampleId, "-", source), patientId = sampleId)
multipleBiopsyStructuralVariants = bind_rows(multipleBiopsyStructuralVariantsProd, multipleBiopsyStructuralVariantsPilot) %>% mutate(sampleId = paste0(sampleId, "-", source), patientId = sampleId)

### Mutational Signatures
multipleBiopsyMutationalSignature = list()
multipleBiopsyIndelSignature = list()
multipleBiopsySVSignature = list()

multipleBiopsyPatientsId = unique(multipleBiopsyCohort$patientId)
for (patientId in multipleBiopsyPatientsId) {
  cat("Processing", patientId, "\n")
  
  # Extract patient data
  patientSampleIds = unique(multipleBiopsySomaticVariants[multipleBiopsyCohort$patientId == patientId, c("sampleId")])
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
  p9 = p8 # Note if samples are identical the following line won't work so we dummy p9 with p8
  p9 <- plot_cosmic_signature(patientMutationalSignature)

    #pdf(file = paste(patientId, "MultipleBiopsies.pdf", sep = "_"), height = 14, width = 20)
  pdf(file = "~/JonMultipleBiopsies.pdf", height = 14, width = 20)
  multiplot(somaticPloidyPlots[[1]], somaticPloidyPlots[[2]], somaticPloidyPlots[[3]],
            structuralPloidyPlots[[1]], structuralPloidyPlots[[2]], structuralPloidyPlots[[3]], 
            p7, p8, p9,
            cols = 3)
  dev.off()
}


