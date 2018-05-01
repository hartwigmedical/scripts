#!/usr/bin/env Rscript
library(purple)
library(RMySQL)
library(data.table)
library(IRanges)
library(grid)
library(MutationalPatterns)
library(ggplot2)
library(dplyr)

ID = commandArgs(trailingOnly = TRUE)
dbProd = dbConnect(MySQL(), user = db_user, password = db_password, dbname = db_name, groups = "RAnalysis")

## Get cohort
multipleBiopsyCohort = purple::query_multiple_biopsy_cohort(dbProd)

##### USE FOLLOWING LINE TO FILTER A PARTICULAR PATIENT ID
multipleBiopsyCohort = multipleBiopsyCohort %>% filter(patientId == ID)

## Get somatic and structural variants
multipleBiopsySomaticVariants = purple::query_somatic_variants(dbProd, multipleBiopsyCohort)
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

