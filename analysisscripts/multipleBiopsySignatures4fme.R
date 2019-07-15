#!/usr/bin/Rscript
library(purple)
library(RMySQL)
library(data.table)
library(IRanges)
library(grid)
library(MutationalPatterns)
library(ggplot2)
library(dplyr)

##### ##### ##### ##### ##### ##### ##### ##### GATHER DATA

##### PILOT DATA
dbPilot = dbConnect(MySQL(), dbname='4fme', groups="RAnalysisWrite", host = "127.0.0.1")
multipleBiopsyCohort = dbGetQuery(dbPilot, "SELECT sampleId, sampleId as patientId from purity")
multipleBiopsySomaticVariantsPilot = purple::query_somatic_variants(dbPilot, multipleBiopsyCohort) %>% mutate(source = "4fme")
multipleBiopsyStructuralVariantsPilot = purple::query_structural_variants(dbPilot, multipleBiopsyCohort) %>% mutate(source = "4fme")
dbDisconnect(dbPilot)
rm(dbPilot)

##### PROD DATA
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis", host = "127.0.0.1")
multipleBiopsySomaticVariantsProd = purple::query_somatic_variants(dbProd, multipleBiopsyCohort) %>% mutate(source = "Prod")
multipleBiopsyStructuralVariantsProd = purple::query_structural_variants(dbProd, multipleBiopsyCohort)  %>% mutate(source = "Prod")
dbDisconnect(dbProd)
rm(dbProd)

##### COMBINE
multipleBiopsySomaticVariants = bind_rows(multipleBiopsySomaticVariantsProd, multipleBiopsySomaticVariantsPilot) %>% mutate(patientId = sampleId, sampleId = paste0(sampleId, "-", source))
multipleBiopsyStructuralVariants = bind_rows(multipleBiopsyStructuralVariantsProd, multipleBiopsyStructuralVariantsPilot) %>% mutate(patientId = sampleId, sampleId = paste0(sampleId, "-", source))
save(multipleBiopsyCohort, multipleBiopsySomaticVariants, multipleBiopsyStructuralVariants, file = "~/hmf/analysis/4fme/data.RData")

##### ##### ##### ##### ##### ##### ##### ##### RUN FROM HERE KORNEEL

load( file = "~/hmf/analysis/4fme/data.RData")
multipleBiopsyMutationalSignature = list()
multipleBiopsyIndelSignature = list()
multipleBiopsySVSignature = list()

multipleBiopsyPatientsId = multipleBiopsyCohort %>% filter(sampleId != 'WIDE01010005T') %>% pull(sampleId)
for (selectedPatientId in multipleBiopsyPatientsId) {
  cat("Processing", selectedPatientId, "\n")
  
  # Extract patient data
  patientSampleIds = multipleBiopsySomaticVariants %>% filter(patientId == selectedPatientId) %>% distinct(sampleId) %>% pull()
  patientSomaticVariants = scope(multipleBiopsySomaticVariants[multipleBiopsySomaticVariants$sample %in% patientSampleIds,])
  patientStructuralVariants = scope(multipleBiopsyStructuralVariants[multipleBiopsyStructuralVariants$sample %in% patientSampleIds,])
  
  # Signatures
  patientMutationalSignature = purple::mutational_signature_by_scope(patientSomaticVariants)
  patientIndelSignature = purple::indel_signature_by_scope(patientSomaticVariants)
  patientSVSignature = purple::sv_signature_by_scope(patientStructuralVariants)
  
  # Save Signatures
  multipleBiopsyMutationalSignature[[selectedPatientId]] <- patientMutationalSignature
  multipleBiopsyIndelSignature[[selectedPatientId]] <- patientIndelSignature
  multipleBiopsySVSignature[[selectedPatientId]] <- patientSVSignature  

  # Plots
  somaticPloidyPlots = somatic_ploidy_plots(patientSomaticVariants)
  structuralPloidyPlots = structural_ploidy_plots(patientStructuralVariants)
  p7 <- plot_indel_signature(patientIndelSignature)
  p8 <- plot_sv_signature(patientSVSignature)
  p9 = p8 # Note if samples are identical the following line won't work so we dummy p9 with p8
  p9 <- plot_cosmic_signature(patientMutationalSignature)

  pdf(file = paste0("~/hmf/analysis/4fme/", selectedPatientId, "_MultipleBiopsies.pdf"), height = 14, width = 20)
  #pdf(file = "~/TestMultipleBiopsies.pdf", height = 14, width = 20)
  multiplot(somaticPloidyPlots[[1]] + ggtitle(patientSampleIds[1]), 
            somaticPloidyPlots[[2]] + ggtitle(patientSampleIds[2]),
            somaticPloidyPlots[[3]],
            structuralPloidyPlots[[1]]  + ggtitle(patientSampleIds[1]), 
            structuralPloidyPlots[[2]] + ggtitle(patientSampleIds[2]), 
            structuralPloidyPlots[[3]], 
            p7, p8, p9,
            cols = 3)
  dev.off()
}


