detach("package:purple", unload=TRUE); library(purple);
library(RMySQL)
library(data.table)
library(IRanges)
library(devtools)
library("NMF")
library(grid)
library(gridExtra)
library(MutationalPatterns)
library(ggplot2)

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")

## Get cohort
cohort = purple::query_purity(dbProd)
patientIdLookups = query_patient_id_lookup(dbProd)
patientIds = purple::apply_to_cohort(cohort, function(x) {purple::sample_to_patient_id(x$sampleId, patientIdLookups)})
cohort$patientId <- patientIds$V1
multipleBiopsyCohort = purple::multiple_biopsy(cohort)
multipleBiopsyPatientsId = unique(multipleBiopsyCohort$patientId)
save(cohort, multipleBiopsyCohort, multipleBiopsyPatientsId, file="~/hmf/RData/multipleBiopsyCohort.RData")

## Get somatic variants
multipleBiopsySomaticVariants = purple::query_somatic_variants(dbProd, multipleBiopsyCohort)
save(multipleBiopsySomaticVariants, file="~/hmf/RData/multipleBiopsySomaticVariants.RData")

## Get structural variants
multipleBiopsyStructuralVariants = purple::query_structural_variants(dbProd, multipleBiopsyCohort)
save(multipleBiopsyStructuralVariants, file="~/hmf/RData/multipleBiopsyStructuralVariants.RData")

### Clean up
dbDisconnect(dbProd)
rm(dbProd)

#### LOAD DATA
load(file="~/hmf/multipleBiopysyCohort.RData")
load(file="~/hmf/multipleBiopsySomaticVariants.RData")
load(file="~/hmf/multipleBiopsyStructuralVariants.RData")

### Mutational Signatures
multipleBiopsyMutationalSignature = list()
multipleBiopsyIndelSignature = list()
multipleBiopsySVSignature = list()

#patientId = multipleBiopsyPatientsId[71]
#patientId = "CPCT02050172"

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

  pdf(file=paste("~/mb/",patientId, "MultipleBiopsies.pdf", sep = ""), height = 14, width = 20)    
  multiplot(somaticPloidyPlots[[1]], somaticPloidyPlots[[2]], somaticPloidyPlots[[3]],
            structuralPloidyPlots[[1]], structuralPloidyPlots[[2]], structuralPloidyPlots[[3]], 
            p7, p8, p9,
            cols = 3)
  dev.off()
}

save(multipleBiopsyMutationalSignature, multipleBiopsyIndelSignature, multipleBiopsySVSignature, file="~/hmf/multipleBiopsySignatures.RData")

