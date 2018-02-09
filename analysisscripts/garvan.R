detach("package:purple", unload=TRUE); 
library(purple);
library(RMySQL)
library(data.table)
library(MutationalPatterns)
library(ggplot2)

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
#dbProd = dbConnect(MySQL(), dbname='hmfpatients', user="build", password="build", port=3307, host = "127.0.0.1")
purity = purple::query_purity(dbProd)
cohort = purity[1:5,  c("sampleId","gender", "purity", "ploidy")]

cat("Querying Clinical Data", "\n")
clinical = purple::query_clinical_data(dbProd)[, c("sampleId", "cancerType")]
cohort =  merge(cohort, clinical, all.x=T, by="sampleId")
rm(clinical)

cat("Querying microsatellite instability", "\n")
msi = purple::apply_to_cohort(cohort, function(x) {purple::query_msi_sample(dbProd, x$sampleId)})
cohort = merge(cohort, msi, all.x=T, by="sampleId")
rm(msi)

cat("Querying whole genome duplicationg", "\n")
wgd = purple::apply_to_cohort(cohort, function(x) {purple::query_whole_genome_duplication(dbProd, x$sampleId)})
cohort = merge(cohort, wgd, all.x=T, by="sampleId")
cohort$duplicatedAutosomes <- ifelse(is.na(cohort$duplicatedAutosomes), 0, cohort$duplicatedAutosomes)
cohort$WGD <- ifelse(is.na(cohort$WGD), FALSE, cohort$WGD)
rm(wgd)

cat("Querying structural variants", "\n")
svTypes = purple::apply_to_cohort(cohort, function(x) {purple::query_structural_variant_overview(dbProd, x$sampleId)})
cohort = merge(cohort, svTypes, all.x=T, by="sampleId")
rm(svTypes)

cat("Querying somatic variants", "\n")
somaticOverview = purple::apply_to_cohort(cohort, function(x) {purple::query_somatic_overview(dbProd, x$sampleId)})
cohort = merge(cohort, somaticOverview, all.x=T, by="sampleId")
rm(somaticOverview)

# Single sample Characteristics
sampleSomatics = purple::query_somatic_variants(dbProd, cohort[1, ])
sampleCosmicSignature = purple::mutational_signature_by_clonality(sampleSomatics)
plot_cosmic_signature(sampleCosmicSignature)

sampleIndelSignature = purple::indel_signature_by_clonality(sampleSomatics)
sampleIndelSignature[, c("NA")] <- NULL
plot_indel_signature(sampleIndelSignature)

sampleSVs = purple::query_structural_variants(dbProd, cohort[1, ])
sampleSVSignature = purple::sv_signature_by_clonality(sampleSVs)
plot_sv_signature(sampleSVSignature)

# Mutliple Biopsy Characteristics
multipleSamples = c("CPCT02080001T", "CPCT02080001TII")
multipleCohort = purity[purity$sampleId %in% multipleSamples,  c("sampleId","gender", "purity", "ploidy")]
patientSomaticVariants = purple::query_somatic_variants(dbProd, multipleCohort)
patientSomaticVariants = purple::add_scope_to_variants(patientSomaticVariants, c("chromosome", "position", "alt", "ref", "type"))
patientStructuralVariants = purple::query_structural_variants(dbProd, multipleCohort)
patientStructuralVariants = purple::add_scope_to_variants(patientStructuralVariants, c("startChromosome", "endChromosome", "startPosition", "endPosition", "startOrientation", "endOrientation", "type"))

patientMutationalSignature = purple::mutational_signature_by_scope(patientSomaticVariants)
plot_cosmic_signature(patientMutationalSignature)

patientIndelSignature = purple::indel_signature_by_scope(patientSomaticVariants)
plot_indel_signature(patientIndelSignature)

patientSVSignature = purple::sv_signature_by_scope(patientStructuralVariants)
plot_sv_signature(patientSVSignature)





dbDisconnect(dbProd)
rm(dbProd)

