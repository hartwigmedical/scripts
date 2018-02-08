detach("package:purple", unload=TRUE);
library(purple);
library(RMySQL)
library(data.table)
library(MutationalPatterns)
library(ggplot2)

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
purity = purple::query_purity(dbProd)
clinical = purple::query_clinical_data(dbProd)[, c("sampleId", "cancerType")]
sampleCohort = merge(purity, clinical, by=c("sampleId"), all.x = TRUE)

patientIdLookups = query_patient_id_lookup(dbProd)
patientIds = purple::apply_to_cohort(sampleCohort, function(x) {purple::sample_to_patient_id(x$sampleId, patientIdLookups)})
colnames(patientIds) <- c("sampleId", "patientId")

sampleCohort = merge(sampleCohort, patientIds, by=c("sampleId"), all.x = TRUE)
patientCohort = highest_purity_patients(sampleCohort)
rm(purity)
rm(patientIds)
rm(patientIdLookups)
rm(clinical)

patientVariants = purple::query_somatic_variants(dbProd, patientCohort)
save(sampleCohort, patientCohort, patientVariants, file = "~/hmf/RData/patientVariants.RData")
load(file = "~/hmf/RData/patientVariants.RData")

dbDisconnect(dbProd)
rm(dbProd)

################ Recurrent Somatics ###################
snps = patientVariants[patientVariants$type == 'SNP', ]
snps = purple::add_scope_to_variants(snps, by = c("chromosome", "position", "alt", "ref", "type"))
snps$mut = purple::standard_mutation(paste(snps$ref, snps$alt, sep = ">"))

recurrentSnps = snps[, .(total = .N, recurrent = sum(.SD$scope == "Shared")), by=c("sampleId", "mut")]
recurrentSnps$percentRecurrent = recurrentSnps$recurrent / recurrentSnps$total

################ Recurrent Single Deletes ###################
indels = patientVariants[patientVariants$type == 'INDEL', ]
singleDels = indels[nchar(indels$ref) == 2 & nchar(indels$alt) == 1, ]
singleDels = purple::add_scope_to_variants(singleDels, by = c("chromosome", "position", "alt", "ref", "type"))
singleDels$mut <- paste(substr(singleDels$ref, 2, 2), "deletion", sep = " ")
singleDels$mut <- sub("T deletion", "A/T deletion", singleDels$mut, fixed = T)
singleDels$mut <- sub("A deletion", "A/T deletion", singleDels$mut, fixed = T)
singleDels$mut <- sub("G deletion", "C/G deletion", singleDels$mut, fixed = T)
singleDels$mut <- sub("C deletion", "C/G deletion", singleDels$mut, fixed = T)

recurrentSingleDels = singleDels[, .(total = .N, recurrent = sum(.SD$scope == "Shared")), by=c("sampleId", "mut")]
recurrentSingleDels$percentRecurrent = recurrentSingleDels$recurrent / recurrentSingleDels$total

################ Recurrent Single Insertions ###################
singleInserts = indels[nchar(indels$ref) == 1 & nchar(indels$alt) == 2, ]
singleInserts = purple::add_scope_to_variants(singleInserts, by = c("chromosome", "position", "alt", "ref", "type"))
singleInserts$mut <- paste(substr(singleInserts$alt, 2, 2), "insertion", sep = " ")
singleInserts$mut <- sub("T insertion", "A/T insertion", singleInserts$mut, fixed = T)
singleInserts$mut <- sub("A insertion", "A/T insertion", singleInserts$mut, fixed = T)
singleInserts$mut <- sub("G insertion", "C/G insertion", singleInserts$mut, fixed = T)
singleInserts$mut <- sub("C insertion", "C/G insertion", singleInserts$mut, fixed = T)

recurrentSingleInserts = singleInserts[, .(total = .N, recurrent = sum(.SD$scope == "Shared")), by=c("sampleId", "mut")]
recurrentSingleInserts$percentRecurrent = recurrentSingleInserts$recurrent / recurrentSingleInserts$total 

################ SAVE DATA ###################
save(sampleCohort, patientCohort, recurrentSnps, recurrentSingleDels, recurrentSingleInserts, file = "~/hmf/RData/recurrentVariants.RData")

############### PLOTS #########################
recurrentPlot<-function(mutation, recurrentSnps, cohort) {
  data = recurrentSnps[recurrentSnps$mut == mutation, ]
  data = merge(data, cohort[,c("sampleId", "cancerType")], by = "sampleId", all.x = T)
  ggplot()+labs(title=mutation, x = "% Recurrent SSMs", y = "# Recurrent SSMs")+
    geom_point(data=data,aes(percentRecurrent, recurrent, color=cancerType))+
    scale_x_continuous(labels = scales::percent, limits=c(0,1))+
    scale_y_log10(
      limits = c(10^0,10^5),
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    annotation_logticks(sides = "l") + theme(legend.position = "none") 
}

purple::multiplot(
  recurrentPlot("C>A", recurrentSnps, patientCohort), 
  recurrentPlot("T>A", recurrentSnps, patientCohort), 
  recurrentPlot("C>G", recurrentSnps, patientCohort), 
  recurrentPlot("T>C", recurrentSnps, patientCohort), 
  recurrentPlot("C>T", recurrentSnps, patientCohort),  
  recurrentPlot("T>G", recurrentSnps, patientCohort), 
cols = 3)

purple::multiplot(
  recurrentPlot("A/T deletion", recurrentSingleDels, patientCohort),
  recurrentPlot("A/T insertion", recurrentSingleInserts, patientCohort),
  recurrentPlot("C/G deletion", recurrentSingleDels, patientCohort),
  recurrentPlot("C/G insertion", recurrentSingleInserts, patientCohort),
  cols = 2)

  