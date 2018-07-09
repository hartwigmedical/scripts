############################################
# How to run analysis for paper:
# 01. cohortClinical.R
# 02. cohort.R
# 03. cohortFusions.R
# 04. dnds.R
# 05. genePanel.R
# 06. dndsClassification.R
# 07. dndsDrivers.R
# 08. ampsDels.R
# 09. ampsDelsTarget.R
# 10. driversByGene.R
# 11. multipleBiopsyDrivers.R


# 12. cohortVisualisation.R
# 13. copyNumberOverview.R
# 14. driversByGeneVisualisation.R
# 15. driversByGeneHeatmap.R

detach("package:purple", unload=TRUE)
library(purple)
library(RMySQL)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(multidplyr)

load(file = '~/hmf/RData/reference/allGeneDeletes.RData')
load(file = "~/hmf/RData/reference/allClinicalData.RData")

### DATABASE
dbProd = dbConnect(MySQL(), dbname='hmfpatients_20180418', groups="RAnalysis")
dbDisconnect(dbProd)
rm(dbProd)

############################################ ENTIRE COHORT QUERIES
cat("Querying gene deletes")
allGeneDeletes = purple::query_gene_deletes(dbProd)
save(allGeneDeletes, file = '~/hmf/RData/reference/allGeneDeletes.RData')

cat("Querying purity")
query_entire_cohort <- function(dbConnect, purityCutoff = 0.2) {
  query = paste(
    "SELECT p.*",
    " FROM purity p",
    sep = " ")
  return (dbGetQuery(dbConnect, query))
}
allPurity = query_entire_cohort(dbProd) %>%
  left_join(allGeneDeletes, by = "sampleId") %>%
  mutate(genesDeleted = ifelse(is.na(genesDeleted), 0, genesDeleted)) %>%
  mutate(qcStatus = ifelse(genesDeleted > 280, "FAIL_DELETED_GENES", qcStatus)) %>%
  filter(sampleId != 'CPCT02050303T')
patientIdLookups = query_patient_id_lookup(dbProd)
allPurity$patientId <- sapply(allPurity$sampleId, function(x) {purple::sample_to_patient_id(x, patientIdLookups)})
allPurity = left_join(allPurity, allClinicalData %>% select(sampleId, cancerType), by = "sampleId")
allPurity$gender = ifelse(substr(allPurity$gender, 1, 4) == "MALE", "MALE", allPurity$gender)
save(allPurity, file = '~/hmf/RData/reference/allPurity.RData' )
write.table(allPurity %>% select(patientId), file = "~/hmf/resources/cpctPatientIds.csv", row.names = F, quote = F, col.names = F)
rm(allGeneDeletes, patientIdLookups)

cat("Querying metrics")
allMetrics = purple::query_metrics(dbProd, allPurity)
save(allMetrics, file = "~/hmf/RData/Reference/allMetrics.RData")

cat("Querying somatics")
allSomatics_p1 = purple::query_somatic_variants(dbProd, allPurity[1:1500, ])
save(allSomatics_p1, file = "~/hmf/RData/reference/allSomatics_p1.RData")
allSomatics_p2 = purple::query_somatic_variants(dbProd, allPurity[1501:nrow(allPurity), ])
save(allSomatics_p2, file = "~/hmf/RData/reference/allSomatics_p2.RData")

allIndelsInPONQuery = "select * from somaticVariant where type = 'INDEL' and filter <> 'PASS' and abs(length(ref) - length(alt)) >= 3"
allIndelsInPON = dbGetQuery(dbProd, allIndelsInPONQuery)
save(allIndelsInPON, file = "~/hmf/RData/reference/allIndelsInPON.RData")

nearPon <- function(samplePon, sampleIndels, distance = 10) {
  hrange <- GRanges(samplePon$chromosome, IRanges(samplePon$position, samplePon$position + nchar(samplePon$ref) + distance - 1))
  mrange <- GRanges(sampleIndels$chromosome, IRanges(sampleIndels$position, sampleIndels$position + nchar(sampleIndels$ref) + distance - 1))
  
  ol = as.matrix(findOverlaps(hrange, mrange, type="any", select="all"))
  sampleIndels$nearPon <- FALSE
  sampleIndels[ol[,2], c("nearPon")] <- TRUE
  return (sampleIndels)
}

nearPonIndel <- function(pon, indels) {
  result = data.frame()
  for (selectedSample in unique(indels$sampleId)) {
    cat("Processing ", selectedSample, "\n")
    samplePon = pon %>% filter(sampleId == selectedSample)
    sampleIndels = indels %>% filter(sampleId == selectedSample)
    result = bind_rows(nearPon(samplePon, sampleIndels), result)
  }
  return (result %>% filter(nearPon))
}

allIndelsNearPON_p1 = nearPonIndel(allIndelsInPON, allSomatics_p1 %>% filter(type == 'INDEL', abs(nchar(ref) - nchar(alt)) >= 3))
allIndelsNearPON_p2 = nearPonIndel(allIndelsInPON, allSomatics_p2 %>% filter(type == 'INDEL', abs(nchar(ref) - nchar(alt)) >= 3))
allIndelsNearPON = bind_rows(allIndelsNearPON_p1, allIndelsNearPON_p2)
rm(allIndelsNearPON_p1, allIndelsNearPON_p2, allIndelsInPON)
save(allIndelsNearPON, file = "~/hmf/RData/reference/allIndelsNearPON.RData")

allSomatics_p1 = allSomatics_p1 %>% left_join(allIndelsNearPON %>% select(id, nearPon), by = "id") %>% filter(is.na(nearPon))
allSomatics_p2 = allSomatics_p2 %>% left_join(allIndelsNearPON %>% select(id, nearPon), by = "id") %>% filter(is.na(nearPon))

fix_splice_effect <- function(input) {
  result = input %>%
    mutate(canonicalCodingEffect=ifelse(type == 'INDEL' & canonicalCodingEffect == 'SPLICE' & canonicalEffect == "splice region variant; intron variant","NONE",canonicalCodingEffect)) %>%
    mutate(canonicalCodingEffect=ifelse(type == 'INDEL' & canonicalCodingEffect == 'SPLICE' & canonicalEffect == "splice region variant; inframe deletion","MISSENSE",canonicalCodingEffect)) %>%
    mutate(canonicalCodingEffect=ifelse(type == 'INDEL' & canonicalCodingEffect == 'SPLICE' & canonicalEffect == "splice region variant","NONE",canonicalCodingEffect) ) %>%
    mutate(canonicalCodingEffect=ifelse(type == 'INDEL' & canonicalCodingEffect == 'SPLICE' & canonicalEffect == "splice region variant; non coding transcript variant","NONE",canonicalCodingEffect)) %>%
    mutate(canonicalCodingEffect=ifelse(type == 'INDEL' & canonicalCodingEffect == 'SPLICE' & canonicalEffect == "start lost; splice region variant; inframe deletion","MISSENSE",canonicalCodingEffect)) %>%
    mutate(canonicalCodingEffect=ifelse(type == 'INDEL' & canonicalCodingEffect == 'SPLICE' & canonicalEffect == "splice region variant; inframe insertion","MISSENSE",canonicalCodingEffect)) %>%
    mutate(canonicalCodingEffect=ifelse(type == 'INDEL' & canonicalCodingEffect == 'SPLICE' & canonicalEffect == "missense variant; splice region variant","MISSENSE",canonicalCodingEffect)) %>%
    mutate(canonicalCodingEffect=ifelse(type == 'INDEL' & canonicalCodingEffect == 'SPLICE' & canonicalEffect == "splice region variant; synonymous variant; inframe deletion","MISSENSE",canonicalCodingEffect))
}

allSomatics_p1 = fix_splice_effect(allSomatics_p1)
allSomatics_p2 = fix_splice_effect(allSomatics_p2)

save(allSomatics_p1, file = "~/hmf/RData/reference/allSomatics_p1.RData")
save(allSomatics_p2, file = "~/hmf/RData/reference/allSomatics_p2.RData")

somatics_summary_p1 = cohort_somatic_summary(allSomatics_p1)
somatics_summary_p2 = cohort_somatic_summary(allSomatics_p2)
allSomaticsSummary = rbind(somatics_summary_p1, somatics_summary_p2)
save(allSomaticsSummary, file = "~/hmf/RData/reference/allSomaticsSummary.RData")

allStructuralVariantSummary = query_structural_variant_summary(dbProd, allPurity)
save(allStructuralVariantSummary, file = "~/hmf/RData/reference/allStructuralVariantSummary.RData")

allWgd = purple::query_whole_genome_duplication(dbProd, allPurity)
save(allWgd, file = "~/hmf/RData/reference/allWgd.RData")

allSampleData = purple::query_sample_data(dbProd)
save(allSampleData, file = "~/hmf/RData/reference/allSampleData.RData")

cat("Querying canonical transcripts")
canonicalTranscripts = purple::query_canonical_transcript(dbProd)
save(canonicalTranscripts, file = "~/hmf/RData/reference/canonicalTranscripts.RData")

allSNPSummary_p1 = allSomatics_p1 %>% filter(filter == 'PASS', type == 'SNP') %>% group_by(sampleId, ref, alt) %>% summarise(n = n())
allSNPSummary_p2 = allSomatics_p2 %>% filter(filter == 'PASS', type == 'SNP') %>% group_by(sampleId, ref, alt) %>% summarise(n = n())
allSNPSummary = bind_rows(allSNPSummary_p1, allSNPSummary_p2)
save(allSNPSummary, file = "~/hmf/RData/Reference/allSNPSummary.RData")

allMNPSummary_p1 = allSomatics_p1 %>% filter(filter == 'PASS', type == 'MNP') %>% group_by(sampleId, ref, alt) %>% summarise(n = n())
allMNPSummary_p2 = allSomatics_p2 %>% filter(filter == 'PASS', type == 'MNP') %>% group_by(sampleId, ref, alt) %>% summarise(n = n())
allMNPSummary = bind_rows(allMNPSummary_p1, allMNPSummary_p2)
save(allMNPSummary, file = "~/hmf/RData/Reference/allMNPSummary.RData")

indel_summary <- function(somatics) {
  result = somatics %>% filter(filter == 'PASS', type == 'INDEL') %>%
    mutate(isRepeat = repeatCount > 3, isMicrohomology = microhomology != "") %>% 
    mutate(
      category1 = ifelse(nchar(ref) > nchar(alt), "DEL", "INS"),  
      category2 = ifelse(category1 == "DEL" & isMicrohomology, "MH", "other"),
      category2 = ifelse(isRepeat, "repeat", category2),
      category = paste(category1, category2, sep = "-")
    ) %>%
    group_by(sampleId, category)  %>% 
    summarise(n = n()) 
}

allIndelSummary_p1 = indel_summary(allSomatics_p1)
allIndelSummary_p2 = indel_summary(allSomatics_p2)
allIndelSummary = bind_rows(allIndelSummary_p1, allIndelSummary_p2)
save(allIndelSummary, file = "~/hmf/RData/Reference/allIndelSummary.RData")


#### COMBINE
load(file = "~/hmf/RData/reference/allPurity.RData")
load(file = "~/hmf/RData/reference/allWgd.RData")
load(file = "~/hmf/RData/reference/allClinicalData.RData")
load(file = "~/hmf/RData/reference/allSampleData.RData")
load(file = "~/hmf/RData/reference/allSomaticsSummary.RData")
load(file = "~/hmf/RData/reference/allStructuralVariantSummary.RData")
load(file = "~/hmf/RData/reference/allMetrics.RData")

clinicalSummary = allClinicalData %>% select(sampleId, primaryTumorLocation, cancerSubtype, biopsyDate, biopsySite, biopsyType, biopsyLocation, birthYear) %>%
  mutate(ageAtBiopsy = as.numeric(substr(biopsyDate, 1, 4)) - birthYear)

acceptableStatus = c('NORMAL','SOMATIC','HIGHLY_DIPLOID')

cohortSummary = allPurity %>%
  select(sampleId, patientId, gender, status, qcStatus, purity, ploidy, genesDeleted, cancerType, score) %>%
  left_join(clinicalSummary, by = "sampleId") %>%
  left_join(allSampleData, by = "sampleId") %>%
  left_join(allSomaticsSummary, by = "sampleId") %>%
  left_join(allStructuralVariantSummary, by = "sampleId") %>%
  left_join(allWgd, by = "sampleId") %>%  mutate(duplicatedAutosomes = ifelse(is.na(duplicatedAutosomes), 0, duplicatedAutosomes), WGD = ifelse(is.na(WGD), F, WGD)) %>%
  left_join(allMetrics %>% select(sampleId, refMeanCoverage, tumorMeanCoverage), by = "sampleId") %>%
  mutate(
    purity = ifelse(status == 'NO_TUMOR', 0, purity),
    status = ifelse(status == 'NO_TUMOR', 'FAIL_TUMOR', status),
    status = ifelse(status %in% acceptableStatus & qcStatus != 'PASS', qcStatus, status),
    status = ifelse(status %in% acceptableStatus & purity < 0.2, 'FAIL_PURITY', status),
    status = ifelse(status %in% acceptableStatus, 'PASS',status)) %>%
  group_by(patientId) %>% 
  arrange(arrivalDate) %>% mutate(patientSample = row_number())

maxPassPuritySampleIds = cohortSummary %>% 
  filter(qcStatus == 'PASS') %>% 
  group_by(patientId, purity) %>% top_n(1, -score) %>% ungroup() %>% group_by(patientId) %>% top_n(1, purity) %>% ungroup() %>% pull(sampleId)

cohortSummary = cohortSummary %>%
  ungroup() %>%
  mutate(patientHighestPurityPassingSample = sampleId %in% maxPassPuritySampleIds) %>%
  select(-qcStatus, -score, -genesDeleted)

cohortSummary %>% group_by(status) %>% count()
cohortSummary %>% filter(status == 'PASS') %>% group_by(patientId) %>% count() %>% filter(n  > 1) %>% arrange(-n)

cohortSummary %>% filter(sample  > 1)

write.csv(cohortSummary, file = "~/hmf/RData/CohortSummary.csv", row.names = F) 

rm(cohortSummary)

##### GENERATE PATIENT/SAMPLEID
#load(file = "~/hmf/RData/reference/allClinicalData.RData")
#load(file = "~/hmf/RData/reference/allPurity.RData")
#cpctPatientIds =  unique(allPurity$patientId)
#hmfPatientIds = paste0("HMF", formatC(c(1:2781), width = 5, format = "d", flag = "0"))
#patientIdMap = data.frame(cpctPatientId = cpctPatientIds, hmfPatientId = hmfPatientIds, stringsAsFactors = F)
#rm(cpctPatientIds, hmfPatientIds)

#sampleIdMap = allPurity %>% select(sampleId, patientId) %>% 
#  left_join(patientIdMap, by = c("patientId" = "cpctPatientId")) %>%
#  left_join(allClinicalData %>% select(sampleId, sampleArrivalDate), by = "sampleId") %>%
#  arrange(hmfPatientId, sampleArrivalDate) %>%
#  group_by(hmfPatientId) %>%
#  mutate(
#    sample = row_number(),
#    sample = chartr("12345", "ABCDE", sample)) %>%
#  mutate(hmfSampleId = paste0(hmfPatientId, sample))

#sampleIdMap$insert = paste0("(\'", sampleIdMap$hmfSampleId, "\',\'", sampleIdMap$patientId, "\',\'", sampleIdMap$sampleId, "\')")
#cat(paste0(sampleIdMap$insert, collapse = ","))

############################################ HIGHEST PURITY
load(file = '~/hmf/RData/reference/allGeneDeletes.RData')
load(file = "~/hmf/RData/reference/allClinicalData.RData")

cat("Querying purity")
highestPurityCohort = purple::query_highest_purity_cohort(dbProd, allGeneDeletes)
highestPurityCohort = left_join(highestPurityCohort, allClinicalData %>% select(sampleId, cancerType), by = "sampleId") %>% filter(!is.na(cancerType))
highestPurityCohort$gender = ifelse(substr(highestPurityCohort$gender, 1, 4) == "MALE", "MALE", highestPurityCohort$gender)
save(highestPurityCohort, file = "~/hmf/RData/reference/highestPurityCohort.RData")

cat("Copy Numbers")
hpcCopyNumbers = purple::query_copy_number(dbProd, highestPurityCohort)
save(hpcCopyNumbers, file = "~/hmf/RData/reference/hpcCopyNumbers.RData")

cat("Querying gene copy number deletes")
hpcGeneCopyNumberDeletes = purple::query_gene_copy_number_deletes(dbProd, highestPurityCohort)
hpcGeneCopyNumberDeletes = left_join(hpcGeneCopyNumberDeletes, highestPurityCohort %>% select(sampleId, cancerType), by = "sampleId")
save(hpcGeneCopyNumberDeletes, file = "~/hmf/RData/reference/hpcGeneCopyNumberDeletes.RData")

cat("Querying gene copy number amplifications")
hpcGeneCopyNumberAmplifications = purple::query_gene_copy_number_amplifications(dbProd, highestPurityCohort)
hpcGeneCopyNumberAmplifications = left_join(hpcGeneCopyNumberAmplifications, highestPurityCohort %>% select(sampleId, cancerType), by = "sampleId")
save(hpcGeneCopyNumberAmplifications, file = "~/hmf/RData/reference/hpcGeneCopyNumberAmplifications.RData")

cat("Determing exonic somatics")
load(file = "~/hmf/RData/reference/HmfRefCDS.RData")
load(file = "~/hmf/RData/reference/allSomatics_p1.RData")
load(file = "~/hmf/RData/reference/allSomatics_p2.RData")

exonic_somatics <- function(somatics, gr_genes) {
  gr_muts = GRanges(somatics$chromosome, IRanges(somatics$position,somatics$position + nchar(somatics$ref) - 1))
  ol = as.matrix(findOverlaps(gr_muts, gr_genes, type="any", select="all"))
  return (somatics[unique(ol[, 1]), ])
}
exonic_p1 = exonic_somatics(allSomatics_p1 %>% filter(sampleId %in% highestPurityCohort$sampleId), gr_genes)
exonic_p2 = exonic_somatics(allSomatics_p2 %>% filter(sampleId %in% highestPurityCohort$sampleId), gr_genes)
hpcExonicSomatics = rbind(exonic_p1, exonic_p2)
save(hpcExonicSomatics, file = "~/hmf/RData/reference/hpcExonicSomatics.RData")
rm(exonic_p1, exonic_p2)

cat("Tert promoters")
hpcTertPromoters = purple::query_tert_promoters(dbProd, highestPurityCohort)
save(hpcTertPromoters, file = "~/hmf/RData/reference/hpcTertPromoters.RData")


#### COMBINE
load(file = "~/hmf/RData/Reference/highestPurityCohort.RData")
load(file = "~/hmf/RData/Reference/allWgd.RData")
load(file = "~/hmf/RData/Reference/allClinicalData.RData")
load(file = "~/hmf/RData/Reference/allSampleData.RData")
load(file = "~/hmf/RData/Reference/allSomaticsSummary.RData")
load(file = "~/hmf/RData/Reference/allMetrics.RData")
load(file = "~/hmf/RData/Reference/allStructuralVariantSummary.RData")

clinicalSummary = allClinicalData %>% select(sampleId, primaryTumorLocation, cancerSubtype, biopsyDate, biopsySite, biopsyType, biopsyLocation, treatment, treatmentType, birthYear) %>%
  mutate(ageAtBiopsy = as.numeric(substr(biopsyDate, 1, 4)) - birthYear)
highestPurityCohortSummary = highestPurityCohort %>%
  select(sampleId, patientId, gender, status, qcStatus, purity, ploidy, genesDeleted, cancerType) %>%
  left_join(clinicalSummary, by = "sampleId") %>%
  left_join(allSampleData, by = "sampleId") %>%
  left_join(allSomaticsSummary, by = "sampleId") %>%
  left_join(allStructuralVariantSummary, by = "sampleId") %>%
  left_join(allWgd, by = "sampleId") %>%  mutate(duplicatedAutosomes = ifelse(is.na(duplicatedAutosomes), 0, duplicatedAutosomes), WGD = ifelse(is.na(WGD), F, WGD)) %>%
  left_join(allMetrics %>% select(sampleId, refMeanCoverage, tumorMeanCoverage), by = "sampleId")

save(highestPurityCohortSummary, file = "~/hmf/RData/Processed/highestPurityCohortSummary.RData")  
write.csv(highestPurityCohortSummary, file = "~/hmf/RData/HighestPurityCohortSummary.csv", row.names = F) 


############## MULTIPLE BIOPSY COHORT
load('~/hmf/RData/reference/allClinicalData.RData')
load('~/hmf/RData/reference/allGeneDeletes.RData')

multipleBiopsyCohort = purple::query_multiple_biopsy_cohort(dbProd, allGeneDeletes)
multipleBiopsyCohort = left_join(multipleBiopsyCohort, allClinicalData %>% select(sampleId, cancerType), by = "sampleId")
save(multipleBiopsyCohort, file = "~/hmf/RData/Reference/multipleBiopsyCohort.RData")

multipleBiopsyScope = multipleBiopsyCohort %>% 
  left_join(allClinicalData %>% select(sampleId, sampleArrivalDate), by = "sampleId") %>%
  arrange(patientId, sampleArrivalDate) %>%
  select(patientId, sampleId) %>% 
  group_by(patientId) %>% 
  mutate(scope = paste0("Sample",row_number())) %>%
  ungroup()
save(multipleBiopsyScope, file = "~/hmf/RData/reference/multipleBiopsyScope.RData")

multipleBiopsyStructuralVariants = query_structural_variants(dbProd, multipleBiopsyCohort)
save(multipleBiopsyStructuralVariants, file = "~/hmf/RData/reference/multipleBiopsyStructuralVariants.RData")

load(file = "~/hmf/RData/reference/allSomatics_p1.RData")
load(file = "~/hmf/RData/reference/allSomatics_p2.RData")
multipleBiopsySomatics = bind_rows(allSomatics_p1 %>% filter(sampleId %in% multipleBiopsyCohort$sampleId), allSomatics_p2 %>% filter(sampleId %in% multipleBiopsyCohort$sampleId)) 
save(multipleBiopsySomatics, file = "~/hmf/RData/reference/multipleBiopsySomatics.RData")

multipleBiopsyStructuralVariantsWithScope = left_join(multipleBiopsyStructuralVariants, multipleBiopsyScope, by = "sampleId") %>%
  group_by(patientId, startChromosome, endChromosome, startPosition,endPosition,startOrientation,endOrientation,type) %>%
  mutate(scope = ifelse(n_distinct(sampleId) > 1, "Shared", scope))
save(multipleBiopsyStructuralVariantsWithScope, file = "~/hmf/RData/reference/multipleBiopsyStructuralVariantsWithScope.RData")

multipleBiopsySomaticsWithScope = multipleBiopsySomatics %>% 
  filter(filter == 'PASS') %>%
  left_join(multipleBiopsyScope, by = "sampleId") %>%
  partition(patientId, chromosome, position, ref, alt, type) %>%
  mutate(scope = ifelse(n_distinct(sampleId) > 1, "Shared", scope)) %>%
  collect() %>% 
  as_tibble() %>%
  ungroup()
save(multipleBiopsySomaticsWithScope, file = "~/hmf/RData/reference/multipleBiopsySomaticsWithScope.Rdata")

mbExonicSomatics = exonic_somatics(multipleBiopsySomaticsWithScope, gr_genes)
save(mbExonicSomatics, file = "~/hmf/RData/reference/mbExonicSomatics.RData")

multipleBiopsyMSI = purple::cohort_msi(multipleBiopsySomaticsWithScope %>% ungroup())
save(multipleBiopsyMSI, file = "~/hmf/RData/reference/multipleBiopsyMSI.RData")

cat("Tert promoters")
mbTertPromoters = purple::query_tert_promoters(dbProd, multipleBiopsyCohort)
mbTertPromoters = mbTertPromoters %>%
  left_join(multipleBiopsyScope, by = "sampleId") %>%
  group_by(patientId, gene) %>%
  mutate(scope = ifelse(n_distinct(sampleId) > 1, "Shared", scope)) %>%
  ungroup()

save(mbTertPromoters, file = "~/hmf/RData/reference/mbTertPromoters.RData")

cat("Querying gene copy number amplifications")
mbGeneCopyNumberAmplifications = purple::query_gene_copy_number_amplifications(dbProd, multipleBiopsyCohort)
mbGeneCopyNumberAmplifications = left_join(mbGeneCopyNumberAmplifications, multipleBiopsyCohort %>% select(sampleId, cancerType), by = "sampleId") %>%
  left_join(multipleBiopsyScope, by = "sampleId") %>%
  group_by(patientId, gene) %>%
  mutate(scope = ifelse(n_distinct(sampleId) > 1, "Shared", scope)) %>%
  ungroup()
save(mbGeneCopyNumberAmplifications, file = "~/hmf/RData/reference/mbGeneCopyNumberAmplifications.RData")

cat("Querying gene copy number deletes")
mbGeneCopyNumberDeletes = purple::query_gene_copy_number_deletes(dbProd, multipleBiopsyCohort)
mbGeneCopyNumberDeletes = left_join(mbGeneCopyNumberDeletes, multipleBiopsyCohort %>% select(sampleId, cancerType), by = "sampleId") %>%
  left_join(multipleBiopsyScope, by = "sampleId") %>%
  group_by(patientId, gene) %>%
  mutate(scope = ifelse(n_distinct(sampleId) > 1, "Shared", scope)) %>%
  ungroup()
save(mbGeneCopyNumberDeletes, file = "~/hmf/RData/reference/mbGeneCopyNumberDeletes.RData")


load(file = "~/hmf/RData/reference/multipleBiopsyCohort.RData")
load(file = "~/hmf/RData/reference/multipleBiopsyMSI.RData")
load(file = "~/hmf/RData/reference/multipleBiopsyScope.RData")
load(file = "~/hmf/RData/reference/multipleBiopsySomaticsWithScope.Rdata")
load(file = "~/hmf/RData/reference/multipleBiopsyStructuralVariantsWithScope.RData")
load(file = "~/hmf/RData/reference/allMetrics.RData")

multipleBiopsyMetrics = allMetrics %>% 
  filter(sampleId %in% multipleBiopsyScope$sampleId) %>%
  select(sampleId, refMeanCoverage, tumorMeanCoverage) %>%
  left_join(multipleBiopsyScope, by = "sampleId") %>%
  gather(type, value, refMeanCoverage, tumorMeanCoverage) %>%
  unite(type, scope, type) %>%
  select(-sampleId) %>%
  spread(type, value)

multipleBiopsyStructuralVariantSummary = multipleBiopsyStructuralVariantsWithScope %>%
  ungroup() %>% 
  distinct(patientId, scope, startChromosome, endChromosome, startPosition, endPosition, startOrientation, endOrientation, type) %>% 
  group_by(patientId, type, scope) %>% summarise(count = n()) %>% 
  unite(type, scope, type, sep = "_") %>% 
  spread(type, count)

multipleBiopsySomaticVariantSummary = multipleBiopsySomaticsWithScope %>%
  ungroup() %>% 
  distinct(patientId, scope, chromosome, position, ref, alt, type, clonality) %>% 
  mutate(
    type = ifelse(type == 'SNP', 'SNV', type),
    type = ifelse(type == 'MNP', 'MNV', type),
    clonality = ifelse(clonality == 'INCONSISTENT', 'CLONAL', clonality)) %>%
  group_by(patientId, scope, type, clonality) %>%
  summarise(count = n()) %>%
  unite(type, scope, type, clonality, sep = "_") %>% 
  spread(type, count, fill = 0)

multipleBiopsyPatientMsi = multipleBiopsyMSI %>% left_join(multipleBiopsyScope, by = "sampleId") %>%
  gather(type, value, msiScore, msiStatus) %>%
  unite(type, scope, type) %>% 
  select(-sampleId) %>%
  spread(type, value)

multipleBiopsySampleIds = multipleBiopsyScope %>% spread(scope, sampleId)

multipleBiopsyPurity = multipleBiopsyCohort %>% select(patientId, sampleId, purity, ploidy) %>% 
  left_join(multipleBiopsyScope %>% select(sampleId, scope), by = "sampleId") %>% 
  gather(type, value, purity, ploidy) %>%  
  unite(type, scope, type) %>%
  select(-sampleId) %>%
  spread(type, value)

multipleBiopsyClinicalSummary = allClinicalData %>% filter(sampleId %in% multipleBiopsyCohort$sampleId) %>% 
  select(sampleId, biopsySite, biopsyLocation, treatment) %>%
  left_join(multipleBiopsyScope, by = "sampleId")  %>%
  gather(type, value, biopsySite, biopsyLocation, treatment)  %>%
  unite(type, scope, type) %>%
  select(-sampleId) %>%
  spread(type, value)
  
multipleBiopsyCohortSummary = multipleBiopsyCohort %>% left_join(allClinicalData %>% select(sampleId, primaryTumorLocation), by = "sampleId") %>%
  select(patientId, sampleId, gender, primaryTumorLocation, cancerType) %>% 
  left_join(multipleBiopsyScope %>% select(sampleId, scope), by = "sampleId") %>% 
  filter(scope == "Sample1") %>% select(patientId, gender, primaryTumorLocation, cancerType) %>%
  left_join(multipleBiopsySampleIds, by = "patientId") %>%
  left_join(multipleBiopsyClinicalSummary, by = "patientId") %>%
  left_join(multipleBiopsyPurity, by = "patientId") %>%
  left_join(multipleBiopsyPatientMsi, by = "patientId") %>%
  left_join(multipleBiopsyStructuralVariantSummary, by = "patientId") %>%
  left_join(multipleBiopsySomaticVariantSummary, by = "patientId") %>%
  left_join(multipleBiopsyMetrics, by = "patientId") 
  
save(multipleBiopsyCohortSummary, file = "~/hmf/RData/processed/multipleBiopsyCohortSummary.RData")
write.csv(multipleBiopsyCohortSummary, file = "~/hmf/RData/MultipleBiopsyCohortSummary.csv", row.names = F) 


################### SANITY CHECKS

load(file = "~/hmf/RData/reference/multipleBiopsyCohort.RData")
load(file = "~/hmf/RData/processed/multipleBiopsyCohortSummary.RData")

inconsistentData = multipleBiopsyCohort %>% select(patientId, sampleId, gender, cancerType) %>% group_by(patientId) %>%
  mutate(n_gender = n_distinct(gender), n_cancerType = n_distinct(cancerType)) %>% 
  filter(n_gender > 1 |  n_cancerType > 1)

jon = allClinicalData %>% filter(patientId %in% inconsistentData$patientId)

clinicalDatabaseSummary = purple::query_clinical_data(dbProd)
clinicalDatabaseSummary = clinicalDatabaseSummary %>% filter(sampleId %in% inconsistentData$sampleId) %>% select(sampleId, primaryTumorLocation)
inconsistentData = left_join(inconsistentData, clinicalDatabaseSummary, by = "sampleId", suffix = c(".curated", ".database"))
