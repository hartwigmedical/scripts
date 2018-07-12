library(RMySQL)
library(dplyr)
library(tidyr)

query_purity<-function(dbConnect, purityCutoff = 0.2) {
  query = paste(
    "SELECT p.*",
    " FROM purity p",
    "WHERE qcStatus = 'PASS'",
    "  AND status <> 'NO_TUMOR'",
    "  AND p.purity >= ", purityCutoff,
    sep = " ")
  return (dbGetQuery(dbConnect, query))
}

query_copy_number <- function(dbConnect, cohort) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT g.*",
    "  FROM copyNumber g",
    " WHERE g.sampleId in (",sampleIdString, ")",
    sep = "")
  return (dbGetQuery(dbConnect, query))
}


query_gene_copy_number_deletes<-function(dbConnect, cohort, cutoff = 0.5) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT g.sampleId, g.chromosome, g.start, g.end, g.gene, g.chromosomeBand, g.minCopyNumber, 1 as score, g.minRegionStartSupport, g.minRegionEndSupport, g.somaticRegions, g.germlineHetRegions, g.germlineHomRegions",
    "  FROM geneCopyNumber g, purity p",
    " WHERE g.sampleId = p.sampleId",
    "   AND g.minCopyNumber < ", cutoff,
    "   AND g.chromosome <> 'Y'",
    "   AND p.sampleId in (",sampleIdString, ")",
    sep = " ")
  return (dbGetQuery(dbConnect, query))
}

query_gene_copy_number_amplifications<-function(dbConnect, cohort, cutoff = 3) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT g.sampleId, g.chromosome, g.start, g.end, g.gene, g.chromosomeBand, p.ploidy, g.minCopyNumber, log2(2 * g.minCopyNumber / p.ploidy / ", cutoff, ") as score, g.minRegionStartSupport, g.minRegionEndSupport, g.somaticRegions, g.germlineHetRegions, g.germlineHomRegions",
    "  FROM geneCopyNumber g, purity p",
    " WHERE g.sampleId = p.sampleId",
    "   AND g.minCopyNumber / p.ploidy > ", cutoff,
    "   AND p.sampleId in (",sampleIdString, ")",
    sep = " ")
  return (dbGetQuery(dbConnect, query))
}

query_somatic_variants <- function(dbConnect, cohort) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT *",
    "FROM somaticVariant",
    "WHERE sampleId in (",sampleIdString, ")",
    " AND filter = 'PASS' ",
    sep = " ")
  
  return (dbGetQuery(dbConnect, query))
}

query_structural_variant_summary<-function(dbConnect, cohort) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT sampleId, type, count(*) as count",
    " FROM structuralVariant ",
    "WHERE sampleId in (",sampleIdString, ")",
    " GROUP BY sampleId, type",
    sep = "")
  result = dbGetQuery(dbConnect, query)
  result = result  %>% group_by(sampleId) %>% spread(type, count)
  result[is.na(result)] <- 0
  return (result)
}

query_whole_genome_duplication<-function(dbConnect, cohort) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT sampleId, count(*) as duplicatedAutosomes",
    "  FROM ",
    " (select sampleId, chromosome, round(sum(bafCount*actualBaf*copyNumber)/sum(bafCount),1) as lwMajorAlleleAvg ",
    " from copyNumber ",
    " where sampleId in (", sampleIdString, ") and chromosome not in ('X', 'Y') ",
    " group by sampleId, chromosome ",
    " HAVING lwMajorAlleleAvg > 1.5) a ",
    " GROUP BY sampleId",
    sep = "")
  
  result = dbGetQuery(dbConnect, query)
  result$WGD <- ifelse(result$duplicatedAutosomes > 10, TRUE, FALSE)
  
  return (result)
}

### VALUES TO CHANGE IF YOU WANT
purityCutoff = 0.2
outputDir = "~/garvan/RData/"

### START OF REFERENCE DATA COLLECTION
referenceDir = paste0(outputDir, "reference/")
processedDir = paste0(outputDir, "processed/")
dbProd = dbConnect(MySQL(), dbname='hmfpatients_20180418', groups="RAnalysis")

cohort = query_purity(dbProd, purityCutoff)
cohortCopyNumbers = query_copy_number(dbProd, cohort) 
cohortGeneCopyNumberAmplifications = query_gene_copy_number_amplifications(dbProd, cohort)
cohortGeneCopyNumberDeletes = query_gene_copy_number_deletes(dbProd, cohort)
cohortSomatics = query_somatic_variants(dbProd, cohort)
cohortStructuralVariantSummary = query_structural_variant_summary(dbProd, cohort)
cohortWGD = query_whole_genome_duplication(dbProd, cohort)

save(cohort, file = paste0(referenceDir, "cohort.RData"))
save(cohortCopyNumbers, file = paste0(referenceDir, "cohortCopyNumbers.RData"))
save(cohortGeneCopyNumberAmplifications, file = paste0(referenceDir, "cohortGeneCopyNumberAmplifications.RData"))
save(cohortGeneCopyNumberDeletes, file = paste0(referenceDir, "cohortGeneCopyNumberDeletes.RData"))
save(cohortSomatics, file = paste0(referenceDir, "cohortSomatics.RData"))
save(cohortStructuralVariantSummary, file = paste0(referenceDir, "cohortStructuralVariantSummary.RData"))
save(cohortWGD, file = paste0(referenceDir, "cohortWGD.RData"))

dbDisconnect(dbProd)
rm(dbProd)

### START OF DATA PROCESSING
cohort_somatic_summary <- function(somatics) {
  result = cohort_somatics_by_type(somatics) %>%
    left_join(cohort_mutational_load(somatics), by = "sampleId") %>%
    left_join(cohort_msi(somatics), by = "sampleId")
  return (result)
}

cohort_mutational_load <- function(somatics) {
  result = somatics %>%
    filter(worstCodingEffect == "MISSENSE") %>%
    group_by(sampleId) %>%
    summarise(mutationalLoad = n())
  return (result)
}

cohort_somatics_by_type <- function(somatics) {
  result = somatics %>%
    filter(filter == 'PASS') %>%
    mutate(
      type = ifelse(type == 'SNP', 'SNV', type),
      type = ifelse(type == 'MNP', 'MNV', type),
      clonality = ifelse(clonality == 'INCONSISTENT', 'CLONAL', clonality)) %>%
    group_by(sampleId, type, clonality) %>%
    summarise(count = n()) %>%
    unite(type, clonality, type) %>%
    spread(type, count, fill = 0) %>%
    mutate(
      TOTAL_INDEL = CLONAL_INDEL + SUBCLONAL_INDEL,
      TOTAL_SNV = CLONAL_SNV +  SUBCLONAL_SNV,
      TOTAL_MNV = CLONAL_MNV +  SUBCLONAL_MNV) %>%
    select( -starts_with("CLONAL"))
  return (result)
}

cohort_msi <- function(somatics) {
  result = somatics %>%
    filter(filter == 'PASS', type == 'INDEL', nchar(alt) <= 50, nchar(ref) <= 50) %>%
    filter((nchar(repeatSequence) %in% c(2:4) & repeatCount >= 4) | (nchar(repeatSequence) == 1 & repeatCount >= 5 )) %>%
    group_by(sampleId) %>%
    summarise(msiScore = n() / 3095) %>%
    mutate(msiStatus = ifelse(msiScore > 4.0, "MSI","MSS"))
  return (result)
}

cohortSomaticsSummary = cohort_somatic_summary(cohortSomatics)
save(cohortSomaticsSummary, file = paste0(processedDir, "cohortSomaticsSummary.RData"))

cohortSNPSummary = cohortSomatics %>% filter(filter == 'PASS', type == 'SNP') %>% group_by(sampleId, ref, alt) %>% summarise(n = n())
save(cohortSomaticsSummary, file = paste0(processedDir, "cohortSNPSummary.RData"))

cohortMNPSummary = cohortSomatics %>% filter(filter == 'PASS', type == 'MNP') %>% group_by(sampleId, ref, alt) %>% summarise(n = n())
save(cohortMNPSummary, file = paste0(processedDir, "cohortMNPSummary.RData"))

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

cohortIndelSummary = indel_summary(cohortSomatics)
save(cohortIndelSummary, file = paste0(processedDir, "cohortIndelSummary.RData"))
