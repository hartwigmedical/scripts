library(RMySQL)
library(dplyr)
library(tidyr)
library(configr)

### Queries
query_sample<-function(dbConnect) {
  query=paste("select patient_id,sample_id,cancer_type,age_at_sample from lkcgp_sample")
  return (dbGetQuery(dbConnect, query))
}

query_purity<-function(dbConnect, purityCutoff = 0.2) {
  query = paste(
  "SELECT p.*",
  " FROM purity p",
  "WHERE qcStatus = 'PASS'",
  #    "  AND status <> 'NO_TUMOR'",
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

### Parameters to change
purityCutoff = 0.1
outputDir = "~/Documents/LKCGP_projects/RData/"
resourceDir = "/Users/marwo2/Documents/LKCGP_projects/RData/Resources/"
dbConf = read.config("~/.mysql/credentials")

### START OF REFERENCE DATA COLLECTION
referenceDir = paste0(outputDir, "reference/")
processedDir = paste0(outputDir, "processed/")
dbProd = dbConnect(MySQL(), user=dbConf$lkcgp_select$user, password=dbConf$lkcgp_select$password, dbname=dbConf$lkcgp_select$database, host=dbConf$lkcgp_select$host)

samples = query_sample(dbProd)
save(samples, file = paste0(referenceDir, "samples.RData"))

cohortTertPromoters = purple::query_tert_promoters(dbProd, cohort)
save(cohortTertPromoters, file = paste0(referenceDir, "cohortTertPromoters.RData"))

cohort = query_purity(dbProd, purityCutoff)
save(cohort, file = paste0(referenceDir, "cohort.RData"))

cohortCopyNumbers = query_copy_number(dbProd, cohort)
save(cohortCopyNumbers, file = paste0(referenceDir, "cohortCopyNumbers.RData"))

cohortGeneCopyNumberAmplifications = query_gene_copy_number_amplifications(dbProd, cohort)
save(cohortGeneCopyNumberAmplifications, file = paste0(referenceDir, "cohortGeneCopyNumberAmplifications.RData"))

cohortGeneCopyNumberDeletes = query_gene_copy_number_deletes(dbProd, cohort)
save(cohortGeneCopyNumberDeletes, file = paste0(referenceDir, "cohortGeneCopyNumberDeletes.RData"))

cohortSomatics = query_somatic_variants(dbProd, cohort)
save(cohortSomatics, file = paste0(referenceDir, "cohortSomatics.RData"))

cohortStructuralVariantSummary = query_structural_variant_summary(dbProd, cohort)
save(cohortStructuralVariantSummary, file = paste0(referenceDir, "cohortStructuralVariantSummary.RData"))

cohortWGD = query_whole_genome_duplication(dbProd, cohort)
save(cohortWGD, file = paste0(referenceDir, "cohortWGD.RData"))

cohortFusions = purple::query_fusions(dbProd, cohort)
save(cohortFusions, file = paste0(referenceDir, "cohortFusions.RData"))

dbEnsemble = dbConnect(MySQL(), dbname='homo_sapiens_core_89_37', groups="RAnalysis")
cohortFusionCodingRegions = purple::query_coding_regions(dbEnsemble, unique(c(cohortFusions$`5pTranscript`, cohortFusions$`3pTranscript`)))
save(cohortFusionCodingRegions, file = paste0(referenceDir, "cohortFusionCodingRegions.RData"))
dbDisconnect(dbEnsemble)
rm(dbEnsemble)

canonicalTranscripts = purple::query_canonical_transcript(dbProd)
save(canonicalTranscripts, file = paste0(referenceDir, "canonicalTranscripts.RData"))

dbDisconnect(dbProd)
rm(dbProd)

### START OF DATA PROCESSING
exonic_somatics <- function(somatics, gr_genes) {
  gr_muts = GRanges(somatics$chromosome, IRanges(somatics$position,somatics$position + nchar(somatics$ref) - 1))
  ol = as.matrix(findOverlaps(gr_muts, gr_genes, type="any", select="all"))
  return (somatics[unique(ol[, 1]), ])
}

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
save(cohortSNPSummary, file = paste0(processedDir, "cohortSNPSummary.RData"))

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

load(paste0(resourceDir, "HmfRefCDS.RData"))
cohortExonicSomatics = exonic_somatics(cohortSomatics, gr_genes)
save(cohortExonicSomatics, file = paste0(processedDir, "cohortExonicSomatics.RData"))
