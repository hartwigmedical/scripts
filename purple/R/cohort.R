highest_purity_cohort<-function(cohort) {
  cohort %>% group_by(patientId, purity) %>% top_n(1, -score) %>% ungroup() %>% group_by(patientId) %>% top_n(1, purity) %>% ungroup()
}

multiple_biopsy_cohort<-function(cohort) {
  cohort %>% group_by(patientId) %>% filter(n() > 1)
}


cohort_deleted_genes <- function(geneCopyNumberDeletes) {
  geneCopyNumberDeletes %>%
    filter(chromosome != 'Y', germlineHomRegions == 0, germlineHetRegions == 0) %>%
    group_by(sampleId) %>%
    summarise(deletedGenes = n())
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
      type = ifelse(type == 'MNP', 'MNV', type)) %>%
    group_by(sampleId, type) %>%
    summarize(
      SUBCLONAL = round(sum(subclonalLikelihood), 0),
      TOTAL = n()) %>%
    gather(clonality, value, SUBCLONAL, TOTAL) %>%
    unite(type, clonality, type) %>%
    spread(type, value, fill = 0)

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
