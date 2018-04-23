highest_purity_patients<-function(cohort) {
  dt = data.table(cohort)
  return (cohort[dt[, .I[which.max(purity)], by=patientId]$V1, ])
}

multiple_biopsy<-function(cohort) {
  dt = data.table(cohort)
  return (cohort[dt[, .I[.N > 1], by=patientId]$V1, ])
}

cohort <-function(cohortRawData) {
  cohort = cohortRawData$purity[, c("sampleId","gender", "version", "purity", "ploidy")]

  # TODO: Investigate this sample for somatics
  cohort = cohort[cohort$sampleId != "CPCT02010063T", ]

  cohort = left_join(cohort, cohortRawData$clinical[, c("sampleId", "cancerType", "cancerSubtype", "biopsySite", "biopsyLocation", "treatment", "treatmentType", "ageAtBiopsyDate")])
  cohort = left_join(cohort, cohortRawData$sample[, c("sampleId", "pathologyPurity")])
  cohort = left_join(cohort, cohortRawData$msi)
  cohort = left_join(cohort, cohortRawData$svTypes)
  cohort = left_join(cohort, cohortRawData$somatics)
  cohort = left_join(cohort, cohortRawData$patientIds)
  cohort = left_join(cohort, cohortRawData$unsupportedSegments)
  cohort = left_join(cohort, cohortRawData$wholeGenomeDuplications)

  # Enrich
  cohort$qcScore <- round(cohort$unsupportedSegments / cohort$ploidy,0)
  cohort$unsupportedSegments <- NULL
  cohort$duplicatedAutosomes <- ifelse(is.na(cohort$duplicatedAutosomes), 0, cohort$duplicatedAutosomes)
  cohort$WGD <- ifelse(is.na(cohort$WGD), FALSE, cohort$WGD)

  # Clean
  cohort["_CLONAL"] <- NULL
  cohort["_SUBCLONAL"] <- NULL
  cohort["_UNKNOWN"] <- NULL

  return (cohort)
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
    group_by(sampleId, type) %>%
    summarise(count = n()) %>%
    spread(type, count)

  return (result)
}

cohort_msi <- function(somatics) {
  result = somatics %>%
    filter(filter == 'PASS', type == 'INDEL', nchar(alt) <= 50, nchar(ref) <= 50, repeatCount >= 5) %>%
    filter(nchar(repeatSequence) %in% c(2:4) | (nchar(repeatSequence) == 1 & repeatCount >= 5 )) %>%
    group_by(sampleId) %>%
    summarise(msiScore = n() / 3095) %>%
    mutate(mscStatus = ifelse(msiScore > 0.909, "MSI","MSS"))

  return (result)
}



cohort_raw_data<-function(dbConnect, limit = 0) {
  # Query Cohort Data
  purity = query_purity(dbConnect)
  sample = query_sample_data(dbConnect)
  clinical = query_clinical_data(dbConnect)
  clinical$ageAtBiopsyDate = year(clinical$biopsyDate) - clinical$birthYear

  cohort = purity[, c("sampleId","gender", "purity", "ploidy")]
  if (limit != 0) {
    cohort = cohort[1:limit, ]
  }

  cat("Querying microsatellite instability", "\n")
  msi = purple::apply_to_cohort(cohort, function(x) {purple::query_msi_sample(dbConnect, x$sampleId)})

  cat("Querying structural variants", "\n")
  svTypes = purple::apply_to_cohort(cohort, function(x) {purple::query_structural_variant_overview(dbConnect, x$sampleId)})

  cat("Querying unsupported copy number segments", "\n")
  unsupportedSegments = purple::apply_to_cohort(cohort, function(x) {purple::query_unsupported_segments(dbConnect, x$sampleId)})

  cat("Querying somatic variants", "\n")
  somatics = purple::apply_to_cohort(cohort, function(x) {purple::query_somatic_overview(dbConnect, x$sampleId)})

  cat("Querying patientIds", "\n")
  patientIdLookups = query_patient_id_lookup(dbConnect)
  patientIds = purple::apply_to_cohort(cohort, function(x) {purple::sample_to_patient_id(x$sampleId, patientIdLookups)})
  colnames(patientIds) <- c("sampleId", "patientId")

  cat("Querying duplicated autosomes", "\n")
  wholeGenomeDuplications = purple::apply_to_cohort(cohort, function(x) {purple::query_whole_genome_duplication(dbConnect, x$sampleId)})

  result = list()
  result$purity <- purity
  result$sample <- sample
  result$clinical <-clinical
  result$msi <- msi
  result$svTypes <- svTypes
  result$unsupportedSegments <- unsupportedSegments
  result$somatics <- somatics
  result$patientIds <- patientIds
  result$wholeGenomeDuplications <- wholeGenomeDuplications

  return(result)
}
