scope <- function(variants) {
  if ("startChromosome" %in% colnames(variants)) {
    by = c("startChromosome", "endChromosome", "startPosition", "endPosition", "startOrientation", "endOrientation", "type")
  } else {
    by = c("chromosome", "position", "alt", "ref", "type")
  }
  return (add_scope_to_variants(variants, by))
}

multiple_biopsy_somatic_variants <- function(somaticVariants) {
  by = c("chromosome", "position", "alt", "ref", "type")
  return (add_scope_to_variants(somaticVariants, by))
}

multiple_biopsy_structural_variants <- function(structuralVariants) {
  by = c("startChromosome", "endChromosome", "startPosition", "endPosition", "startOrientation", "endOrientation", "type")
  return (add_scope_to_variants(structuralVariants, by))
}

add_scope_to_variants <- function(variants, by) {
  DT = data.table(variants)

  DT$scope <- DT$sampleId
  DT[DT[, .I[.N > 1], by=by]$V1, ]$scope <- "Shared"

  return (DT)
}

somatic_ploidy_plots<-function(somaticVariants) {
  by = c("chromosome", "position", "alt", "ref")
  return (ploidy_plots("Somatic Variant", somaticVariants, by))
}

structural_ploidy_plots<-function(structuralVariants) {
  by = c("startChromosome", "endChromosome", "startPosition", "endPosition", "startOrientation", "endOrientation", "type")
  return (ploidy_plots("Structural Variant", structuralVariants, by, size = 1))
}

ploidy_plots<-function(label, variants, by, size = 0.2) {
  patientSampleIds = unique(variants$sampleId)

  patient1Variants = variants[variants$scope == patientSampleIds[1], ]
  patient2Variants = variants[variants$scope == patientSampleIds[2], ]
  sharedVariants = variants[variants$scope == "Shared", ]

  p1 <- ploidy_histogram(label, patient1Variants, sharedVariants)
  p2 <- ploidy_histogram(label, patient2Variants, sharedVariants)

  p3Title = paste("Shared", label, "Ploidy", sep = " ")
  p3 <- ggplot()+labs(title=p3Title, x = patientSampleIds[1], y = patientSampleIds[2])+xlim(0,3)+ylim(0,3)
  if (nrow(sharedVariants) > 0 && length(patientSampleIds) == 2) {
    sharedPloidy = sharedVariants[, lapply(.SD$ploidy, c)[1:2], by=by,.SDcols = c("sampleId","ploidy")]
    p3 <- p3+geom_point(data=sharedPloidy,aes(V1,V2), size = size)
  }

  return (list(p1, p2, p3))
}

ploidy_histogram<-function(label, patientData, sharedData) {
  binsize<-0.05
  title = unique(patientData$sampleId)

  label = paste(label, "Ploidy", sep = " ")
  p1 <- ggplot() + xlim(0,3.5)+ labs(title = title,x=label) +
    geom_histogram(aes(x=ploidy),data=patientData,fill = "red", alpha = 0.6, binwidth = binsize) +
    geom_histogram(aes(x=ploidy),data=sharedData,fill = "blue", alpha = 0.2, binwidth = binsize)

  return (p1)
}
