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
