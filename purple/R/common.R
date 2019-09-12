redGradient = c("#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15")
greenGradient = c("#edf8e9", "#bae4b3", "#74c476", "#31a354", "#006d2c")
blackGradient = c("#f7f7f7", "#cccccc", "#969696", "#636363", "#252525")

slookup<-function(df, lookup, on="sampleId") {
  return (sapply(df[[on]], function(x) {lookup[match(x, lookup[[on]]), -c(1)] }))
}

slookup2<-function(v, lookup, on="sampleId") {
  return (sapply(v, function(x) {lookup[match(x, lookup[[on]]), -c(1)] }))
}

apply_to_cohort<-function(cohort, sampleFunction) {
  result = ddply(cohort, 1, sampleFunction, .drop = FALSE, .progress = "text", .parallel = FALSE)
  return (result)
}

sample_to_patient_id<-function(sampleId, lookup) {
  colnames(lookup) <- c("truncatedSampleIds", "patientIds")

  # TODO: THIS IS OUT OF DATE

  index = match(substr(sampleId, 1, 12) , substr(lookup[[1]], 1, 12))
  if (is.na(index)) {
    substr(sampleId, 1, 12)
  } else {
    lookup[index, c(2)]
  }
}
