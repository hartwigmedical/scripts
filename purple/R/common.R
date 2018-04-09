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

  lookup = rbind(manual_patient_id(), lookup)

  index = match(substr(sampleId, 1, 12) , substr(lookup[[1]], 1, 12))
  if (is.na(index)) {
    substr(sampleId, 1, 12)
  } else {
    lookup[index, c(2)]
  }
}

manual_patient_id<-function() {
  truncatedSampleIds  = c("CPCT02020192", "CPCT02030224", "DRUP01010007", "DRUP01070024", "DRUP01050008",
                          "DRUP01010065", "DRUP01330002", "DRUP01340004", "DRUP01340003", "DRUP01340002", "DRUP01070008")
  patientIds = c("CPCT02020438", "CPCT02030292", "DRUP01010044", "CPCT02070110", "CPCT02050116",
                 "CPCT02010639", "CPCT02330049", "CPCT02340029", "CPCT02340014", "CPCT02340026", "CPCT02070023")
  return (data.frame(truncatedSampleIds, patientIds, stringsAsFactors = FALSE))
}

#this is a slow function - don't use!
dflookup<-function(df, lookup, on="sampleId") {
  return(as.data.frame(t(sapply(df[[on]], function(x) {lookup[match(x, lookup[[on]]), -c(1)] }))))
}

