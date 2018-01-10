left_join<-function(left, right, by="sampleId") {
  tmp = merge(x = left, y = right, by=by, all.x=TRUE)
  return (tmp)
}

slookup<-function(df, lookup, on="sampleId") {
  return (sapply(df[[on]], function(x) {lookup[match(x, lookup[[on]]), -c(1)] }))
}

apply_to_cohort<-function(cohort, sampleFunction) {
  result = ddply(cohort, 1, sampleFunction, .drop = FALSE, .progress = "text", .parallel = FALSE)
  return (result)
}

sample_to_patient_id<-function(sampleId, lookup) {
  index = match(sampleId, lookup[[1]])
  if (is.na(index)) {
    substr(sampleId, 1, 12)
  } else {
    lookup[index, c(2)]
  }
}

manual_patient_id<-function() {
  sampleId  = c("CPCT02020192T", "CPCT02030224T", "DRUP01010044T", "DRUP01070024T", "DRUP01050008T")
  patientId = c("CPCT02020438",  "CPCT02030292",  "DRUP01010044",  "CPCT02070110",  "CPCT02050116")
  return (data.frame(sampleId, patientId, stringsAsFactors = FALSE))
}


#this is a slow function - don't use!
dflookup<-function(df, lookup, on="sampleId") {
  return(as.data.frame(t(sapply(df[[on]], function(x) {lookup[match(x, lookup[[on]]), -c(1)] }))))
}

