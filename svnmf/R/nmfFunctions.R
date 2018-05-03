convert_summary_counts_to_nmf<-function(summaryCounts) {
  matrixData1 <- summaryCounts %>% gather(CountField, CountVal, -SampleId)
  return (matrixData1 %>% spread(SampleId, CountVal))
}

get_bucket_names<-function(nmfMatrixData) {
  return (nmfMatrixData$CountField)
}

remove_bucket_names<-function(nmfMatrixData) {
  matrixData = within(nmfMatrixData, rm(CountField))
  return (matrixData)
}
