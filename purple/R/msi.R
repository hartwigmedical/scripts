msi<-function(somaticVariants) {
  somatics = data.table(somaticVariants)
  result = somatics[type =='INDEL' & repeatCount > 0 & nchar(alt) <= 50 & nchar(ref) <= 50, .N, by=sampleId]$N/3095
}





