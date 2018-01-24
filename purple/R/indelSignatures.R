indel_signature_by_scope <- function(variants) {

  DT = data.table(variants)
  sampleIds = unique(DT$sampleId)

  ## scope
  DT$scope <- DT$sampleId
  DT[DT[, .I[.N > 1], by=.(chr, pos, alt, ref)]$V1, ]$scope <- "Shared"

  ## Length
  DT$type <- nchar(DT$alt) - nchar(DT$ref)
  DT$type <- ifelse(DT$type > 5, 5, DT$type)
  DT$type <- ifelse(DT$type < -5, -5, DT$type)


  empty = create_empty_indel_signature()
  result = dcast(DT, type ~ scope, value.var = "sampleId", fun.aggregate = length)
  result = merge(create_empty_indel_signature(), result, all.x = TRUE)
  result[is.na(result)] <- 0
  result$type <- as.character(result$type)
  return (result)
}

create_empty_indel_signature<-function() {
  return (data.frame(type =  c(-5:-1, 1:5), stringsAsFactors = FALSE))
}
