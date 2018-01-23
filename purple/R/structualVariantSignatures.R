
sv_signature_by_scope <- function(variants) {
  DT = data.table(variants)
  sampleIds = unique(DT$sampleId)

  ## scope
  DT$scope <- DT$sampleId
  DT[DT[, .I[.N > 1], by=.(startChromosome, endChromosome, startPosition, endPosition, startOrientation, endOrientation, type)]$V1, ]$scope <- "Shared"

  ## Length
  DT$length <- pmax(DT$startPosition, DT$endPosition) - pmin(DT$startPosition, DT$endPosition)
  DT$length <- ifelse(DT$type == "BND", 0, DT$length)
  DT$length <- ifelse(DT$type == "INS", 0, DT$length)
  DT$length <- cut(DT$length, right = FALSE, breaks = c(-Inf, 1, 1000, 10000, 1e+05, 1e+06, Inf), labels = length_buckets())

  result = dcast(DT, type + length ~ scope, value.var = "sampleId", fun.aggregate = length)
  result$type = paste(result$type, result$length, sep="_")

  result = result[, -c("length")]
  result = merge(create_empty_sv_signature(), result, all.x = TRUE)
  result[is.na(result)] <- 0
  return (result)
}

length_buckets<-function() {
  return (c("0", "<1k", "1k-10k", "10k-100k", "100k-1M", ">1M"))
}

create_empty_sv_signature<-function() {
  buckets = length_buckets()[-1]
  labels = c("BND_0", paste("DEL", buckets, sep = "_"), paste("DUP", buckets, sep = "_"), "INS_0", paste("INV", buckets, sep = "_"))
  return (data.frame(type = labels, stringsAsFactors = FALSE))
}
