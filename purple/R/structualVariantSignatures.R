
sv_signature_by_scope <- function(variants) {
  DT = data.table(variants)

  ## Length
  DT$length <- pmax(DT$startPosition, DT$endPosition) - pmin(DT$startPosition, DT$endPosition)
  DT$length <- ifelse(DT$type == "BND", 0, DT$length)
  DT$length <- ifelse(DT$type == "INS", 0, DT$length)
  DT$length <- cut(DT$length, right = FALSE, breaks = c(-Inf, 1, 1000, 10000, 1e+05, 1e+06, Inf), labels = sv_length_buckets())
  DT$type <- factor(paste(DT$type, DT$length, sep="_"), levels=sv_type_length_buckets(), ordered = TRUE)

  result = dcast(DT, type ~ scope, value.var = "sampleId", fun.aggregate = length)
  result = merge(create_empty_sv_signature(), result, all.x = TRUE)
  result[is.na(result)] <- 0
  return (result)
}

sv_length_buckets<-function() {
  return (c("0", "<1k", "1k-10k", "10k-100k", "100k-1M", ">1M"))
}

sv_type_length_buckets<-function() {
  buckets = length_buckets()[-1]
  labels = c("INS_0","BND_0", paste("DEL", buckets, sep = "_"), paste("DUP", buckets, sep = "_"), paste("INV", buckets, sep = "_"))
  return (labels)
}

create_empty_sv_signature<-function() {
  buckets = sv_type_length_buckets()
  type = factor(buckets,levels=buckets,ordered=TRUE)
  return (data.frame(type))
}

