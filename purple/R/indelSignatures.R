indel_signature_by_scope <- function(somaticVariants) {

  ## scope
  DT = data.table(somaticVariants[somaticVariants$type == "INDEL", ])

  ## Length
  DT$type <- nchar(DT$alt) - nchar(DT$ref)
  DT$type <- ifelse(DT$type > 5, 5, DT$type)
  DT$type <- ifelse(DT$type < -5, -5, DT$type)
  DT$type <- factor(DT$type, levels=indel_length_buckets(), ordered = TRUE)

  empty = create_empty_indel_signature()
  result = dcast(DT, type ~ scope, value.var = "sampleId", fun.aggregate = length)
  result = merge(create_empty_indel_signature(), result, all.x = TRUE)
  result[is.na(result)] <- 0
  return (result)
}

indel_length_buckets<-function() {
  buckets = c("-5","-4","-3","-2","-1", "1", "2", "3", "4", "5")
  return (buckets)
}

create_empty_indel_signature<-function() {
  buckets = indel_length_buckets()
  type = factor(buckets,levels=buckets,ordered=TRUE)
  return (data.frame(type))
}
