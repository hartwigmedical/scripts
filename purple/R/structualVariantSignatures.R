plot_sv_signature<-function(svSignature) {
  require(MutationalPatterns)
  require(ggplot2)

  insYellow = "#e3d200"
  bndBlue = "#6baed6"
  svSignatureColours = c(insYellow, bndBlue, redGradient, greenGradient, blackGradient)
  p1 <- plot_absolute_contribution(svSignature) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10),legend.text=element_text(size=5),axis.title.y = element_text(size=10))+
    ggtitle("Structural Variant Signatures")+
    scale_fill_manual(name = "", values = svSignatureColours, drop = FALSE)

  return (p1)
}

sv_signature_by_scope <- function(variants) {
  DT = sv_signature_data(variants)

  result = dcast(DT, type ~ scope, value.var = "sampleId", fun.aggregate = length)
  result = merge(create_empty_sv_signature(), result, all.x = TRUE)
  result[is.na(result)] <- 0
  return (result)
}

sv_signature_by_clonality <- function(variants) {
  DT = sv_signature_data(variants)

  result = DT[, .(sample=.N), by=type]
  result = merge(create_empty_sv_signature(), result, all.x = TRUE)
  result[is.na(result)] <- 0
  return (result)
}

sv_signature_data <- function(variants) {
  DT = data.table(variants)

  ## Length
  DT$length <- pmax(DT$startPosition, DT$endPosition) - pmin(DT$startPosition, DT$endPosition)
  DT$length <- ifelse(DT$type == "BND", 0, DT$length)
  DT$length <- ifelse(DT$type == "INS", 0, DT$length)
  DT$length <- cut(DT$length, right = FALSE, breaks = c(-Inf, 1, 1000, 10000, 1e+05, 1e+06, Inf), labels = sv_length_buckets())
  DT$type <- factor(paste(DT$type, DT$length, sep="_"), levels=sv_type_length_buckets(), ordered = TRUE)

  return (DT)
}

sv_length_buckets<-function() {
  return (c("0", "<1k", "1k-10k", "10k-100k", "100k-1M", ">1M"))
}

sv_type_length_buckets<-function() {
  buckets = sv_length_buckets()[-1]
  labels = c("INS_0","BND_0", paste("DEL", buckets, sep = "_"), paste("DUP", buckets, sep = "_"), paste("INV", buckets, sep = "_"))
  return (labels)
}

create_empty_sv_signature<-function() {
  buckets = sv_type_length_buckets()
  type = factor(buckets,levels=buckets,ordered=TRUE)
  return (data.frame(type))
}

