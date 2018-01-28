
plot_indel_signature<-function(indelSignature) {
  require(MutationalPatterns)
  require(ggplot2)

  indelSignatureColours = c(rev(redGradient), greenGradient);
  p1 <- plot_absolute_contribution(indelSignature)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10),legend.text=element_text(size=5),axis.title.y = element_text(size=10))+
    ggtitle("Indel Signatures")+
    scale_fill_manual(name = "", values =indelSignatureColours)
    #scale_fill_manual(name = "", values =indelSignatureColours, guide=guide_legend(ncol = 2))
  return (p1)
}


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
