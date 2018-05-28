plot_cosmic_signature<-function(mutationalSignature, mode = "absolute") {
  require(MutationalPatterns)
  require(ggplot2)

  contribution = fit_to_signatures(mutationalSignature[, -c(1, 2)], cosmicSignatures)$contribution
  p1 <- plot_contribution(contribution, cosmicSignatures, mode = mode)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10),legend.text=element_text(size=5),axis.title.y = element_text(size=10))+
    scale_fill_manual( values= cosmicSignatureColours)+labs(fill="")+ggtitle("Mutational Signatures") +
    labs(x = "", y = "Absolute contribution")

  return (p1)
}

getCOSMICSignatures <- function() {
  sp_url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
  cancer_signatures = read.table(sp_url, sep = "\t", header = T)
  # reorder (to make the order of the trinucleotide changes the same)
  cancer_signatures = cancer_signatures[order(cancer_signatures[, 1]),]
  # only signatures in matrix
  cancer_signatures = as.matrix(cancer_signatures[, 4:33])
}

standard_double_mutation <- function(types) {
  result = ifelse(types %in% c("CC>AA", "GG>TT"), "CC>AA", NA)
  result = ifelse(is.na(result) & types %in% c("CC>TT", "GG>AA"), "CC>TT", result)
  result = ifelse(is.na(result) & substr(types, 1, 2) %in% c("CC", "GGT"), "CC>Other", result)

  result = ifelse(is.na(result) & substr(types, 1, 2) %in% c("TC", "AG"), "TC", result)
  result = ifelse(is.na(result) & substr(types, 1, 2) %in% c("TT", "AA"), "TT", result)
  result = ifelse(is.na(result) & substr(types, 1, 2) %in% c("CT", "GA"), "CT", result)
  result = ifelse(is.na(result) & substr(types, 1, 2) == "TG", "TG", result)
  result = ifelse(is.na(result) & substr(types, 1, 2) == "AC", "AC", result)
  result = ifelse(is.na(result) & substr(types, 1, 2) %in% c("GC", "CG"), "GC", result)
  result = ifelse(is.na(result), "Other", result)

  return(result)
}

standard_mutation <- function(types) {
  types = gsub("G>T", "C>A", types)
  types = gsub("G>C", "C>G", types)
  types = gsub("G>A", "C>T", types)
  types = gsub("A>T", "T>A", types)
  types = gsub("A>G", "T>C", types)
  types = gsub("A>C", "T>G", types)
  return(types)
}

standard_context <- function(raw_type, standard_type, context) {
  x = which(raw_type != standard_type)
  context[x] = reverse(chartr("ATGC", "TACG", context[x]))
  return(context)
}


create_empty_mutational_signature <- function() {
  DF <- data.frame(type = character(), context = character(), stringsAsFactors = FALSE)
  ref_bases = c("C", "T")
  bases = c("A", "C", "G", "T")
  for (ref in ref_bases) {
    for (alt in bases) {
      if (alt != ref) {
        type = paste(ref, alt, sep = ">")
        for (before in bases) {
          for (after in bases) {
            context = paste(before, after, sep = ref)
            DF = rbind(DF, data.frame(type, context, stringsAsFactors = FALSE))
          }
        }
      }
    }
  }
  return(DF)
}

extract_mutational_signature_data <- function(somaticVariants, feature="scope") {
  raw_data = somaticVariants[somaticVariants$type == 'SNP', ]
  raw_data = raw_data[!(raw_data$trinucleotideContext %in% 'N'), ]

  raw_types <- paste(raw_data$ref, raw_data$alt, sep = ">")
  standard_types = standard_mutation(raw_types)

  raw_context = raw_data$trinucleotideContext
  context = standard_context(raw_types, standard_types, raw_context)

  DT = data.table(
    sampleId = raw_data$sampleId,
    type = standard_types,
    context = context,
    ploidy = raw_data$ploidy,
    clonality = raw_data$clonality,
    chromosome = raw_data$chromosome,
    position = raw_data$position,
    scope = raw_data[[feature]])

  return (DT)
}

mutational_signature_by_clonality <- function(somaticVariants) {
  empty = create_empty_mutational_signature()
  DT = extract_mutational_signature_data(somaticVariants, feature="clonality")

  result = dcast(DT, type + context ~ clonality, value.var = "clonality", fun.aggregate = length)
  result = merge(empty, result, all.x=TRUE)
  result[is.na(result)] <- 0

  return (result)
}

mutational_signature_by_scope <- function(somaticVariants) {
  empty = create_empty_mutational_signature()
  DT = extract_mutational_signature_data(somaticVariants, feature="scope")

  result = dcast(DT, type + context ~ scope, value.var = "scope", fun.aggregate = length)
  result = merge(empty, result, all.x=TRUE)
  result[is.na(result)] <- 0

  return (result)
}

cosmicSignatures = getCOSMICSignatures()
cosmicSignatureColours = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
                           "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a",
                           "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
                           "#dea185","#a0729d","#8a392f")




