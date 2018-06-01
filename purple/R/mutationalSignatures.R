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



base_complements <- function(bases) {
  complements = setNames(c("A", "C", "G","T"), c("T", "G", "C","A"))

  point_complement <- function(base) {
    paste(rev(sapply(strsplit(base, split = ""), function (x) {complements[x]})), collapse = "")
  }

  sapply(bases, point_complement)
}


standard_double_mutation <- function(types) {

  single_standard_mutation <- function(type) {
    if (type %in% c("CC>AA","GG>TT")) {
      return ("CC>AA")
    }

    if (type %in% c("CC>TT","GG>AA")) {
      return ("CC>TT")
    }

    if (substr(type, 3, 3) != '>' | nchar(type) != 5) {
      return ("Other")
    }

    if (substr(type, 1, 2) %in% c("AC","GT")) {return ("AC>NN")}
    if (substr(type, 1, 2) %in% c("AT","AT")) {return ("AT>NN")}
    if (substr(type, 1, 2) %in% c("CC","GG")) {return ("CC>NN")}
    if (substr(type, 1, 2) %in% c("CG","CG")) {return ("CG>NN")}
    if (substr(type, 1, 2) %in% c("CT","AG")) {return ("CT>NN")}
    if (substr(type, 1, 2) %in% c("GC","GC")) {return ("GC>NN")}
    if (substr(type, 1, 2) %in% c("TA","TA")) {return ("TA>NN")}
    if (substr(type, 1, 2) %in% c("TC","GA")) {return ("TC>NN")}
    if (substr(type, 1, 2) %in% c("TG","CA")) {return ("TG>NN")}
    if (substr(type, 1, 2) %in% c("TT","AA")) {return ("TT>NN")}
  }

  sapply(types, single_standard_mutation)
}

#single_standard_mutation("AC>AC")
#standard_double_mutation(c("CC>TA", "AT>AC"))

standard_double_mutation3 <- function(types) {

  lookup = data.frame(base = c("AC", "AT", "CC", "CG", "CT", "GC", "TA","TC", "TG", "TT"), stringsAsFactors = F)
  lookup$complement = base_complements(lookup$base)

  df = data.frame(type = types, stringsAsFactors = F)
  df = df %>% mutate(ref = substr(type, 1, 2), refComplement = base_complements(ref)) %>%
    left_join(lookup %>% select(ref = base, refMatch = complement), by = "ref") %>%
    left_join(lookup %>% select(refComplement = base, refComplementMatch = complement), by = "refComplement") %>%
    mutate(
      match = coalesce(refMatch, refComplementMatch),
      match = ifelse(is.na(match), "Other", paste0(match, ">NN")))

  df$match = ifelse(nchar(types) != 5, "Other", df$match)
  df$match = ifelse(types %in% c("CC>AA", "GG>TT"), "CC>AA", df$match)
  df$match = ifelse(types %in% c("CC>TT", "GG>AA"), "CC>TT", df$match)

  return(df$match)
}



standard_double_mutation2 <- function(types) {

  isMatch <- function(type, base) {
    return (type == base | type == base_complements(base))
  }

  result = ifelse(substr(types, 3, 3) != '>', "Other", NA)

  result = ifelse(is.na(result) & isMatch(substr(types, 1, 2), "CC"), "CC>Other", NA)
  result = ifelse(is.na(result) & types %in% c("CC>AA", "GG>TT"), "CC>AA", NA)
  result = ifelse(is.na(result) & types %in% c("CC>TT", "GG>AA"), "CC>TT", NA)

  result = ifelse(is.na(result) & isMatch(substr(types, 1, 2), "AC"), "AC>NN", NA)
  result = ifelse(is.na(result) & isMatch(substr(types, 1, 2), "AT"), "AT>NN", NA)

  result = ifelse(is.na(result) & isMatch(substr(types, 1, 2), "CC"), "CC>NN", NA)
  result = ifelse(is.na(result) & isMatch(substr(types, 1, 2), "CG"), "CG>NN", NA)
  result = ifelse(is.na(result) & isMatch(substr(types, 1, 2), "CT"), "CT>NN", NA)

  result = ifelse(is.na(result) & isMatch(substr(types, 1, 2), "GC"), "GC>NN", NA)

  result = ifelse(is.na(result) & isMatch(substr(types, 1, 2), "TA"), "TA>NN", NA)
  result = ifelse(is.na(result) & isMatch(substr(types, 1, 2), "TC"), "TC>NN", NA)
  result = ifelse(is.na(result) & isMatch(substr(types, 1, 2), "TG"), "TG>NN", NA)
  result = ifelse(is.na(result) & isMatch(substr(types, 1, 2), "TT"), "TT>NN", NA)

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




