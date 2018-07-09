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
  sp_url = "https://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
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

standard_mutation <- function(types) {
  types = gsub("G>T", "C>A", types)
  types = gsub("G>C", "C>G", types)
  types = gsub("G>A", "C>T", types)
  types = gsub("A>T", "T>A", types)
  types = gsub("A>G", "T>C", types)
  types = gsub("A>C", "T>G", types)
  return(types)
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

extract_mutational_signature_data <- function(somatics) {
  result = somatics %>%
    filter(type == 'SNP', !grepl('N', trinucleotideContext)) %>%
    mutate(
      raw_types = paste(ref, alt, sep = ">"),
      type = standard_mutation(raw_types),
      context = standard_context(raw_types, type, trinucleotideContext)) %>%
    select(-raw_types)
  return (result)
}

mutational_signature <- function(somaticVariants) {
  signature_data = extract_mutational_signature_data(somaticVariants)
  signature = signature_data %>% group_by(type, context) %>% count()
  return (ensure_96_types(signature))
}

mutational_signature_by_scope <- function(somaticVariants) {
  signature_data = extract_mutational_signature_data(somaticVariants)
  signature = signature_data %>% group_by(type, context, scope) %>% count() %>% spread(scope, n, fill = 0)
  return (ensure_96_types(signature))
}

mutational_signature_by_clonality <- function(somaticVariants) {
  signature_data = extract_mutational_signature_data(somaticVariants)
  signature = signature_data %>% group_by(type, context, clonality) %>% count() %>% spread(scope, n, fill = 0)
  return (ensure_96_types(signature))
}

ensure_96_types <- function(signature) {
  empty = create_empty_mutational_signature()
  result = merge(empty, signature, all.x=TRUE)
  result[is.na(result)] <- 0
  return (result)
}

cosmicSignatures = getCOSMICSignatures()
cosmicSignatureColours = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
                           "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a",
                           "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
                           "#dea185","#a0729d","#8a392f")




