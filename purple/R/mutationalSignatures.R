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


create_empty_signature <- function() {
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

signature_matrix_by_scope <- function(variants) {

  empty = create_empty_signature()
  DT = data.table(variants)

  sampleIds = unique(DT$sample)
  DT$scope <- paste("Private", match(DT$sample, sampleIds), sep ="")

  DT[DT[, .I[.N > 1], by=.(chromosome, position, type)]$V1, ]$scope <- "Shared"
  variantsByScope = dcast(DT, type + context ~ scope, value.var = "scope", fun.aggregate = length)

  result = merge(empty, variantsByScope, all.x=TRUE)
  return (result)
}
