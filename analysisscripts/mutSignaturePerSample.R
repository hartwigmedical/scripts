#!/usr/bin/Rscript

library(devtools) #; install_github("im3sanger/dndscv")
library(MutationalPatterns)
library(RMySQL)
library(data.table)
library("NMF")

# plotting
library(grid)
library(gridExtra)
library(ggplot2)

myCOLORS = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
             "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a",
             "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
             "#dea185","#a0729d","#8a392f")


getCOSMICSignatures <- function() {
  sp_url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
  cancer_signatures = read.table(sp_url, sep = "\t", header = T)
  # reorder (to make the order of the trinucleotide changes the same)
  cancer_signatures = cancer_signatures[order(cancer_signatures[, 1]),]
  # only signatures in matrix
  cancer_signatures = as.matrix(cancer_signatures[, 4:33])
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

# warning: this takes ages!
query_sample_variants <- function(dbConnect,sampleId) {
  query = paste(
  "SELECT s.sampleId, trinucleotideContext as context, concat(ref,'>', alt) as snv, adjustedVaf * copyNumber as ploidy ",
    " FROM somaticVariant s inner join purity p on p.sampleId = s.sampleId ",
    " WHERE filter = 'PASS' and length(alt) = length(ref) and length(alt) = 1 and trinucleotideContext not like '%N%' ",
    #" and qcStatus = 'PASS' and status <> 'NO_TUMOR' and p.sampleId ='",
    " and p.sampleId ='",
    sampleId,
    "'",
    sep = "")
  print(query)
  raw_data = dbGetQuery(dbConnect, query)
  raw_types = raw_data$snv
  standard_types = standard_mutation(raw_types)
  raw_context = raw_data$context
  context = standard_context(raw_types, standard_types, raw_context)

  DT = data.table(
    sample = raw_data$sampleId,
    type = standard_types,
    context = context,
    ploidy = raw_data$ploidy,
    clonality = raw_data$clonality)

  return(DT)
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

process_variants <- function(variants) {
  samples = unique(variants$sample)
  empty = create_empty_signature()

  result = list()
  for (s in samples) {

    # slice for our variants
    sample_variants = variants[sample == s]

    # TODO: do we want to ignore unknown clonality??
    total = sample_variants[clonality != 'UNKNOWN', .(total = .N), keyby = .(type, context)]
    subclonal = sample_variants[clonality == 'SUBCLONAL', .(subclonal = .N), keyby = .(type, context)]
    clonal = sample_variants[clonality == 'CLONAL', .(clonal = .N), keyby = .(type, context)]
    clonalA = sample_variants[clonality == 'CLONAL' & ploidy < 1.5, .(clonalLowPloidy = .N), keyby = .(type, context)]
    clonalB = sample_variants[clonality == 'CLONAL' & ploidy >= 1.5, .(clonalHighPloidy = .N), keyby = .(type, context)]

    # cleanup
    rm(sample_variants)

    tmp = merge(empty, total, all=TRUE)
    tmp = merge(tmp, subclonal, all=TRUE)
    tmp = merge(tmp, clonal, all=TRUE)
    tmp = merge(tmp, clonalA, all=TRUE)
    tmp = merge(tmp, clonalB, all=TRUE)
    tmp[is.na(tmp)] <- 0 # TODO check this works
    stopifnot(nrow(tmp) == 96)

    result[[s]] = tmp
  }

  return(result)
}

##### START PROCESSING
sampleIds = commandArgs(trailingOnly=TRUE)
cancer_signatures = getCOSMICSignatures()

## the DB params are set from /usr/lib64/R/etc/Rprofile.site
dbConnect = dbConnect(MySQL(), user=db_user, password=db_password, dbname = db_name, groups = "RAnalysis")

## run analysis per sampleId
for ( sampleId in sampleIds ){
    
    cat( paste( "[INFO] Creating mutational signature for: ", sampleId, sep=""), "\n" )
    variants = query_sample_variants(dbConnect,sampleId) # returns a DT
    mutation_vectors = process_variants(variants)
    # we need to slice out only the mutation count columns (delete col 1 and 2)
    signatures = fit_to_signatures(mutation_vectors[[sampleId]][, -c(1, 2)], cancer_signatures)$contribution    

    ## table output
    txtFileName = paste( sampleId, "_mutSig.tsv", sep="" )
    write.table(data.frame("Signature" = rownames(signatures), signatures), txtFileName, row.names=FALSE, quote=FALSE, sep="\t")

    ## plotting    
    pdfFileName = paste( sampleId, "_mutSig.pdf", sep="" )
    pdf( pdfFileName )
    print( 
      plot_contribution(signatures,cancerSignatures)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10),legend.text=element_text(size=5),axis.title.y = element_text(size=10))+
      scale_fill_manual( values= myCOLORS)+
      labs(fill="")+
      ggtitle(paste("Mutational signatures by clonality for",sampleId)) 
    )
    dev.off()

    cat( paste( "[INFO] Output in: ", txtFileName, " and ", pdfFileName, sep=""), "\n" )
}

## finish up
dbDisconnect(dbConnect)

