library(devtools) #; install_github("im3sanger/dndscv")
library(MutationalPatterns)
library(RMySQL)
library(data.table)
library("NMF")

# plotting
library(grid)
library(gridExtra)
library(ggplot2)

getCOSMICSignatures <- function() {
    sp_url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
    cancer_signatures = read.table(sp_url, sep = "\t", header = T)
    # reorder (to make the order of the trinucleotide changes the same)
    cancer_signatures = cancer_signatures[order(cancer_signatures[, 1]),]
    # only signatures in matrix
    cancer_signatures = as.matrix(cancer_signatures[, 4:33])
}

select_cohort <- function(dbConnect) {
    query = paste(
        "SELECT s.sampleId, biopsySite, cancerType, cancerSubtype",
        "FROM hmfpatients.clinical c, sample s, purity p",
        "WHERE s.sampleId = c.sampleId and s.sampleId = p.sampleId and qcStatus = 'PASS' and status <> 'NO_TUMOR'",
        sep = " ")
    DF = dbGetQuery(dbConnect, query)
    setDT(DF);
    return(DF)
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
query_variants <- function(dbConnect) {
    query = paste(
        "SELECT sampleId, trinucleotideContext as context, concat(ref,'>', alt) as snv, adjustedVaf * adjustedCopyNumber as ploidy, clonality",
        "FROM somaticVariant",
        "WHERE filter = 'PASS' and length(alt) = length(ref) and length(alt) = 1 and trinucleotideContext not like '%N%'",
        sep = " ")

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

### START -> DATA SETUP

dataFile = "~/hmf/mutSignature2.RData"

cancer_signatures = getCOSMICSignatures()

dbConnect = dbConnect(MySQL(), dbname = 'hmfpatients', groups = "RAnalysis")
cohort = select_cohort(dbConnect) # returns a DT
variants = query_variants(dbConnect) # returns a DT
dbDisconnect(dbConnect)

# list of patients -> data.table of mutation counts
mutation_vectors = process_variants(variants)

signatures = list()
for (s in cohort$sampleId) {
    if (!is.null(mutation_vectors[[s]])) {
        # we need to slice out only the mutation count columns (delete col 1 and 2)
        res = fit_to_signatures(mutation_vectors[[s]][, -c(1, 2)], cancer_signatures)
        signatures[[s]] <- res$contribution
    }
}

save(cancer_signatures,
     cohort,
     mutation_vectors,
     signatures,
     file = dataFile)
