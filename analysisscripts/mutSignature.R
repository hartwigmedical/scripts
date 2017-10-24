library(devtools)#; install_github("im3sanger/dndscv")
library(MutationalPatterns)
library(RMySQL)
library(data.table)
library("NMF")

standard_mutation<-function(types)
{
  types = gsub("G>T", "C>A", types)
  types = gsub("G>C", "C>G", types)
  types = gsub("G>A", "C>T", types)
  types = gsub("A>T", "T>A", types)
  types = gsub("A>G", "T>C", types)
  types = gsub("A>C", "T>G", types)
  return(types)
}

standard_context<-function(raw_type, standard_type, context) {
  x = which(raw_type != standard_type)
  y = context[x]
  y = reverse(chartr("ATGC", "TACG", y))
  context[x] = y
  return(context)
}

create_empty_signature2<-function() {
  DF <- data.frame(type=rep("", 96), context=rep("   ", 96), count=rep(0L, 96), stringsAsFactors=FALSE)
  ref_bases = c("C", "T")
  alt_bases = list(c("A", "G", "T"), c("A", "C", "G"))
  bases = c("A", "C", "G", "T")
  for (i in  0:1) {
    ref = ref_bases[i+1]
    for (j in 0:2) {
      alt = alt_bases[[i+1]][j+1]
      type = paste(ref, alt, sep = ">")
      for (x in 0:3) {
        for (y in 0:3) {
          before = bases[x+1]
          after = bases[y+1]
          index = i * 48 + j * 16 + x * 4 + y + 1
          context = paste(before, after, sep=ref)
          DF[index, ] <- list(type, context, 0)
        }
      }
    }
  }
  return (DF)
}



create_empty_signature<-function() {
  DF <- data.frame(type=character(), context=character(), count=double(), stringsAsFactors=FALSE)
  ref_bases = c("C", "T")
  bases = c("A", "C", "G", "T")
  for (ref in ref_bases) {
    for (alt in bases) {
      if (alt != ref) {
        type = paste(ref, alt, sep = ">")
        for (before in bases) {
          for (after in bases) {
            context = paste(before, after, sep=ref)
            DF = rbind(DF, data.frame(type, context, count=0, stringsAsFactors=FALSE))
          }
        }
      }
    }
  }
  return(DF)
}

select_cohort<-function(dbConnect, type)
{
  query = paste(
    "select s.sampleId from patient p, sample s where s.patientId = p.id and upper(primaryTumorLocation) like '%",
    type,
    "%'",
    sep = "")
  return (dbGetQuery(dbConnect, query))
}

sample_signature<-function(dbConnect, sample, empty_signature)
{
  query = paste(
    "select sampleId, trinucleotideContext as context, concat(ref,'>', alt) as snv, count(*) as count from somaticVariant where sampleId = '",
    sample,
    "' and length(alt) = length(ref) and length(alt) = 1 group by 1, 2, 3",
    sep = "")

  raw_data = dbGetQuery(dbConnect, query)
  raw_types = raw_data$snv
  standard_types = standard_mutation(raw_types)
  raw_context = raw_data$context
  standard_context = standard_context(raw_types, standard_types, raw_context)

  DT=data.table(type=standard_types, context=standard_context, count=raw_data$count)
  DT=rbind(DT, empty_signature)

  grouped=DT[, .(count=sum(count)), keyby=.(type,context)]
  return(grouped)
}

cohort_signature<-function(dbConnect, cohort)
{
  empty_signature = create_empty_signature()
  result = matrix(, nrow = 96, ncol = 0)
  for (sample in cohort) {
    signature = sample_signature(dbConnect, sample, empty_signature)
    total_count = sum(signature$count)
    if (total_count > 0)
    {
      print(paste("Processing sample", sample))
      result = cbind(result, signature$count)
      rownames(result) <- signature$context
      colnames(result)[ncol(result)] <- sample
    }
    else
    {
      print(paste("No mutations available in sample", sample))
    }
  }
  return(result)
}

getCOSMICSignatures<-function(){
  sp_url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
  cancer_signatures = read.table(sp_url, sep = "\t", header = T)
  # reorder (to make the order of the trinucleotide changes the same)
  cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]
  # only signatures in matrix
  cancer_signatures = as.matrix(cancer_signatures[,4:33])
}

# Create mutational matrix
dbConnect = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
cohort = select_cohort(dbConnect, "BREAST")
mutation_matrix = cohort_signature(dbConnect, cohort[1:40,])
dbDisconnect(dbConnect)
plot_96_profile(mutation_matrix)

# Fit to cosmic signatures
cancer_signatures = getCOSMICSignatures()
fit_res = fit_to_signatures(mutation_matrix, cancer_signatures)
fit_res$contribution<-fit_res$contribution[, order(colSums(fit_res$contribution),decreasing=F),drop=FALSE]
select = which(rowSums(fit_res$contribution) > 0.01 * sum(fit_res$contribution))
plot_contribution(fit_res$contribution[select, ], cancer_signatures[,select], mode = "relative")

# Non-negative matrix factorization
mutation_matrix = mutation_matrix + 0.0001
estimate = nmf(mutation_matrix, rank=2:5, method="brunet", nrun=100, seed=123456)
plot(estimate)
nmf_res <- extract_signatures(mutation_matrix, rank = 2)



######## RANDOM TESTING UNDER HERE
library(MASS)
library(dndscv)

sample_mutations<-function(dbConnect, cohort)
{
  query = paste(
    "select sv.sampleId, chromosome as chr, position as pos, ref,alt from somaticVariant sv,sample s,patient p ",
    "where sv.sampleId = s.sampleId and s.patientId = p.id and p.primaryTumorLocation='",
    cohort,
    "' and length(alt) = length(ref) and filter = 'PASS' and gene <> '' limit 500000",
    sep = "")
  
  raw_data = dbGetQuery(dbConnect, query)
}

dbConnect = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
mutations<-sample_mutations(dbConnect,'Melanoma')
dbDisconnect(dbConnect)
dndscv(mutations)

paste((cohort),collapse=" ")
c(cohort)
