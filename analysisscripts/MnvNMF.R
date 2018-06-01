library(devtools) 
library(MutationalPatterns)
library(RMySQL)
library(data.table)
library(NMF)

# plotting
library(grid)
library(gridExtra)
library(ggplot2)

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

# INDELS DB count = 12M
# MNVs count = 1M


query_indels_mnvs<-function(dbConnect, sampleIdStr, typeStr) {
  query = paste(
    "SELECT SampleId, Type, Ref, Alt, Microhomology, RepeatCount, Clonality",
    "FROM somaticVariant",
    "WHERE filter = 'PASS' and type <> 'SNV'",
    sep = " ")
  
  if(sampleIdStr != "")
  {
    query = paste(query, " and sampleId in (", sampleIdStr, ")", sep='')
  }
  
  if(typeStr != "")
  {
    query = paste(query, " and type = '", typeStr, "'", sep='')
  }
  
  print(query)

  queryResults = dbGetQuery(dbConnect, query)
  # raw_types = raw_data$snv
  # standard_types = standard_mutation(raw_types)
  # raw_context = raw_data$context
  # context = standard_context(raw_types, standard_types, raw_context)
  # 
  # DT = data.table(
  #   sample = raw_data$sampleId,
  #   type = standard_types,
  #   context = context,
  #   ploidy = raw_data$ploidy,
  #   clonality = raw_data$clonality)
  
  return(queryResults) # was DT
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


# load("~/data/highestPurityCohortSummary.RData")
#View(highestPurityCohortSummary)

dbProd = dbConnect(MySQL(), user='hmf', password='HMFhmf@1', dbname='hmfpatients', groups = "RAnalysis")
# dbDisconnect(dbProd)

cancerTypes = highestPurityCohortSummary %>% group_by(cancerType) %>% count()
View(cancerTypes)

# download MNVs and INDELs 
mnvData = query_indels_mnvs(dbProd, "", "MNP")
View(mnvData)
nrow(mnvData)
write.csv(mnvData, "~/logs/r_output/mnv_prod_all.csv")

indelData = query_indels_mnvs(dbProd, "", "INDEL")
nrow(indelData)
write.csv(indelData, "~/logs/r_output/indel_prod_all.csv")

# limit to samples in the high-confidence set
mnvData = mnvData %>% filter(SampleId %in% highestPurityCohortSummary$sampleId)
nrow(mnvData)

indelData = indelData %>% filter(SampleId %in% highestPurityCohortSummary$sampleId)
nrow(indelData)


# MNV Buckets
View(head(mnvData,10))

mnvData$Length = stringi::stri_length(mnvData$Alt)
mnvData$Mutation = ifelse(mnvData$Length==2,paste(mnvData$Ref, mnvData$Alt, sep='>'), ifelse(mnvData$Length==3,"3Bases","4+Bases"))

View(mnvData %>% group_by(Length) %>% count())
View(mnvData %>% group_by(Mutation) %>% count())

mnvSampleCounts = mnvData %>% group_by(SampleId,Mutation) %>% count()
# mnvSampleCounts = mnvSampleCounts %>% filter(grepl("CPCT020101", SampleId)) # temp
View(mnvSampleCounts)
n_distinct(mnvSampleCounts$SampleId)

mnvMatrixData = mnvSampleCounts %>% spread(SampleId, n)
mnvMatrixData[is.na(mnvMatrixData)] = 0
View(mnvMatrixData[,1:10])

bucketNames = mnvMatrixData$Mutation
View(bucketNames)
mnvMatrixData = within(mnvMatrixData, rm(Mutation))

mnvNmfEstimate <- nmf(mnvMatrixData, rank=6:15, method="brunet", nrun=4, seed=123456, .opt='vp6')
plot(mnvNmfEstimate)
save(mnvNmfEstimate, file="~/logs/r_output/mnvNmfEstimate.RData")


# generate the actual NMF results
mnvSigCount = 11
mnvNmfResult <- nmf(mnvMatrixData, rank=mnvSigCount, method="brunet", nrun=5, seed=123456, .opt='vp6')
save(mnvNmfResult, file="~/logs/r_output/mnvNmfResult_sig11.RData")




# INDEL Buckets






