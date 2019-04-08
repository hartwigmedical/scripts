library(devtools) #; install_github("im3sanger/dndscv")
library(purple);
library(data.table)
library(dplyr)
library(tidyr)
library(stringi)

library(MutationalPatterns)
library(RMySQL)
library(data.table)
library("NMF")
library("pracma")

# plotting
library(grid)
library(gridExtra)
library(ggplot2)
library(cowplot)

getCOSMICSignatures <- function()
{
  sp_url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
  cancer_signatures = read.table(sp_url, sep = "\t", header = T)
  # reorder (to make the order of the trinucleotide changes the same)
  cancer_signatures = cancer_signatures[order(cancer_signatures[, 1]),]
  # only signatures in matrix
  cancer_signatures = as.matrix(cancer_signatures[, 4:33])

  return (cancer_signatures)
}

getCOSMICSigBuckets <- function()
{
  sp_url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
  cancer_signatures = read.table(sp_url, sep = "\t", header = T)
  # reorder (to make the order of the trinucleotide changes the same)
  cancer_signatures = cancer_signatures[order(cancer_signatures[, 1]),]
  # only signatures in matrix
  sig_buckets = as.matrix(cancer_signatures[, 1])

  return (sig_buckets)
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

query_variants <- function(dbConnect, sampleIdStr) {
    query = paste(
        "SELECT sampleId, trinucleotideContext as context, concat(ref,'>', alt) as snv, adjustedVaf * adjustedCopyNumber as ploidy, clonality",
        "FROM somaticVariant",
        "WHERE filter = 'PASS' and length(alt) = length(ref) and length(alt) = 1 and trinucleotideContext not like '%N%'",
        sep = " ")

    if(sampleIdStr != "")
    {
      query = paste(query, " and sampleId in (", sampleIdStr, ")", sep='')
    }

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

getSampleIdsStr<-function(samples)
{
  sampleIdsStr = ""

  for(i in 1:nrow(samples))
  {
    sample <- samples[i,]

    if(i > 1)
      sampleIdStr = paste(",'", sample$sampleId, "'", sep="")
    else
      sampleIdStr = paste("'", sample$sampleId, "'", sep="")

    # sampleIdStr = stringi::stri_replace_all_fixed(sampleIdStr, " ", "")
    sampleIdsStr = paste(sampleIdsStr, sampleIdStr)
  }

  return (sampleIdsStr)
}


# loading frozen data files and prep for sample counts

# packages required
#library(devtools) #; install_github("im3sanger/dndscv")
#library(purple);
library(data.table)
library(dplyr)
#library(tidyr)
#library(stringi)

library(MutationalPatterns)
# library(RMySQL)
# library(data.table)



load('/data/data_archive/180921_paper_db_files/Reference/allSomatics_p1.RData')
load('/data/data_archive/180921_paper_db_files/Reference/allSomatics_p2.RData')

sampleSNVs = read.csv('~/logs/CPCT02020323T_SNVs.csv')
nrow(sampleSNVs)
View(sampleSNVs)





# filter for SNV types

# "SELECT sampleId, trinucleotideContext as context, concat(ref,'>', alt) as snv, adjustedVaf * adjustedCopyNumber as ploidy, clonality",
# "FROM somaticVariant",
# "WHERE filter = 'PASS' and length(alt) = length(ref) and length(alt) = 1 and trinucleotideContext not like '%N%'",

# View(sampleSNVs %>% filter(type=='SNP'&!grepl('N',trinucleotideContext)&filter=='PASS'))

allSomatics = (allSomatics_p2 %>% filter(type=='SNP'&!grepl('N',trinucleotideContext)&filter=='PASS') 
                  %>% mutate(context=trinucleotideContext,
                             ploidy=adjustedVaf*adjustedCopyNumber,
                             snv=paste(ref,alt,sep='>'))
                  %>% select(sampleId,context,snv,ploidy,clonality))

View(allSomatics)

View(allSomatics %>% group_by(context) %>% count())

raw_types = allSomatics$snv
View(head(raw_types,100))

# load functions: standard_mutation and standard_context

# prepare bucket data
standard_types = standard_mutation(raw_types)
View(head(standard_types,100))
raw_context = allSomatics$context
View(head(raw_context,100))


context = standard_context(raw_types, standard_types, raw_context)
View(context)

allSomaticsDT = data.table(
  sample = allSomatics$sampleId,
  type = standard_types,
  context = context,
  ploidy = allSomatics$ploidy,
  clonality = allSomatics$clonality)

View(allSomaticsDT)

# filter for high-purity samples (2405 samples)
load('~/data/r_data/highestPurityCohortSummary.RData')
nrow(highestPurityCohortSummary)
highestPurityCohort = highestPurityCohortSummary 
allSomaticsDT = allSomaticsDT %>% filter(sample %in% highestPurityCohort$sampleId)
nrow(allSomaticsDT)

# prepare sample count bucket data

View(allSomaticsDT %>% group_by(sample,context,type) %>% summarise(Count=n()))

sampleCounts = allSomaticsDT %>% group_by(sample,context,type) %>% summarise(Count=n())

emptySigs = create_empty_signature()
View(emptySigs)

sampleResult = merge(emptySigs, sampleCounts, all=T)
sampleResult[is.na(sampleResult)] <- 0
View(sampleResult)

sampleResult = sampleResult %>% mutate(Bucket=paste(type,context,sep='_'))


colnames(allSomaticsDT) = c("SampleId", "Context", "SNV", "Ploidy", "Clonality")

snvSampleNames = unique(allSomaticsDT$sample)

View(snvSampleNames)

allSampleCounts = data.frame()

index = 1
for(s in snvSampleNames)
{
  index = index + 1
  
  emptySigs = create_empty_signature()
  
  sampleCounts = allSomaticsDT %>% filter(sample==s) %>% group_by(sample,context,type) %>% summarise(Count=n())

  sampleResult = merge(emptySigs, sampleCounts, all=T)
  sampleResult[is.na(sampleResult)] <- 0
  
  allSampleCounts = rbind(allSampleCounts,sampleResult)
}

View(allSampleCounts)
allSampleCounts = allSampleCounts %>% mutate(Bucket=paste(type,context,sep='_')) %>% select(sample,Bucket,Count)

allSampleCounts = allSampleCounts %>% select(sample,Bucket,Count)
colnames(allSampleCounts) = c('SampleId', "Bucket",'Count')








## NMF functions

# restrict to samples in high-purity and DR022 set
load("~/data/r_data/highestPurityCohortSummary.RData")
dr022Samples  = read.csv("~/data/DR-022_metadata.tsv", sep='\t')
View(dr022Samples)
hpcSamples = highestPurityCohortSummary %>% filter(sampleId %in% dr022Samples$sampleId)
# hpcSamples = highestPurityCohortSummary
View(hpcSamples)
nrow(hpcSamples)

write.csv(highestPurityCohortSummary %>% select(sampleId,patientId), "~/dev/bachelor/dr22_samples.csv", quote=F,row.names = F)
nonDr22Samples = highestPurityCohortSummary %>% filter(!(sampleId %in% dr022Samples$sampleId))
nrow(nonDr22Samples)
write.csv(nonDr22Samples %>% select(sampleId,patientId), "~/dev/bachelor/non_dr22_samples.csv", quote=F,row.names = F)


load("~/data/r_data/highestPurityCohortSummary.RData")

write.csv(highestPurityCohortSummary %>% select(sampleId,patientId), "~/dev/bachelor/driver_paper_samples.csv", row.names = F, quote=F)


View(highestPurityCohortSummary)

# set MSI status
baSampleData = read.csv("~/dev/nmf/sample_ext_data3.csv")
View(baSampleData)

baSampleData = merge(baSampleData, hpcSamples %>% select (sampleId, msiStatus), by.x='SampleId', by.y='sampleId', all.x=T)
baSampleData$msiStatus = ifelse(baSampleData$msiStatus=='MSI','MSI','')
baSampleData$MSI = (baSampleData$msiStatus=='MSI')
baSampleData = within(baSampleData, rm(msiStatus))
write.csv(baSampleData, "~/dev/nmf/sample_ext_data4.csv", quote=F, row.names=F)

highestPurityCohortSummary$IndelByMsiScore = highestPurityCohortSummary$TOTAL_INDEL/highestPurityCohortSummary$msiScore

cosmicSignatures = getCOSMICSignatures()
View(cosmicSignatures)
write.csv(cosmicSignatures, "~/data/r_data/snv_cosmic_sigs.csv", quote=F, row.names=F)
cosmicSignatures = as.matrix(read.csv("~/data/r_data/snv_cosmic_sigs.csv"))

# dbDisconnect(dbProd)

cancerTypes = hpcSamples %>% group_by(cancerType) %>% count()
View(cancerTypes)

# download SNVs by cancer type (to keep the queries smaller)
for(i in 1:nrow(cancerTypes))
{
  cancerTypeRow = cancerTypes[i,]
  cancerTypeStr = cancerTypeRow$cancerType

  print(paste(i, ": retrieving SNVs for cancer=", cancerTypeStr, sep=''))

  sampleIdsStr = getSampleIdsStr(hpcSamples %>% filter(cancerType==cancerTypeStr))
  # print(sampleIdsStr)
  snvVariants = query_variants(dbProd, sampleIdsStr)
  print(paste(i, ":cancer=", cancerTypeStr, ", SNV count=", nrow(snvVariants),  sep=''))

  ctStr = stringi::stri_replace_all_fixed(cancerTypeStr, '/', '')
  snvFilename = paste("~/logs/r_output/snv_", ctStr, ".csv", sep='')

  write.csv(snvVariants, snvFilename, row.names=F, quote=F)
}

rm(snvVariants)

# now reload these same variants from file into a single massive set
allSnvVariants = data.frame(matrix(ncol = 5, nrow = 0))
allSnvVariants = setNames(allSnvVariants, c("SampleId", "Context", "SNV", "Ploidy", "Clonality"))

# View(cancerTypes)
for(i in 1:nrow(cancerTypes))
{
  cancerTypeRow = cancerTypes[i,]
  cancerTypeStr = cancerTypeRow$CancerType
  print(paste(i, ": loading SNVs for cancer=", cancerTypeStr, sep=''))

  ctStr = stringi::stri_replace_all_fixed(cancerTypeStr, '/', '')
  snvFilename = paste("~/logs/r_output/snv_", ctStr, ".csv", sep='')

  allSnvVariants = rbind(allSnvVariants, read.csv(snvFilename))
}

nrow(allSnvVariants)
# rm(allSnvVariants)

svnSampleStats = allSnvVariants %>% group_by(sample) %>% summarise(Count=n())
View(svnSampleStats)
# load("~/logs/r_output/snvSampleCounts.RData")
# View(result)

# sampleCountResults = result
index = 1
for(s in snvSampleNames)
{
  print(paste(index, ": collating data for ", s, sep=''))
  index = index + 1

  emptySigs = create_empty_signature()

  # allSnvVariants
  sampleCounts = (allSnvVariants %>% filter(sample==s) %>% group_by(type,context)
                  %>% summarise(total=n(),
                                subclonal=sum(clonality=='SUBCLONAL'),
                                clonal=sum(clonality=='CLONAL')))

  sampleResult = merge(emptySigs, sampleCounts, all=T)
  sampleResult[is.na(sampleResult)] <- 0

  sampleCountResults[[s]] = sampleResult
}

save(sampleCountResults, file="~/logs/r_output/snvSampleCounts2.RData")
load("~/logs/r_output/snvSampleCounts2.RData")

View(sampleCountResults)

# release the mammoth
rm(allSnvVariants)

# fit the samples counts to the cosmic signatures
snvSampleSigs = list()
for(s in snvSampleNames)
# for(s in sampleCountResults)
  {
  print(paste("fitting signatures for ", s, sep=''))

  if (!is.null(sampleCountResults[[s]]))
  {
    # we need to slice out only the mutation count columns (delete col 1 and 2)
    res = fit_to_signatures(sampleCountResults[[s]][, -c(1, 2)], cosmicSignatures)
    snvSampleSigs[[s]] <- res$contribution
  }
}

View(snvSampleSigs)
View(cosmicSignatures)

snvContribsTotal = data.frame(matrix(ncol = 0, nrow = snvSigCount))
snvContribsSubclonal = data.frame(matrix(ncol = 0, nrow = snvSigCount))

for(s in snvSampleNames)
{
  # print(paste("extracting sig data for ", s, sep=''))

  if (!is.null(snvSampleSigs[[s]]))
  {
    sampleSigData = snvSampleSigs[[s]]
    totalValues = sampleSigData[,1]
    snvContribsTotal = cbind(snvContribsTotal, totalValues)

    subclonalValues = sampleSigData[,2]
    snvContribsSubclonal = cbind(snvContribsSubclonal, subclonalValues)

    # clonalValues = sampleSigData[,3]
    # contributionClonal = cbind(contributionClonal, clonalValues)
  }
}

snvContribsTotal = setNames(snvContribsTotal, snvSampleNames)
snvContribsSubclonal = setNames(snvContribsSubclonal, snvSampleNames)
# contributionClonal = setNames(contributionClonal, snvSampleNames)
View(snvContribsTotal)
View(snvContribsSubclonal)


# SNV NMF Standard Routine

# Signature Creation (ie independently from COSMIC)

# create bucket names
snvBucketNames = create_empty_signature()
snvBucketNames = unite(snvBucketNames, "Bucket", type, context, sep='_')
View(snvBucketNames)

# create sample counts per bucket
snvBucketCount = nrow(snvBucketNames)
print(snvBucketCount)

snvMatrixData = data.frame(matrix(ncol = 0, nrow = snvBucketCount))
snvSampleTotals = data.frame(matrix(ncol = 2, nrow = 0))
snvSampleTotals = setNames(snvSampleTotals, c("SampleId", "SampleCount"))
snvSampleCounts = data.frame(matrix(ncol = 3, nrow = 0))

totalColIndex = 3 # (could check on the fly)
subclonalColIndex = 4

snvSampleNames = hpcSamples$sampleId
View(snvSampleNames)

colIndex = totalColIndex
# colIndex = subclonalColIndex
for(s in snvSampleNames)
{
  # print(paste("extracting sig data for ", s, sep=''))
  sampleCounts = sampleCountResults[[s]]

  snvMatrixData = cbind(snvMatrixData, sampleCounts[,colIndex])

  sampleCountData = data.frame(matrix(ncol = 0, nrow=snvBucketCount))
  sampleCountData$SampleId = s
  sampleCountData = cbind(sampleCountData, snvBucketNames)
  sampleCountData = cbind(sampleCountData, sampleCounts[,colIndex])
  snvSampleCounts = rbind(snvSampleCounts, sampleCountData)

  rowIndex = nrow(snvSampleTotals)+1
  snvSampleTotals[rowIndex,1] = s
  snvSampleTotals[rowIndex,2] = sum(sampleCounts[,colIndex])
  # index = index + 1
  # if(index > 10)
  #   break
}

snvMatrixData = setNames(snvMatrixData, snvSampleNames)
snvSampleCounts = setNames(snvSampleCounts, c("SampleId", "Bucket", "Count"))
snvSampleCounts = snvSampleCounts %>% filter(Count>0) # since this DF does not need to contain a row per bucket if zero
snvSampleTotals = snvSampleTotals %>% arrange(-SampleCount) # for use with high and low mutational load analysis
View(snvSampleTotals)
View(snvSampleCounts)

# View(sampleBucketCounts)
View(snvMatrixData[,1:10])
nrow(snvMatrixData)
ncol(snvMatrixData)

# repeat high and low mutational loads
highMutLoadSamples = head(snvSampleTotals,nrow(snvSampleTotals)*0.1)
View(highMutLoadSamples)

snvHighMatrixData = data.frame(matrix(ncol = 0, nrow = 96))

for(s in highMutLoadSamples$SampleId)
{
  sampleCounts = sampleCountResults[[s]]
  snvHighMatrixData = cbind(snvHighMatrixData, sampleCounts[,totalColIndex])
}

colnames(snvHighMatrixData) <- highMutLoadSamples$SampleId
View(snvHighMatrixData[,1:20])

lowMutLoadSamples = tail(snvSampleTotals,nrow(snvSampleTotals)*0.5)
View(lowMutLoadSamples)

snvLowMatrixData = data.frame(matrix(ncol = 0, nrow = 96))

for(s in lowMutLoadSamples$SampleId)
{
  sampleCounts = sampleCountResults[[s]]
  snvLowMatrixData = cbind(snvLowMatrixData, sampleCounts[,totalColIndex])
}

snvLowMatrixData = setNames(snvLowMatrixData, highMutLoadSamples$SampleId)

View(snvLowMatrixData[,1:20])

write.csv(snvHighMatrixData, file="~/logs/r_output/snv_nmf_high_counts.csv", row.names=F, quote=F)
write.csv(snvLowMatrixData, file="~/logs/r_output/snv_nmf_low_counts.csv", row.names=F, quote=F)

load("~/data/snvNmfResult_sig30.RData")
View(snvNmfResult)

snvSigCount = 30
snvSigNamesNum = get_signame_list(snvSigCount, F)
snvSigNamesStr = get_signame_list(snvSigCount, T)

snvSignatures = cosmicSignatures
colnames(snvSignatures) <- NULL
rownames(snvContribsTotal) <- NULL
rownames(snvContribsSubclonal) <- NULL
View(snvSignatures)

View(snvSampleCounts)

write.csv(snvMatrixData, file="~/data/r_data/snv_nmf_matrix_data.csv", row.names=F, quote=F)
write.csv(snvSampleCounts, file="~/data/r_data/snv_nmf_sample_counts.csv", row.names=F, quote=F)

snvMatrixData = as.matrix(read.csv("~/data/r_data/snv_nmf_matrix_data.csv", stringsAsFactors=F))
snvSampleCounts = as.data.frame(read.csv("~/data/r_data/snv_nmf_sample_counts.csv", stringsAsFactors=F))
sum(snvMatrixData) # check sum..
sum(snvSampleCounts$Count)
ncol(snvMatrixData)


# write.csv(snvMatrixData, file="~/data/r_data/snv_nmf_subc_matrix_data.csv", row.names=F, quote=F)
# write.csv(snvContribsTotal, file="~/data/r_data/snv_nmf_contribs.csv", row.names=F, quote=F)
# write.csv(snvSampleCounts, file="~/data/r_data/snv_nmf_subc_sample_counts.csv", row.names=F, quote=F)
# write.csv(snvContribsSubclonal, file="~/data/r_data/snv_nmf_subc_contribs.csv", row.names=F, quote=F)

snvContribsTotal = as.matrix(read.csv(file="~/data/r_data/snv_nmf_contribs.csv", stringsAsFactors=F))
is.matrix(snvContribsTotal)
snvSignatures = as.matrix(snvSignatures)
is.matrix(snvSignatures)

# evaluation using cosmic signatures
evaluate_nmf_data("SNV", "sig30_cosmic_2", snvSigCount, snvSignatures, snvContribsTotal, snvMatrixData, snvSampleCounts,
                 sampleCancerTypes, snvBucketNames, snvSigNamesNum, snvSigNamesStr, FALSE, FALSE, TRUE)



# subclonal evaluation
evaluate_nmf_data("SNV_Subclonal", "sig30_cosmic", snvSigCount, snvSignatures, snvContribsSubclonal, snvMatrixData, snvSampleCounts,
                  sampleCancerTypes, snvBucketNames, snvSigNamesNum, snvSigNamesStr, TRUE, FALSE, TRUE)


snvSampleSigData = get_sig_data(snvSignatures, snvContribsTotal, snvSigNamesStr, colnames(snvContribsTotal))
snvSampleSigData$Count = round(snvSampleSigData$Count,0)
View(snvSampleSigData)
write.csv(snvSampleSigData, "~/logs/r_output/snvSampleSigData.csv", row.names=F, quote=F)
snvSampleSigData = read.csv("~/logs/r_output/snvSampleSigData.csv")
View(snvSampleSigData)

cancerSampleCounts = highestPurityCohortSummary %>% group_by(cancerType) %>% summarise(CancerSampleCount=n())
View(cancerSampleCounts)

snvSampleSigData = merge(snvSampleSigData, highestPurityCohortSummary %>% select(sampleId,cancerType), by.x="SampleId", by.y="sampleId", all.x=T)
View(snvSampleSigData)

snvSigCancerTotals = (snvSampleSigData %>% filter(Count>0) %>% group_by(SigName,cancerType)
                      %>% summarise(SampleCount=n(),
                                    Min=min(Count),
                                    Max=max(Count),
                                    TotalCount=sum(Count),
                                    Median=round(median(Count),0),
                                    Mean=round(mean(Count),0)))

snvSigCancerTotals = merge(snvSigCancerTotals, cancerSampleCounts, all.x=T)
View(snvSigCancerTotals)
write.csv(snvSigCancerTotals, "~/logs/r_output/snvSigCancerTotals.csv", quote=F, row.names=F)


# evaluation simulated counts and cosmic signatures
snvSigSimCount = 13
snvSimSigNamesNum = get_signame_list(snvSigSimCount, F)
snvSimSigNamesStr = get_signame_list(snvSigSimCount, T)



# all SNVs actuals fit with Cosmic sigs
sum(snvSampleCounts$Count)
n_distinct(snvSampleCounts$SampleId)

View(snvSampleCounts)
snvMatrixData = snvSampleCounts %>% select(SampleId,Bucket,Count) %>% spread(SampleId,Count)
snvMatrixData[is.na(snvMatrixData)] = 0
snvMatrixData = within(snvMatrixData, rm(Bucket))

snvCosmicContribs = apply_signatures(snvMatrixData, cosmicSignatures)
snvSig30Count = 30
snvSig30NamesNum = get_signame_list(snvSig30Count, F)
snvSig30NamesStr = get_signame_list(snvSig30Count, T)
colnames(cosmicSignatures) = snvSig30NamesNum

evaluate_nmf_data("SNV", "cosmic30_fit_r_all", cosmicSignatures, snvCosmicContribs, snvMatrixData, snvSampleCounts,
                  sampleCancerTypes, snvBucketNames, snvSig30NamesNum, snvSig30NamesStr, T, F)

# all SNVs actuals fit with PCAWG sigs
pcawgSigs = read.csv("~/dev/nmf/snv_pcawg_sigs.csv")
# pcawgSigs = pcawgSigs %>% select(-Type,-SubType)
pcawgSigs = as.matrix(pcawgSigs)
write.csv(pcawgSigs, "~/dev/nmf/snv_pcawg_sigs.csv", quote=F, row.names=F)
View(pcawgSigs)
pcawgSigCount = ncol(pcawgSigs)

snvPcawgContribs = apply_signatures(snvMatrixData, as.matrix(pcawgSigs))
# snvPcawgSigsNamesStr = get_signame_list(pcawgSigCount, T)

evaluate_nmf_data("SNV", "pcawg_fit_r_all", pcawgSigs, snvPcawgContribs, snvMatrixData, snvSampleCounts,
                  sampleCancerTypes, snvBucketNames, snvPcawgSigsNamesStr, T, F, T)

snvPcawgReducedContribs = as.matrix(read.csv("~/dev/nmf/logs/snv_pcawg_fit_nmf_contribs.csv"))

tmpT = t(snvPcawgReducedContribs)
View(tmpT)

snvPcawgSigsNamesStr = c("01", "02", "03", "04", "05", "06", "07a", "07b", "07c", "07d", "08", "09", "10a", "10b",
                         "11", "12", "13", "14", "15", "16", "17a", "17b", "18", "19", "20",
                         "21", "22", "23", "24", "25", "26", "27", "28", "29", "30",
                         "31", "32", "33", "34", "35", "36", "37", "38", "39", "40",
                         "41", "42", "43", "44", "45", "46", "47", "48", "49", "50",
                         "51", "52", "53", "54", "55", "56", "57", "58", "59", "60")

length(snvPcawgSigsNamesStr)

evaluate_nmf_data("SNV", "pcawg_reduced_fit_v4", pcawgSigs, snvPcawgReducedContribs, snvMatrixData, snvSampleCounts,
                  sampleCancerTypes, snvBucketNames, snvPcawgSigsNamesStr, T, F)


sampleSigData = get_sig_data(pcawgSigs, snvPcawgReducedContribs, snvPcawgSigsNamesStr, colnames(snvPcawgReducedContribs))
sigStats = get_sig_stats(sampleSigData)
View(sigStats)

# just on Kidney
snvPcawgReducedContribs = as.matrix(read.csv("~/dev/nmf/logs/snv_kidney_nmf_contribs.csv"))
ncol(snvPcawgReducedContribs)
View(snvPcawgReducedContribs)

evaluate_nmf_data("SNV", "pcawg_reduced_fit_kidney", pcawgSigs, snvPcawgReducedContribs, snvSimMatrixData, snvSimSampleCounts,
                  sampleCancerTypes, snvBucketNames, snvPcawgSigsNamesStr, F, F)



# all SNVs with de-novo fit (30 sigs)
snv30Signatures = as.matrix(read.csv("~/dev/nmf/logs/snv_all_denovo_nmf_sigs.csv"))
colnames(snv30Signatures) = snvSig30NamesNum
snv30Contribs = as.matrix(read.csv("~/dev/nmf/logs/snv_all_denovo_nmf_contribs.csv"))

write.csv(snv30Contribs, "~/dev/nmf/logs/snv_all_denovo_nmf_contribs_4proposed.csv", row.names=F, quote=F)
write.csv(snv30Signatures, "~/dev/nmf/logs/snv_all_denovo_nmf_sigs_4proposed.csv", row.names=F, quote=F)


evaluate_nmf_data("SNV", "denovo_sig30_all_30prop", snvSig30Count, snv30Signatures, snv30Contribs, snvMatrixData, snvSampleCounts,
                  sampleCancerTypes, snvBucketNames, snvSig30NamesNum, snvSig30NamesStr, T, F, T)


# Bucket-Analyser sigs and contributions
snvBaContribs = as.matrix(read.csv(file="~/dev/nmf/logs/snv_rr_ba_contribs.csv", stringsAsFactors=F))
View(snvBaContribs[,50:70])
nrow(snvBaContribs)
snvBaSigs = as.matrix(read.csv(file="~/dev/nmf/logs/snv_rr_ba_sigs.csv", stringsAsFactors=F))
View(snvBaSigs)
snvBaSigNames = colnames(snvBaSigs)
print(svnBaSigNames)
ncol(snvBaContribs)
snvBaSigCount = ncol(snvBaSigs)
print(snvBaSigCount)

snvBaSigNames = colnames(snvBaSigs)
snvBaSigNames = trim_ba_sig_names(snvBaSigNames)
print(snvBaSigNames)


evaluate_nmf_data("SNV", "ba_denovo_all_RR", snvBaSigs, snvBaContribs, snvMatrixData, snvSampleCounts,
                  sampleCancerTypes2, snvBucketNames, snvBaSigNames, T, F, bgSigCount, F)

produce_signature_report("SNV", "ba_denovo_all_RR_20", "~/dev/nmf/logs/snv_rr_ba_sigs.csv", "~/dev/nmf/logs/snv_rr_ba_contribs.csv",
                         "~/dev/nmf/logs/snv_rr_ba_group_data.csv", "~/dev/nmf/logs/snv_rr_ba_sample_sig_allocs.csv",
                         snvMatrixData, snvSampleCounts, sampleCancerTypes2, snvBucketNames, T, F, F)

produce_signature_report("SNV", "ba_denovo_no_MSI", "~/dev/nmf/logs/snv_no_msi_ba_sigs.csv", "~/dev/nmf/logs/snv_no_msi_ba_contribs.csv",
                         "~/dev/nmf/logs/snv_no_msi_ba_group_data.csv", snvMatrixData, snvSampleCounts,
                         sampleCancerTypes2, snvBucketNames, T, F, F)

produce_signature_report("SNV", "ba_PCAWG_fit", "~/dev/nmf/logs/snv_pcwg_ba_sigs.csv", "~/dev/nmf/logs/snv_pcwg_ba_contribs.csv",
                         "~/dev/nmf/logs/snv_pcwg_ba_group_data.csv", snvMatrixData, snvSampleCounts,
                         sampleCancerTypes2, snvBucketNames, T, F, F)

# baSigInfo = as.data.frame(read.csv(file="~/dev/nmf/logs/snv_ba_group_data_test.csv", stringsAsFactors=F))
baSigInfo = as.data.frame(read.csv(file="~/dev/nmf/logs/snv_rr_ba_group_data.csv", stringsAsFactors=F))
View(baSigInfo)
sigCount = nrow(baSigInfo)




# DEBUG ONLY
contribution = snvBaContribs
signatures = snvBaSigs
summaryCounts = snvSampleCounts
bucketNames = snvBucketNames
sigNamesNamed = snvBaSigNames
matrixData = snvMatrixData

sigCount = nrow(contribution)

sampleNames = colnames(contribution)

origSampleCounts = summaryCounts %>% group_by(SampleId) %>% summarise(OrigSampleCount=sum(Count))

sigNamesUnamed = get_signame_list(sigCount, F)
colnames(signatures) = sigNamesUnamed

sigNamesCombined = cbind(sigNamesUnamed, sigNamesNamed)
colnames(sigNamesCombined) <- c("Signature", "SigName")
print(sigNamesCombined)

bucketIndex = data.frame(as.numeric(as.character(rownames(bucketNames))))
colnames(bucketIndex) <- c("BucketIndex")
bucketNamesIndexed = cbind(bucketNames, bucketIndex)
bucketNamesIndexed$BucketIndex = bucketNamesIndexed$BucketIndex-1

sigBucketData = get_bucket_data(signatures, contribution, bucketNames)
sigBucketData = merge(sigBucketData,bucketNamesIndexed,by="Bucket",all.x=T)
sigBucketData = merge(sigBucketData,sigNamesCombined,by="Signature",all.x=T)
sigBucketStats = get_sig_bucket_stats(sigBucketData)
sigBucketTopN = get_top_buckets(sigBucketData)
bucketSummaryData = get_bucket_stats(sigBucketData)
sampleBucketData = get_sample_bucket_data(matrixData, origSampleCounts, bucketNames)
sampleSigData = get_sig_data(signatures, contribution, sigNamesNamed, sampleNames)
View(sampleSigData)



sigStats = get_sig_stats(sampleSigData)

# calculate and factor in residuals
sampleSigData = append_residuals(contribution, signatures, matrixData, bucketNames, sampleSigData)

sampleSigData = merge(sampleSigData, sampleCancerTypes,by.x="SampleId",by.y="SampleId",all.x=TRUE)
sampleSigData$CancerType = ifelse(is.na(sampleSigData$CancerType), 'N/A', paste(sampleSigData$CancerType, sep=""))

View(sampleSigData)

nrow(sampleSigData %>% filter(Count > 0))
View(sampleSigData %>% filter(Count > 0))


sampleSigData = merge(sampleSigData, origSampleCounts, by.x="SampleId",by.y="SampleId",all.x=TRUE)
View(sampleSigData)
allSamplesBySampleCount = sampleSigData %>% group_by(SampleId) %>% summarise(SampleCount=first(OrigSampleCount)) %>% arrange(-SampleCount)
View(allSamplesBySampleCount)

sampleSigData = merge(sampleSigData, sampleCancerTypes,by.x="SampleId",by.y="SampleId",all.x=TRUE)

names(sampleSigData)[names(sampleSigData) == 'OrigSampleCount.x'] <- 'SampleCount'

sigColours = get_sig_colours(40)
plot_sig_samples(sampleSigData, "Skin", sigColours, "SV", 0, T)

View(sampleSigData %>% group_by(SigName) %>% count())

topNSamplesBySig = head(sampleSigData %>% filter(SigName=="01_BG_Stoma"&Count>0) %>% arrange(-Count),50)
View(topNSamplesBySig)
topNSampleSigData = sampleSigData %>% filter(Count>0&SampleId %in% topNSamplesBySig$SampleId)
View(topNSampleSigData)

bgSigList = topNSampleSigData %>% filter(grepl("BG_",SigName)) %>% group_by(SigName) %>% count()
View(bgSigList)
bgSigCount = nrow(bgSigList)

print(sigColours)
sampleSigDataNoResiduals = sampleSigData %>% filter(SigName!="Unalloc"&SigName!="Excess")
plot_top_n_samples_by_sig(sampleSigDataNoResiduals, sigNamesNamed, 50, sigColours, "SV", T)

# try to match sampleSigData using external sig allocs data
sigAllocs = as.data.frame(read.csv(file="~/dev/nmf/logs/snv_rr_ba_sample_sig_allocs.csv", stringsAsFactors=F))
colnames(sigAllocs) = c("Signature", "BgId", "SampleId", "Count", "SigPercent")
sigAllocs$Signature = ifelse(sigAllocs$BgId=='Excess',42,sigAllocs$Signature)
View(sigAllocs)
nrow(sigAllocs %>% filter(Count > 0))
sampleSigData2 = sigAllocs
sigCount = length(sigNamesCombined)
sigNamesCombined = rbind(sigNamesCombined, c(sigCount+1,"Unalloc"))
sigNamesCombined = rbind(sigNamesCombined, c(sigCount+2,"Excess"))
View(sigNamesCombined)
sampleSigData2$PercBucket = round(sampleSigData2$SigPercent/0.1)*0.1
sampleSigData2 = merge(sampleSigData2, sigNamesCombined, by.x="Signature", by.y="Signature", all=T)
sampleSigData2 = merge(sampleSigData2, sampleSigFullSet, by=c("SampleId", "Signature"), all=T)
sampleSigData2[is.na(sampleSigData2)] = 0
View(sampleSigData2)

sampleSigData2 = merge(sigAllocs, sampleSigFullSet, by=c("SampleId", "Signature"), all=T)
sampleSigData2[is.na(sampleSigData2)] = 0


sampleSigFullSet = merge(sampleNames, sigNamesCombined)
colnames(sampleSigFullSet) = c("SampleId", "Signature", "SigName")
View(sampleSigFullSet)


nrow(sampleSigData2 %>% filter(Count > 0))

View(origSampleCounts)





# ratio range data
ratioRangeData = as.data.frame(read.csv(file="~/dev/nmf/logs/snv_rr_ba_grp_ratio_ranges.csv", stringsAsFactors=F))
View(ratioRangeData)


print(plot_ratio_range_distribution(ratioRangeData, 443, T, 0.1, 0.02))

tmpRR = ratioRangeData %>% filter(BgId == -1)
bucketList = tmpRR %>% group_by(Bucket) %>% count()
View(bucketList)
bucketCount = nrow(bucketList)
print(bucketCount)

bucketPlotList = bucketList[1:20,]
View(bucketPlotList)

print(paste("pageCount=", pageCount, ", bucketCount=", bucketCount, ", start=", startBucketIndex, ", end=", endBucketIndex, sep=''))

plotData = tmpRR %>% filter(Bucket %in% bucketPlotList$Bucket)
View(plotData)

# create PDF of ratio range distributions
runType = "SNV"
runId = "ratio_ranges"
outputFile = paste("~/logs/r_output/pdfs/", runType, "_", runId, ".pdf", sep = "")
pdf(file=outputFile, height = 14, width = 20)
par(mar=c(1,1,1,1))

#ratioRangeData = as.data.frame(read.csv(file="~/dev/nmf/logs/snv_rr_ba_grp_ratio_ranges.csv", stringsAsFactors=F))
#ratioRangeData = within(ratioRangeData, rm(Type.1))

bgRRSummary = ratioRangeData %>% group_by(BgId) %>% summarise(SampleCount=first(SampleCount))
View(bgRRSummary)
bgList = unique(ratioRangeData$BgId)

for(bgId in bgList)
{
  minRatioSeg = 0.005
  if(bgId == -1)
    minRatioSeg = 0.02

  print(plot_ratio_range_distribution(ratioRangeData, bgId, T, 0.05, minRatioSeg))
  print(plot_ratio_range_distribution(ratioRangeData, bgId, F, 0.05, minRatioSeg))
}

dev.off()



print(plot_ratio_range_distribution(ratioRangeData, 443, F, 0.05))
print(plot_ratio_range_distribution(ratioRangeData, 443, T, 0.05))

print(plot_ratio_range_distribution(ratioRangeData, -1, T, 100, 0.1, 0.02))

colCount = ncol(ratioRangeData)
print(colCount)
freqColStart = colCount - 100 + 1
weightColEnd = freqColStart -1
weightColStart = weightColEnd - 100 + 1

nrow(ratioRangeData)


tmpRRData = ratioRangeData %>% filter(BgId==443)
tmpRRData = ratioRangeData %>% filter(BgId==-1)
maxRatio = max(tmpRRData$SigRatio)
tmpRRData = tmpRRData %>% filter(SigRatio >= maxRatio * 0.1)
sampleFreqData = tmpRRData[,freqColStart:colCount]
View(sampleFreqData)


sampleWeightData = tmpRRData[,weightColStart:weightColEnd]
xAxisValues = colnames(sampleWeightData)
print(xAxisValues)
ratioSegmentNames = stri_replace_all_fixed(xAxisValues,'X','')
print(ratioSegmentNames)

rsTmp = as.numeric(ratioSegmentNames)
View(rsTmp)
colnames(sampleWeightData) = ratioSegmentNames
View(sampleWeightData)
# sampleWeightData = cbind(tmpRRData$Bucket, sampleWeightData)

swTrans = as.data.frame(t(sampleWeightData))
colnames(swTrans) = tmpRRData$Bucket
swTrans = cbind(ratioSegmentNames, swTrans)
rownames(swTrans) = NULL
View(swTrans)

# sampleSigPercData = sigPercentsBySample %>% gather(SigName, SigPercent, -SampleId) %>% arrange(SampleId, SigName)
# sampleBucketCounts = cbind(bucketNames,matrixData)
# gatherIndex = ncol(sampleBucketCounts)
# sampleBucketCounts2 = gather(sampleBucketCounts, "SampleId", "ActualCount", 2:gatherIndex)

gatherIndex = ncol(swTrans)
swTrans2 = gather(swTrans, "Bucket", "Weight", 2:gatherIndex)
colnames(swTrans2) = c("RatioSegment", "Bucket", "Weight")
View(swTrans2)

weightPlot = (ggplot(data = swTrans2, aes(x=RatioSegment, y=Weight))
              + geom_line(aes(group=Bucket, colour=Bucket))
              + facet_wrap( ~ Bucket, ncol=3)
              + theme(axis.text.x = element_text(angle = 90, hjust = 1))
              + labs(title="Sample Ratios by Weight", x="Ratio Segment", fill="Bucket"))

print(weightPlot)

# theme_set(theme_classic())

weightPlot = (ggplot(data = swTrans2, aes(x=RatioSegment, y=Weight))
              + stat_density(aes(group = Bucket, color = Bucket), position="identity", geom="line")
              + theme(axis.text.x = element_text(angle = 90, hjust = 1))
              + labs(title="Sample Ratios by Weight", x="Ratio Segment", fill="Bucket"))

print(weightPlot)

weightPlot = (ggplot(data = swTrans2, aes(Weight))
              + geom_density(aes(fill=factor(Bucket)), alpha=0.99)
              + theme(axis.text.x = element_text(angle = 90, hjust = 1))
              + labs(title="Sample Ratios by Weight", x="Ratio Segment", fill="Bucket"))

print(weightPlot)






singleBgColours = strip_multi_bg_colours(newSigColours, 20)
print(singleBgColours)
length(singleBgColours)

View(hpcSamples)
sampleCancerTypes2 = hpcSamples %>% select(sampleId, cancerType)
colnames(sampleCancerTypes2) <- c("SampleId", "CancerType")
View(sampleCancerTypes2)
View(sampleCancerTypes)

topBGType = head(sampleCancerTypes,1)$SampleId
print(topBGType)

View(sampleSigData %>% filter(Count>0))
View(sampleSigData %>% filter((!grepl("BG_",SigName)|(grepl("BG_",SigName)&Count>0))))
View(sampleSigData %>% filter(CancerType=='Ovary'&((!grepl("BG_",SigName)|(grepl("BG_",SigName)&Count>0)))))

View(sigColours)
plot_sig_samples(sampleSigData, "Ovary", sigColours, "SV", 4)

baExtData = read.csv("~/dev/nmf/sample_ext_data2.csv")
View(baExtData)

View(hpcSamples)

baExtData2 = merge(baExtData, hpcSamples %>% select(sampleId, cancerType), by.x="SampleId", by.y="sampleId", all.x=T)
View(baExtData2)
baExtData2$Cancer=baExtData2$cancerType
baExtData2$Cancer = ifelse(baExtData2$Cancer=="Other","Minors")
baExtData2 = within(baExtData2, rm(cancerType))
write.csv(baExtData2, "~/dev/nmf/sample_ext_data3.csv", quote=F, row.names=F)



snvBaSigsR25 = as.matrix(read.csv(file="~/dev/nmf/logs/snv_rr_ba_sigs.csv", stringsAsFactors=F))

css = cosine_sim(snvBaSigsR25[,23], snvBaSigsR25[,34])
css = cosine_sim(snvBaSigsR25[,23], snvBaSigsR25[,37])
print(css)
















# close comparison of sig 4 (smoking)
css = cosine_sim(cosmicSignatures[,4], snv30Signatures[,4])
print(css)
sig4Comp = cbind(cosmicSignatures[,4], snv30Signatures[,4])
View(sig4Comp)

View(snvSignatures)
colnames(snvSignatures) = snvSigNamesNum
snvSigSubset = as.matrix(snvSignatures %>% select(1,2,3,5,6,7,8,9,12,13,16,17,18))
View(snvSigSubset)
is.numeric(snvSigSubset)
colnames(snvSigSubset) = snvSimSigNamesNum
write.csv(snvSigSubset, "~/dev/nmf/cosmic_skin_sigs13.csv", row.names=F, quote=F)

# evaluate the java NMF run against the simulated data input

snvSimSignatures = as.matrix(read.csv("~/dev/nmf/logs/snv_nmf_sigs.csv"))
colnames(snvSimSignatures) = snvSimSigNamesNum
snvSimContribs = as.matrix(read.csv("~/dev/nmf/logs/snv_nmf_contribs.csv"))
View(snvSimContribs[,1:50])
View(snvSimSignatures)

View(snvSimSampleCounts %>% group_by(Bucket) %>% summarise(Count=sum(Count)))

evaluate_nmf_data("SNV", "sim_denovo_sig13_skin", snvSigSimCount, snvSimSignatures, snvSimContribs, snvSimMatrixData, snvSimSampleCounts,
                  snvSimCancerTypes, snvBucketNames, snvSimSigNamesNum, snvSimSigNamesStr, FALSE, FALSE, TRUE)

# ref then floating
evaluate_nmf_data("SNV", "sim_denovo_sig13_skin_sig_ref_float", snvSigSimCount, snvSimSignatures, snvSimContribs, snvSimMatrixData, snvSimSampleCounts,
                  snvSimCancerTypes, snvBucketNames, snvSimSigNamesNum, snvSimSigNamesStr, FALSE, FALSE, TRUE)

evaluate_nmf_data("SNV", "sim_nmf_sig13_skin_cosmic_sig_ref_fixed", snvSigSimCount, snvSimSignatures, snvSimContribs, snvSimMatrixData, snvSimSampleCounts,
                  snvSimCancerTypes, snvBucketNames, snvSimSigNamesNum, snvSimSigNamesStr, FALSE, FALSE, TRUE)

# actual skin samples with de-novo

snvSkinSampleCounts = gather(skinMatrixData2, "SampleId", "Count", -Bucket)
View(snvSkinSampleCounts)

snvSkinCancerTypes = as.data.frame(colnames(snvSimContribs))
snvSkinCancerTypes$CancerType = "Skin"
colnames(snvSkinCancerTypes) = c("SampleId", "CancerType")
View(snvSkinCancerTypes)

evaluate_nmf_data("SNV", "actual_denovo_sig13_skin", snvSigSimCount, snvSimSignatures, snvSimContribs, skinMatrixData, snvSkinSampleCounts,
                  snvSkinCancerTypes, snvBucketNames, snvSimSigNamesNum, snvSimSigNamesStr, FALSE, FALSE, TRUE)

# kidney
snvSimMatrixData = as.matrix(read.csv("~/data/r_data/kidney_matrix_data.csv", stringsAsFactors=F))
snvSimMatrixData = kidneyMatrixData
View(snvBucketNames)
snvSimMatrixData2 = cbind(snvBucketNames,snvSimMatrixData)
snvSimSampleCounts = gather(snvSimMatrixData2, "SampleId", "Count", -Bucket)

snvSig10Count = 10
snvSig10NamesNum = get_signame_list(snvSig10Count, F)
snvSig10NamesStr = get_signame_list(snvSig10Count, T)

sigColours = get_sig_colours(5)
View(sigColours)

snvSimSignatures = as.matrix(read.csv("~/dev/nmf/logs/snv_kidney_nmf_sigs.csv"))
# colnames(snvSimSignatures) = snvSig10NamesNum
snvSimContribs = as.matrix(read.csv("~/dev/nmf/logs/snv_kidney_nmf_contribs.csv"))

snvSimCancerTypes = as.data.frame(colnames(snvSimContribs))
snvSimCancerTypes$CancerType = "Kidney"
colnames(snvSimCancerTypes) = c("SampleId", "CancerType")


evaluate_nmf_data("SNV", "actual_denovo_sig10_kidney", snvSimSignatures, snvSimContribs, snvSimMatrixData, snvSimSampleCounts,
                  sampleCancerTypes, snvBucketNames, snvSig10NamesStr, FALSE, FALSE, TRUE)

evaluate_nmf_data("SNV", "actual_cos30_minfit_kidney", snvSimSignatures, snvSimContribs, snvSimMatrixData, snvSimSampleCounts,
                  sampleCancerTypes, snvBucketNames, snvSig30NamesStr, FALSE, FALSE, TRUE)


View(highestPurityCohortSummary)

View(snvSimCancerTypes)





# MMR/MSI actuals fit with Cosmic sigs
snvCosmicMmrContribs = apply_signatures(mmrMatrixData, cosmicSignatures)
snvSig30Count = 30
snvSig30NamesNum = get_signame_list(snvSig30Count, F)
snvSig30NamesStr = get_signame_list(snvSig30Count, T)

colnames(cosmicSignatures) = snvSig30NamesNum

evaluate_nmf_data("SNV", "cosmic30_fit_mmr", snvSig30Count, cosmicSignatures, snvCosmicMmrContribs, mmrMatrixData, mmrSampleCounts,
                  sampleCancerTypes, snvBucketNames, snvSig30NamesNum, snvSig30NamesStr, F, F, T)

snvMmrSigCount = 10
snvMmrSigNum = get_signame_list(snvMmrSigCount, F)
snvMmrSigStr = get_signame_list(snvMmrSigCount, T)
snvSimSignatures = as.matrix(read.csv("~/dev/nmf/logs/snv_mmr_nmf_sigs.csv"))
colnames(snvSimSignatures) = snvMmrSigNum
snvSimContribs = as.matrix(read.csv("~/dev/nmf/logs/snv_mmr_nmf_contribs.csv"))
View(mmrMatrixData)

evaluate_nmf_data("SNV", "denovo_sig10_mmr", snvMmrSigCount, snvSimSignatures, snvSimContribs, mmrMatrixData, mmrSampleCounts,
                  sampleCancerTypes, snvBucketNames, snvMmrSigNum, snvMmrSigStr, F, F, T)


mmrBucketData = get_bucket_data(cosmicSignatures, snvCosmicMmrContribs, snvBucketNames)
mmrBucketData = merge(mmrBucketData,snvBucketNames,by="Bucket",all.x=T)
View(mmrBucketData)



# prostate
snvSimMatrixData = as.matrix(read.csv("~/dev/nmf/logs/snv_prostate_sim_sc.csv", stringsAsFactors=F))
View(snvSimMatrixData[,1:30])
View(snvBucketNames)
snvSimMatrixData2 = cbind(snvBucketNames,snvSimMatrixData)
snvSimSampleCounts = gather(snvSimMatrixData2, "SampleId", "Count", -Bucket)

snvSimCancerTypes = as.data.frame(colnames(snvSimContribs))
snvSimCancerTypes$CancerType = "Prostate"
colnames(snvSimCancerTypes) = c("SampleId", "CancerType")

snvSimSignatures = as.matrix(read.csv("~/dev/nmf/logs/snv_nmf_sigs.csv"))
colnames(snvSimSignatures) = snvSimSigNamesNum
snvSimContribs = as.matrix(read.csv("~/dev/nmf/logs/snv_nmf_contribs.csv"))

evaluate_nmf_data("SNV", "sim_denovo_nosf_sig13_prostate", snvSigSimCount, snvSimSignatures, snvSimContribs, snvSimMatrixData, snvSimSampleCounts,
                  snvSimCancerTypes, snvBucketNames, snvSimSigNamesNum, snvSimSigNamesStr, FALSE, FALSE, TRUE)

# using R's NMF
load("~/data/r_data/snv_sim_prostate_sig13.RData")
View(nmfResult)
snvSimProstateNmfResult = nmfResult
evaluate_nmf_run("SNV", "sim_rnmf_sig13_prostate", snvSigSimCount, snvSimProstateNmfResult, snvSimMatrixData, snvSimSampleCounts,
                  snvSimCancerTypes, snvBucketNames, snvSimSigNamesNum, snvSimSigNamesStr, FALSE, FALSE, TRUE)

snvCosmicProstateContribs = apply_signatures(snvSimMatrixData, snvSigSubset)
colnames(snvSigSubset) = snvSimSigNamesNum

evaluate_nmf_data("SNV", "sim_cosmic_fit_sig13_prostate", snvSigSimCount, snvSigSubset, snvCosmicProstateContribs, snvSimMatrixData, snvSimSampleCounts,
                  snvSimCancerTypes, snvBucketNames, snvSimSigNamesNum, snvSimSigNamesStr, FALSE, FALSE, TRUE)

evaluate_nmf_data("SNV", "sim_sig13_prostate_using_refs", snvSigSimCount, snvSigSubset, snvCosmicProstateContribs, snvSimMatrixData, snvSimSampleCounts,
                  snvSimCancerTypes, snvBucketNames, snvSimSigNamesNum, snvSimSigNamesStr, FALSE, FALSE, TRUE)




# the actual allocation per sig from the simulation run
snvSimSigContribs = as.matrix(read.csv("~/dev/nmf/logs/snv_sim_contributions.csv", stringsAsFactors=F))

snvSimMatrixData = as.matrix(read.csv("~/dev/nmf/logs/snv_sim_sample_counts.csv", stringsAsFactors=F))
View(snvSimMatrixData[,1:30])
View(snvBucketNames)
snvSimMatrixData2 = cbind(snvBucketNames,snvSimMatrixData)
snvSimSampleCounts = gather(snvSimMatrixData2, "SampleId", "Count", -Bucket)
View(snvSimSampleCounts)




snvSimCancerTypes = as.data.frame(colnames(snvSimContribs))
snvSimCancerTypes$CancerType = "Skin"
colnames(snvSimCancerTypes) = c("SampleId", "CancerType")
View(snvSimCancerTypes)

snvSimSampleTotals = snvSimSampleCounts %>% group_by(SampleId) %>% summarise(Count=sum(Count))
View(snvSimSampleTotals)

# de-novo sigs, less than input (ie 13 vs 5)
snvSigSimCount5 = 5
snvSimSigNamesNum5 = get_signame_list(snvSigSimCount5, F)
snvSimSigNamesStr5 = get_signame_list(snvSigSimCount5, T)
colnames(snvSimSignatures) = snvSimSigNamesNum5
evaluate_nmf_data("SNV", "sim_nmf_sig5_skin", snvSigSimCount5, snvSimSignatures, snvSimContribs, snvSimMatrixData, snvSimSampleCounts,
                  snvSimCancerTypes, snvBucketNames, snvSimSigNamesNum5, snvSimSigNamesStr5, FALSE, FALSE, TRUE)


# fit sim counts to cosmic 13
snvFitContribs = apply_signatures(snvSimMatrixData, snvSigSubset)
View(snvFitContribs[,1:20])
nrow(snvFitContribs)
colnames(snvSigSubset) = snvSimSigNamesNum
View(snvSigSubset)
write.csv(snvFitContribs, "~/dev/nmf/cosmic_skin_contribs13.csv", row.names=F, quote=F)
max(snvSigSubset)


# report using cosmic sigs and cosmic fit on simulated data
evaluate_nmf_data("SNV", "sim_sig13_skin_cosmic_ref", snvSigSimCount, snvSigSubset, snvFitContribs, snvSimMatrixData, snvSimSampleCounts,
                  snvSimCancerTypes, snvBucketNames, snvSimSigNamesNum, snvSimSigNamesStr, FALSE, FALSE, TRUE)



# fit to cosmic signatures

snvSignatures = as.data.frame(cosmicSignatures)
View(snvSignatures)
colnames(snvSignatures) = snvSigNamesNum
snvSigSubset = as.matrix(snvSignatures %>% select(1,2,3,5,6,7,8,9,12,13,16,17,18))
View(snvSigSubset)
is.numeric(snvSigSubset)
colnames(snvSigSubset) = snvSimSigNamesNum


# Signature to Age Correlation

# first get sig data
sampleNames = colnames(snvCosmicContribs)

sampleSigData = get_sig_data(cosmicSignatures, snvCosmicContribs, snvSig30NamesStr, sampleNames)
# sampleSigData = get_sig_data(snv30Signatures, snv30Contribs, snvSig30NamesStr, sampleNames)
View(sampleSigData)

# and on the PCAWG reduced fit
sampleSigData = get_sig_data(pcawgSigs, snvPcawgReducedContribs, snvPcawgSigsNamesStr, sampleNames)


View(sampleSigSuppData %>% filter(SigName=="06"&Count>0))


nrow(highestPurityCohortSummary)
# hpcSampleData = highestPurityCohortSummary %>% select(sampleId, patientId, cancerType, cancerSubtype, ageAtBiopsy, msiStatus)
hpcSampleData = hpcSamples %>% select(sampleId, patientId, cancerType, cancerSubtype, ageAtBiopsy, msiStatus)
colnames(hpcSampleData) = c("SampleId", "PatientId", "CancerType", "CancerSubtype", "Age", "MsiStatus")
View(hpcSampleData)

sampleSigSuppData = merge(sampleSigData, hpcSampleData, by="SampleId", all.x=T)
cancerTypesList = unique(sampleSigSuppData$CancerType)
View(cancerTypesList)

View()
View(sampleSigSuppData %>% filter(SigPercent > 0 & grepl("BG_",SigName)))

testSigName = "05"

outputFile = paste("~/logs/r_output/pdfs/snv_", "denovo_sigs", "_age_correlations.pdf", sep='')
# outputFile = paste("~/logs/r_output/pdfs/snv_", "denovo_sig", testSigName, "new_bgs_age_correlations.pdf", sep='')
pdf(file=outputFile, height = 14, width = 20)
par(mar=c(1,1,1,1))

plotList = list()
plotCount = 1
index = 1
for(cancerType in cancerTypesList)
{
  # sigAgeData = sampleSigSuppData %>% filter(SigName==testSigName&CancerType==cancerType&Count>0)
  sigAgeData = sampleSigSuppData %>% filter(SigPercent>0&grepl("BG_",SigName)&CancerType==cancerType&Count>0)

  # cap any sample with a count-age ratio more than 2x the median ratio
  sigAgeData$Gradient = sigAgeData$Count / sigAgeData$Age
  medGradient = median(sigAgeData$Gradient)

  sigAgeData = sigAgeData %>% filter(sigAgeData$Gradient <= 2*medGradient)

  spearman = calc_spearman(sigAgeData)

  print(paste(cancerType, ": sampleCount=", nrow(sigAgeData), ", spearman=", round(spearman,4), sep=''))
  # sigAgeData$AdjCount = ifelse(sigAgeData$Gradient>2*medGradient,2*medGradient*sigAgeData$Age,sigAgeData$Count)

  plot = (ggplot(data=sigAgeData, aes(x=Age, y=Count))
          + geom_point()
          + xlim(0,100)
          + ggtitle(cancerType))

  plotList[[plotCount]] = plot

  if(plotCount >= 12 || length(cancerTypesList) == index)
  {
    multiplot(plotlist = plotList, cols=3)
    plotList = list()
    plotCount = 1
  }
  else
  {
    plotCount = plotCount + 1
  }

  index = index + 1
}

dev.off()


print(plotList[[7]])

multiplot(plotlist = plotList, cols=3)

sigAgeData = sampleSigSuppData %>% filter(SigName=="01"&CancerType=="Colon/Rectum"&Count>0)

calc_spearman<-function(sigAgeData)
{
  # calc spearman's correlation
  sigAgeData = sigAgeData %>% arrange(Count)
  rank = data.frame(as.numeric(as.character(rownames(sigAgeData))))
  colnames(rank) <- c("CountRank")
  sigAgeData = cbind(sigAgeData, rank)

  sigAgeData = sigAgeData %>% arrange(Age)
  rank = data.frame(as.numeric(as.character(rownames(sigAgeData))))
  colnames(rank) <- c("AgeRank")
  sigAgeData = cbind(sigAgeData, rank)

  sigAgeData$dSqrd = (sigAgeData$CountRank - sigAgeData$AgeRank)^2
  sumSqrd = sum(sigAgeData$dSqrd)
  sampleCount = nrow(sigAgeData)
  sprearman = 1 - (6 *sumSqrd / (sampleCount*(sampleCount^2-1)))

  return (sprearman)
}

sigAgeData = sampleSigSuppData %>% filter(SigName=="01"&CancerType=="Skin")
sigAgeData$Gradient = sigAgeData$Count / sigAgeData$Age
medGradient = median(sigAgeData$Gradient)
sigAgeData$Gradient = ifelse(sigAgeData$Gradient>2*medGradient,2*medGradient,sigAgeData$Gradient)
sigAgeData$AdjCount = ifelse(sigAgeData$Gradient>2*medGradient,2*medGradient*sigAgeData$Age,sigAgeData$Count)

View(sigAgeData)

plot = (ggplot(data=sigAgeData, aes(x=Age, y=AdjCount))
        + geom_point())

print(plot)

+
  geom_segment(aes(x = 1e2, xend=1e6, y = 12400, yend=12380), linetype = "dashed") + annotate("text", x = 1e2, y = 15000, label = "MSI Threshold", size = 3, hjust = 0) +
  geom_segment(aes(y = 1e2, yend=1e6, x = 30950, xend=30950), linetype = "dashed") + annotate("text", x = 32000, y = 1.1e2, label = "TMB High Threshold", size = 3, hjust = 0) +
  geom_point(aes(color = cancerType)) +
  scale_color_manual(values = cancerTypeColours) +
  scale_x_continuous(trans="log10", limits = c(1e2, 1e6)) +
  scale_y_continuous(trans="log10", limits = c(1e2, 1e6)) +
  theme(legend.position = "right", legend.title = element_blank()) + guides(colour = guide_legend(ncol = 1)) +
  xlab("SNVs") + ylab("Indels") + ggtitle("")









# estimates
load("~/data/snvNmfEstimate_highML_10_15.RData")
View(snvNmfHighEstimate)
plot(snvNmfHighEstimate)

load("~/data/snvNmfResult_highML_sig13.RData")
#load("~/data/snvNmfResult_highML_sig30.RData")
#load("~/data/snvNmfResult_lowML_sig30.RData")
View(snvNmfHighResult)

snvHighSigCount = 13
sig13NamesNums = get_signame_list(13, F)
sig13NamesStr = get_signame_list(13, T)





highMLSampleCounts = snvSampleCounts %>% filter(SampleId %in% highMutLoadSamples$SampleId)
highMLSampleCancerTypes = sampleCancerTypes %>% filter(SampleId %in% highMutLoadSamples$SampleId)
write.csv(highMLSampleCounts, "~/logs/r_output/snvHighMLCounts.csv", row.names=F, quote=F)
View(highMLSampleCounts)

evaluate_nmf_run(
  "SNV", "sig13_highML", snvHighSigCount, snvNmfHighResult, highMLSampleCounts, highMLSampleCancerTypes,
  snvBucketNames, sig13NamesNums, sig13NamesStr, TRUE, FALSE)



# least squares fit

sig1 = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)
sig2 = c(0.3, 0.2, 0.1, 0.2, 0.3, 0.4)
sig3 = c(0.4, 0.5, 0.6, 0.6, 0.4, 0.2)
model = matrix(sig1, nrow=6, ncol=1)
model = cbind(model, sig2)
model = cbind(model, sig3)
colnames(model) = c(1,2,3)
View(model)
is.matrix(model)

data = matrix(c(16, 18, 22, 28, 30, 32), nrow=6, ncol=1)

lsr = calcLeastSquares(model, data[,1])
print(lsr$x)
print(lsr$resid.norm)

result = c()

C = model
d = data[,1]

is.matrix(C)
is.vector(d)
length(d)

output = model %*% lsr$x
View(output)

print(eps())

normC = norm(model, type="2")
print(normC)


cols = nrow(C)
rows = ncol(C)
print(cols)
print(rows)
tol = 10 * eps() * norm(C, type = "2") * (max(cols, rows) + 1)
print(tol)

x = rep(0, rows)  # empty result vector, 6 rows
P = logical(rows); # vector of FALSE values, 6 rows
Z = !P   # vector of TRUE values, 6 rows
View(x)
View(P)
View(Z)

# calc initial residuals with zeros for result
resid = d - C %*% x
print(resid) # residuals per original bucket, in this case the exact inputs since the result ie 'x' is all zeros
w <- t(C) %*% resid # purpose of multiplying transpose of C by residuals? some cost function?
print(w) # 3 rows, 1 column, the current
wz <- numeric(cols)
print(wz) # empty vector of 6 x zero


# iteration parameters
outeriter = 0;
it = 0
itmax = 5 * cols # is 30 times
exitflag <- 1
print(itmax)


# single iteration
z <- numeric(cols) # a vector of 6 zeros
print(z)
wz <- rep(-Inf, cols) # a vector of 6 -ve infinities
print(wz)
print(w[Z])
print(wz[Z])
wz[Z] <- w[Z] # allocates the contents of w, which is the set of 3 residual-types values, twice into the 6 size vector??
print(wz)

im <- which.max(wz) # find the index of the largest value
print(im)

# the index of this max value in the P vector (originally all falses) to TRUE, and reverse to Z
P[im] <- TRUE
Z[im] <- FALSE
print(P)
print(Z)

print(z[P])
print(z[Z])
print(z)

z[P] <- qr.solve(C[, P], d)


# internal iteration, controlled by only the while loop condition



while (any(Z) && any(w[Z] > tol)) # while any of the Z values remain true and any of the w vector values are greater than tol (some min, in this case 2e-14)
{
  z <- numeric(cols) # a vector of 6 zeros
  wz <- rep(-Inf, cols)
  wz[Z] <- w[Z]
  im <- which.max(wz)
  P[im] <- TRUE; Z[im] <- FALSE
  z[P] <- qr.solve(C[, P], d)

  while (any(z[P] <= 0))
  {
    it <- it + 1
    if(it > itmax)
    {
      print(paste("iteration count exceed at i=", i, sep=''))
      break;
      # stop("Iteration count exceeded")
    }

    Q <- (z <= 0) & P
    alpha <- min(x[Q] / (x[Q] - z[Q]))
    x <- x + alpha*(z - x)
    Z <- ((abs(x) < tol) & P) | Z
    P <- !Z
    z <- numeric(cols)
    z[P] <- qr.solve(C[, P], d)
  }
  x <- z
  resid <- d - C %*% x
  w <- t(C) %*% resid
}



calcLeastSquares<-function(C, d)
{
  rows <- nrow(C);
  cols <- ncol(C)

  tol = 10 * eps() * norm(C, type = "2") * (max(cols, rows) + 1)

  x  <- rep(0, cols)             # initial point
  P  <- logical(cols); Z <- !P   # non-active / active columns

  resid <- d - C %*% x
  w <- t(C) %*% resid
  wz <- numeric(cols)

  # iteration parameters
  outeriter <- 0; it <- 0
  itmax <- 5 * cols; exitflag <- 1

  while (any(Z) && any(w[Z] > tol))
  {
    outeriter <- outeriter + 1
    z <- numeric(cols)
    wz <- rep(-Inf, cols)
    wz[Z] <- w[Z]
    im <- which.max(wz)
    P[im] <- TRUE; Z[im] <- FALSE
    z[P] <- qr.solve(C[, P], d)

    while (any(z[P] <= 0))
    {
      it <- it + 1
      if(it > itmax)
      {
        print(paste("iteration count exceed at i=", i, sep=''))
        break;
        # stop("Iteration count exceeded")
      }

      Q <- (z <= 0) & P
      alpha <- min(x[Q] / (x[Q] - z[Q]))
      x <- x + alpha*(z - x)
      Z <- ((abs(x) < tol) & P) | Z
      P <- !Z
      z <- numeric(cols)
      z[P] <- qr.solve(C[, P], d)
    }
    x <- z
    resid <- d - C %*% x
    w <- t(C) %*% resid
  }
  return(list(x = x, resid.norm = sum(resid*resid)))
}












