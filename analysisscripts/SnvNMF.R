library(devtools) #; install_github("im3sanger/dndscv")
library(MutationalPatterns)
library(RMySQL)
library(data.table)
library("NMF")

# plotting
library(grid)
library(gridExtra)
library(ggplot2)

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

## NMF functions

# restrict to samples in high-purity and DR022 set
load("~/data/highestPurityCohortSummary.RData")
dr022Samples = read.csv("~/data/DR-022_metadata.tsv", sep='\t')
hpcSamples = highestPurityCohortSummary %>% filter(sampleId %in% dr022Samples$sampleId)
# hpcSamples = highestPurityCohortSummary
nrow(hpcSamples)

cosmicSignatures = getCOSMICSignatures()
View(cosmicSignatures)

dbProd = dbConnect(MySQL(), user='hmf', password='HMFhmf@1', dbname='hmfpatients', groups = "RAnalysis")
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

# colIndex = totalColIndex
colIndex = subclonalColIndex
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

sum(snvHighMatrixData$CPCT02010503TII)
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

snvSigNamesNum = get_signame_list(30, F)
snvSigNamesStr = get_signame_list(30, T)
snvSigCount = 30
print(snvSigNamesNum)
print(snvSigNamesStr)

snvSignatures = cosmicSignatures
colnames(snvSignatures) <- NULL
rownames(snvContribsTotal) <- NULL
rownames(snvContribsSubclonal) <- NULL

# evaluate_nmf_run("SNV", "sig30", snvSigCount, snvNmfResult, snvMatrixData, snvSampleCounts,
#                  sampleCancerTypes, snvBucketNames, snvSigNamesNum, snvSigNamesStr, FALSE, FALSE)

View(snvSampleCounts)

# write.csv(snvMatrixData, file="~/logs/r_output/snv_nmf_matrix_data.csv", row.names=F, quote=F)
# write.csv(snvMatrixData, file="~/logs/r_output/snv_nmf_subc_matrix_data.csv", row.names=F, quote=F)
# write.csv(snvSampleCounts, file="~/logs/r_output/snv_nmf_sample_counts.csv", row.names=F, quote=F)
# write.csv(snvContribsTotal, file="~/logs/r_output/snv_nmf_contribs.csv", row.names=F, quote=F)
# write.csv(snvSampleCounts, file="~/logs/r_output/snv_nmf_subc_sample_counts.csv", row.names=F, quote=F)
# write.csv(snvContribsSubclonal, file="~/logs/r_output/snv_nmf_subc_contribs.csv", row.names=F, quote=F)


# evaluation using cosmic signatures
evaluate_nmf_data("SNV", "sig30_cosmic", snvSigCount, snvSignatures, snvContribsTotal, snvMatrixData, snvSampleCounts,
                 sampleCancerTypes, snvBucketNames, snvSigNamesNum, snvSigNamesStr, TRUE, FALSE, TRUE)



# subclonal evaluation
evaluate_nmf_data("SNV_Subclonal", "sig30_cosmic", snvSigCount, snvSignatures, snvContribsSubclonal, snvMatrixData, snvSampleCounts,
                  sampleCancerTypes, snvBucketNames, snvSigNamesNum, snvSigNamesStr, TRUE, FALSE, TRUE)


snvSampleSigData = get_sig_data(snvSignatures, snvContribsTotal, snvSigNamesStr, colnames(snvContribsTotal))
snvSampleSigData$Count = round(snvSampleSigData$Count,0)
View(snvSampleSigData)
write.csv(snvSampleSigData, "~/logs/r_output/snvSampleSigData.csv", row.names=F, quote=F)


View(snvContribsSubclonal[,1:10])

sampleNames = colnames(snvContribsSubclonal)
View(sampleNames)
origSampleCounts = snvSampleCounts %>% group_by(SampleId) %>% summarise(OrigSampleCount=sum(Count))
View(origSampleCounts)
sampleBucketData = get_sample_bucket_data(snvMatrixData, origSampleCounts, sampleCancerTypes, snvBucketNames)
View(sampleBucketData)

plot_sample_bucket_contrib(sampleBucketData, "", snvBucketCount, "SNV", 0)



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

signatures = NMF::basis(snvNmfHighResult)
contribution = NMF::coef(snvNmfHighResult)

# TMP: bucket count plot:

# snvBucketData = get_bucket_data(signatures, contribution, snvBucketNames)
# View(snvBucketData)
#
# origSampleCounts = highMLSampleCounts %>% group_by(SampleId) %>% summarise(OrigSampleCount=sum(Count))
# View(origSampleCounts)
# View(highMLSampleCancerTypes)
#
# sampleBucketTopN = get_top_buckets_by_sample(highMLSampleCounts, origSampleCounts, highMLSampleCancerTypes, 0)
# View(sampleBucketTopN)
# tmp1 = sampleBucketTopN %>% filter(SampleId=="CPCT02030213T")
# sum(tmp1$Count)

plot_sample_bucket_contrib(sampleBucketTopN, "", 96, "SNV", 2)

print(plot_contribution(snvHMLContribution, snvNmfHighResult$signature, mode = "relative"))
print(plot_contribution(snvHMLContribution, signatures, mode = "absolute"))

snvHMLContribution = contribution
row.names(snvHMLContribution) <- sig13NamesStr
View(snvHMLContribution)

plot_contribution_heatmap(snvHMLContribution, cluster_samples = TRUE, method = "complete")

plot_compare_profiles(snvHighMatrixData[,1], snvNmfHighResult$reconstructed[,1], profile_names = c("Original", "Reconstructed"), condensed = TRUE)

View(snvNmfHighResult)
View(snvNmfHighResult$reconstructed)

View(snvBucketNames)

snvHighMLResiduals = calc_sample_residuals(contribution, signatures, snvBucketNames, highMLSampleCounts)
write.csv(snvHighMLResiduals, "~/logs/r_output/snvHighMLResiduals.csv", row.names=F, quote=F)
View(snvHighMLResiduals)
sum(snvHighMLResiduals$Count)
sum(snvHighMLResiduals$ResidualTotal)
print(residuals(snvNmfHighResult))


###################
# Cosine Similarity

cosine_sim<-function(vec1, vec2)
{
  cosineSim = (vec1 %*% vec2) / (sqrt(vec1 %*% vec1)*sqrt(vec2 %*% vec2))
  return (cosineSim)
}

calc_cosine_sims<-function(dataMatrix, cssCutoff, itemName, logMatches=F)
{
  numItems = ncol(dataMatrix)
  itemNames = colnames(dataMatrix)

  cssResults = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(cssResults) <- c(paste(itemName,'1',sep=''), paste(itemName,'2',sep=''), "CSS")

  for(i in 1:numItems)
  {
    data1 = dataMatrix[,i]

    for(j in i+1:numItems)
    {
      if(j>numItems)
        break

      data2 = dataMatrix[,j]

      css = cosine_sim(data1,data2)

      if(css >= cssCutoff)
      {
        if(logMatches)
        {
          print(paste(i, "=", itemNames[i], ", count=", sum(data1), ", ", j, "=", itemNames[j], ", count=", sum(data2), ", CSS=", round(css,4), sep=''))
        }

        rowIndex = nrow(cssResults)+1
        cssResults[rowIndex,1] = itemNames[i]
        cssResults[rowIndex,2] = itemNames[j]
        cssResults[rowIndex,3] = round(css,6)
      }
    }

    if(i %% 100 == 0)
    {
      print(paste("processed ", i, sep=''))
    }
  }

  return (cssResults)
}

# all samples, all CSS values - note this will be 2400^2/2 entries
snvAllCssResults = calc_cosine_sims(snvMatrixData, 0.5, "Sample", F)
View(snvAllCssResults)

snvHighMLCssResults = calc_cosine_sims(snvHighMatrixData, 0.5, "Sample", F)
View(snvHighMLCssResults)

snvHighMLCssResults = merge(snvHighMLCssResults, sampleCancerTypes, by.x="Sample1", by.y="SampleId", ,all.x=TRUE)
snvHighMLCssResults = merge(snvHighMLCssResults, sampleCancerTypes, by.x="Sample2", by.y="SampleId", ,all.x=TRUE)
snvHighMLCssResults = snvHighMLCssResults %>% arrange(-CSS)
colnames(snvHighMLCssResults) <- c("Sample1", "Sample2", "CSS", "CancerType1", "CancerType2")
write.csv(snvHighMLCssResults, "~/logs/r_output/snvHighMLCssResults.csv")

numSamples = ncol(snvHighMatrixData)
snvHighMLSampleNames = colnames(snvHighMatrixData)
View(snvHighMLSampleNames)

View(snvBucketNames)

cssRCount = 100
cssRInv = 1/cssRCount
cssResultGroups = data.frame(matrix(ncol = 2, nrow = cssRCount))
colnames(cssResultGroups) <- c("CSSBand", "Count")

for(i in 1:nrow(cssResultGroups))
{
  cssResultGroups[i,1] = i*cssRInv
  # cssResultGroups[i,2] = 0
}



View(cssResultGroups)


print(cosine_sim(snvHighMatrixData$CPCT02060041T, snvHighMatrixData$DRUP01230001T))

View(cosineSimResults)

View(snvHighMatrixData %>% select(CPCT02060041T,DRUP01230001T))
View(highMLSampleCancerTypes %>% filter(SampleId=="CPCT02060041T"|SampleId=="DRUP01230001T"))
View(sampleSigData %>% filter(SampleId=="CPCT02060041T"|SampleId=="DRUP01230001T"))

View(sampleSigData)

View(snvHighMatrixData)

# same again but for buckets

cssBucketResults = data.frame(matrix(ncol = 3, nrow = 0))
colnames(cssBucketResults) <- c("Bucket1", "Bucket2", "CSS")

cssRCount = 100
cssRInv = 1/cssRCount
cssBucketResultGroups = data.frame(matrix(ncol = 2, nrow = cssRCount))
colnames(cssBucketResultGroups) <- c("CSSBand", "Count")

for(i in 1:nrow(cssBucketResultGroups))
{
  cssBucketResultGroups[i,1] = i*cssRInv
  cssBucketResultGroups[i,2] = 0
}

numBuckets = nrow(snvBucketNames)
snvHighMatrixDataTrans = t(snvHighMatrixData)
# ncol(snvHighMatrixDataTrans)
# nrow(snvHighMatrixDataTrans)
# View(snvHighMatrixDataTrans[,1:10])
rownames(snvHighMatrixDataTrans) <- NULL

View(snvBucketNames)
snvBucketNamesVec = snvBucketNames$Bucket

for(i in 1:numBuckets)
{
  data1 = snvHighMatrixDataTrans[,i]

  for(j in i+1:numBuckets)
  {
    if(j>numBuckets)
      break

    data2 = snvHighMatrixDataTrans[,j]

    css = cosine_sim(data1,data2)

    # print(paste("s1=", snvHighMLSampleNames[i], ", count=", sum(sam1Data), ", s2=", snvHighMLSampleNames[j], ", count=", sum(sam2Data), sep=''))
    # print(paste(i, j, css, sep=', '))

    if(css >= 0.8)
    {
      print(paste("b1=", snvBucketNamesVec[i], ", count=", sum(data1), ", b2=", snvBucketNamesVec[j], ", count=", sum(data2), sep=''))
      print(paste(i, j, css, sep=', '))
      rowIndex = nrow(cssBucketResults)+1
      cssBucketResults[rowIndex,1] = snvBucketNamesVec[i]
      cssBucketResults[rowIndex,2] = snvBucketNamesVec[j]
      cssBucketResults[rowIndex,3] = round(css,6)
    }

    cssRounded = round(css/cssRInv)
    cssBucketResultGroups[cssRounded,2] = cssBucketResultGroups[cssRounded,2]+1
  }
}

View(cssBucketResultGroups)
View(cssBucketResults)


