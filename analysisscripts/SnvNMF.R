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

load("~/data/highestPurityCohortSummary.RData")
View(highestPurityCohortSummary)

cosmicSignatures = getCOSMICSignatures()
View(cosmicSignatures)

dbProd = dbConnect(MySQL(), user='hmf', password='HMFhmf@1', dbname='hmfpatients', groups = "RAnalysis")
# dbDisconnect(dbProd)

cancerTypes = highestPurityCohortSummary %>% group_by(cancerType) %>% count()
View(cancerTypes)

# download SNVs by cancer type (to keep the queries smaller)
for(i in 1:nrow(cancerTypes))
{
  cancerTypeRow = cancerTypes[i,]
  cancerTypeStr = cancerTypeRow$cancerType

  print(paste(i, ": retrieving SNVs for cancer=", cancerTypeStr, sep=''))

  sampleIdsStr = getSampleIdsStr(highestPurityCohortSummary %>% filter(cancerType==cancerTypeStr))
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

samplesList = unique(highestPurityCohortSummary$sampleId)
View(samplesList)

# load("~/logs/r_output/snvSampleCounts.RData")
# View(result)

# sampleCountResults = result
index = 1
for(s in samplesList)
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
sampleSigs = list()
for(s in sampleCountResults)
{
  print(paste("fitting signatures for ", s, sep=''))

  if (!is.null(sampleCountResults[[s]]))
  {
    # we need to slice out only the mutation count columns (delete col 1 and 2)
    res = fit_to_signatures(sampleCountResults[[s]][, -c(1, 2)], cosmicSignatures)
    sampleSigs[[s]] <- res$contribution
  }
}

View(sampleSigs)
View(cosmicSignatures)

contributionTotal = data.frame(matrix(ncol = 0, nrow = 30))
contributionClonal = data.frame(matrix(ncol = 0, nrow = 30))
contributionSubclonal = data.frame(matrix(ncol = 0, nrow = 30))

for(s in samplesList)
{
  # print(paste("extracting sig data for ", s, sep=''))

  if (!is.null(sampleSigs[[s]]))
  {
    sampleSigData = sampleSigs[[s]]
    totalValues = sampleSigData[,1]
    contributionTotal = cbind(contributionTotal, totalValues)

    subclonalValues = sampleSigData[,2]
    contributionSubclonal = cbind(contributionSubclonal, subclonalValues)

    clonalValues = sampleSigData[,3]
    contributionClonal = cbind(contributionClonal, clonalValues)
  }
}

contributionTotal = setNames(contributionTotal, samplesList)
contributionSubclonal = setNames(contributionSubclonal, samplesList)
contributionClonal = setNames(contributionClonal, samplesList)
# View(contributionTotal)


# SNV NMF Standard Routine

# Signature Creation (ie independently from COSMIC)

# create bucket names
snvBucketNames = create_empty_signature()
snvBucketNames = unite(snvBucketNames, "Bucket", type, context, sep='_')
View(snvBucketNames)

# create sample counts per bucket
bucketCount = nrow(snvBucketNames)
# print(bucketCount)

snvMatrixData = data.frame(matrix(ncol = 0, nrow = 96))
snvSampleTotals = data.frame(matrix(ncol = 2, nrow = 0))
snvSampleTotals = setNames(snvSampleTotals, c("SampleId", "SampleCount"))
snvSampleCounts = data.frame(matrix(ncol = 3, nrow = 0))

totalColIndex = 3 # (could check on the fly)

sampleNames = samplesList

for(s in sampleNames)
{
  # print(paste("extracting sig data for ", s, sep=''))
  sampleCounts = sampleCountResults[[s]]

  snvMatrixData = cbind(snvMatrixData, sampleCounts[,totalColIndex])

  sampleCountData = data.frame(matrix(ncol = 0, nrow=bucketCount))
  sampleCountData$SampleId = s
  sampleCountData = cbind(sampleCountData, snvBucketNames)
  sampleCountData = cbind(sampleCountData, sampleCounts[,totalColIndex])
  snvSampleCounts = rbind(snvSampleCounts, sampleCountData)

  rowIndex = nrow(snvSampleTotals)+1
  snvSampleTotals[rowIndex,1] = s
  snvSampleTotals[rowIndex,2] = sum(sampleCounts[,totalColIndex])
  # index = index + 1
  # if(index > 10)
  #   break
}

snvMatrixData = setNames(snvMatrixData, sampleNames)
snvSampleCounts = setNames(snvSampleCounts, c("SampleId", "Bucket", "Count"))

snvSampleTotals = snvSampleTotals %>% arrange(-SampleCount)
View(snvSampleTotals)
View(snvSampleCounts)

# View(sampleBucketCounts)
View(snvMatrixData[,1:10])
nrow(snvMatrixData)
ncol(snvMatrixData)

write.csv(snvMatrixData, file="~/logs/r_output/snv_nmf_counts.csv", row.names=F, quote=F)

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

sigNamesUnamed = get_signame_list(30, F)
sigNamesNamed = get_signame_list(30, T)

snvSigCount = 30
bucketCounts = snvSampleCounts %>% group_by(Bucket) %>% count()

evaluate_nmf_run("SNV", "sig30", snvSigCount, snvNmfResult, snvSampleCounts, sampleCancerTypes, snvBucketNames,
                 sigNamesUnamed, sigNamesNamed, TRUE, FALSE)

print(get_signame_list(10, F))

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

evaluate_nmf_run(
  "SNV", "sig13_highML", snvHighSigCount, snvNmfHighResult, highMLSampleCounts, highMLSampleCancerTypes,
  snvBucketNames, sig13NamesNums, sig13NamesStr, TRUE, FALSE)

View(sampleNames)

signatures = NMF::basis(snvNmfHighResult)
contribution = NMF::coef(snvNmfHighResult)
sampleNames = colnames(contribution)

View(snvBucketNames)

snvHighMLResiduals = calc_sample_residuals(contribution, signatures, snvBucketNames, highMLSampleCounts)
write.csv(snvHighMLResiduals, "~/logs/r_output/snvHighMLResiduals.csv", row.names=F, quote=F)
View(snvHighMLResiduals)
sum(snvHighMLResiduals$Count)
sum(snvHighMLResiduals$ResidualTotal)
print(residuals(snvNmfHighResult))


# TMP: testing plotting residuals
sampleSigData = get_sig_data(signatures, contribution, sig13NamesStr, sampleNames)
sampleResiduals = snvHighMLResiduals %>% select(SampleId,ResidualTotal)
colnames(sampleResiduals) <- c("SampleId","Count")
sampleResiduals$Count = round(sampleResiduals$Count,0)
sampleResiduals$SigName = "Residual"
sampleResiduals$SigPercent = 0
sampleResiduals$PercBucket = 0
View(sampleResiduals)
sampleSigData2 = rbind(sampleSigData, sampleResiduals %>% select(SampleId,SigName,SigPercent,PercBucket,Count))
sampleSigData2 = merge(sampleSigData2, highMLSampleCancerTypes,by.x="SampleId",by.y="SampleId",all.x=TRUE)

sampleSigCounts = sampleSigData2 %>% group_by(SampleId) %>% summarise(SampleCount=sum(Count))
sampleSigData2 = merge(sampleSigData2, sampleSigCounts, by.x="SampleId",by.y="SampleId",all.x=TRUE)
View(sampleSigData2)

plot_sig_samples(sampleSigData2, "", get_sig_colours(snvHighSigCount), 'SNV') # all samples




# Cosine Simarity

numSamples = ncol(snvHighMatrixData)
snvHighMLSampleNames = colnames(snvHighMatrixData)
View(snvHighMLSampleNames)

View(snvBucketNames)

cosineSimResults = data.frame(matrix(ncol = 3, nrow = 0))
colnames(cosineSimResults) <- c("Sample1", "Sample2", "CosineSim")

cssRCount = 100
cssRInv = 1/cssRCount
cssResultGroups = data.frame(matrix(ncol = 2, nrow = cssRCount))
colnames(cssResultGroups) <- c("CSSBand", "Count")

for(i in 1:nrow(cssResultGroups))
{
  cssResultGroups[i,1] = i*cssRInv
  # cssResultGroups[i,2] = 0
}

for(i in 1:numSamples)
{
  sam1Data = snvHighMatrixData[,i]

  for(j in i+1:numSamples)
  {
    if(j>numSamples)
      break

    sam2Data = snvHighMatrixData[,j]

    css = cosine_sim(sam1Data,sam2Data)

    # print(paste("s1=", snvHighMLSampleNames[i], ", count=", sum(sam1Data), ", s2=", snvHighMLSampleNames[j], ", count=", sum(sam2Data), sep=''))
    # print(paste(i, j, css, sep=', '))

    if(css >= 0.8)
    {
      # print(paste("s1=", snvHighMLSampleNames[i], ", count=", sum(sam1Data), ", s2=", snvHighMLSampleNames[j], ", count=", sum(sam2Data), sep=''))
      # print(paste(i, j, css, sep=', '))
      rowIndex = nrow(cosineSimResults)+1
      cosineSimResults[rowIndex,1] = snvHighMLSampleNames[i]
      cosineSimResults[rowIndex,2] = snvHighMLSampleNames[j]
      cosineSimResults[rowIndex,3] = round(css,6)
    }

    cssRounded = round(css/cssRInv)
    cssResultGroups[cssRounded,2] = cssResultGroups[cssRounded,2]+1
  }
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

View(t(data1))
View(data2)
td1 = t(data1)
rownames(td1) <- NULL
colnames(td1) <- NULL
td2 = t(data2)
rownames(td2) <- NULL
colnames(td2) <- NULL
tmp1 = td1 %*% td2

tmpTrans = t(snvHighMatrixData)
View(tmpTrans)

View(td2)

View(sam1Data)


View(cssBucketResults)


View(snvHighMatrixData)
tmpSamples = snvHighMatrixData %>% select(CPCT02060041T,DRUP01230001T)
colnames(tmpSamples) <- c("A", "B")
View(tmpSamples)

sam1Values = tmpSamples[,1]
sam2Values = tmpSamples[,2]


sam3Values = round(sam1Values * 0.1,0)
sam4Values = round(sam2Values * 0.1,0)


print(cosine_sim(sam1Values,sam2Values))
print(cosine_sim(sam3Values,sam4Values))
print(cosine_sim(sam4Values,sam3Values))

abDP = sam1Values %*% sam2Values
aaDP = sam1Values %*% sam1Values
bbDP = sam2Values %*% sam2Values
cosSim2 = abDP/(sqrt(aaDP)*sqrt(bbDP))
print(cosSim2)
print(aaDP)

cosine_sim<-function(vec1, vec2)
{
  cosineSim = (vec1 %*% vec2) / (sqrt(vec1 %*% vec1)*sqrt(vec2 %*% vec2))
  return (cosineSim)
}

cosSim3 = (tmpSamples$A %*% tmpSamples$B) / (sqrt(tmpSamples$A %*% tmpSamples$A)*sqrt(tmpSamples$B %*% tmpSamples$B))
print(cosSim3)

