library(devtools)
library(RMySQL)
library(data.table)

library(purple);
library(IRanges)
library(dplyr)
library(tidyr)
library(stringi)

# NMF related
library(NMF)
library(MutationalPatterns)
# detach("package:svnmf", unload=TRUE);
library(svnmf)

# plotting
library(grid)
library(gridExtra)
library(ggplot2)
library(ggpubr)

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

  return(queryResults)
}

baseConvert<-function(base)
{
  if(base=="A")
    newBase = "T"
  else if(base=="T")
    newBase = "A"
  else if(base=="C")
    newBase = "G"
  else if(base=="G")
    newBase = "C"

  return (newBase)
}

# convert the pairs so as to match the PCAWG set
setMutationAdjusted<-function(ref, alt)
{
  # if(!is.integer(length))
  #   return ("Invalid")
  length = stri_length(ref)

  if(length==2)
  {
    if(ref=="CA"|ref=="AA"|ref=="AG"|ref=="GA"|ref=="GG"|ref=="GT")
    {
      mutation = paste(baseConvert(substr(ref,2,2)), baseConvert(substr(ref,1,1)), '>', baseConvert(substr(alt,2,2)), baseConvert(substr(alt,1,1)), sep='')
      return (mutation)
    }

    if(ref=="AT"|ref=="TA"|ref=="CG"|ref=="GC")
    {
      # these refs will all remain the same with a reversal, so then convert the alts
      if(alt=="CA"|alt=="AA"|alt=="AG"|alt=="GA"|alt=="GG"|alt=="GT")
      {
        mutation = paste(ref, '>', baseConvert(substr(alt,2,2)), baseConvert(substr(alt,1,1)), sep='')
        return (mutation)
      }
    }

    mutation = paste(ref, '>', alt, sep='')
    return (mutation)
  }

  if(length==3)
  {
    return("3Bases")
  }
  else
  {
    return("4+Bases")
  }
}

create_mnv_fields<-function(mnvData)
{
  mnvData$Length = stringi::stri_length(mnvData$Alt)
  mnvData$Mutation = ifelse(mnvData$Length==2,paste(mnvData$Ref, mnvData$Alt, sep='>'), ifelse(mnvData$Length==3,"3Bases","4+Bases"))
  mnvData$Bucket = apply(mnvData[,c('Ref','Alt')], 1, function(x) setMutationAdjusted(x[1],x[2]))

  return (mnvData)
}


load("~/data/highestPurityCohortSummary.RData")
View(highestPurityCohortSummary)

dr022Samples = read.csv("~/data/DR-022_metadata.tsv", sep='\t')
nrow(dr022Samples)
View(dr022Samples)

hpcSamples = highestPurityCohortSummary %>% filter(sampleId %in% dr022Samples$sampleId)
# hpcSamples = highestPurityCohortSummary
nrow(hpcSamples)

# dbDisconnect(dbProd)

cancerTypes = hpcSamples %>% group_by(cancerType) %>% count()
cancerTypes = setNames(cancerTypes, c("CancerType", "Count"))
sampleCancerTypes = hpcSamples %>% select(sampleId, cancerType)
sampleCancerTypes = setNames(sampleCancerTypes, c("SampleId", "CancerType"))
View(sampleCancerTypes)
View(cancerTypes)


# MNV Handling

# download from prod around 1.1M (Jun 18)
mnvData = query_indels_mnvs(dbProd, "", "MNP")
View(mnvData)
nrow(mnvData)
write.csv(mnvData, "~/logs/r_output/mnv_prod_all_20180418.csv", row.names=F, quote=F)
#mnvData = read.csv("~/logs/r_output/mnv_prod_all.csv")
mnvData = read.csv("~/logs/r_output/mnv_prod_all_20180418.csv") # with 3-base PON fix
#mnvData = within(mnvData, rm(X))
nrow(mnvData)

# limit to samples in the high-confidence set
nrow(hpcSamples)
mnvData = mnvData %>% filter(SampleId %in% hpcSamples$sampleId)
mnvData = mnvData %>% filter(SampleId %in% highestPurityCohortSummary$sampleId)

nrow(mnvData) # 864K at last count 796K within HPC & DR022
View(mnvData)



mnvData = create_mnv_fields(mnvData)
View(mnvData)

# create and store Buckets
mnvBuckets = mnvData %>% group_by(Bucket) %>% count() %>% arrange(Bucket)
nrow(mnvBuckets)
View(mnvBuckets)
write.csv(mnvBuckets, "~/data/r_data/mnv_buckets.csv", row.names=F, quote=F)

View(mnvData %>% group_by(Length) %>% count())
View(mnvData %>% group_by(Bucket) %>% count())

allMnvData = mnvData # last set to full cohort (HPC, 2405 samples)
mnvData = allMnvData
nrow(mnvData %>% filter(Clonality=="SUBCLONAL"))
mnvData = mnvData %>% filter(Clonality=="SUBCLONAL")
# mnvData = mnvMBData # use multiple biopsy samples only
nrow(mnvData)

mnvSampleCounts = mnvData %>% group_by(SampleId,Bucket) %>% summarise(Count=n()) %>% arrange(SampleId,Bucket)
# colnames(mnvSampleCounts) <- c("SampleId","Bucket","Count")
View(mnvSampleCounts)
n_distinct(mnvSampleCounts$SampleId)

# mnvBucketCounts = mnvData %>% group_by(Bucket) %>% summarise(Count=n()) %>% arrange(Bucket)
# View(mnvBucketCounts)


mnvMatrixData = mnvSampleCounts %>% spread(SampleId, Count)
mnvMatrixData[is.na(mnvMatrixData)] = 0
View(mnvMatrixData[,1:10]) # check a subset

mnvBucketNames = mnvMatrixData$Bucket
mnvBucketNames = data.frame(mnvBucketNames)
colnames(mnvBucketNames) <- c("Bucket")
View(mnvBucketNames)
mnvMatrixData = within(mnvMatrixData, rm(Bucket))

write.csv(mnvMatrixData, file="~/data/r_data/mnv_nmf_matrix_data_dr22.csv", row.names=F, quote=F)
# write.csv(mnvMatrixData, file="~/data/r_data/mnv_nmf_subc_matrix_data_dr22.csv", row.names=F, quote=F)

# run estimations
mnvNmfEstimate <- nmf(mnvMatrixData, rank=8:14, method="brunet", nrun=4, seed=123456, .opt='vp4')
save(mnvNmfEstimate, file="~/logs/r_output/mnvNmfEstimate.RData")
load("~/data/mnvNmfEstimate_6_10.RData")
plot(mnvNmfEstimate)

# generate the actual NMF results
mnvSigCount = 8
mnvNmfResult <- nmf(mnvMatrixData, rank=mnvSigCount, method="brunet", nrun=10, seed=123456, .opt='vp1')
# save(mnvNmfResult, file="~/logs/r_output/mnvNmfResult_sig8_dr22.RData")

rm(nmfResult)
rm(mnvNmfResult)
load("~/data/r_data/mnvNmfResult_sig8_dr22.RData")
# load("~/data/r_data/mnvNmfResult_sig8.RData")
View(nmfResult)
mnvNmfResult = nmfResult

mnvSigCount = 8
mnvSigNamesNum = get_signame_list(mnvSigCount, F)
mnvSigNamesStr = get_signame_list(mnvSigCount, T)

View(mnvBucketNames)

evaluate_nmf_run("MNV", "sig08_DR22", mnvSigCount, mnvNmfResult, mnvMatrixData, mnvSampleCounts,
                 sampleCancerTypes, mnvBucketNames, mnvSigNamesNum, mnvSigNamesStr, TRUE, FALSE, TRUE)


View(mnvSampleCounts %>% group_by(Bucket) %>% summarise(Count=sum(Count)))

# using own NMF
mnvSimSignatures = as.matrix(read.csv("~/dev/nmf/logs/mnv_nmf_sigs.csv"))
colnames(mnvSimSignatures) = mnvSigNamesNum
mnvSimContribs = as.matrix(read.csv("~/dev/nmf/logs/mnv_nmf_contribs.csv"))

evaluate_nmf_data("MNV", "nmf_sf_sig08_dr22", mnvSigCount, mnvSimSignatures, mnvSimContribs, mnvMatrixData, mnvSampleCounts,
                 sampleCancerTypes, mnvBucketNames, mnvSigNamesNum, mnvSigNamesStr, FALSE, FALSE, TRUE)



# Subclonal Evalaution
rm(nmfResult)
load("~/data/r_data/mnvNmfResult_subc_sig8_dr22.RData")
View(nmfResult)
mnvNmfResult = nmfResult

evaluate_nmf_run("MNV_Subclonal", "sig08_DR22", mnvSigCount, mnvNmfResult, mnvMatrixData, mnvSampleCounts,
                 sampleCancerTypes, mnvBucketNames, mnvSigNamesNum, mnvSigNamesStr, TRUE, FALSE, TRUE)


# using Clonal signatures
load("~/data/mnvNmfResult_sig8_dr22.RData")

mnvSignatures = NMF::basis(mnvNmfResult)
View(mnvSignatures)
mnvBucketCount = nrow(mnvBucketNames)
mnvContributions = apply_signatures(mnvMatrixData, mnvSignatures)

# sum(mnvSampleCounts$Count)

evaluate_nmf_data("MNV_Subclonal_CS", "sig08_DR22", mnvSigCount, mnvSignatures, mnvContributions, mnvMatrixData, mnvSampleCounts,
                 sampleCancerTypes, mnvBucketNames, mnvSigNamesNum, mnvSigNamesStr, TRUE, FALSE, TRUE)


mnvSampleSigData = get_sig_data(NMF::basis(mnvNmfResult), NMF::coef(mnvNmfResult), mnvSigNamesStr, colnames(NMF::coef(mnvNmfResult)))
mnvSampleSigData$Count = round(mnvSampleSigData$Count,0)
View(mnvSampleSigData)
write.csv(mnvSampleSigData, "~/logs/r_output/mnvSampleSigData.csv", row.names=F, quote=F)


# Bucket-Analyser sigs and contributions
mnvMatrixData = read.csv("~/data/r_data/mnv_nmf_matrix_data_dr22.csv", stringsAsFactors=F)

mnvBaContribs = as.matrix(read.csv(file="~/dev/nmf/logs/mnv_ba_contribs.csv", stringsAsFactors=F))
View(mnvBaContribs[,50:70])
nrow(mnvBaContribs)
mnvBaSigs = as.matrix(read.csv(file="~/dev/nmf/logs/mnv_ba_sigs.csv", stringsAsFactors=F))
View(mnvBaSigs)
mnvBaSigNames = colnames(mnvBaSigs)
print(mnvBaSigNames)
ncol(mnvBaContribs)
mnvBaSigCount = ncol(mnvBaSigs)
print(mnvBaSigCount)

svnBaSigNumList = get_signame_list(mnvBaSigCount, F)
mnvBaSigStrList = get_signame_list(mnvBaSigCount, T)
colnames(mnvBaSigs) <- svnBaSigNumList
View(mnvBaSigs)

# trim sig names
for(i in 1:mnvBaSigCount)
{
  mnvBaSigNames[i] = paste(mnvBaSigStrList[i], mnvBaSigNames[i], sep="_")
  sigLen = stringi::stri_length(mnvBaSigNames[i])

  if(sigLen > 8)
  {
    catIndex = stri_locate_first_fixed(mnvBaSigNames[i], "_cat")
    if(!is.na(catIndex[1]))
    {
      mnvBaSigNames[i] = substring(mnvBaSigNames[i], 1, catIndex[1]-1)
      # print(paste("sig name shorted: ", mnvBaSigNames[i], sep=''))
    }

    mnvBaSigNames[i] = stri_replace_all_fixed(mnvBaSigNames[i], ".", "")
  }
}

print(mnvBaSigNames)

bgSigCount = 20

evaluate_nmf_data("MNV", "ba_denovo_all_CT", mnvBaSigs, mnvBaContribs, mnvMatrixData, mnvSampleCounts,
                  sampleCancerTypes2, mnvBucketNames, mnvBaSigNames, T, F, bgSigCount, F)

mnvSampleTotals = mnvSampleCounts %>% group_by(SampleId) %>% summarise(SampleCount=sum(Count)) %>% arrange(-SampleCount)
View(mnvSampleTotals)

# DEBUG
contribution = mnvBaContribs
signatures = mnvBaSigs
summaryCounts = mnvSampleCounts
bucketNames = mnvBucketNames
sigNamesNamed = mnvBaSigNames
matrixData = mnvMatrixData

sigCount = nrow(contribution)

sampleNames = colnames(contribution)

origSampleCounts = summaryCounts %>% group_by(SampleId) %>% summarise(OrigSampleCount=sum(Count))

sigNamesUnamed = get_signame_list(sigCount, F)
colnames(signatures) = sigNamesUnamed

sigNamesCombined = cbind(sigNamesUnamed, sigNamesNamed)
colnames(sigNamesCombined) <- c("Signature", "SigName")

bucketIndex = data.frame(as.numeric(as.character(rownames(bucketNames))))
colnames(bucketIndex) <- c("BucketIndex")
bucketNamesIndexed = cbind(bucketNames, bucketIndex)
bucketNamesIndexed$BucketIndex = bucketNamesIndexed$BucketIndex-1

sigBucketData = get_bucket_data(signatures, contribution, bucketNames)
sigBucketData = merge(sigBucketData,bucketNamesIndexed,by="Bucket",all.x=T)
sigBucketData = merge(sigBucketData,sigNamesCombined,by="Signature",all.x=T)

sigBucketStats = get_sig_bucket_stats(sigBucketData)
sigBucketTopN = get_top_buckets(sigBucketData)
View(sigBucketTopN)
bucketSummaryData = get_bucket_stats(sigBucketData)

plot_bucket_summary_data(bucketSummaryData, sigBucketTopN, "title", 40)





# INDEL Handling

# download from prod - around 12M (Jun 18)
# indelData = query_indels_mnvs(dbProd, "", "INDEL")
# write.csv(indelData, "~/logs/r_output/indel_prod_all.csv", row.names=F, quote=F)
# indelData = read.csv("~/logs/r_output/indel_prod_all.csv")
indelData = read.csv("~/data/indel_prod_all_v3.csv")
nrow(indelData)

indelData = indelData %>% filter(SampleId %in% hpcSamples$sampleId)
nrow(indelData)

View(head(indelData,100))

getFirstAlt<-function(alt, patternChar)
{
  if(grepl(patternChar, alt))
  {
    commaIndex = stri_locate(pattern=patternChar, alt, fixed = TRUE)[1,1]
    altFirst = substr(alt,1,commaIndex-1)
    return (altFirst)
  }
  else
  {
    return (alt)
  }
}

calcSplitAltLength<-function(alt, patternChar)
{
  if(grepl(',', alt))
  {
    commaIndex = stri_locate(pattern=patternChar, alt, fixed = TRUE)[1,1]
    altFirst = substr(alt,1,commaIndex-1)
    altLen = stri_length(altFirst)
  }
  else
  {
    altLen = stri_length(alt)
  }
  return (altLen)
}

create_indel_fields<-function(indelData)
{
  indelData$SubType = ifelse(indelData$AltLength > stringi::stri_length(indelData$Ref),'INS','DEL')

  indelData$Length = ifelse(indelData$SubType=='INS',indelData$AltLength-stringi::stri_length(indelData$Ref),
                            stringi::stri_length(indelData$Ref)-indelData$AltLength)

  indelData$LengthGrp = ifelse(indelData$Length<=3,indelData$Length,4)

  indelData$RepeatCountHigh = ifelse(indelData$RepeatCount>=4,T,F)
  indelData$RepeatCountType = ifelse(indelData$RepeatCountHigh,ifelse(indelData$RepeatCount>8,'High','Med'),'Low')

  indelData$MH = ifelse(indelData$SubType=='DEL'&!is.na(indelData$Microhomology)&stringi::stri_length(indelData$Microhomology)>0,T,F)

  return (indelData)
}

get_base_ins_or_del<-function(ref,alt,length)
{
  if(length != 1)
    return ('')

  longStr = ifelse(stri_length(ref) > stri_length(alt), ref, alt)
  shortStr = ifelse(stri_length(ref) < stri_length(alt), ref, alt)
  shortLen = stri_length(shortStr)

  curLetter = ''
  for(i in 1:stri_length(longStr))
  {
    curLetter = stri_sub(longStr, i, i)
    if(i > shortLen)
      break

    if(curLetter != stri_sub(shortStr, i, i))
      break
  }

  if(curLetter == 'G' || curLetter == 'C')
    return ('C')
  else
    return ('T')
}

# make a combined bucket name: DEL vs INS, LengthGroup, RepeatCount high and type, C or T for 1-base ins and dels, and for DELs whether has MH or not
create_indel_bucket<-function(indels)
{
  indels$Bucket = paste(indels$SubType, ifelse(indels$LengthGrp>3,3,indels$LengthGrp), sep="_")
  indels$Bucket = ifelse(indels$RepeatCountHigh, paste(indels$Bucket, "REP", sep="_"), indels$Bucket)
  indels$Bucket = ifelse(indels$LengthGrp==1&(indels$SubType=="DEL"|indels$RepeatCountHigh), paste(indels$Bucket, indels$BaseInsDel, sep="_"), indels$Bucket)
  indels$Bucket = ifelse(!indels$RepeatCountHigh&indels$MH, paste(indels$Bucket, "MH", sep="_"), indels$Bucket)
  indels$Bucket = ifelse(indels$LengthGrp==1&indels$RepeatCountHigh, paste(indels$Bucket, indels$RepeatCountType, sep="_"), indels$Bucket)

  return (indels)
}

create_indel_bucket_old<-function(indels)
{
  indels$Bucket = paste(indels$SubType, indels$LengthGrp, sep="_")
  indels$Bucket = ifelse(indels$RepeatCountHigh, paste(indels$Bucket, "REP", sep="_"), indels$Bucket)
  indels$Bucket = ifelse(!indels$RepeatCountHigh&indels$MH, paste(indels$Bucket, "MH", sep="_"), indels$Bucket)

  return (indels)
}

# write.csv(indelData3, "~/data/indel_prod_all_v3.csv", row.names=F, quote=F)

# View(indelData %>% group_by(Clonality) %>% count())

#allIndelData = indelData
#indelData = allIndelData
#indelData = indelData %>% filter(Clonality=="SUBCLONAL")
# nrow(indelData %>% filter(Clonality=="SUBCLONAL"))
nrow(indelData)
View(head(indelData,1000))
nrow(indelData %>% group_by(SampleId) %>% count())


# note: allIndelData already has the additional bucket fields
indelData$AltLength = apply(indelData[,c('Alt')], 1, function(x) calcSplitAltLength(x[1], ';'))
indelData$AltFirst = apply(indelData[,c('Alt'),drop=F], 1, function(x) getFirstAlt(x[1], ';'))
indelData$BaseInsDel = apply(indelData[,c('Ref','AltFirst'),drop=F], 1, function(x) getBaseInsertedDeleted(x[1], x[2]))


# print(getBaseInsertedDeleted("ACGT","ACG"))

indelData = create_indel_fields(indelData)

View(head(indelData,1000))

indelSampleCounts = indelData %>% group_by(SampleId,SubType,LengthGrp,MH,RepeatCountHigh,RepeatCountType,BaseInsDel) %>% summarise(Count=n())
indelSampleCountsOld = indelData %>% group_by(SampleId,SubType,LengthGrp,MH,RepeatCountHigh) %>% summarise(Count=n())

indelSampleCountsOld = create_indel_bucket(indelSampleCountsOld)
write.csv(indelSampleCounts, "~/data/r_data/indel_sample_counts.csv", row.names=F, quote=F)

# create and store Buckets
indelBuckets = indelSampleCounts %>% group_by(Bucket) %>% summarise(VarCount=sum(Count)) %>% arrange(Bucket)
indelBucketsNew = indelBuckets
sum(indelBuckets$VarCount)

indelSampleCountsOld = create_indel_bucket_old(indelSampleCounts)
View(indelSampleCountsOld)
indelBucketsOld = indelSampleCountsOld %>% group_by(Bucket) %>% summarise(VarCount=sum(Count)) %>% arrange(Bucket)
sum(indelBucketsOld$VarCount)
View(indelBucketsOld)
View(indelBucketsNew)
# write.csv(indelBuckets, "~/data/r_data/indel_buckets.csv", row.names=F, quote=F)

# filter for DR22 samples only
indelSampleCounts = indelSampleCounts %>% filter(SampleId %in% hpcSamples$sampleId)

# convert to sample and bucket counts for signature analysis
indelSampleBucketCounts = indelSampleCounts %>% group_by(SampleId,Bucket) %>% summarise(Count=sum(Count))
View(indelSampleBucketCounts)

indelMatrixData = indelSampleBucketCounts %>% spread(SampleId, Count) # %>% select(SampleId,BucketName,Count)
indelMatrixData[is.na(indelMatrixData)] = 0
View(indelMatrixData[,1:10]) # check a subset
nrow(indelMatrixData)
ncol(indelMatrixData)

indelBucketNames = data.frame(indelMatrixData$Bucket)
colnames(indelBucketNames) <- c("Bucket")
View(indelBucketNames)

# ordering is critical for when signatures and then applying to new samples
indelMatrixData = indelMatrixData %>% arrange(Bucket)
indelMatrixData = within(indelMatrixData, rm(Bucket))

write.csv(indelMatrixData, file="~/data/r_data/indel_nmf_matrix_data.csv", row.names=F, quote=F)
write.csv(indelMatrixData, file="~/data/r_data/indel_matrix_data_dr22_23bkts.csv", row.names=F, quote=F)
#write.csv(indelMatrixData, file="~/logs/r_output/indel_nmf_subc_counts_dr22.csv", row.names=F, quote=F)

# indelSampleCounts = read.csv("~/data/r_data/indel_sample_counts.csv")
# nrow(indelSampleCounts)
# ncol(indelSampleCounts)
#
# indelMatrixData = read.csv("~/data/r_data/indel_nmf_counts.csv")
# nrow(indelMatrixData)
# ncol(indelMatrixData)

# using bucket analyser

indBaContribs = as.matrix(read.csv(file="~/dev/nmf/logs/indel_ba_contribs.csv", stringsAsFactors=F))
View(indBaContribs[,50:70])
nrow(indBaContribs)
indBaSigs = as.matrix(read.csv(file="~/dev/nmf/logs/indel_ba_sigs.csv", stringsAsFactors=F))
View(indBaSigs)
indBaSigNames = colnames(indBaSigs)
print(indBaSigNames)
ncol(indBaContribs)
indBaSigCount = ncol(indBaSigs)
print(indBaSigCount)

indBaSigNumList = get_signame_list(indBaSigCount, F)
indBaSigStrList = get_signame_list(indBaSigCount, T)
colnames(indBaSigs) <- indBaSigNumList
View(indBaSigs)

# trim sig names
for(i in 1:indBaSigCount)
{
  indBaSigNames[i] = paste(indBaSigStrList[i], indBaSigNames[i], sep="_")
  sigLen = stringi::stri_length(indBaSigNames[i])

  if(sigLen > 8)
  {
    catIndex = stri_locate_first_fixed(indBaSigNames[i], "_cat")
    if(!is.na(catIndex[1]))
    {
      indBaSigNames[i] = substring(indBaSigNames[i], 1, catIndex[1]-1)
      # print(paste("sig name shorted: ", indBaSigNames[i], sep=''))
    }

    indBaSigNames[i] = stri_replace_all_fixed(indBaSigNames[i], ".", "")
  }
}

print(indBaSigNames)

bgSigCount = 20

evaluate_nmf_data("INDEL", "ba_denovo_all_CT_v2", indBaSigs, indBaContribs, indelMatrixData, indelSampleBucketCounts,
                  sampleCancerTypes2, indelBucketNames, indBaSigNames, T, F, bgSigCount, F)


origSampleCounts = indelSampleBucketCounts %>% group_by(SampleId) %>% summarise(SampleCount=sum(Count))
View(origSampleCounts)




# run estimations
indelNmfEstimate <- nmf(indelMatrixData, rank=8:17, method="brunet", nrun=4, seed=123456, .opt='vp4')
# save(indelNmfEstimate, file="~/logs/r_output/indelNmfEstimate.RData")
load(file="~/data/indelNmfEstimate.RData")
load(file="~/data/indelNmfEstimate_4_8.RData")
plot(indelNmfEstimate)

# generate the actual NMF results
indelSigCount = 5
indelNmfResult <- nmf(indelMatrixData, rank=indelSigCount, method="brunet", nrun=5, seed=123456, .opt='vp5')
save(indelNmfResult, file="~/logs/r_output/indelNmfResult_sig12.RData")
load(file="~/data/r_data/indelNmfResult_sig5.RData")
load(file="~/data/indelNmfResult_sig5_dr22.RData")
load(file="~/data/indelNmfResult_sig5_subc_dr22.RData")
View(indelNmfResult)


indelSigNamesNum = get_signame_list(indelSigCount, F)
indelSigNamesStr = get_signame_list(indelSigCount, T)
print(indelSigNamesNum)
print(indelSigNamesStr)
View(indelSampleCounts)

evaluate_nmf_run("INDEL", "sig5_DR22", indelSigCount, indelNmfResult, indelMatrixData, indelSampleBucketCounts,
                 sampleCancerTypes, indelBucketNames, indelSigNamesNum, indelSigNamesStr, TRUE, FALSE, TRUE)

 # using external fitting, yeah baby
indelExtSignatures = as.matrix(read.csv("~/logs/nmf_signatures.csv"))
indelExtContributions = read.csv("~/logs/nmf_contributions.csv")
colnames(indelExtSignatures) <- NULL
View(indelExtSignatures)
View(indelExtContributions[,1:10])
indelExtTotals = indelExtSignatures %*% as.matrix(indelExtContributions)
nrow(indelExtTotals)
ncol(indelExtTotals)
sum(indelExtTotals)
sum(indelMatrixData)
sum(indelMatrixData[,2])

indelSignatures = NMF::basis(indelNmfResult)
indelContributions = NMF::coef(indelNmfResult)

View(indelSignatures)
View(indelContributions[,1:10])

evaluate_nmf_data("INDEL", "sig5_ext2", indelSigCount, indelExtSignatures, indelExtContributions, indelMatrixData, indelSampleBucketCounts,
                 sampleCancerTypes, indelBucketNames, indelSigNamesNum, indelSigNamesStr, FALSE, FALSE, TRUE)



# evaluate_nmf_run("INDEL_Subclonal", "sig5_DR22", indelSigCount, indelNmfResult, indelMatrixData, indelSampleCounts,
#                  sampleCancerTypes, indelBucketNames, indelSigNamesNum, indelSigNamesStr, TRUE, FALSE)


# Subclonal evaluation using Clonal signatures

indelMatrixData = read.csv("~/logs/r_output/indel_nmf_subc_counts_dr22.csv")

indelSignatures = NMF::basis(indelNmfResult)
indelBucketCount = nrow(indelBucketNames)
indelContributions = apply_signatures(indelMatrixData, indelSignatures, indelBucketCount)

evaluate_nmf_data("INDEL_Subclonal_CS", "sig5_DR22", indelSigCount, indelSignatures, indelContributions, indelMatrixData, indelSampleBucketCounts,
                 sampleCancerTypes, indelBucketNames, indelSigNamesNum, indelSigNamesStr, TRUE, FALSE, TRUE)


svSignatures = NMF::basis(indelNmfResult)
indelContribution = NMF::coef(indelNmfResult)
svSampleNames = colnames(svContribution)

View(indelNmfResult)
sum(indelSampleBucketCounts$Count)

indelSampleSigData = get_sig_data(NMF::basis(indelNmfResult), NMF::coef(indelNmfResult), indelSigNamesStr, colnames(NMF::coef(indelNmfResult)))
indelSampleSigData$Count = round(indelSampleSigData$Count,0)
View(indelSampleSigData)
write.csv(indelSampleSigData, "~/logs/r_output/indelSampleSigData.csv", row.names=F, quote=F)




# indelData2 = read.csv("~/data/indel_prod_all_v2.csv")
# indelData2$HasSplit = ifelse(indelData2$Clonality!="CLONAL"&indelData2$Clonality!="SUBCLONAL"&indelData2$Clonality!="INCONSISTENT",1,0)
# indelSplitData = indelData2 %>% filter(HasSplit==1)
# indelNonSplitData = indelData2 %>% filter(!HasSplit==0)
#
# View(indelSplitData)
# indelSplitData$AltTest = paste(indelSplitData$Alt,indelSplitData$Microhomology, sep=';')
# indelSplitData$Alt = indelSplitData$AltTest
# indelSplitData$Microhomology = indelSplitData$RepeatCount
# indelSplitData$RepeatCount = indelSplitData$Clonality
# nrow(indelSplitData)
#
# indelSplitData2 = indelSplitData %>% filter(RepeatCount!=""&!grepl("A",RepeatCount)&!grepl("C",RepeatCount)&!grepl("G",RepeatCount)&!grepl("T",RepeatCount))
# View(indelSplitData %>% filter(RepeatCount==""|grepl("A",RepeatCount)))
# nrow(indelSplitData2)
# indelSplitData2 = within(indelSplitData2, rm(AltTest))
# View(indelSplitData2)




# New Bucket Analysis
origSampleCounts = indelSampleBucketCounts %>% group_by(SampleId) %>% summarise(SampleCount=sum(Count))

View(indelSampleBucketCounts)

sampleBucketData = get_sample_bucket_data(indelMatrixData, origSampleCounts, sampleCancerTypes, indelBucketNames)
sampleBucketData = sampleBucketData %>% arrange(-SampleCount,Bucket)
nrow(sampleBucketData)
View(sampleBucketData)

rowIndex = data.frame(as.numeric(as.character(rownames(sampleBucketData))))
colnames(rowIndex) <- c("RowIndex")
View(rowIndex)
sampleBucketData = cbind(sampleBucketData, rowIndex)

sampleBucketData$RowGroup = floor(sampleBucketData$RowIndex/groupSize)
View(sampleBucketData)

sampleBucketStats = (sampleBucketData %>% group_by(Bucket,RowGroup)
                    %>% summarise(Count=sum(Count))
                    %>% arrange(Bucket,RowGroup))

View(sampleBucketStats)

bucketStatsSpread = sampleBucketStats %>% spread(Bucket,Count)
View(bucketStatsSpread)

indelBucketData = data.frame(matrix(ncol = 5, nrow = 0))
colnames(indelBucketData) <- c("Group", "AvgSampleCount", "BucketCount", "BucketPercent")
sampleCount = n_distinct(indelCounts$SampleId)
samplesPerGroup = 100
bucketCount = nrow(indelBucketNames)
groupCount = round(sampleCount/samPerGroup,0)
groupSize = samplesPerGroup*bucketCount


rowStart = 1
for(i in 1:groupCount)
{
  rowStart = i*samPerGroup*bucketCount
  rowEnd = rowStart + samPerGroup*bucketCount - 1

  print(paste("start=", rowStart, ", end=", rowEnd, sep=''))
  bucketData = slice(sampleBucketData, rowStart:rowEnd)


}



