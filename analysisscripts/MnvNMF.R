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

# limit to samples in the high-confidence set
nrow(hpcSamples)
mnvData = mnvData %>% filter(SampleId %in% hpcSamples$sampleId)
mnvData = mnvData %>% filter(SampleId %in% highestPurityCohortSummary$sampleId)

nrow(mnvData) # 864K at last count 796K within HPC & DR022


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

# write.csv(mnvMatrixData, file="~/logs/r_output/mnv_nmf_subc_counts_dr22.csv", row.names=F, quote=F)
write.csv(mnvMatrixData, file="~/logs/r_output/mnv_nmf_matrix_data.csv", row.names=F, quote=F)

# run estimations
mnvNmfEstimate <- nmf(mnvMatrixData, rank=8:14, method="brunet", nrun=4, seed=123456, .opt='vp4')
save(mnvNmfEstimate, file="~/logs/r_output/mnvNmfEstimate.RData")
load("~/data/mnvNmfEstimate_6_10.RData")
plot(mnvNmfEstimate)

# generate the actual NMF results
mnvSigCount = 8
mnvNmfResult <- nmf(mnvMatrixData, rank=mnvSigCount, method="brunet", nrun=10, seed=123456, .opt='vp1')
save(mnvNmfResult, file="~/logs/r_output/mnvNmfResult_sig8_dr22.RData")
# load("~/data/mnvNmfResult_sig8_dr22.RData")
load("~/data/r_data/mnvNmfResult_sig8.RData")
View(mnvNmfResult)

mnvSigNamesNum = get_signame_list(8, F)
mnvSigNamesStr = get_signame_list(8, T)

View(mnvBucketNames)

evaluate_nmf_run("MNV", "sig08", mnvSigCount, mnvNmfResult, mnvMatrixData, mnvSampleCounts,
                 sampleCancerTypes, mnvBucketNames, mnvSigNamesNum, mnvSigNamesStr, TRUE, FALSE, TRUE)


View(mnvSampleCounts %>% group_by(Bucket) %>% summarise(Count=sum(Count)))

mnvSignatures = NMF::basis(mnvNmfResult)
mnvContributions = NMF::coef(mnvNmfResult)
View(mnvSignatures)
View(mnvContributions[,1:10])
nrow(mnvContributions)
ncol(mnvContributions)
nrow(mnvSignatures)
ncol(mnvSignatures)
View(mnvBucketNames)
View(mnvMatrixData[,1:10])
ncol(mnvMatrixData)

View(mnvSampleCounts)

mnvResiduals2 = calc_sample_residuals_v2(mnvContributions, mnvSignatures, mnvMatrixData, mnvBucketNames)
View(mnvResiduals2)

sum(mnvResiduals2$Count)
sum(mnvResiduals2$ResidualTotal)

sampleSigData = get_sig_data(mnvSignatures, mnvContributions, mnvSigNamesStr, colnames(mnvContributions))
sampleSigData = append_residuals(mnvContributions, mnvSignatures, mnvBucketNames, mnvSampleCounts, sampleSigData)




# Subclonal Evalaution
load("~/data/mnvNmfResult_sig8_subc_dr22.RData")

evaluate_nmf_run("MNV_Subclonal", "sig08_DR22", mnvSigCount, mnvNmfResult, mnvMatrixData, mnvSampleCounts,
                 sampleCancerTypes, mnvBucketNames, mnvSigNamesNum, mnvSigNamesStr, FALSE, FALSE, TRUE)


# using Clonal signatures
load("~/data/mnvNmfResult_sig8_dr22.RData")

mnvSignatures = NMF::basis(mnvNmfResult)
View(mnvSignatures)
mnvBucketCount = nrow(mnvBucketNames)
mnvContributions = apply_signatures(mnvMatrixData, mnvSignatures, mnvBucketCount)
# View(mnvMatrixData[,1:10])
# nrow(mnvMatrixData)
# ncol(mnvMatrixData)
# View(mnvContributions[,1:10])

# sum(mnvSampleCounts$Count)

evaluate_nmf_data("MNV_Subclonal_CS", "sig08_DR22", mnvSigCount, mnvSignatures, mnvContributions, mnvMatrixData, mnvSampleCounts,
                 sampleCancerTypes, mnvBucketNames, mnvSigNamesNum, mnvSigNamesStr, TRUE, FALSE, TRUE)


mnvSampleSigData = get_sig_data(NMF::basis(mnvNmfResult), NMF::coef(mnvNmfResult), mnvSigNamesStr, colnames(NMF::coef(mnvNmfResult)))
mnvSampleSigData$Count = round(mnvSampleSigData$Count,0)
View(mnvSampleSigData)
write.csv(mnvSampleSigData, "~/logs/r_output/mnvSampleSigData.csv", row.names=F, quote=F)





# View(mnvSampleCounts)
# bucketNames = mnvBucketNames
# View(bucketNames)
# View(mnvSampleCounts %>% select(Bucket))
# mnvResiduals = calc_sample_residuals(contribution, signatures, mnvBucketNames, mnvSampleCounts)
# View(mnvResiduals)

sum(mnvResiduals$Count)
sum(mnvResiduals$ResidualTotal)
mnvResiduals$ResidualPerc = round(mnvResiduals$ResidualTotal/mnvResiduals$Count,2)
print(residuals(mnvNmfResult))



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

  indelData$MH = ifelse(indelData$SubType=='DEL'&!is.na(indelData$Microhomology)&stringi::stri_length(indelData$Microhomology)>0,T,F)

  return (indelData)
}

# make a combined bucket name: DEL vs INS, LengthGroup, RepeatCount high and for DELs whether has MH or not
create_indel_bucket<-function(indelSampleCounts)
{
  indelSampleCounts$Bucket = paste(indelSampleCounts$SubType, indelSampleCounts$LengthGrp, sep="_")
  indelSampleCounts$Bucket = ifelse(indelSampleCounts$RepeatCountHigh, paste(indelSampleCounts$Bucket, "REP", sep="_"), indelSampleCounts$Bucket)
  indelSampleCounts$Bucket = ifelse(!indelSampleCounts$RepeatCountHigh&indelSampleCounts$MH, paste(indelSampleCounts$Bucket, "MH", sep="_"), indelSampleCounts$Bucket)

  return (indelSampleCounts)
}

# write.csv(indelData3, "~/data/indel_prod_all_v3.csv", row.names=F, quote=F)

nrow(indelData)

View(indelData %>% group_by(Clonality) %>% count())

#allIndelData = indelData
indelData = allIndelData
indelData = indelData %>% filter(Clonality=="SUBCLONAL")
# nrow(indelData %>% filter(Clonality=="SUBCLONAL"))
nrow(indelData)
View(head(indelData,1000))
nrow(indelData %>% group_by(SampleId) %>% count())


# note: allIndelData already has the additional bucket fields
indelData$AltLength = apply(indelData[,c('Alt')], 1, function(x) calcSplitAltLength(x[1], ';'))
indelData = create_indel_fields(indelData)

# View(head(indelData,1000))

indelSampleCounts = indelData %>% group_by(SampleId,SubType,LengthGrp,MH,RepeatCountHigh) %>% summarise(Count=n())

indelSampleCounts = create_indel_bucket(indelSampleCounts)
View(indelSampleCounts)
write.csv(indelSampleCounts, "~/data/r_data/indel_sample_counts.csv", row.names=F, quote=F)


# create and store Buckets
indelBuckets = indelSampleCounts %>% group_by(Bucket) %>% summarise(VarCount=sum(Count)) %>% arrange(Bucket)
View(indelBuckets)
write.csv(indelBuckets, "~/data/r_data/indel_buckets.csv", row.names=F, quote=F)


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
write.csv(indelMatrixData, file="~/logs/r_output/indel_nmf_counts_dr22.csv", row.names=F, quote=F)
#write.csv(indelMatrixData, file="~/logs/r_output/indel_nmf_subc_counts_dr22.csv", row.names=F, quote=F)


indelMatrixData = read.csv("~/data/r_data/indel_nmf_counts.csv")


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



