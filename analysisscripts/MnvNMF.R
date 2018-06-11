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

dbProd = dbConnect(MySQL(), user='hmf', password='HMFhmf@1', dbname='hmfpatients', groups = "RAnalysis")
# dbDisconnect(dbProd)

cancerTypes = highestPurityCohortSummary %>% group_by(cancerType) %>% count()
cancerTypes = setNames(cancerTypes, c("CancerType", "Count"))
sampleCancerTypes = highestPurityCohortSummary %>% select(sampleId, cancerType)
sampleCancerTypes = setNames(sampleCancerTypes, c("SampleId", "CancerType"))
View(sampleCancerTypes)
View(cancerTypes)


# MNV Handling

# download from prod around 1.1M (Jun 18)
mnvData = query_indels_mnvs(dbProd, "", "MNP")
View(mnvData)
nrow(mnvData)
write.csv(mnvData, "~/logs/r_output/mnv_prod_all.csv", row.names=F, quote=F)
#mnvData = read.csv("~/logs/r_output/mnv_prod_all.csv")
#mnvData = within(mnvData, rm(X))

# limit to samples in the high-confidence set
mnvData = mnvData %>% filter(SampleId %in% highestPurityCohortSummary$sampleId)
nrow(mnvData) # 864K at last count


# create Buckets
View(head(mnvData,10))

mnvData$Length = stringi::stri_length(mnvData$Alt)
mnvData$Mutation = ifelse(mnvData$Length==2,paste(mnvData$Ref, mnvData$Alt, sep='>'), ifelse(mnvData$Length==3,"3Bases","4+Bases"))

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

# print(setMutationAdjusted("AT", "GA"))

mnvData$MutationAdj = apply(mnvData[,c('Ref','Alt')], 1, function(x) setMutationAdjusted(x[1],x[2]))

print(n_distinct((mnvData %>% filter(Length==2))$Mutation))
print(n_distinct(mnvData$MutationAdj))

View(mnvData %>% group_by(Length) %>% count())
View(mnvData %>% group_by(Mutation) %>% count())

mnvSampleCounts = mnvData %>% group_by(SampleId,MutationAdj) %>% summarise(Count=n())
colnames(mnvSampleCounts) <- c("SampleId", "Bucket", "Count")
View(mnvSampleCounts)
n_distinct(mnvSampleCounts$SampleId)
sum(mnvSampleCounts$n)

mnvMatrixData = mnvSampleCounts %>% spread(SampleId, Count)
mnvMatrixData[is.na(mnvMatrixData)] = 0
View(mnvMatrixData[,1:10]) # check a subset

mnvBucketNames = mnvMatrixData$Bucket
View(mnvBucketNames)
mnvMatrixData = within(mnvMatrixData, rm(Bucket))

write.csv(mnvMatrixData, file="~/logs/r_output/mnv_nmf_counts.csv", row.names=F, quote=F)

# run estimations
mnvNmfEstimate <- nmf(mnvMatrixData, rank=8:14, method="brunet", nrun=4, seed=123456, .opt='vp4')
save(mnvNmfEstimate, file="~/logs/r_output/mnvNmfEstimate.RData")
load("~/data/mnvNmfEstimate_6_10.RData")
plot(mnvNmfEstimate)

# generate the actual NMF results
mnvSigCount = 8
mnvNmfResult <- nmf(mnvMatrixData, rank=mnvSigCount, method="brunet", nrun=5, seed=123456, .opt='vp5')
save(mnvNmfResult, file="~/logs/r_output/mnvNmfResult_sig8.RData")
# load("~/logs/r_output/mnvNmfResult_sig8.RData")
load("~/data/mnvNmfResult_sig8_buck80.RData")
View(mnvNmfResult)

sigNamesUnamed = c(1, 2, 3, 4, 5, 6, 7, 8)
sigNamesNamed = c("01", "02", "03", "04", "05", "06", "07", "08")

View(mnvBucketNames)
evaluate_nmf_run("MNV", "sig08", mnvSigCount, mnvNmfResult, mnvSampleCounts, sampleCancerTypes, mnvBucketNames,
                 sigNamesUnamed, sigNamesNamed, TRUE, FALSE)

View(mnvNmfResult)

# TMP testing
signatures = NMF::basis(mnvNmfResult)
contribution = NMF::coef(mnvNmfResult)
sampleNames = colnames(mnvNmfResult)
View(contribution[,1:2])
View(signatures)
View(sampleNames)

View(mnvSampleCounts)
bucketNames = mnvBucketNames
View(bucketNames)
View(mnvSampleCounts %>% select(Bucket))
mnvResiduals = calc_sample_residuals(contribution, signatures, mnvBucketNames, mnvSampleCounts)
View(mnvResiduals)

sum(mnvResiduals$Count)
sum(mnvResiduals$ResidualTotal)
mnvResiduals$ResidualPerc = round(mnvResiduals$ResidualTotal/mnvResiduals$Count,2)
print(residuals(mnvNmfResult))

sigCount = ncol(signatures)
contribTrans = t(contribution) %>% as.data.frame()
sampleNames = colnames(contribution)
sampleSigContribs = cbind(sampleNames, contribTrans)
rownames(sampleSigContribs) <- NULL
View(sampleSigContribs)

# colnames(sampleSigContribs) <- c("SampleId", "SS_1", "SS_2", "SS_3", "SS_4", "SS_5", "SS_6", "SS_7", "SS_8")

# merge with bucket-sig data
sigContribs = cbind(signatures, bucketNames)
colnames(sigContribs) <- sigNamesTmp
View(sigContribs)
# colnames(sigContribs2) <- c("SB_1", "SB_2", "SB_3", "SB_4", "SB_5", "SB_6", "SB_7", "SB_8", "Bucket")

sigNames = c()
for(i in 1:sigCount)
{
  sigNamesTmp[i] = i
}
sigNamesTmp[sigCount+1] = "Bucket"
View(sigNamesTmp)

# sampleSigContribs2 = merge(sampleSigContribs, signatures, by.x="SampleId",by.y="SampleId",all.x=TRUE)

# initially don't merge on any common fields
samSigContribs2 = merge(sampleSigContribs, sigContribs, all.x=TRUE)
View(samSigContribs2)

names(samSigContribs2)[names(samSigContribs2) == 'sampleNames'] <- 'SampleId'

# colnames(mnvSampleCounts) <- c("SampleId", "Bucket", "Count")
sampleSigContribs3 = merge(samSigContribs2, mnvSampleCounts, by.x=c("SampleId","Bucket"),by.y=c("SampleId","Bucket"),all.x=TRUE)
sampleSigContribs3[is.na(sampleSigContribs3)] <- 0

View(mnvSampleCounts %>% filter(SampleId=="CPCT02010003T"))
View(mnvSampleCounts)
View(sampleSigContribs3)

sampleSigContribs3[, c(3:18)] <- sapply(sampleSigContribs3[, c(3:18)], as.character)
sampleSigContribs3[, c(3:18)] <- sapply(sampleSigContribs3[, c(3:18)], as.numeric)

sampleSigContribs3$SigAlloc = 0
for(i in 1:mnvSigCount)
{
  sampleSigContribs3$SigAlloc = sampleSigContribs3$SigAlloc + apply(sampleSigContribs3[,c(i+2,i+2+mnvSigCount)], 1, function(x) x[1]*x[2])
}

View(mnvNmfResult)

tmpSample = sampleSigContribs3 %>% filter(SampleId=="CPCT02010003T")
View(tmpSample)

tmpSample$SigAlloc = 0
for(i in 1:mnvSigCount)
{
  sigAlloc = apply(tmpSample[,c(i+2,i+2+mnvSigCount)], 1, function(x) x[1]*x[2])
  tmpSample$SigAlloc = tmpSample$SigAlloc + sigAlloc
}

View(tmpSample)
sum(tmpSample$ResidualDiff)
sum(tmpSample$SigAlloc)
sum(tmpSample$Count)

tmpSample$SigAlloc2 = (tmpSample$V1*tmpSample$`1`
                     + tmpSample$V2*tmpSample$`2`
                     + tmpSample$V3*tmpSample$`3`
                     + tmpSample$V4*tmpSample$`4`
                     + tmpSample$V5*tmpSample$`5`
                     + tmpSample$V6*tmpSample$`6`
                     + tmpSample$V7*tmpSample$`7`
                     + tmpSample$V8*tmpSample$`8`)

tmpSample$ResidualDiff = abs(tmpSample$Count-tmpSample$SigAlloc)

View(tmpSample)
sum(tmpSample$ResidualDiff)
sum(tmpSample$SigAlloc)
sum(tmpSample$Count)

View(contribution[,1])
View(signatures)

sampleSigContribs3$ResidualDiff = abs(sampleSigContribs3$Count-sampleSigContribs3$SigAlloc)
View(sampleSigContribs3)

sampleResidualData = (sampleSigContribs3 %>% group_by(sampleNames)
                      %>% summarise(Count=sum(Count),
                                    ResidualTotal=round(sum(ResidualDiff),1)))
View(sampleResidualData)



mnvSampleBucketTopN = get_top_buckets_by_sample(mnvSampleCounts, origSampleCounts, sampleCancerTypes)

plot_sample_bucket_contrib(mnvSampleBucketTopN, "Breast", 10, "MNV", 4)



# INDEL Handling

# download from prod - around 12M (Jun 18)
indelData = query_indels_mnvs(dbProd, "", "INDEL")
nrow(indelData)
write.csv(indelData, "~/logs/r_output/indel_prod_all.csv", row.names=F, quote=F)
indelData = read.csv("~/logs/r_output/indel_prod_all.csv")
indelData = within(indelData, rm(X))

indelData = indelData %>% filter(SampleId %in% highestPurityCohortSummary$sampleId)
nrow(indelData)

# create buckets
View(head(indelData,10))

# if an alt has different insert or del strings, split these and take the first (arbitrarily)
indelData$HasSplit = grepl(",",indelData$Alt)
nrow(indelData %>% filter(HasSplit))
indelData$AltAdjusted = indelData$Alt
indelData = within(indelData, rm(AltAdjusted))

indelData$AltLength = stri_length(indelData$Alt)
indelSplitData = indelData %>% filter(HasSplit)
indelNonSplitData = indelData %>% filter(!HasSplit)

calcSplitAltLength<-function(ref, alt)
{
  commaIndex = stri_locate(pattern=',', alt, fixed = TRUE)[1,1]
  altFirst = substr(alt,1,commaIndex-1)
  altLen = stri_length(altFirst)
  # print(paste("alt: ", alt, ", altFirst: ", altFirst, ", commaIndex: ", commaIndex, ", altLen: ", altLen, sep=''))
  return (altLen)
}

indelSplitData$AltLength = apply(indelSplitData[,c('Ref', 'Alt')], 1, function(x) calcSplitAltLength(x[1], x[2]))

View(indelSplitData)

# join back together
indelData = rbind(indelSplitData, indelNonSplitData)

indelData$SubType = ifelse(indelData$AltLength > stringi::stri_length(indelData$Ref),'INS','DEL')

indelData$Length = ifelse(indelData$SubType=='INS',indelData$AltLength-stringi::stri_length(indelData$Ref),
                          stringi::stri_length(indelData$Ref)-indelData$AltLength)

indelData$LengthGrp = ifelse(indelData$Length<=3,indelData$Length,4)

indelData$RepeatCountHigh = ifelse(indelData$RepeatCount>=4,T,F)
indelData$MH = ifelse(indelData$SubType=='DEL'&indelData$RepeatCountHigh&!is.na(indelData$Microhomology)&stringi::stri_length(indelData$Microhomology)>0,T,F)

View(indelData)
View(indelData %>% filter(Length==0))

indelSampleCounts = indelData %>% group_by(SampleId,SubType,LengthGrp,MH,RepeatCountHigh) %>% count()
View(indelSampleCounts)

indelMatrixData = indelSampleCounts %>% spread(SampleId, n)
indelMatrixData[is.na(indelMatrixData)] = 0
View(indelMatrixData[,1:10]) # check a subset

# make a combined bucket name
indelMatrixData$BucketName = paste(indelMatrixData$SubType, indelMatrixData$LengthGrp,
                                   ifelse(indelMatrixData$MH,"MH","NoMH"), ifelse(indelMatrixData$RepeatCountHigh,"RCH", "RCL"), sep="_")

indelBucketNames = indelMatrixData$BucketName
View(indelBucketNames)
indelMatrixData = within(indelMatrixData, rm(BucketName))
indelMatrixData = within(indelMatrixData, rm(SubType))
indelMatrixData = within(indelMatrixData, rm(MH))
indelMatrixData = within(indelMatrixData, rm(RepeatCountHigh))
indelMatrixData = within(indelMatrixData, rm(LengthGrp))

write.csv(indelMatrixData, file="~/logs/r_output/indel_nmf_counts.csv", row.names=F, quote=F)


# run estimations
indelNmfEstimate <- nmf(indelMatrixData, rank=8:17, method="brunet", nrun=4, seed=123456, .opt='vp4')
# save(indelNmfEstimate, file="~/logs/r_output/indelNmfEstimate.RData")
load(file="~/data/indelNmfEstimate.RData")
load(file="~/data/indelNmfEstimate_4_8.RData")
plot(indelNmfEstimate)

# generate the actual NMF results
indelSigCount = 12
indelNmfResult <- nmf(indelMatrixData, rank=indelSigCount, method="brunet", nrun=5, seed=123456, .opt='vp5')
save(indelNmfResult, file="~/logs/r_output/indelNmfResult_sig12.RData")
load(file="~/data/indelNmfResult_sig5.RData")

sigNamesUnamed = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
sigNamesNamed = c("01", "02", "03", "04", "05", "06", "07", "09", "10", "11", "12")

evaluate_nmf_run("INDEL", "sig12", indelSigCount, indelNmfResult, indelSampleCounts,
                 sampleCancerTypes, indelBucketNames, sigNamesUnamed, sigNamesNamed, TRUE, FALSE)






