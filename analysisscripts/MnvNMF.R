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

# download from prod
mnvData = query_indels_mnvs(dbProd, "", "MNP")
View(mnvData)
nrow(mnvData)
write.csv(mnvData, "~/logs/r_output/mnv_prod_all.csv", row.names=F, quote=F)
mnvData = read.csv("~/logs/r_output/mnv_prod_all.csv")
mnvData = within(mnvData, rm(X))

# limit to samples in the high-confidence set
mnvData = mnvData %>% filter(SampleId %in% highestPurityCohortSummary$sampleId)
nrow(mnvData)


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
mnvData = within(mnvData, rm(NeedsConverting))
View(mnvData)

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

print(setMutationAdjusted("AT", "GA", 2))

mnvData$MutationAdj2 = apply(mnvData[,c('Ref','Alt')], 1, function(x) setMutationAdjusted(x[1],x[2]))

mnvData = within(mnvData, rm(MutationAdj2))
# mnvData$MutationAdjTest = ifelse(mnvData$Length==2,mnvData$MutationAdj, ifelse(mnvData$Length==3,"3Bases","4+Bases"))

print(n_distinct((mnvData %>% filter(Length==2))$Mutation))
print(n_distinct(mnvData$MutationAdj2))
print(n_distinct(mnvData$MutationAdj))
View(mnvData %>% group_by(MutationAdj2) %>% count())

View(mnvData %>% group_by(Length) %>% count())
View(mnvData %>% group_by(Mutation) %>% count())

mnvSampleCounts = mnvData %>% group_by(SampleId,MutationAdj) %>% count()
# mnvSampleCounts = mnvSampleCounts %>% filter(grepl("CPCT020101", SampleId)) # temp
View(mnvSampleCounts)
n_distinct(mnvSampleCounts$SampleId)

mnvMatrixData = mnvSampleCounts %>% spread(SampleId, n)
mnvMatrixData[is.na(mnvMatrixData)] = 0
View(mnvMatrixData[,1:10]) # check a subset

mnvBucketNames = mnvMatrixData$MutationAdj
View(mnvBucketNames)
mnvMatrixData = within(mnvMatrixData, rm(MutationAdj))

write.csv(mnvMatrixData, file="~/logs/r_output/mnv_nmf_counts.csv", row.names=F, quote=F)

# run estimations
mnvNmfEstimate <- nmf(mnvMatrixData, rank=8:14, method="brunet", nrun=4, seed=123456, .opt='vp4')
save(mnvNmfEstimate, file="~/logs/r_output/mnvNmfEstimate.RData")
plot(mnvNmfEstimate)

# generate the actual NMF results
mnvSigCount = 8
mnvNmfResult <- nmf(mnvMatrixData, rank=mnvSigCount, method="brunet", nrun=5, seed=123456, .opt='vp5')
save(mnvNmfResult, file="~/logs/r_output/mnvNmfResult_sig8.RData")

sigNamesUnamed = c("1", "2", "3", "4", "5", "6", "7", "8")
sigNamesNamed = c("01", "02", "03", "04", "05", "06", "07", "08")

View(cancerTypes)

evaluate_nmf_run("MNV", "sig08", mnvSigCount, mnvNmfResult, mnvSampleCounts, sampleCancerTypes, mnvBucketNames, sigNamesUnamed, sigNamesNamed, FALSE) 



# INDEL Handling

# download from prod
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

View(head(indelData,1000))

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
indelSigCount = 5
indelNmfResult <- nmf(indelMatrixData, rank=indelSigCount, method="brunet", nrun=5, seed=123456, .opt='vp5')
save(indelNmfResult, file="~/logs/r_output/indelNmfResult_sig12.RData")
load(file="~/data/indelNmfResult_sig5.RData")

sigNamesUnamed = c("1", "2", "3", "4", "5")
sigNamesNamed = c("01", "02", "03", "04", "05")
# sigNamesUnamed = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
# sigNamesNamed = c("01", "02", "03", "04", "05", "06", "07", "09", "10", "11", "12")

evaluate_nmf_run("INDEL", "sig5", indelSigCount, indelNmfResult, indelSampleCounts, 
                 sampleCancerTypes, indelBucketNames, sigNamesUnamed, sigNamesNamed, FALSE) 






