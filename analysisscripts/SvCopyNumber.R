library(purple);
library(RMySQL)
library(data.table)
library(IRanges)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(stringi)
library(devtools)


# to get arm lengths
detach("package:svnmf", unload=TRUE);
library(svnmf)

get_copy_number_range<-function(cnData, sample, chr, posStart = 0, posEnd = 2e8)
{
  cnSection = cnData %>% filter(SampleId==sample & Chr==chr & PosEnd < posEnd & PosStart > posStart)
  return (cnSection)
}

get_copy_number_plot<-function(cnData, sampleId, chr, posStart = 0, posEnd = 2e8) {

  cnSection = get_copy_number_range(cnData, sampleId, chr, posStart, posEnd)
  plot = (ggplot(data = cnSection, aes(x=PosStart, y=CopyNumber))
          + geom_line())

  return (plot)
}

plot_copy_number_section<-function(cnSection) {

  plot = (ggplot(data = cnSection, aes(x=start, y=copyNumber))
          + geom_line())

  print(plot)
}


# DB methods for accessing CN data
getSampleIdsStr<-function(samples)
{
  sampleIdsStr = ""

  for(i in 1:nrow(samples))
  {
    sample <- samples[i,]

    if(i > 1)
      sampleIdStr = paste(",'", sample$SampleId, "'", sep="")
    else
      sampleIdStr = paste("'", sample$SampleId, "'", sep="")

    # sampleIdStr = stringi::stri_replace_all_fixed(sampleIdStr, " ", "")
    sampleIdsStr = paste(sampleIdsStr, sampleIdStr)
  }

  return (sampleIdsStr)
}

getCopyNumber<-function(dbConnect,sampleId,chromosome)
{
  sampleIdStr = paste("SampleId='", sampleId, "'")
  sampleIdStr = stringi::stri_replace_all_fixed(sampleIdStr, " ", "")
  sql = paste("select SampleId, SegmentStartSupport, SegmentEndSupport, CopyNumber from copyNumber where ", sampleIdStr, " and chromosome=",chromosome)
  sql = paste(sql, " and (segmentStartSupport = 'CENTROMERE' or segmentStartSupport = 'TELOMERE' or segmentEndSupport = 'CENTROMERE' or segmentEndSupport = 'TELOMERE')")
  # View(sql)
  return ((dbGetQuery(dbConnect, sql)))
}

getCopyNumbers<-function(dbConnect,sampleIds)
{
  sql = paste("select SampleId, Chromosome, SegmentStartSupport, SegmentEndSupport, CopyNumber from copyNumber where SampleId in(", sampleIds, ")")
  sql = paste(sql, " and (segmentStartSupport = 'CENTROMERE' or segmentStartSupport = 'TELOMERE' or segmentEndSupport = 'CENTROMERE' or segmentEndSupport = 'TELOMERE')")
  return ((dbGetQuery(dbConnect, sql)))
}

annotateWithCopyNumber<-function(dbConnect,samples)
{
  samples$CopyNumberStart = 2
  samples$CopyNumberEnd = 2

  sampleIdsStr = getSampleIdsStr(samples)

  cnResults = getCopyNumbers(dbConnect, sampleIdsStr)

  for(i in 1:nrow(samples))
  {
    sample <- samples[i,]

    cnSampleData = cnResults %>% filter(SampleId==sample$SampleId,Chromosome==sample$Chr)

    cnStart = head(cnSampleData %>% filter(SegmentStartSupport=='TELOMERE'),1)

    if(nrow(cnStart) == 1)
    {
      samples[i,]$CopyNumberStart = cnStart$CopyNumber
    }

    cnEnd = head(cnSampleData %>% filter(SegmentEndSupport=='TELOMERE'),1)

    if(nrow(cnEnd) == 1)
    {
      samples[i,]$CopyNumberEnd = cnEnd$CopyNumber
    }
  }

  return (samples)

}


## CODE STARTS HERE
###################

cnArmData = read.csv('~/logs/CN_ANALYSIS_V4.csv')
View(cnArmData)

# cnRawData = read.csv("~/data/prod_cn_data.csv")

patientCancerTypes = read.csv('~/data/patient_cancertypes.csv')
View(patientCancerTypes)
cnArmData = (merge(cnArmData, patientCancerTypes, by.x="SampleId", by.y="SampleId", all.x=TRUE))
cnArmData$CancerType = ifelse(is.na(cnArmData$CancerType), 'N/A', paste(cnArmData$CancerType, sep="")) # set 'N/A' for unknowns

cancerTypes = cnArmData %>% group_by(CancerType,SampleId) %>% summarise(Count=n())
cancerTypes = cancerTypes %>% group_by(CancerType) %>% summarise(SampleCount=n(), SvCount=sum(Count))
View(cancerTypes)

# basic stats about flips and deletion bridges across samples
# cnArmData$FlipCountBucket = 2**round(log(cnArmData$FlipCount,2),0)
cnArmData$FlipCountBucket = (ifelse(cnArmData$FlipCount==0,0,ifelse(cnArmData$FlipCount<=2,2,ifelse(cnArmData$FlipCount<=5,5,
                            ifelse(cnArmData$FlipCount<=50,round(cnArmData$FlipCount/10)*10,round(cnArmData$FlipCount/50)*50)))))

cnArmData$FlipPercent = round(2 * (cnArmData$FlipCount/cnArmData$SegCount),2)
cnArmData$DBPercent = ifelse(cnArmData$FlipCount>0,round(cnArmData$DBCount/cnArmData$FlipCount,2),0)
cnArmData$DelPercent = round(cnArmData$DelCount/cnArmData$SegCount,2)
cnArmData$DupPercent = round(cnArmData$DupCount/cnArmData$SegCount,2)
cnArmData$InvPercent = round(cnArmData$InvCount/cnArmData$SegCount,2)
cnArmData$BndPercent = round(cnArmData$BndCount/cnArmData$SegCount,2)
cnArmData$ConsistentCNLossPercent = ifelse(cnArmData$FlipCount>0,round(cnArmData$MaxFlips/cnArmData$FlipCount,2),0)
cnArmData$NFlipCount = round(cnArmData$FlipCount * cnArmData$ArmLenRatio,1)
cnArmData$NSegCount = round(cnArmData$SegCount * cnArmData$ArmLenRatio,1)
View(cnArmData)

# counts by flip count bucket
cnFlipStats = (cnArmData %>% group_by(FlipCountBucket)
               %>% summarise(Count=n(),
                             SegCount=sum(SegCount),
                             AvgFlipPercent=round(sum(FlipPercent)/n(),2),
                             AvgConsistency=round(sum(ConsistentCNLossPercent)/n(),2),
                             AvgDelPercent=round(sum(DelPercent)/n(),2),
                             AvgDupPercent=round(sum(DupPercent)/n(),2),
                             AvgInvPercent=round(sum(InvPercent)/n(),2),
                             AvgBndPercent=round(sum(BndPercent)/n(),2),
                             AvgDBPercent=round(sum(DBPercent)/n(),3))
               %>% arrange(FlipCountBucket))

View(cnFlipStats)

# findings:
# FlipPercent increases with FlipCount
# consistency (ie a single CN Change) increases with FlipCount
# BND count drops with increases in FlipCount
# 9% of arms have >= 6 flips and 45% of segments are in these arms
# INV:DUP:DEL ratio falls into 2:1:2


cnArmData$FlipPercentBucket = round(cnArmData$FlipPercent/0.2)*0.2
cnArmData$DBPercentBucket = round(cnArmData$DBPercent/0.2)*0.2

# flips and DBs across cohort
cnFlipDBStats = (cnArmData %>% filter(FlipPercent > 0) %>% group_by(FlipPercentBucket,DBPercentBucket)
               %>% summarise(ArmCount=n(),
                             SegCount=sum(SegCount),
                             AvgConsistency=round(sum(ConsistentCNLossPercent)/n(),2),
                             AvgDelPercent=round(sum(DelPercent)/n(),2),
                             AvgDupPercent=round(sum(DupPercent)/n(),2),
                             AvgInvPercent=round(sum(InvPercent)/n(),2),
                             AvgBndPercent=round(sum(BndPercent)/n(),2))
               %>% arrange(FlipPercentBucket,DBPercentBucket))

View(cnFlipDBStats)

cnArmData$DBPercentBucket = round(cnArmData$DBPercent/0.1)*0.1

View(cnArmData %>% filter(DBPercentBucket>= 0.9))

cnDBStats = (cnArmData %>% filter(FlipPercent > 0 & FlipCount >= 5) %>% group_by(DBPercentBucket)
                 %>% summarise(ArmCount=n(),
                               BucketPerc=round(n()/nrow(cnArmData %>% filter(FlipPercent > 0 & FlipCount >= 5)),2),
                               SegCount=sum(SegCount),
                               AvgFlipPercent=round(sum(FlipPercent)/n(),2),
                               AvgConsistency=round(sum(ConsistentCNLossPercent)/n(),2),
                               AvgDelPercent=round(sum(DelPercent)/n(),2),
                               AvgDupPercent=round(sum(DupPercent)/n(),2),
                               AvgInvPercent=round(sum(InvPercent)/n(),2),
                               AvgBndPercent=round(sum(BndPercent)/n(),2))
                 %>% arrange(DBPercentBucket))

View(cnDBStats)

write.csv(cnDBStats, "~/logs/r_output/db_summary.csv")


# DB Length buckets
dbLengths = data_frame()
for(i in 1: nrow(cnArmData))
{
  cnArmDataRow = cnArmData[i,]
  dbLenStr = cnArmDataRow$DelLens

  if(stri_length(dbLenStr) > 0)
  {
    dbLenList = stri_split_fixed(dbLenStr, ";") %>% as.data.frame()
    colnames(dbLenList) <- c("DBLen")

    # print(paste("DB Lengths=", dbLenStr, ", listLen=", nrow(dbLenList), sep=''))
    dbLengths = rbind(dbLengths,dbLenList)
  }

  # if(i > 200)
  #   break
}

nrow(cnArmData %>% filter(DBCount>0))
sum(cnArmData$DBCount)

nrow(dbLengths)
View(dbLengths)

dbLengths$Length = as.integer(as.character(dbLengths$DBLen))
dbLengths$LengthBucket = 2**round(log(dbLengths$Length,2),0)

dbLenBucketed = (dbLengths %>% group_by(LengthBucket)
                 %>% summarise(Count=n()))

View(dbLenBucketed)


dbLenPlot = (ggplot(data = dbLenBucketed, aes(x = LengthBucket))
                      + geom_line(aes(y=Count, colour='LengthBucket'))
                      + scale_x_log10()
                      # + facet_wrap(as.formula(paste("~", facetWrap)))
                      + ylab("DB Count") + labs(title = "Deletion Bridge Lengths")
)

print(dbLenPlot)

dbLenPlot = (ggplot(data = dbLenBucketed, aes(x=LengthBucket, y=Count), fill=LengthBucket)
                + geom_bar(stat = "identity", colour = "black", size=2)
                + ylab("Count") + xlab("DB Length")
                + theme(legend.position="none"))

print(dbLenPlot)


# by chromosomal arm
cnChrArmFlipStats = (cnArmData %>% group_by(Chromosome,Arm,FlipCountBucket)
               %>% summarise(Count=n(),
                             SegCount=sum(SegCount),
                             NSegCount=sum(NSegCount),
                             AvgFlipPerc=round(sum(NFlipCount/NSegCount)/n(),2),
                             AvgDBCount=round(sum(DBCount)/n(),1),
                             MaxDBCount=max(DBCount))
               %>% arrange(FlipCountBucket))

# viewing only significant flip values
View(cnChrArmFlipStats %>% filter(FlipCountBucket>=4))

cnChrArmFlipMedStats = (cnArmData %>% group_by(Chromosome,Arm)
                     %>% summarise(Count=n(),
                                   SegCount=sum(SegCount),
                                   MedFlipCount=median(FlipCount),
                                   AvgFlipCount=round(sum(FlipCount)/n(),1),
                                   MaxFlipCount=max(FlipCount),
                                   FlipLT5=sum(FlipCount<5),
                                   Flip5to10=sum(FlipCount>=5&FlipCount<10),
                                   Flip10to25=sum(FlipCount>=10&FlipCount<25),
                                   Flip25to50=sum(FlipCount>=25&FlipCount<50),
                                   Flip50to100=sum(FlipCount>=50&FlipCount<100),
                                   Flip100Plus=sum(FlipCount>=100),
                                   AvgDBCount=round(sum(DBCount)/n(),1),
                                   MaxDBCount=max(DBCount))
                     %>% arrange(-MedFlipCount))

View(cnChrArmFlipMedStats)


# samples with lots of flips / DBs
topNCount = 200
topNSamples = head(cnArmData %>% arrange(-FlipCount), topNCount)
View(topNSamples)


# Chromothripsis Criteria
# FlipCount >= 10
# major event on the (ie relates to Flip Percent)
# and mostly related to a single fluctuating CN Change (ie MaxFlips vs FlipCount)
# not dominated by BNDs (a little imprecise since no check for overlaps)
# some presence of DBs
# start CN and max CN? not critical since may be other events, although max CN > 8 may be suspicious
# ratio of INVs, DUPs and DELs
ctSamples = cnArmData %>% filter(FlipCount >= 6&FlipPercent >= 0.25)
# nrow(cnArmData)
nrow(ctSamples)
View(ctSamples)

ctSamples = ctSamples %>% filter(InvPercent > DelPercent & InvPercent > DupPercent)
nrow(ctSamples)
View(ctSamples)

ctSamples = ctSamples %>% filter(DBPercent >= 0.15)
nrow(ctSamples)
View(ctSamples)

ctSamples = ctSamples %>% filter(ConsistentCNLossPercent >= 0.7)
nrow(ctSamples)
View(ctSamples)


View(cnArmData %>% filter(FlipCount>=50&FlipPercent <= 0.25))


View(head(ctSamples %>% arrange(-FlipPercent), 400))

View(cnArmData %>% filter(SampleId=="CPCT02020258T"))



# by cancer
cancerCNData = cnArmData %>% filter(CancerType=='Kidney')
View(cancerCNData %>% filter(FlipCount>=4 & Chromosome == '3') %>% arrange(-FlipCount))

cancerChrArmStats = (cancerCNData %>% group_by(Chromosome,Arm)
                        %>% summarise(Count=n(),
                                      SegCount=sum(SegCount),
                                      MedFlipCount=median(FlipCount),
                                      AvgFlipCount=round(sum(FlipCount)/n(),1),
                                      MaxFlipCount=max(FlipCount),
                                      AvgDBCount=round(sum(DBCount)/n(),1),
                                      MaxDBCount=max(DBCount))
                        %>% arrange(-MedFlipCount))

View(cancerChrArmStats)


# access to raw CN data
cnSection = get_copy_number_range(cnRawData, "CPCT02110023T", "21", 15450001)
View(cnSection)

cnRawPlot = get_copy_number_plot(cnRawData, "CPCT02110023T", "21", 15450001)
print(cnRawPlot)




# copy number plotting

cnData = read.csv('~/data/cn_stressed_arms.csv')
colnames(cnData)

print(get_copy_number_plot(cnData, 'CPCT02020258T', '13', 0, 200e6))

cnSection = get_copy_number_range(cnData, 'CPCT02020258T', '13', 0, 200e6)
cnSection = cnSection %>% filter(copyNumber <= 20)
View(cnSection)
plot_copy_number_section(cnSection)

sampleArmSvs = cnData %>% filter(sampleId=='CPCT02020258T' & chromosome=='13')
View(sampleArmSvs)

plot = (ggplot(data = sampleArmSvs, aes(x=start, y=copyNumber))
  + geom_line())

print(plot)

+ geom_bar(stat = "identity", colour = "black", size = 0.2)
  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  + ylab("Copy Number") + xlab("Position")
  + theme(legend.position="none"))

sigStatsPlot = (ggplot(data = cancerSigStats, aes(x = SigName, y = SvCount, group = 1), fill = SigName)
                + geom_bar(stat = "identity", colour = "black", size = 0.2)
                + theme(axis.text.x = element_text(angle = 90, hjust = 1))
                + ylab("SV Count") + xlab("Signature") + ggtitle(title)
                + theme(legend.position="none"))
