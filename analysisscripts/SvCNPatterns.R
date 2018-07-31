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



svData = read.csv('~/logs/CLUSTER_V23.csv')
nrow(svData)

svData = svData %>% filter(svData$Type=='BND'|svData$ArmEnd==svData$ArmStart)

svData$IsLINE = ifelse(svData$LEStart!='false'|svData$LEEnd!='false',1,0)
svData$IsFS = ifelse(svData$FSStart!='false'|svData$FSEnd!='false',1,0)
svData$IsSpan = ifelse(svData$TransType=='SPAN',1,0)
svData$DoubleDupBE = ifelse(svData$DupBEStart=='true'&svData$DupBEEnd=='true',1,0)
svData$SingleDupBE = ifelse(svData$DoubleDupBE==0&(svData$DupBEStart=='true'|svData$DupBEEnd=='true'),1,0)

svData$Length = ifelse(svData$Type=='BND'|svData$Type=='INS'|svData$ArmEnd!=svData$ArmStart, -1, svData$PosEnd-svData$PosStart)


# first group SVs by distinct chromosomal arm
svArmData = svData
svArmData$Chr = svArmData$ChrStart
svArmData$Arm = svArmData$ArmStart
svArmData$ArmCount = svArmData$ArmCountStart
svArmData$ArmExpected = svArmData$ArmExpStart
svArmData$Position = svArmData$PosStart

svArmBndData = svData %>% filter(Type=='BND')
svArmBndData$Chr = svArmBndData$ChrEnd
svArmBndData$Arm = svArmBndData$ArmEnd
svArmBndData$ArmCount = svArmBndData$ArmCountEnd
svArmBndData$ArmExpected = svArmBndData$ArmExpEnd
svArmBndData$Position = svArmBndData$PosEnd

# merge rows prior to arm grouping
combinedArmData = rbind(svArmData, svArmBndData)

# determine IsStressed using poisson distribution using expected vs actual SV counts per arm
combinedArmData$ArmExpected = ifelse(combinedArmData$ArmExpected>0,combinedArmData$ArmExpected,0.1)

combinedArmData$StressedPoissonProb = round(1 - ppois(combinedArmData$ArmCount - 1, combinedArmData$ArmExpected),4)
combinedArmData$IsStressed = ifelse(combinedArmData$StressedPoissonProb <= 0.001 & combinedArmData$ArmCount >= 10,1,0)

# prepare arm stats
armStats = (combinedArmData %>% group_by(SampleId, Chr, Arm)
            %>% summarise(SvCount=n(),
                          ClusteredPerc=round(sum(ClusterCount>1)/n(),2),
                          DelCount=sum(Type=='DEL'),
                          DupCount=sum(Type=='DUP'),
                          InsCount=sum(Type=='INS'),
                          InvCount=sum(Type=='INV'),
                          BndCount=sum(Type=='BND'),
                          DupBECount=sum(SingleDupBE|DoubleDupBE),
                          MinPosition=min(Position),
                          MaxPosition=max(Position),
                          MaxClusterCount=max(ClusterCount),
                          AvgClusterCount=round(sum(ClusterCount)/n(),0),
                          LEPerc=round(sum(IsLINE==1)/n(),3),
                          FSPerc=round(sum(IsFS==1)/n(),3),
                          IsStressed=max(IsStressed))
            %>% arrange(SampleId, Chr, Arm))

View(armStats)
nrow(armStats)

ncArmData = armStats %>% filter(MaxClusterCount==1&BndCount==0)
nrow(ncArmData)
sum(ncArmData$SvCount)
View(ncArmData)

simpleArms = armStats %>% filter(SvCount==2&DupCount==2&DupBECount==0)
View(simpleArms)

View(svData %>% filter(SampleId=='CPCT02010037T'&ChrStart==5))
View(combinedArmData %>% filter(SampleId=='CPCT02020667T'&Chr==17&Arm=='Q'))



View(armStats %>% filter(SvCount==1&InvCount==1&DupBECount==0))



rm(cnRawData)
cnRawData = read.csv("~/data/prod_cn_data3.csv")

cnRawData$MinorAllele = (1-cnRawData$ActualBaf)*cnRawData$CopyNumber
cnRawData$MajorAllele = cnRawData$ActualBaf*cnRawData$CopyNumber

View(cnRawData %>% filter(SampleId=='CPCT02020667T'&Chromosome==17))
View(cnRawData %>% filter(SampleId=='CPCT02010037T'&Chromosome==5))
