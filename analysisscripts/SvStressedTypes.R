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

detach("package:svnmf", unload=TRUE);
library(svnmf)


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

####################
### CODE STARTS HERE
####################

rm(svData)

# SV data file
svData = read.csv('~/logs/CLUSTER_V16.csv')

# filter out multiple biopsy (approximately)
svData = svData %>% filter(!grepl("DRUP", SampleId)&!grepl("TIII", SampleId)&!grepl("TII", SampleId))

# FILTER FOR PONCount <2 for all subsequent analyses - no longer required since done already
svData = svData %>% filter(PONCount<2)

svData$IsLINE = ifelse(svData$LEStart!='false'|svData$LEEnd!='false',1,0)
svData$IsFS = ifelse(svData$FSStart!='false'|svData$FSEnd!='false',1,0)
svData$Length = ifelse(svData$Type=='BND'|svData$Type=='INS'|svData$ArmEnd!=svData$ArmStart, -1, svData$PosEnd-svData$PosStart)
svData$IsDB=ifelse(svData$NearestDBLen>-1&(svData$NearestDBLen<svData$NearestTILen|svData$NearestTILen<30),1,0)
svData$IsTI=ifelse(svData$NearestTILen>=30&svData$IsDB==0,1,0)


svData$ClusterNone=ifelse(svData$ClusterCount==1,1,0)
svData$ClusterSmall=ifelse(svData$ClusterCount>1&svData$ClusterCount<=3,1,0)
svData$ClusterLarge=ifelse(svData$ClusterCount>3,1,0)
View(svData)


# templated insertions
svData$ClusterSize = ifelse(svData$ClusterCount==1, 'None', ifelse(svData$ClusterCount < 4, 'Small', 'Large'))
svData$ArmExpStart = ifelse(svData$ArmExpStart>0,svData$ArmExpStart,0.1)
svData$StressedPoissonProb = round(1 - ppois(svData$ArmCountStart - 1, svData$ArmExpStart),4)
View(svData)

svData$IsStressed = ifelse(svData$ArmCountStart>=10 & svData$StressedPoissonProb <= 0.001,1,0)
nrow(svData %>% filter(IsStressed==0))

svData$TILengthBucket=ifelse(svData$IsTI==1&svData$NearestTILen>0,2**round(log(svData$NearestTILen,2),0),0)

tiBucketed = (svData %>% filter(IsLINE==0 & IsStressed==0 & TILengthBucket > 0 & TILengthBucket <= 1e7) %>% group_by(TILengthBucket)
              %>% summarise(Count=n(),
                            ClusterNone=sum(ClusterCount==1),
                            Cluster2_3=sum(ClusterCount>1&ClusterCount>4),
                            Cluster4_10=sum(ClusterCount>=4&ClusterCount<10),
                            Cluster10_25=sum(ClusterCount>=10&ClusterCount<25),
                            Cluster25_plus=sum(ClusterCount>=25))
              %>% arrange(TILengthBucket))

View(tiBucketed)

tiStressedBucketed = (svData %>% filter(IsLINE==0 & IsStressed==1 & TILengthBucket > 0 & TILengthBucket <= 1e7) %>% group_by(TILengthBucket)
              %>% summarise(Count=n(),
                            ClusterNone=sum(ClusterCount==1),
                            Cluster2_3=sum(ClusterCount>1&ClusterCount>4),
                            Cluster4_10=sum(ClusterCount>=4&ClusterCount<10),
                            Cluster10_25=sum(ClusterCount>=10&ClusterCount<25),
                            Cluster25_plus=sum(ClusterCount>=25))
              %>% arrange(TILengthBucket))

View(tiStressedBucketed)

tiLenClusteredPlot = (ggplot(data = tiStressedBucketed %>% filter(TILengthBucket <= 1e4), aes(x = TILengthBucket))
                 + geom_line(aes(y=Cluster2_3, colour='Cluster2_3'))
                 + geom_line(aes(y=Cluster4_10, colour='Cluster4_10'))
                 + geom_line(aes(y=Cluster10_25, colour='Cluster10_25'))
                 + geom_line(aes(y=Cluster25_plus, colour='Cluster25_plus'))
                 + scale_x_log10()
                 # + facet_wrap(as.formula(paste("~", facetWrap)))
                 + ylab("SV Count") + labs(title = "Stressed TI Length by Cluster Count")
)

print(tiLenClusteredPlot)


# templated insertions from cluster-count 2, back to same location
cc2Data = svData %>% filter(ClusterCount==2)
nrow(cc2Data)

# group by clusterId and then look at the BE pairs
cc2Clusters = (cc2Data %>% group_by(SampleId,ClusterId)
               %>% summarise(Id1=first(Id),
                             Id2=last(Id),
                             ChrStart1=first(ChrStart),ChrEnd1=first(ChrEnd),
                             ChrStart2=last(ChrStart),ChrEnd2=last(ChrEnd),
                             CrossChr=ifelse(first(Type)=='BND'|last(Type)=='BND',1,0),
                             IsStressed=ifelse(first(IsStressed)==1,1,0),
                             SS_IsTI=ifelse(first(ChrStart)==last(ChrStart),ifelse((first(PosStart)<last(PosStart) & first(OrientStart)==-1 & last(OrientStart)==1)|(first(PosStart)>last(PosStart) & first(OrientStart)==1 & last(OrientStart)==-1),1,0),0),
                             EE_IsTI=ifelse(first(ChrEnd)==last(ChrEnd),ifelse((first(PosEnd)<last(PosEnd) & first(OrientEnd)==-1 & last(OrientEnd)==1)|(first(PosEnd)>last(PosEnd) & first(OrientEnd)==1 & last(OrientEnd)==-1),1,0),0),
                             SE_IsTI=ifelse(first(ChrStart)==last(ChrEnd),ifelse((first(PosStart)<last(PosEnd) & first(OrientStart)==-1 & last(OrientEnd)==1)|(first(PosStart)>last(PosEnd) & first(OrientStart)==1 & last(OrientEnd)==-1),1,0),0),
                             ES_IsTI=ifelse(first(ChrEnd)==last(ChrStart),ifelse((first(PosEnd)<last(PosStart) & first(OrientEnd)==-1 & last(OrientStart)==1)|(first(PosEnd)>last(PosStart) & first(OrientEnd)==1 & last(OrientStart)==-1),1,0),0),
                             SS_IsDB=ifelse(first(ChrStart)==last(ChrStart),ifelse((first(PosStart)<last(PosStart) & first(OrientStart)==1 & last(OrientStart)==-1)|(first(PosStart)>last(PosStart) & first(OrientStart)==-1 & last(OrientStart)==1),1,0),0),
                             EE_IsDB=ifelse(first(ChrEnd)==last(ChrEnd),ifelse((first(PosEnd)<last(PosEnd) & first(OrientEnd)==1 & last(OrientEnd)==-1)|(first(PosEnd)>last(PosEnd) & first(OrientEnd)==-1 & last(OrientEnd)==1),1,0),0),
                             SE_IsDB=ifelse(first(ChrStart)==last(ChrEnd),ifelse((first(PosStart)<last(PosEnd) & first(OrientStart)==1 & last(OrientEnd)==-1)|(first(PosStart)>last(PosEnd) & first(OrientStart)==-1 & last(OrientEnd)==1),1,0),0),
                             ES_IsDB=ifelse(first(ChrEnd)==last(ChrStart),ifelse((first(PosEnd)<last(PosStart) & first(OrientEnd)==1 & last(OrientStart)==-1)|(first(PosEnd)>last(PosStart) & first(OrientEnd)==-1 & last(OrientStart)==1),1,0),0),
                             SS_Len=ifelse(first(ChrStart)==last(ChrStart),abs(first(PosStart)-last(PosStart)),1e7),
                             EE_Len=ifelse(first(ChrEnd)==last(ChrEnd),abs(first(PosEnd)-last(PosEnd)),1e7),
                             SE_Len=ifelse(first(ChrStart)==last(ChrEnd),abs(first(PosStart)-last(PosEnd)),1e7),
                             ES_Len=ifelse(first(ChrEnd)==last(ChrStart),abs(first(PosEnd)-last(PosStart)),1e7)
                             )
               %>% arrange(SampleId,ClusterId))

View(cc2Clusters)
nrow(cc2Clusters)

cc2Clusters$SS_Len = ifelse(cc2Clusters$SS_IsTI==1|cc2Clusters$SS_IsDB==1,cc2Clusters$SS_Len,1e7)
cc2Clusters$EE_Len = ifelse(cc2Clusters$EE_IsTI==1|cc2Clusters$EE_IsDB==1,cc2Clusters$EE_Len,1e7)
cc2Clusters$SE_Len = ifelse(cc2Clusters$SE_IsTI==1|cc2Clusters$SE_IsDB==1,cc2Clusters$SE_Len,1e7)
cc2Clusters$ES_Len = ifelse(cc2Clusters$ES_IsTI==1|cc2Clusters$ES_IsDB==1,cc2Clusters$ES_Len,1e7)

# both ends must be near each other
cc2Tmp2 = cc2Clusters %>% filter((SS_Len<=1e4&EE_Len<=1e4)|(SE_Len<=1e4&ES_Len<=1e4))
nrow(cc2Tmp2)
View(cc2Tmp2)

cc2Tmp2$TICount = cc2Tmp2$SS_IsTI+cc2Tmp2$EE_IsTI+cc2Tmp2$SE_IsTI+cc2Tmp2$ES_IsTI
cc2Tmp2$DBCount = cc2Tmp2$SS_IsDB+cc2Tmp2$EE_IsDB+cc2Tmp2$SE_IsDB+cc2Tmp2$ES_IsDB

# must have a TI in one of the pairs and either a TI or DB in the other
# cc2Clusters3 = cc2Clusters2 %>% filter(SS_IsTI==1|EE_IsTI==1|SE_IsTI==1|ES_IsTI==1)
cc2Tmp3 = cc2Tmp2 %>% filter(TICount>1|(TICount==1&DBCount>0))
nrow(cc2Tmp3)
View(cc2Tmp3)

# tangent to show DSBs (which cannot also be TIs)
cc2DSBs = cc2Tmp2 %>% filter(DBCount==4)
nrow(cc2DSBs)
View(cc2DSBs)

# use the shortest length to infer what sort of link is made (ie a TI or DB)

# whether linked start-to-start or start-to-end positions
cc2Tmp3$SS_EE = ifelse((cc2Tmp3$SS_IsTI==1|cc2Tmp3$SS_IsDB==1)&(cc2Tmp3$EE_IsTI==1|cc2Tmp3$EE_IsDB==1)&(cc2Tmp3$SS_Len+cc2Tmp3$EE_Len<cc2Tmp3$SE_Len+cc2Tmp3$ES_Len),1,0)

# is this 2x templated insertions, or a TI and a DB
cc2Tmp3$IsDoubleTI = ifelse(cc2Tmp3$SS_EE==1,ifelse(cc2Tmp3$SS_IsTI==1&cc2Tmp3$EE_IsTI==1,1,0),ifelse(cc2Tmp3$SE_IsTI==1&cc2Tmp3$ES_IsTI==1,1,0))

# get the length of the templated insertion
cc2Tmp3$TILen = ifelse(cc2Tmp3$IsDoubleTI==1,ifelse(cc2Tmp3$SS_EE==1,cc2Tmp3$SS_Len,cc2Tmp3$SE_Len),
                       ifelse(cc2Tmp3$SS_EE==1,ifelse(cc2Tmp3$SS_IsTI==1,cc2Tmp3$SS_Len,cc2Tmp3$EE_Len),ifelse(cc2Tmp3$SE_IsTI==1,cc2Tmp3$SE_Len,cc2Tmp3$ES_Len)))

# get the other length (either a DB or another TI)
cc2Tmp3$OtherLen = ifelse(cc2Tmp3$IsDoubleTI==1,ifelse(cc2Tmp3$SS_EE==1,cc2Tmp3$EE_Len,cc2Tmp3$ES_Len),
                          ifelse(cc2Tmp3$SS_EE==1,ifelse(cc2Tmp3$SS_IsTI==1,cc2Tmp3$EE_Len,cc2Tmp3$SS_Len),ifelse(cc2Tmp3$SE_IsTI==1,cc2Tmp3$ES_Len,cc2Tmp3$SE_Len)))

cc2Tmp3$TILenBucket = 2**round(log(cc2Tmp3$TILen,2),0)
cc2Tmp3$OtherLenBucket = 2**round(log(cc2Tmp3$OtherLen,2),0)

View(cc2Tmp3)

# overall stats
cc2Stats = (cc2Tmp3 %>% group_by(IsStressed,CrossChr,IsDoubleTI)
                      %>% summarise(Count=n(),
                                    AvgTILen=round(sum(TILen)/n(),0),
                                    TL_LT200=round(sum(TILen<=200)/n(),2),
                                    TL_200_1K=round(sum(TILen>200&TILen<=1e3)/n(),2),
                                    TL_GT1K=round(sum(TILen>1e3)/n(),2),
                                    AvgOtherLen=round(sum(OtherLen)/n(),0),
                                    OT_LT200=round(sum(OtherLen<=200)/n(),2),
                                    OT_200_1K=round(sum(OtherLen>200&OtherLen<=1e3)/n(),2),
                                    OT_GT1K=round(sum(OtherLen>1e3)/n(),2))
                      %>% arrange(IsStressed,CrossChr,IsDoubleTI))

View(cc2Stats)

# distribution of TI lengths, by whether cross arm or not
cc2StatsByTILength = (cc2Tmp3 %>% group_by(CrossChr,TILenBucket)
            %>% summarise(Count=n(),
                          AvgOtherLen=round(sum(OtherLen)/n(),0))
            %>% arrange(CrossChr,TILenBucket))

View(cc2StatsByTILength)

# distribution of TI lengths, by whether stressed or not
cc2StressedByTILength = (cc2Tmp3 %>% group_by(IsStressed,TILenBucket)
                      %>% summarise(Count=n(),
                                    AvgOtherLen=round(sum(OtherLen)/n(),0))
                      %>% arrange(IsStressed,TILenBucket))

View(cc2StressedByTILength)

# distribution of DB lengths
cc2StatsByDBLength = (cc2Tmp3 %>% filter(IsDoubleTI==0) %>% group_by(OtherLenBucket)
                      %>% summarise(Count=n(),
                                    AvgITLen=round(sum(TILen)/n(),0))
                      %>% arrange(OtherLenBucket))

View(cc2StatsByDBLength)


cc2Stats = (cc2Tmp3 %>% group_by(IsDoubleTI,TILenBucket,OtherLenBucket)
            %>% summarise(Count=n(),
                        AvgLen1=round(sum(BridgeLen1)/n(),0),
                        AvgLen2=round(sum(BridgeLen2)/n(),0))
            %>% arrange(IsDoubleTI,TILenBucket,OtherLenBucket))

View(cc2Stats)



View(ssTmp1)


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
combinedArmData$IsStressedOld = ifelse(combinedArmData$ArmCount>=10 & combinedArmData$ArmCount >= 2.5 * combinedArmData$ArmExpected,1,0)
combinedArmData$IsStressed = ifelse(combinedArmData$StressedPoissonProb <= 0.001 & combinedArmData$ArmCount >= 10,1,0)


# prepare arm stats
armStats = (combinedArmData %>% group_by(SampleId, Chr, Arm)
            %>% summarise(SvCount=n(),
                          MaxCN=round(max((AdjCNStart+AdjCNEnd)*0.5),2),
                          AvgCN=round(sum((AdjCNStart+AdjCNEnd)*0.5)/n(),2),
                          AvgPloidy=round(sum(Ploidy)/n(),2),
                          MaxInvCN=round(max(ifelse(Type=='INV',(AdjCNStart+AdjCNEnd)*0.5,0)),2),
                          InvMinPosStart=min(ifelse(Type=='INV'&AdjCNStart>=50&AdjCNEnd>=50&NearestTILen>=0,PosStart,3e8)),
                          InvMaxPosEnd=max(ifelse(Type=='INV'&AdjCNStart>=50&AdjCNEnd>=50&NearestTILen>=0,PosEnd,-1)),
                          ClusteredPerc=round(sum(ClusterCount>1)/n(),2),
                          MaxClusterCount=max(ClusterCount),
                          AvgClusterCount=round(sum(ClusterCount)/n(),0),
                          LEPerc=round(sum(IsLINE=='true')/n(),3),
                          FSPerc=round(sum(IsFS=='true')/n(),3),
                          BndCount=sum(Type=='BND'),
                          BndPerc=round(sum(Type=='BND')/n(),2),
                          IsStressed=max(IsStressed))
            %>% arrange(SampleId, Chr, Arm))

View(armStats)



combinedArmData$IsStressed = ifelse(combinedArmData$ArmCount>=10 & combinedArmData$StressedPoissonProb<=0.001,1,0)
nrow(combinedArmData %>% filter(IsStressed==1))
# combinedArmData$IsStressed = ifelse(combinedArmData$ArmCount>=combinedArmData$ArmExpected*2.5&combinedArmData$ArmCount>=10,1,0)


# prepare arm stats
armStats = (combinedArmData %>% group_by(SampleId, Chr, Arm)
            %>% summarise(SvCount=n(),
                      MaxCN=round(max((AdjCNStart+AdjCNEnd)*0.5),2),
                      AvgCN=round(sum((AdjCNStart+AdjCNEnd)*0.5)/n(),2),
                      AvgPloidy=round(sum(Ploidy)/n(),2),
                      MaxInvCN=round(max(ifelse(Type=='INV',(AdjCNStart+AdjCNEnd)*0.5,0)),2),
                      InvMinPosStart=min(ifelse(Type=='INV'&AdjCNStart>=50&AdjCNEnd>=50&NearestTILen>=0,PosStart,3e8)),
                      InvMaxPosEnd=max(ifelse(Type=='INV'&AdjCNStart>=50&AdjCNEnd>=50&NearestTILen>=0,PosEnd,-1)),
                      ClusteredPerc=round(sum(ClusterCount>1)/n(),2),
                      MaxClusterCount=max(ClusterCount),
                      AvgClusterCount=round(sum(ClusterCount)/n(),0),
                      LEPerc=round(sum(IsLINE=='true')/n(),3),
                      FSPerc=round(sum(IsFS=='true')/n(),3),
                      BndCount=sum(Type=='BND'),
                      BndPerc=round(sum(Type=='BND')/n(),2),
                      IsStressed=max(IsStressed))
            %>% arrange(SampleId, Chr, Arm))

View(armStats)

# now filter for high CN vs avg ploidy
allDmCandidates = armStats %>% filter(InvMinPosStart < 3e8 & InvMaxPosEnd > -1 & MaxInvCN >= 6 * AvgPloidy) %>% arrange(BndPerc)
View(allDmCandidates)

# and only arms with few BNDs to keep things simple
dmCandidates = allDmCandidates %>% filter(BndPerc <= 0.1)
View(dmCandidates)


# supplement with copy number info for start and end of each arm
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")

dmCandidates = annotateWithCopyNumber(dbProd, dmCandidates)
View(dmCandidates)

# view all SVs within the bands of the high CN INVs
sampleId = 'CPCT02010549T'

dmSampleData = dmCandidates %>% filter(SampleId==sampleId)

dmAllArmSVs = combinedArmData %>% filter(SampleId==sampleId,Chr==dmSampleData$Chr,Arm==dmSampleData$Arm) %>% arrange(Id)
View(dmAllArmSVs)

dmContainedArmSVs = (dmAllArmSVs %>% filter((dmAllArmSVs$ChrStart == dmSampleData$Chr & dmAllArmSVs$ArmStart == dmSampleData$Arm
                                           & dmAllArmSVs$PosStart >= dmSampleData$InvMinPosStart & dmAllArmSVs$PosEnd <= dmSampleData$InvMaxPosEnd)
                                          | (dmAllArmSVs$ChrEnd == dmSampleData$Chr & dmAllArmSVs$ArmEnd == dmSampleData$Arm
                                           & dmAllArmSVs$PosEnd >= dmSampleData$InvMinPosStart & dmAllArmSVs$PosEnd <= dmSampleData$InvMaxPosEnd)) %>% arrange(Id))
View(dmContainedArmSVs)


# view most affected chromosomal arms and positions
mostAffectedChrArms = (allDmCandidates %>% group_by(Chr,Arm)
                       %>% summarise(ArmCount=n(),
                                     BndCount=sum(BndCount),
                                     StressedArmCount=sum(IsStressed),
                                     AvgPosStart=round(sum(InvMinPosStart)/n(),0),
                                     AvgPosEnd=round(sum(InvMaxPosEnd)/n(),0),
                                     AvgPloidy=round(sum(AvgPloidy)/n(),2))
                       %>% arrange(-ArmCount))

View(mostAffectedChrArms)

# samples with more BNDs and arms affected could be complicated DMs or BFBs
# probably need to adjust this to show connections to any other arm in the same via a BND from within the potential DM site
multiArmSamples = (allDmCandidates %>% group_by(SampleId)
                       %>% summarise(ArmCount=n(),
                                     BndCount=sum(BndCount),
                                     StressedArmCount=sum(IsStressed),
                                     AvgPloidy=round(sum(AvgPloidy)/n(),2))
                       %>% arrange(-ArmCount))

View(multiArmSamples)


# Possible Chromothripsis - defined as shattering and general CN loss across a single arm
ctAllArms = (combinedArmData %>% group_by(SampleId, Chr, Arm)
            %>% summarise(SvCount=n(),
                          DelCount=sum(Type=='DEL'),
                          DelPerc=round(sum(Type=='DEL')/n(),2),
                          DupCount=sum(Type=='DUP'),
                          DupPerc=round(sum(Type=='DUP')/n(),2),
                          InvCount=sum(Type=='INV'),
                          InvPerc=round(sum(Type=='INV')/n(),2),
                          BndCount=sum(Type=='BND'),
                          BndPerc=round(sum(Type=='BND')/n(),2),
                          LEPerc=round(sum(IsLINE=='true')/n(),3),
                          FSPerc=round(sum(IsFS=='true')/n(),3),
                          MaxCN=round(max((AdjCNStart+AdjCNEnd)*0.5),2),
                          AvgCN=round(sum((AdjCNStart+AdjCNEnd)*0.5)/n(),2),
                          AvgPloidy=round(sum(Ploidy)/n(),2),
                          AvgSvLength=round(sum(ifelse(Length>-1,Length,0))/sum(Length>-1),0),
                          DBPerc=round(sum(IsDB)/n(),2),
                          TIPerc=round(sum(IsTI)/n(),2),
                          AvgDBLen=round(sum(ifelse(IsDB==1,NearestDBLen,0))/sum(IsDB==1),0),
                          ClusteredPerc=round(sum(ClusterCount>1)/n(),2),
                          IsStressed=round(sum(IsStressed)/n(),2))
            %>% arrange(SampleId, Chr, Arm))

View(ctAllArms)

# now filter for the INV:DUP:DEL 2:1:1 ratio, low avg ploidy, few BNDs and predominantly deletion bridges
chrThripCandidates = ctAllArms %>% filter(InvCount >= 10 & abs(InvPerc - DupPerc - DelPerc) <= 0.2 & abs(DupPerc - DelPerc) <= 0.05)
chrThripCandidates = chrThripCandidates %>% filter(AvgPloidy <= 2.5 & BndPerc <= 0.1 & DBPerc >= 0.7) %>% arrange(-InvCount)

# annotate with the arm start & end copy numbers
chrThripCandidates = annotateWithCopyNumber(dbProd, chrThripCandidates)

# stricter check of arm CN
chrThripCandidates = chrThripCandidates %>% filter(CopyNumberStart<=2.2&CopyNumberEnd<=2.2)
View(chrThripCandidates)



# Samples dominated by BNDs in non-stressed arms
bndStats = (svData %>% group_by(SampleId) %>% filter(IsLINE=='false')
             %>% summarise(SvCount=n(),
                           StressedPerc=round(sum(IsStressed==1)/n(),2),
                           NonStrPerc=round(sum(IsStressed==0)/n(),2),
                           BndCount=sum(IsStressed==0&Type=='BND'),
                           BndPerc=round(sum(IsStressed==0&Type=='BND')/sum(IsStressed==0),2),
                           ClusteredPerc=round(sum(IsStressed==0&ClusterCount>1)/sum(IsStressed==0),2),
                           BndClusteredPerc=round(sum(IsStressed==0&Type=='BND'&ClusterCount>1)/sum(IsStressed==0),2),
                           ShortTIPerc=round(sum(IsStressed==0&Type=='BND'&NearestTILen<=1e3&NearestTILen>=30)/sum(IsStressed==0),2),
                           AvgTILen=round(sum(NearestTILen)/sum(IsStressed==0),0),
                           AvgDBLen=round(sum(NearestDBLen)/sum(IsStressed==0),0))
             %>% arrange(-BndPerc))

bndNonStressed = bndStats %>% filter(SvCount>5&BndPerc>=0.5)

bndNonStressed$ShortBndTIRatio = bndNonStressed$ShortTIPerc/bndNonStressed$BndPerc
bndNonStressed = bndNonStressed %>% filter(ShortBndTIRatio>=0.5) %>% arrange(-ShortBndTIRatio)

View(bndNonStressed)


# LINE element analysis
# number of non-LINE SVs part of clusters with LINE SVs
# and number of non-LINE arms affected
svData$IsLINE = ifelse(svData$LEStart!='false'|svData$LEEnd!='false',1,0)
nrow(svData)
nrow(svData %>% filter(IsLINE==1 & ClusterCount > 1))

clustersWithLE = (svData %>% filter(ClusterCount > 1) %>% group_by(SampleId,ClusterId)
                 %>% summarise(Count=n(),
                               LECount=sum(IsLINE==1),
                               LEPerc=round(sum(IsLINE==1)/n(),2))
                 %>% arrange(SampleId,ClusterId))

View(clustersWithLE %>% filter(LEPerc > 0))

clustersWithLE = clustersWithLE %>% filter(LEPerc > 0)
clustersWithLE = unite(clustersWithLE, "SampClustId", SampleId, ClusterId, sep='_')
nrow(clustersWithLE)
View(clustersWithLE)

leData = svData %>% filter(paste(SampleId, ClusterId, sep='_') %in% clustersWithLE$SampClustId)
nrow(leData)
View(leData)

# check filtering process
leTmp1 = (leData %>% filter(ClusterCount > 1) %>% group_by(SampleId,ClusterId)
                  %>% summarise(Count=n(),
                                LECount=sum(IsLINE==1),
                                LEPerc=round(sum(IsLINE==1)/n(),2))
                  %>% arrange(SampleId,ClusterId))

View(leTmp1)

# now group by arm to see non-LINE element arms affected
leData$ClusterSize = 2**round(log(leData$ClusterCount,2),0)

leArmData = (leData %>% group_by(SampleId,ClusterId,ChrStart,ChrEnd)
             %>% summarise(SvCount=n(),
                           ClusterCount=first(ClusterCount),
                           NonLECount=sum(IsLINE==0),
                           LECount=sum(IsLINE==1),
                           NonLEPerc=round(sum(IsLINE==0)/n(),2),
                           LEPerc=round(sum(IsLINE==1)/n(),2))
             %>% arrange(SampleId,ClusterId,ChrStart,ChrEnd))

View(leArmData)
nrow(leArmData)

leArmData$LEStatus = ifelse(leArmData$LEPerc==0,'None',ifelse(leArmData$LEPerc==1,'All','Mixed'))

leTmp2 = (leArmData %>% group_by(LEStatus)
                    %>% summarise(ArmCount=n(),
                                  SvCount=sum(SvCount))
                    %>% arrange(-SvCount))

View(leTmp2)

leTmp3 = (leArmData %>% group_by(SampleId,LEStatus,ChrStart,ChrEnd)
          %>% summarise(ClusterCount=n(),
                        SvCount=sum(SvCount))
          %>% arrange(SampleId,LEStatus,ChrStart,ChrEnd))

View(leTmp3)

leTmp4 = (leTmp3 %>% group_by(SampleId,LEStatus)
          %>% summarise(ArmCount=n(),
                        SvCount=sum(SvCount))
          %>% arrange(SampleId,LEStatus))

View(leTmp4)

leTmp4 = (leTmp3 %>% group_by(LEStatus))

clustersWithLE = (svData %>% group_by(SampleId,ClusterId)
                  %>% summarise(Count=n(),
                                LECount=sum(IsLINE==1),
                                LEPerc=round(sum(IsLINE==1)/n(),2))
                  %>% arrange(SampleId,ChrStart,ChrEnd,ClusterId))


# summary data per sample of remaining clustering and stressed arms

# first get all stats by arm
summaryAllArms = (combinedArmData %>% group_by(SampleId, Chr, Arm)
             %>% summarise(SvCount=n(),
                           DelCount=sum(Type=='DEL'),
                           DupCount=sum(Type=='DUP'),
                           InvCount=sum(Type=='INV'),
                           BndCount=sum(Type=='BND'),
                           LECount=sum(IsLINE=='true'),
                           FSCount=sum(IsFS=='true'),
                           DBCount=sum(IsDB),
                           TICount=sum(IsTI),
                           ClusterNone=sum(ClusterCount==1),
                           ClusterSmall=sum(ClusterCount>1&ClusterCount<4),
                           ClusterLarge=sum(ClusterCount>=4),
                           IsStressed=max(IsStressed))
             %>% arrange(SampleId, Chr, Arm))

View(summaryAllArms)

# then group into sample info
sampleSummaryData = (summaryAllArms %>% group_by(SampleId)
                 %>% summarise(SvCount=sum(SvCount),
                               DelPerc=round(sum(DelCount)/sum(SvCount),2),
                               DupPerc=round(sum(DupCount)/sum(SvCount),2),
                               InvPerc=round(sum(InvCount)/sum(SvCount),2),
                               BndPerc=round(sum(BndCount)/sum(SvCount),2),
                               LEPerc=round(sum(LECount)/sum(SvCount),2),
                               FSPerc=round(sum(FSCount)/sum(SvCount),2),
                               DBPerc=round(sum(DBCount)/sum(SvCount),2),
                               TIPerc=round(sum(TICount)/sum(SvCount),2),
                               ClusterNonePerc=round(sum(ClusterNone)/sum(SvCount),2),
                               ClusterSmallPerc=round(sum(ClusterSmall)/sum(SvCount),2),
                               ClusterLargePerc=round(sum(ClusterLarge)/sum(SvCount),2),
                               StressedArmPerc=round(sum(IsStressed)/n(),2))
                 %>% arrange(SampleId))

View(sampleSummaryData)



# DUPs and DELs - relationship with stressed arms
stdSigData = svData
# View(svData)

stdSigData$LengthBucket = ifelse(stdSigData$Length <= 0, 0, 10**(round(log(stdSigData$Length,10),0)))
stdSigData$ClusterSize = ifelse(stdSigData$ClusterCount==1, 'None', ifelse(stdSigData$ClusterCount < 4, 'Small', 'Large'))
stdSigData$IsDB = ifelse(stdSigData$IsDB==1&stdSigData$NearestDBLen<=1e2,1,0)
# svData$HighMutLoad = ifelse()

dbExamples = stdSigData %>% filter(IsDB==1)
View(dbExamples)

# get sample arm data together
stdSigArmData = (stdSigData %>% group_by(SampleId, ChrStart, ArmStart, IsStressed) # ClusterSize
                 %>% summarise(SvCount=n(),
                               DB=sum(IsDB==1),
                               Del_LT10K=sum(Type=='DEL'&Length<=1e4),
                               Del_10Kto100K=sum(Type=='DEL'&Length>1e4&Length<=1e5),
                               Del_100Kto500K=sum(Type=='DEL'&Length>1e5&Length<=5e5),
                               Del_500Kto5M=sum(Type=='DEL'&Length>5e5&Length<=5e6),
                               Del_GT5M=sum(Type=='DEL'&Length>5e6),
                               Dup_LT10K=sum(Type=='DUP'&Length<=1e4),
                               Dup_10Kto100K=sum(Type=='DUP'&Length>1e4&Length<=1e5),
                               Dup_100Kto500K=sum(Type=='DUP'&Length>1e5&Length<=5e5),
                               Dup_500Kto5M=sum(Type=='DUP'&Length>5e5&Length<=5e6),
                               Dup_GT5M=sum(Type=='DUP'&Length>5e6),
                               InvCount=sum(Type=='INV'),
                               BndCount=sum(Type=='BND'))
                   %>% arrange(-SvCount))

View(stdSigArmData)

stdSigSumData = (stdSigArmData %>% group_by(IsStressed)
                 %>% summarise(ArmCount=n(),
                               SvCount=sum(SvCount),
                               Sv_Avg=round(sum(SvCount)/n(),1),
                               DB_Count=sum(DB),
                               DB_Avg=round(sum(DB)/n(),1),
                               DB_Perc=round(sum(DB)/SvCount,3),
                               Del_LT10K=sum(Del_LT10K),
                               Del_LT10K_Avg=round(sum(Del_LT10K)/n(),1),
                               Del_LT10K_Perc=round(sum(Del_LT10K)/SvCount,3),
                               Del_10Kto100K=sum(Del_10Kto100K),
                               Del_10Kto100K_Avg=round(sum(Del_10Kto100K)/n(),1),
                               Del_10Kto100K_Perc=round(sum(Del_10Kto100K)/SvCount,3),
                               Del_100Kto500K=sum(Del_100Kto500K),
                               Del_100Kto500K_Avg=round(sum(Del_100Kto500K)/n(),1),
                               Del_100Kto500K_Perc=round(sum(Del_100Kto500K)/SvCount,3),
                               Del_500Kto5M=sum(Del_500Kto5M),
                               Del_500Kto5M_Avg=round(sum(Del_500Kto5M)/n(),1),
                               Del_500Kto5M_Perc=round(sum(Del_500Kto5M)/SvCount,3),
                               Del_GT5M=sum(Del_GT5M),
                               Del_GT5M_Avg=round(sum(Del_GT5M)/n(),1),
                               Del_GT5M_Perc=round(sum(Del_GT5M)/SvCount,3),
                               Dup_LT10K=sum(Dup_LT10K),
                               Dup_LT10K_Avg=round(sum(Dup_LT10K)/n(),1),
                               Dup_LT10K_Perc=round(sum(Dup_LT10K)/SvCount,3),
                               Dup_10Kto100K=sum(Dup_10Kto100K),
                               Dup_10Kto100K_Avg=round(sum(Dup_10Kto100K)/n(),1),
                               Dup_10Kto100K_Perc=round(sum(Dup_10Kto100K)/SvCount,3),
                               Dup_100Kto500K=sum(Dup_100Kto500K),
                               Dup_100Kto500K_Avg=round(sum(Dup_100Kto500K)/n(),1),
                               Dup_100Kto500K_Perc=round(sum(Dup_100Kto500K)/SvCount,3),
                               Dup_500Kto5M=sum(Dup_500Kto5M),
                               Dup_500Kto5M_Avg=round(sum(Dup_500Kto5M)/n(),1),
                               Dup_500Kto5M_Perc=round(sum(Dup_500Kto5M)/SvCount,3),
                               Dup_GT5M=sum(Dup_GT5M),
                               Dup_GT5M_Avg=round(sum(Dup_GT5M)/n(),1),
                               Dup_GT5M_Perc=round(sum(Dup_GT5M)/SvCount,3),
                               Inv_Count=sum(InvCount),
                               Inv_Avg=round(sum(InvCount)/n(),1),
                               Inv_Perc=round(sum(InvCount)/SvCount,3),
                               Bnd_Count=sum(BndCount),
                               Bnd_Avg=round(sum(BndCount)/n(),1),
                               Bnd_Perc=round(sum(BndCount)/SvCount,3))
                 %>% arrange(-SvCount))

View(stdSigSumData)





# stressed arms by chr and arm
# get sample arm data togethe

temp = combinedArmData %>% filter(IsStressed==1&IsLINE=='false')
View(temp)


stressedArmsBySample = (combinedArmData %>% filter(IsStressed==1&IsLINE=="false") %>% group_by(SampleId, Chr, Arm)
                %>% summarise(SvCount=n(),
                              DelCount=sum(Type=='DEL'),
                              DelPerc=round(sum(Type=='DEL')/n(),2),
                              DupCount=sum(Type=='DUP'),
                              DupPerc=round(sum(Type=='DUP')/n(),2),
                              InvCount=sum(Type=='INV'),
                              InvPerc=round(sum(Type=='INV')/n(),2),
                              BndCount=sum(Type=='BND'),
                              BndPerc=round(sum(Type=='BND')/n(),2),
                              LEPerc=round(sum(IsLINE=='true')/n(),3),
                              FSPerc=round(sum(IsFS=='true')/n(),3),
                              MaxCN=round(max((AdjCNStart+AdjCNEnd)*0.5),2),
                              AvgCN=round(sum((AdjCNStart+AdjCNEnd)*0.5)/n(),2),
                              AvgPloidy=round(sum(Ploidy)/n(),2),
                              AvgSvLength=round(sum(ifelse(Length>-1,Length,0))/sum(Length>-1),0),
                              DBPerc=round(sum(IsDB)/n(),2),
                              TIPerc=round(sum(IsTI)/n(),2),
                              AvgDBLen=round(sum(ifelse(IsDB==1,NearestDBLen,0))/sum(IsDB==1),0),
                              ClusteredPerc=round(sum(ClusterCount>1)/n(),2))
                %>% arrange(SampleId, Chr, Arm))

View(stressedArmsBySample)


stressedArms = (stressedArmsBySample %>% group_by(Chr, Arm)
             %>% summarise(SampleCount=n(),
                           SvCount=sum(SvCount),
                           DelCount=sum(DelCount),
                           DupCount=sum(DupCount),
                           InvCount=sum(InvCount),
                           BndCount=sum(BndCount))
             %>% arrange(Chr, Arm))

View(stressedArms)








