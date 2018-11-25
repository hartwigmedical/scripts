library(data.table)
library(IRanges)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(stringi)
library(devtools)

sv_set_common_fields<-function(svData)
{
  svData$IsLINE = ifelse(svData$LEStart!='None'|svData$LEEnd!='None',T,F)
  svData$IsFS = ifelse(svData$FSStart!='false'|svData$FSEnd!='false',T,F)
  svData$Length = ifelse(as.character(svData$ChrStart)!=as.character(svData$ChrEnd)|svData$Type=='INS'|svData$ArmEnd!=svData$ArmStart, -1, svData$PosEnd-svData$PosStart)
  svData$DoubleDupBE = ifelse(svData$DupBEStart=='true'&svData$DupBEEnd=='true',T,F)
  svData$SingleDupBE = ifelse(svData$DoubleDupBE==0&(svData$DupBEStart=='true'|svData$DupBEEnd=='true'),T,F)
  svData$TICount = ifelse(svData$LnkTypeStart=='TI',0.5,0)+ifelse(svData$LnkTypeEnd=='TI',0.5,0)
  svData$DBCount = ifelse(svData$DBLenStart>=0,0.5,0)+ifelse(svData$DBLenEnd>=0,0.5,0)
  svData$IsSglTI = ifelse(svData$LnkTypeStart=='SGL',0.5,0)
  svData$AsmbTICount = ifelse(svData$AsmbMatchStart=='MATCH',0.5,0)+ifelse(svData$AsmbMatchEnd=='MATCH',0.5,0)
  svData$InferTICount = svData$TICount - svData$AsmbTICount
  svData$ShortTICount=ifelse(svData$LnkTypeStart=='TI'&svData$LnkLenStart<=1000,0.5,0)+ifelse(svData$LnkTypeEnd=='TI'&svData$LnkLenEnd<=1000,0.5,0)
  svData$ClusterSize = ifelse(svData$ClusterCount==1,'Single',ifelse(svData$ClusterCount<=4,'Small','Large'))
  svData$IsConsistent = ifelse(svData$Consistency==0,T,F)
  svData$IsChained = (svData$ChainCount>=1)
  svData$FoldbackCount = ifelse(svData$FoldbackLenStart>=0,0.5,0)+ifelse(svData$FoldbackLenEnd>=0,0.5,0)
  svData$RepeatedChainLink = (svData$ChainCount>0 & grepl(';',svData$ChainIndex))
  return (svData)
}


svData = read.csv('~/data/sv/CLUSTER.csv')
# svData = read.csv('~/logs/CLUSTER_GRIDSS.csv')
# svData = read.csv('~/logs/COLO829T.csv')
nrow(svData)
View(svData)
View(head(svData,1000))

sampleCancerTypes = read.csv('~/data/sample_cancer_types.csv')
# View(sampleCancerTypes)
svData = merge(svData, sampleCancerTypes, by='SampleId', all.x=T)

# filter out Low-Quality SVs
svData = svData %>% filter(ResolvedType!='LowQual')
nrow(svData)

# simple annotations
svData = sv_set_common_fields(svData)


# OVERALL SAMPLE STATS
sampleClusterSummary = (svData %>% group_by(SampleId,ClusterId)
                        %>% summarise(ClusterCount=first(ClusterCount),
                                      IsResolved=first(IsResolved),
                                      ResolvedType=first(ResolvedType),
                                      ArmCount=first(ArmCount),
                                      DBCount=sum(DBCount),
                                      TICount=sum(TICount),
                                      ShortTICount=sum(ShortTICount),
                                      AsmbLinkCount=sum(AsmbTICount),
                                      InferLinkCount=sum(InferTICount),
                                      ChainedCount=sum(IsChained),
                                      RepeatedChainLinkCount=sum(RepeatedChainLink),
                                      FoldbackCount=sum(FoldbackCount),
                                      Consistent=first(IsConsistent),
                                      DupBECount=sum(DoubleDupBE),
                                      DelCount=sum(Type=='DEL'),
                                      DupCount=sum(Type=='DUP'),
                                      InsCount=sum(Type=='INS'),
                                      InvCount=sum(Type=='INV'),
                                      BndCount=sum(Type=='BND'),
                                      SglCount=sum(Type=='SGL'),
                                      NoneCount=sum(Type=='NONE'),
                                      LineCount=sum(IsLINE),
                                      FragileSiteCount=sum(IsFS)
                        )
                        %>% arrange(SampleId,ClusterId))

View(sampleClusterSummary)

sampleClusterSummary$SimpleSVCluster = (sampleClusterSummary$ClusterCount<=2&sampleClusterSummary$IsResolved=='true'&sampleClusterSummary$ResolvedType!='LowQual')
sampleClusterSummary$LowQual = (sampleClusterSummary$ResolvedType=='LowQual')
sampleClusterSummary$SimpleCluster = (sampleClusterSummary$ClusterCount>2&sampleClusterSummary$FoldbackCount==0&sampleClusterSummary$RepeatedChainLinkCount==0)
sampleClusterSummary$ComplexCluster = (sampleClusterSummary$ClusterCount>2&(sampleClusterSummary$FoldbackCount>0&sampleClusterSummary$RepeatedChainLinkCount>0))
sampleClusterSummary$UnresolvedSmallCluster = (sampleClusterSummary$ClusterCount<=2&sampleClusterSummary$IsResolved=='false')
  
sampleSummary = (sampleClusterSummary %>% group_by(SampleId)
                 %>% summarise(SvCount=sum(ClusterCount),
                               Clusters=n(),
                               SimpleSVClusters=sum(SimpleSVCluster),
                               SmallClusters=sum(ClusterCount>1&ClusterCount<=3),
                               LargeClusters=sum(ClusterCount>3),
                               LowQuals=sum(LowQual),
                               SimpleClusters=sum(SimpleCluster),
                               ComplexClusters=sum(ComplexCluster),
                               UnresolvedSmallClusters=sum(UnresolvedSmallCluster),
                               # DBCount=sum(DBCount),
                               # TICount=sum(TICount),
                               # ShortTICount=sum(ShortTICount),
                               # AsmbLinkCount=sum(AsmbLinkCount),
                               # InferLinkCount=sum(InferLinkCount),
                               DupBECount=sum(DupBECount),
                               DelCount=sum(DelCount),
                               DupCount=sum(DupCount),
                               InsCount=sum(InsCount),
                               InvCount=sum(InvCount),
                               BndCount=sum(BndCount),
                               SglCount=sum(SglCount),
                               NoneCount=sum(NoneCount),
                               ChainedCount=sum(ChainedCount),
                               RepeatedChainLinkCount=sum(RepeatedChainLinkCount),
                               FoldbackCount=sum(FoldbackCount),
                               LineCount=sum(LineCount),
                               FragileSiteCount=sum(FragileSiteCount)
                 )
                 %>% arrange(SampleId))

View(sampleSummary)


# Resolved Types
totalSVCount = nrow(svData)
View(svData %>% group_by(ResolvedType,ClusterSize) 
     %>% summarise(Clusters=n_distinct(paste(SampleId,ClusterId,sep='_')), TotalSVs=n(), AsPerc=round(n()/totalSVCount,2))
     %>% arrange(-AsPerc))









