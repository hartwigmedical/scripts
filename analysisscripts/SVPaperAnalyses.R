library(data.table)
library(IRanges)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(stringi)
library(devtools)


# library(MutationalPatterns)
library(grid)
library(gridExtra)
library(cowplot)


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

sv_load_and_prepare<-function(filename)
{
  svData = read.csv(filename)
  sampleCancerTypes = read.csv('~/data/sample_cancer_types.csv')
  svData = merge(svData, sampleCancerTypes, by='SampleId', all.x=T)
  svData = svData %>% filter(ResolvedType!='LowQual')
  svData = sv_set_common_fields(svData)
  return (svData)  
}

clusters_load<-function(filename)
{
  clusters = read.csv(filename)
  sampleCancerTypes = read.csv('~/data/sample_cancer_types.csv')
  clusters = merge(clusters, sampleCancerTypes, by='SampleId', all.x=T)
  return (clusters)  
}

svData = sv_load_and_prepare('~/data/sv/CLUSTER.csv')
nrow(svData)
View(svData)
View(head(svData,1000))


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
# sampleClusterSummary$LowQual = (sampleClusterSummary$ResolvedType=='LowQual')
sampleClusterSummary$SimpleCluster = (sampleClusterSummary$ClusterCount>2&sampleClusterSummary$FoldbackCount==0&sampleClusterSummary$RepeatedChainLinkCount==0)
sampleClusterSummary$ComplexCluster = (sampleClusterSummary$ClusterCount>2&(sampleClusterSummary$FoldbackCount>0|sampleClusterSummary$RepeatedChainLinkCount>0))
sampleClusterSummary$UnresolvedSmallCluster = (sampleClusterSummary$ClusterCount<=2&sampleClusterSummary$IsResolved=='false')

sampleSummary = (sampleClusterSummary %>% group_by(SampleId)
                 %>% summarise(SvCount=sum(ClusterCount),
                               Clusters=n(),
                               SimpleSVClusters=sum(SimpleSVCluster),
                               SmallClusters=sum(ClusterCount>1&ClusterCount<=3),
                               LargeClusters=sum(ClusterCount>3),
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

totalSVCount = nrow(svData)
View(svData %>% group_by(ResolvedType,ClusterSize) 
     %>% summarise(Clusters=n_distinct(paste(SampleId,ClusterId,sep='_')), TotalSVs=n(), AsPerc=round(n()/totalSVCount,2))
     %>% arrange(-AsPerc))


# CLUSTER ANALYSIS
clusters = clusters_load('~/logs/SVA_CLUSTERS.csv')
nrow(clusters)

View(clusters %>% group_by(ResolvedType) %>% count())


# Simple chained clusters
simpleChainedClusters = clusters %>% filter(ResolvedType=='SimpleChain'|ResolvedType=='SimplePartialChain')
View(simpleChainedClusters)

View(clusters %>% filter(SampleId=='CPCT02020670TII'&ClusterId==3))
colnames(clusters)

View(clusters %>% filter(SampleId=='CPCT02020670TII'&ClusterId==3) 
     %>% select(SampleId,ClusterId,ClusterDesc,ResolvedType,FullyChained,ChainCount,Consistency,ArmCount,AssemblyLinks,ShortTIRemotes,Annotations))

View(clusters %>% filter(SampleId=='CPCT02020670TII'&ClusterId==3) 
     %>% select(Annotations))


simpleChainedClusters = sampleClusterSummary %>% filter(ResolvedType=='SimpleChain'|ResolvedType=='SimplePartialChain')
simpleChainedClusters$ClusterCountBucket = ifelse(simpleChainedClusters$ClusterCount<=5,simpleChainedClusters$ClusterCount,2**round(log(simpleChainedClusters$ClusterCount,2)))
View(simpleChainedClusters %>% group_by(ResolvedType,ClusterCountBucket) %>% count() %>% spread(ResolvedType,n))
View(simpleChainedClusters %>% filter(ResolvedType=='SimpleChain'&SglCount==0&NoneCount==0))


View(svData %>% filter(SampleId=='CPCT02020536T'&ClusterId==526))


nrow(sampleClusterSummary %>% filter(ResolvedType=='Line'&KnownLineCount<ClusterCount&SuspectLineCount==0&PolyAorTCount==0))



# plot counts of simple DELs, DUPs, INSs, synthetic DELs and DUPs, complex, simple and line clusters

sampleClusterSummary$ClusterType = ifelse(sampleClusterSummary$ResolvedType=='ComplexChain'|sampleClusterSummary$ResolvedType=='ComplexPartialChain','ComplexChain',
                                   ifelse(sampleClusterSummary$ResolvedType=='SimpleChain'|sampleClusterSummary$ResolvedType=='SimplePartialChain','SimpleChain',
                                   ifelse(sampleClusterSummary$ResolvedType=='DEL_Ext_TI'|sampleClusterSummary$ResolvedType=='DUP_Ext_TI'
                                          |grepl('Sgl',sampleClusterSummary$ResolvedType)|sampleClusterSummary$ResolvedType=='RecipTrans','SimpleSV',
                                   ifelse(sampleClusterSummary$ResolvedType=='DEL_Int_TI'|sampleClusterSummary$ResolvedType=='DUP_Int_TI','SimpleChain',
                                          as.character(sampleClusterSummary$ResolvedType)))))

sampleCounts = sampleClusterSummary %>% group_by(SampleId) %>% summarise(SampleCount=sum(ClusterCount))

View(sampleClusterSummary)

View(sampleClusterSummary %>% group_by(ClusterType) %>% summarise(Clusters=n(),SvCount=sum(ClusterCount)))

sampleClusterTypeData = sampleClusterSummary %>% group_by(SampleId,ClusterType) %>% summarise(SvCount=sum(ClusterCount)) %>% spread(ClusterType,SvCount)
sampleClusterTypeData[is.na(sampleClusterTypeData)] <- 0
View(sampleClusterTypeData)
gatherIndex = ncol(sampleClusterTypeData)
sampleClusterTypeData2 = gather(sampleClusterTypeData, "ClusterType", "SvCount", 2:gatherIndex)
View(sampleClusterTypeData2)

sampleClusterTypeData2 = merge(sampleClusterTypeData2,sampleCounts,by='SampleId',all.x=T)

# sampleClusterTypeData = sampleClusterSummary %>% group_by(SampleId,ClusterType) %>% summarise(SvCount=sum(ClusterCount),SampleCount=first(SampleCount))

clusterTypeColours = c("yellow", "blue", "green", "red", "orange", "purple", "pink", "brown", "darkgreen", "deepskyblue", "tan")

maxSamples = 50
maxRows = maxSamples * 5

sampleClusterTypePlot = (ggplot(data = head(sampleClusterTypeData2,maxRows), aes(x = reorder(SampleId, -SampleCount), y = SvCount, fill = ClusterType))
                + geom_bar(stat = "identity", colour = "black")
                + scale_fill_manual(values = clusterTypeColours)
                + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15))
                + ylab("SV Count") + xlab("Sample"))

print(sampleClusterTypePlot)


sampleSigPlot <- (ggplot(plotDataSet, aes(x = reorder(SampleId, -SampleCount), y = Count, fill = SigName))
                  + geom_bar(stat = "identity", colour = "black")
                  + labs(x = "", y = paste(varType, " Count by Sample", sep=''))
                  + scale_fill_manual(values = sigColours)
                  + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                  + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                  + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))



# Resolved Types












