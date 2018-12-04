library(data.table)
library(IRanges)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(stringi)
library(devtools)
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

svData =  read.csv('~/Dropbox/HMF Australia team folder/Structural Variant Analysis/CLUSTER_GRIDSS.csv', header = T, stringsAsFactors = F)#gridssCohortVariantsread.csv('~/data/sv/CLUSTER.csv')


##############################################################################
# simple annotations
svData = sv_set_common_fields(svData)

# Summary stats per cluster
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
                                      KnownLineCount=sum(LEStart=='Known'|LEEnd=='Known'|LEStart=='Ident'|LEEnd=='Ident'),
                                      SuspectLineCount=sum(LEStart=='Suspect'|LEEnd=='Suspect'),
                                      PolyAorTCount=sum(grepl('AAAAAAAA',InsertSeq)|grepl('TTTTTTTT',InsertSeq)),
                                      FragileSiteCount=sum(IsFS)
                        )
                        %>% arrange(SampleId,ClusterId))

sampleClusterSummary$SimpleSVCluster = (sampleClusterSummary$ClusterCount<=2&sampleClusterSummary$IsResolved=='true'&sampleClusterSummary$ResolvedType!='LowQual')
sampleClusterSummary$SimpleCluster = (sampleClusterSummary$ClusterCount>2&sampleClusterSummary$FoldbackCount==0&sampleClusterSummary$RepeatedChainLinkCount==0)
sampleClusterSummary$ComplexCluster = (sampleClusterSummary$ClusterCount>2&(sampleClusterSummary$FoldbackCount>0|sampleClusterSummary$RepeatedChainLinkCount>0))
sampleClusterSummary$UnresolvedSmallCluster = (sampleClusterSummary$ClusterCount<=2&sampleClusterSummary$IsResolved=='false')

# Stats per sample
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

sampleClusterSummary$ClusterType = ifelse(sampleClusterSummary$ResolvedType=='ComplexChain'|sampleClusterSummary$ResolvedType=='ComplexPartialChain','ComplexChain',
                                   ifelse(sampleClusterSummary$ResolvedType=='SimpleChain'|sampleClusterSummary$ResolvedType=='SimplePartialChain','SimpleChain',
                                   ifelse(sampleClusterSummary$ResolvedType=='DEL_Ext_TI'|sampleClusterSummary$ResolvedType=='DUP_Ext_TI'
                                          |grepl('Sgl',sampleClusterSummary$ResolvedType)|sampleClusterSummary$ResolvedType=='RecipTrans','SimpleSV',
                                   ifelse(sampleClusterSummary$ResolvedType=='DEL_Int_TI'|sampleClusterSummary$ResolvedType=='DUP_Int_TI','SimpleChain',
                                          as.character(sampleClusterSummary$ResolvedType)))))

sampleCounts = sampleClusterSummary %>% group_by(SampleId) %>% summarise(SampleCount=sum(ClusterCount))

######################################
########### ANALYSES #################
######################################

#1. Counts by cluster type 
View(sampleClusterSummary %>% group_by(ClusterType) %>% summarise(Clusters=n(),SvCount=sum(ClusterCount)))


#2. Basic signature type plot
sampleClusterTypeData = sampleClusterSummary %>% group_by(SampleId,ClusterType) %>% summarise(SvCount=sum(ClusterCount)) %>% spread(ClusterType,SvCount,fill=0)
sampleClusterTypeData2 = gather(sampleClusterTypeData, "ClusterType", "SvCount", 2:ncol(sampleClusterTypeData))
sampleClusterTypeData2 = merge(sampleClusterTypeData2,sampleCounts,by='SampleId',all.x=T)

clusterTypeColours = c("yellow", "blue", "green", "red", "orange", "purple", "pink", "brown", "darkgreen", "deepskyblue", "tan")
maxSamples = 50
maxRows = maxSamples * 5

print(ggplot(data = head(sampleClusterTypeData2,maxRows), aes(x = reorder(SampleId, -SampleCount), y = SvCount, fill = ClusterType))
                + geom_bar(stat = "identity", colour = "black")
                + scale_fill_manual(values = clusterTypeColours)
                + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15))
                + ylab("SV Count") + xlab("Sample"))

####################################










