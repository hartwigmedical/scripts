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


######################
# SV Signatures
######################

svData = sv_load_and_prepare('~/data/sv/CLUSTER.csv')


# first collect simple, unclustered DELs and DUPs

simpleDels = svData %>% filter(ResolvedType=='SimpleSV'&Type=='DEL')
nrow(simpleDels)
simpleDels$Bucket = ifelse(simpleDels$Length<=1e2,'DEL_100',ifelse(simpleDels$Length<=2e3,'DEL_1K','DEL_LONG'))

simpleDups = svData %>% filter(ResolvedType=='SimpleSV'&Type=='DUP')
nrow(simpleDups)
simpleDups$Bucket = ifelse(simpleDups$Length<=1e2,'DUP_100',ifelse(simpleDups$Length<=2e4,'DUP_20K',ifelse(simpleDups$Length<=2e5,'DUP_200K','DUP_LONG')))

inserts = svData %>% filter(ResolvedType=='SimpleSV'&Type=='INS')
inserts$Bucket = 'INS'

# synthetic DELs and DUPs
delExtTIs = svData %>% filter(ResolvedType=='DEL_Ext_TI'&ChainIndex=='0s')
nrow(delExtTIs)
delExtTIs$Bucket = ifelse(delExtTIs$SynDelDupLen<=1e2,'DEL_SE_100',ifelse(delExtTIs$SynDelDupLen<=2e3,'DEL_SE_1K','DEL_SE_LONG'))

#recipInvs = svData %>% filter(ResolvedType=='DEL_Int_TI'&ChainIndex=='0s'&ClusterDesc=='INV=2'&SynDelDupTILen>=0.99*SynDelDupLen)
#nrow(recipInvs)
#delIntTIs = svData %>% filter(ResolvedType=='DEL_Int_TI'&ChainIndex=='0s'&!(ClusterDesc=='INV=2'&SynDelDupTILen>=0.99*SynDelDupLen))
delIntTIs = svData %>% filter(ResolvedType=='DEL_Int_TI'&ChainIndex=='0s')
View(delIntTIs)
nrow(delIntTIs)

delIntTIs$Bucket = ifelse(delIntTIs$SynDelDupLen<=1e2,'DEL_SI_100',ifelse(delIntTIs$SynDelDupLen<=2e3,'DEL_SI_1K','DEL_SI_LONG'))

delIntTIs$DelDupType = ifelse(delIntTIs$ClusterDesc=='INV=2'&delIntTIs$SynDelDupTILen >= 0.99*delIntTIs$SynDelDupLen,'RecipInv',as.character(delIntTIs$ResolvedType))
delIntTIs = svData %>% filter(ResolvedType=='DEL_Int_TI'&ChainCount>0)

dupExtTIs = svData %>% filter(ResolvedType=='DUP_Ext_TI'&ChainIndex=='0s')
dupExtTIs$Bucket = ifelse(dupExtTIs$SynDelDupLen<=1e2,'DUP_SE_100',ifelse(dupExtTIs$SynDelDupLen<=2e4,'DUP_SE_20K',ifelse(dupExtTIs$SynDelDupLen<=2e5,'DUP_SE_200K','DUP_SE_LONG')))

# balance and unbalanced translocations
recipTrans = svData %>% filter(ResolvedType=='RecipTrans')
recipTrans$Bucket = 'RECIP_TRANS'

unbalancedTrans = svData %>% filter(ClusterCount==1&Type=='BND'&IsConsistent)
unbalancedTrans$Bucket = 'UNBAL_TRANS'
View(unbalancedTrans)

# line elements no further splits for now
lineSVs = svData %>% filter(ResolvedType=='Line')
lineSVs$Bucket = 'Line'

# synthetic simple types from SGL pairs
sglPairSimples = (svData %>% filter(ResolvedType=='SglPair_INS'|ResolvedType=='SglPair_DEL'|ResolvedType=='SglPair_DUP')
                  %>% group_by(SampleId,ClusterId,ResolvedType) %>% summarise(Length=abs(first(PosStart) - last(PosStart))))

sglPairSimples$Bucket

View(sglPairSimples)


clusters = clusters_load("~/data/sv/SVA_CLUSTERS.csv")
View(clusters)

# simple chains split by their DB count
simpleChainSVs = svData %>% filter(ResolvedType=='SimpleChain')
View(simpleChainSVs)






nrow(svData %>% filter(IsFS))
nrow(svData %>% filter(FSStart!='false'|FSEnd!='false'))

View(svData)


                           
                           
                           
                           
                           
                           
sampleBucketData = (svData %>% filter(ResolvedType=='Simple') %>% group_by(SampleId)
                        %>% summarise(ClusterCount=first(ClusterCount),
                                      IsResolved=first(IsResolved),
                
sampleBucketData = 


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

View(clusters %>% filter(SampleId=='CPCT02050052T'))
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

print(sampleClusterTypePlot)


sampleSigPlot <- (ggplot(plotDataSet, aes(x = reorder(SampleId, -SampleCount), y = Count, fill = SigName))
                  + geom_bar(stat = "identity", colour = "black")
                  + labs(x = "", y = paste(varType, " Count by Sample", sep=''))
                  + scale_fill_manual(values = sigColours)
                  + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                  + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                  + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))





