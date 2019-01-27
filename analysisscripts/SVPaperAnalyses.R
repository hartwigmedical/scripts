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

# svData = sv_set_common_fields(svData)


DEL_LENGTH_1 = 1e2
DEL_LENGTH_2 = 2e3
DEL_LENGTH_3 = 1e5 # not required
DEL_LENGTH_MAX = 1e8

DUP_LENGTH_1 = 3e2
DUP_LENGTH_2 = 5e4
DUP_LENGTH_3 = 5e5
DUP_LENGTH_4 = 3e6 # not required

length_to_bucket<-function(type, length)
{
  lenBucket = ifelse(length<=DEL_LENGTH_1,DEL_LENGTH_1,ifelse(length>=DEL_LENGTH_MAX,DEL_LENGTH_MAX,10**round(log(length,10))))
  bucket = paste(type, format(lenBucket,scientific=T), sep='_')
  bucketShort = stri_replace_all_fixed(bucket, '+0', '')
  return (bucketShort)
}

# first collect simple, unclustered DELs and DUPs

simpleDels = svData %>% filter(ResolvedType=='SimpleSV'&Type=='DEL')

# filter us suspect DEL with a matching insert sequence
#View(simpleDels %>% filter(stri_length(InsertSeq) == Length - 1) %>% group_by(Length) %>% count())
simpleDels = simpleDels %>% filter(stri_length(InsertSeq) != Length - 1)

#simpleDels$Bucket = ifelse(simpleDels$Length<=DEL_LENGTH_1,'DEL_100',ifelse(simpleDels$Length<=DEL_LENGTH_2,'DEL_1K','DEL_LONG'))

simpleDels$Bucket = apply(simpleDels[,c('Length'),drop=F], 1, function(x) length_to_bucket('DEL',x[1]))

View(simpleDels %>% group_by(Bucket) %>% count())
View(simpleDels %>% group_by(Bucket1) %>% count())
View(simpleDels %>% group_by(LenBucket) %>% count())

simpleDups = svData %>% filter(ResolvedType=='SimpleSV'&Type=='DUP')
simpleDups$Bucket = ifelse(simpleDups$Length<=DUP_LENGTH_1,'DUP_100',ifelse(simpleDups$Length<=DUP_LENGTH_2,'DUP_20K',ifelse(simpleDups$Length<=DUP_LENGTH_3,'DUP_200K','DUP_LONG')))
simpleDups$Bucket = apply(simpleDups[,c('Length'),drop=F], 1, function(x) length_to_bucket('DUP',x[1]))
View(simpleDups %>% group_by(Bucket) %>% count())

inserts = svData %>% filter(ResolvedType=='SimpleSV'&Type=='INS')
nrow(inserts)
inserts$Bucket = 'INS'

# synthetic DELs and DUPs
delExtTIs = svData %>% filter(ResolvedType=='DEL_Ext_TI'&ChainIndex=='0s')
delExtTIs$Bucket = ifelse(delExtTIs$SynDelDupLen<=DEL_LENGTH_1,'DEL_SE_100',ifelse(delExtTIs$SynDelDupLen<=DEL_LENGTH_2,'DEL_SE_1K','DEL_SE_LONG'))
delExtTIs$Bucket = apply(delExtTIs[,c('SynDelDupLen'),drop=F], 1, function(x) length_to_bucket('DEL_SE',x[1]))
View(delExtTIs %>% group_by(Bucket) %>% count())

#recipInvs = svData %>% filter(ResolvedType=='DEL_Int_TI'&ChainIndex=='0s'&ClusterDesc=='INV=2'&SynDelDupTILen>=0.99*SynDelDupLen)
#nrow(recipInvs)
#delIntTIs = svData %>% filter(ResolvedType=='DEL_Int_TI'&ChainIndex=='0s'&!(ClusterDesc=='INV=2'&SynDelDupTILen>=0.99*SynDelDupLen))
delIntTIs = svData %>% filter(ResolvedType=='DEL_Int_TI'&ChainIndex=='0s')
delIntTIs$Bucket = ifelse(delIntTIs$SynDelDupLen<=DEL_LENGTH_1,'DEL_SI_100',ifelse(delIntTIs$SynDelDupLen<=DEL_LENGTH_2,'DEL_SI_1K','DEL_SI_LONG'))
delIntTIs$Bucket = apply(delIntTIs[,c('SynDelDupLen'),drop=F], 1, function(x) length_to_bucket('DEL_SI',x[1]))
View(delIntTIs %>% group_by(Bucket) %>% count())

#delIntTIs$DelDupType = ifelse(delIntTIs$ClusterDesc=='INV=2'&delIntTIs$SynDelDupTILen >= 0.99*delIntTIs$SynDelDupLen,'RecipInv',as.character(delIntTIs$ResolvedType))

dupExtTIs = svData %>% filter(ResolvedType=='DUP_Ext_TI'&ChainIndex=='0s')
dupExtTIs$Bucket = ifelse(dupExtTIs$SynDelDupLen<=DUP_LENGTH_1,'DUP_SE_100',ifelse(dupExtTIs$SynDelDupLen<=DUP_LENGTH_2,'DUP_SE_20K',ifelse(dupExtTIs$SynDelDupLen<=DUP_LENGTH_3,'DUP_SE_200K','DUP_SE_LONG')))
dupExtTIs$Bucket = apply(dupExtTIs[,c('SynDelDupLen'),drop=F], 1, function(x) length_to_bucket('DUP_SE',x[1]))
View(dupExtTIs %>% group_by(Bucket) %>% count())

# balance and unbalanced translocations
recipTrans = svData %>% filter(ResolvedType=='RecipTrans')
recipTrans$Bucket = 'RECIP_TRANS'
nrow(recipTrans)

unbalancedTrans = svData %>% filter(ClusterCount==1&Type=='BND'&IsConsistent)
unbalancedTrans$Bucket = 'UNBAL_TRANS'
nrow(unbalancedTrans)

# line elements no further splits for now
lineSVs = svData %>% filter(ResolvedType=='Line')
lineSVs$Bucket = 'Line'
nrow(lineSVs)

# synthetic simple types from SGL pairs
sglPairSimples = (svData %>% filter(ResolvedType=='SglPair_INS'|ResolvedType=='SglPair_DEL'|ResolvedType=='SglPair_DUP')
                  %>% group_by(SampleId,ClusterId,ResolvedType) %>% summarise(Length=abs(first(PosStart) - last(PosStart))))
View(sglPairSimples)

sglPairSimples$Bucket = ifelse(sglPairSimples$ResolvedType=='SglPair_DEL',
                               ifelse(sglPairSimples$Length<=DEL_LENGTH_1,'DEL_SG_100',ifelse(sglPairSimples$Length<=DEL_LENGTH_2,'DEL_SG_1K','DEL_SG_LONG')),
                        ifelse(sglPairSimples$ResolvedType=='SglPair_DUP',
                               ifelse(sglPairSimples$Length<=DUP_LENGTH_1,'DUP_SG_100',ifelse(sglPairSimples$Length<=DUP_LENGTH_2,'DUP_SG_20K',ifelse(sglPairSimples$Length<=DUP_LENGTH_3,'DUP_SG_200K','DUP_SG_LONG'))),
                               'INS_SG'))

sglPairDels = sglPairSimples %>% filter(ResolvedType=='SglPair_DEL')
sglPairDels$Bucket = apply(sglPairDels[,c('Length'),drop=F], 1, function(x) length_to_bucket('DEL_SG',x[1]))
View(sglPairDels %>% group_by(Bucket) %>% count())

sglPairDups = sglPairSimples %>% filter(ResolvedType=='SglPair_DUP')
sglPairDups$Bucket = apply(sglPairDups[,c('Length'),drop=F], 1, function(x) length_to_bucket('DUP_SG',x[1]))
View(sglPairDups %>% group_by(Bucket) %>% count())

sglPairDels = sglPairSimples %>% filter(ResolvedType=='SglPair_DEL')
sglPairDels$Bucket = apply(sglPairDels[,c('Length'),drop=F], 1, function(x) length_to_bucket('DEL_SG',x[1]))


View(sglPairSimples %>% group_by(Bucket) %>% count())



clusters = clusters_load("~/data/sv/SVA_CLUSTERS.csv")
View(clusters)

# simple chains split by their DB count

# what is likely to be significant for chromothriptic like events - severity?
# number of breaks
# length of chain(s)
# number of origin arm

simpleChains = clusters %>% filter(ResolvedType=='SimpleChain'&HasReplicated=='false')
simpleChains$Bucket = ifelse(simpleChains$DSBs==0,'SIM_CHAIN_DSB_NONE',ifelse(simpleChains$DSBs<=2,'SIM_CHAIN_DB_LOW',ifelse(simpleChains$DSBs<=4,'SIM_CHAIN_DB_MED','SIM_CHAIN_DB_HIGH')))
#View(simpleChains)
View(simpleChains %>% group_by(Bucket) %>% count())

View(clusters %>% filter(BndCount==1&DelCount>10))
View(clusters %>% filter(ClusterCount>1&DSBs==ClusterCount))

View(simpleChains %>% filter(DSBs==0))
View(simpleChains %>% group_by(DSBs) %>% count())


View(simpleChains %>% group_by(DSBs) %>% count())
View(simpleChains %>% group_by(ClusterCount) %>% count())
View(simpleChains %>% group_by(FullyChained) %>% count())

View(complexChains %>% filter(ResolvedType=='SimpleChain'&HasReplicated=='true'))

# other complex chains and clusters

complexChains = clusters %>% filter((ResolvedType=='SimpleChain'&HasReplicated=='true')|ResolvedType=='ComplexChain'|(ClusterCount>=2&ResolvedType=='None'))
View(complexChains %>% group_by(ResolvedType) %>% count())
complexChains$FoldbackGroup = ifelse(complexChains$Foldbacks==0,'FB=0',ifelse(complexChains$Foldbacks==1,'FB=1',ifelse(complexChains$Foldbacks<=3,'FB=2-3','FB=HIGH')))
View(complexChains %>% group_by(FoldbackGroup) %>% count())
complexChains$CopyNumberGroup = ifelse(complexChains$MaxCopyNumber<=1,'CN=NORM',ifelse(complexChains$MaxCopyNumber<=4,'CN=MED','CN=HIGH'))
View(complexChains %>% group_by(CopyNumberGroup) %>% count())
complexChains$ClusterSize = ifelse(complexChains$ClusterCount<=3,'SMALL','LARGE')
View(complexChains %>% group_by(ClusterSize) %>% count())

complexChains$Bucket = paste(complexChains$ClusterSize,complexChains$FoldbackGroup,complexChains$CopyNumberGroup,sep='_')
View(complexChains %>% group_by(Bucket) %>% count())


# now combine all groups together - just SampleId, CancerType and Count
get_sample_count_fields<-function(data)
{
  return (data %>% group_by(SampleId,Bucket) %>% summarise(Count=n()))
}

svSampleCounts = get_sample_count_fields(simpleDels)
svSampleCounts = rbind(svSampleCounts, get_sample_count_fields(simpleDups))
svSampleCounts = rbind(svSampleCounts, get_sample_count_fields(inserts))
svSampleCounts = rbind(svSampleCounts, get_sample_count_fields(delExtTIs))
svSampleCounts = rbind(svSampleCounts, get_sample_count_fields(delIntTIs))
svSampleCounts = rbind(svSampleCounts, get_sample_count_fields(dupExtTIs))
svSampleCounts = rbind(svSampleCounts, get_sample_count_fields(recipTrans))
svSampleCounts = rbind(svSampleCounts, get_sample_count_fields(unbalancedTrans))
svSampleCounts = rbind(svSampleCounts, get_sample_count_fields(lineSVs))
svSampleCounts = rbind(svSampleCounts, get_sample_count_fields(sglPairSimples))
svSampleCounts = rbind(svSampleCounts, get_sample_count_fields(simpleChains))
svSampleCounts = rbind(svSampleCounts, get_sample_count_fields(complexChains))

# only DELs and DUPs
svSampleCounts = get_sample_count_fields(simpleDels)
svSampleCounts = rbind(svSampleCounts, get_sample_count_fields(simpleDups))
svSampleCounts = rbind(svSampleCounts, get_sample_count_fields(delExtTIs))
svSampleCounts = rbind(svSampleCounts, get_sample_count_fields(delIntTIs))
svSampleCounts = rbind(svSampleCounts, get_sample_count_fields(dupExtTIs))
svSampleCounts = rbind(svSampleCounts, get_sample_count_fields(sglPairDels))
svSampleCounts = rbind(svSampleCounts, get_sample_count_fields(sglPairDups))

# merge synthetic DELs and DUPs in with standard
svSampleCounts$Bucket = stri_replace_all_fixed(svSampleCounts$Bucket, '_SI', '')
svSampleCounts$Bucket = stri_replace_all_fixed(svSampleCounts$Bucket, '_SE', '')
svSampleCounts$Bucket = stri_replace_all_fixed(svSampleCounts$Bucket, '_SG', '')
svSampleCounts = svSampleCounts %>% group_by(SampleId,Bucket) %>% summarise(Count=sum(Count))

View(svSampleCounts %>% group_by(Bucket) %>% count())
View(svSampleCounts)

View(svSampleCounts %>% group_by(SampleId,Bucket) %>% count())

# transform into signature matrix
svMatrixData = svSampleCounts %>% spread(SampleId,Count)
svMatrixData[is.na(svMatrixData)] = 0
svBucketNames = data.frame(svMatrixData$Bucket)
colnames(svBucketNames) <- c("Bucket")
View(svBucketNames)
  
svMatrixData = within(svMatrixData, rm(Bucket))

write.csv(svMatrixData, file="~/logs/r_output/sv_sig_matrix_data.csv", row.names=F, quote=F)
write.csv(svSampleCounts, file="~/logs/r_output/sv_sig_sample_counts.csv", row.names=F, quote=F)
write.csv(svBucketNames, file="~/logs/r_output/sv_bucket_names.csv", row.names=F, quote=F)


svSampleTotals = svSampleCounts %>% group_by(SampleId) %>% summarise(SampleTotal=sum(Count)) %>% arrange(-SampleTotal)
View(svSampleTotals)

View(svSampleTotals %>% group_by(SampleId,CancerType) %>% arrange(CancerType,-SampleTotal))

sampleCancerTypes = read.csv('~/data/sample_cancer_types.csv')
svSampleTotals = merge(svSampleTotals, sampleCancerTypes, by='SampleId', all.x=T)




######################
# DRIVER GENE ANALYSIS
######################

svData = svData %>% separate(DriverStart, c('DriverTypeStart','DriverGeneStart','DriverInfoStart'), sep = ';')
svData = svData %>% separate(DriverEnd, c('DriverTypeEnd','DriverGeneEnd','DriverInfoEnd'), sep = ';')

driverGeneSvs = (svData %>% filter(DriverTypeStart!=''|DriverTypeEnd!='') 
     %>% select(SampleId,ClusterId,ClusterDesc,ResolvedType,Id,Type,DriverTypeStart,DriverGeneStart,DriverInfoStart,DriverTypeEnd,DriverGeneEnd,DriverInfoEnd,
                ChrStart,PosStart,OrientStart,ChrEnd,PosEnd,OrientEnd,Ploidy,AdjCNStart,AdjCNChgStart,AdjCNEnd,AdjCNChgEnd))

View(driverGeneSvs)

write.csv(driverGeneSvs, "~/logs/driver_gene_svs.csv", quote = F, row.names = F)

driverGeneStart = driverGeneSvs %>% filter(DriverTypeStart!='')
driverGeneStart$DriverGeneType = driverGeneStart$DriverTypeStart
driverGeneStart$DriverGene = driverGeneStart$DriverGeneStart
driverGeneStart$DriverGeneInfo = driverGeneStart$DriverInfoStart

driverGeneEnd = driverGeneSvs %>% filter(DriverTypeStart==''&DriverTypeEnd!='')
driverGeneEnd$DriverGeneType = driverGeneEnd$DriverTypeEnd
driverGeneEnd$DriverGene = driverGeneEnd$DriverGeneEnd
driverGeneEnd$DriverGeneInfo = driverGeneEnd$DriverInfoEnd
driverGenes = rbind(driverGeneStart,driverGeneEnd)

View(driverGenes)
View(driverGenes %>% filter(DriverGeneType=='AMP'&DriverGeneInfo=='SV') %>% group_by(DriverGeneType,DriverGene,DriverGeneInfo) %>% count())
View(driverGenes %>% filter(DriverGeneType=='AMP'&DriverGeneInfo=='SV') %>% group_by(Type) %>% count())

View(driverGenes %>% group_by(DriverGeneType,DriverGene,DriverGeneInfo) %>% count())


View(svData %>% group_by(DriverTypeStart,DriverGeneStart,DriverInfoStart) %>% count())



###############
# Genic Regions
###############

View(tiDirectData)
nrow(tiDirectData)

tiDirectData$GeneStart = as.character(tiDirectData$GeneStart)
tiDirectData$GeneEnd = as.character(tiDirectData$GeneEnd)
tiDirectData$BothEndsGenic = (tiDirectData$GeneStart!=''&tiDirectData$GeneEnd!='')
View(tiDirectData %>% filter(BothEndsGenic))
nrow(tiDirectData %>% filter(BothEndsGenic&GeneStart==GeneEnd&TILength<300&NextSVDistance>5e3))
nrow(tiDirectData %>% filter(TILength<300&NextSVDistance>5e3))

shortTIData = svData %>% filter(ShortTICount>0)
nrow(shortTIData)

StoSLinks = merge(shortTIData %>% filter(shortTIData$LnkLenStart>0),shortTIData %>% filter(shortTIData$LnkLenStart>0),by.x=c('Id','LnkSvStart'),by.y=c('LnkSvStart','Id'),all.x=T)

View(StoSLinks %>% select(Id,LnkSvStart,SampleId.x,Type.x,ChrStart.x,PosStart.x,OrientStart.x,LnkLenStart.x,GeneStart.x,
                          Type.x,ChrStart.y,PosStart.y,OrientStart.y,LnkLenStart.y,GeneStart.y))

View(StoSLinks)

endTILinks = merge(shortTIData,shortTIData %>% filter(shortTIData$LnkLenEnd>0),by.x='Id',by.y='LnkSvEnd',all.x=T)

View(startTILinks)
colnames(startTILinks)
tmpColnames = colnames(startTILinks)
tmpColnames[1] = 'Id1'
colnames(startTILinks) = tmpColnames
print(tmpColnames)

View(startTILinks %>% select(SampleId.x,Id1,Id,ChrStart.x,PosStart.x,LnkLenStart.x,GeneStart.x,GeneStart.y,GeneEnd.x,GeneEnd.y))


shortTIData$LinkIdStart = ifelse(shortTIData$LnkLenStart>0,ifelse(shortTIData$Id<shortTIData$LnkSvStart,
                                paste(shortTIData$Id,shortTIData$LnkSvStart,sep='_'),paste(shortTIData$LnkSvStart,shortTIData$Id,sep='_')),'')

shortTIData$LinkIdEnd = ifelse(shortTIData$LnkLenEnd>0,ifelse(shortTIData$Id<shortTIData$LnkSvEnd,
                                                                  paste(shortTIData$Id,shortTIData$LnkSvEnd,sep='_'),paste(shortTIData$LnkSvEnd,shortTIData$Id,sep='_')),'')

startTIs = shortTIData %>% filter(LinkIdStart!='') %>% group_by(LinkIdStart) %>% count()
View(startTIs)

shortTIData$ShortTIGenic = ((shortTIData$LnkLenStart>0&shortTIData$LnkLenStart<=1000&shortTIData$GeneStart!="")
                          | (shortTIData$LnkLenEnd>0&shortTIData$LnkLenEnd<=1000&shortTIData$GeneEnd!=""))

View(shortTIData %>% group_by(ShortTIGenic) %>% count())

View(shortTIData %>% group_by(ShortTIGenic) %>% count())


svData$ShortTICount=ifelse(svData$LnkLenStart>0&svData$LnkLenStart<=1000,0.5,0)+ifelse(svData$LnkLenEnd>0&svData$LnkLenEnd<=1000,0.5,0)





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
totalClusterCount = nrow(svData %>% group_by(SampleId,ClusterId) %>% count())

View(svData %>% group_by(SampleId,ClusterId,ResolvedType) 
     %>% summarise(SvCount=n()) %>% group_by(ResolvedType)
     %>% summarise(Clusters=n(), ClusterPerc=round(n()/totalClusterCount,2), SvTotal=sum(SvCount),AsPerc=round(sum(SvCount)/totalSVCount,2))
     %>% arrange(-AsPerc))

View(svData %>% group_by(ResolvedType,ClusterSize) 
     %>% summarise(Clusters=n_distinct(paste(SampleId,ClusterId,sep='_')), AsPerc=round(n()/totalSVCount,2))
     %>% arrange(-AsPerc))


svData484 = svData = sv_load_and_prepare('~/data/sv/CLUSTER_003.csv')

totalSVCount484 = nrow(svData)
totalClusterCount484 = nrow(svData484 %>% group_by(SampleId,ClusterId) %>% count())

View(svData484 %>% group_by(SampleId,ClusterId,ResolvedType) 
     %>% summarise(SvCount=n()) %>% group_by(ResolvedType)
     %>% summarise(Clusters=n(), ClusterPerc=round(n()/totalClusterCount484,2), SvTotal=sum(SvCount),AsPerc=round(sum(SvCount)/totalSVCount484,2))
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





