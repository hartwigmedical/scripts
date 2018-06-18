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


## Cluster and Chain analysis

rm(svData)

# SV data file
svData = read.csv('~/logs/CLUSTER_V25.csv')
nrow(svData)
View(head(svData,100))

# filter out multiple biopsy (approximately)
svData = svData %>% filter(!grepl("DRUP", SampleId)&!grepl("TIII", SampleId)&!grepl("TII", SampleId))
nrow(svData)

# FILTER FOR PONCount <2 for all subsequent analyses - no longer required since done already
View(svData %>% filter(PONCount>=2))
svData = svData %>% filter(PONCount<2)

svData$IsLINE = ifelse(svData$LEStart!='false'|svData$LEEnd!='false',1,0)
svData$IsFS = ifelse(svData$FSStart!='false'|svData$FSEnd!='false',1,0)
svData$Length = ifelse(svData$Type=='BND'|svData$Type=='INS'|svData$ArmEnd!=svData$ArmStart, -1, svData$PosEnd-svData$PosStart)
svData$DoubleDupBE = ifelse(svData$DupBEStart=='true'&svData$DupBEEnd=='true',1,0)
svData$SingleDupBE = ifelse(svData$DoubleDupBE==0&(svData$DupBEStart=='true'|svData$DupBEEnd=='true'),1,0)
svData$MantaPrecise = ifelse(svData$MantaPrecise=='true',1,0)
svData$TICount = ifelse(svData$LnkTypeStart=='TI',0.5,0)+ifelse(svData$LnkTypeEnd=='TI',0.5,0)
svData$DBCount = ifelse(svData$LnkTypeStart=='DB',0.5,0)+ifelse(svData$LnkTypeEnd=='DB',0.5,0)
svData$IsTI = ifelse(svData$TICount>0,1,0)
svData$IsDB = ifelse(svData$DBCount>0,1,0)
svData$ShortTICount=ifelse(svData$LnkTypeStart=='TI'&svData$LnkLenStart<=500,0.5,0)+ifelse(svData$LnkTypeEnd=='TI'&svData$LnkLenEnd<=500,0.5,0)
# svData$IsDB=ifelse(svData$NearestDBLen>-1&(svData$NearestDBLen<svData$NearestTILen|svData$NearestTILen<30),1,0)
# svData$IsTI=ifelse(svData$NearestTILen>=30&svData$IsDB==0,1,0)
svData$DoubleTI = ifelse(svData$TICount==1,1,0)
svData$DoubleDB = ifelse(svData$DBCount==1,1,0)
svData$IsSpan = ifelse(svData$TransType=='SPAN',1,0)
svData$IsTrans = ifelse(svData$TransType=='TRANS',1,0)
svData$ClusterSize = ifelse(svData$ClusterCount==1,'None',ifelse(svData$ClusterCount<=4,'Small','Large'))
svData$IsConsistent = ifelse(svData$Consistency==0,1,0)
svData$ChainCount = ifelse(svData$ChainCount>0,svData$ChainCount,ifelse(svData$IsTI==0&svData$IsDB==0,0,1)) # set ChainCount to 1 for single link

View(svData)


# set stressed state
svData$ArmExpStart = ifelse(svData$ArmExpStart>0,svData$ArmExpStart,0.1)
svData$StressedPPStart = round(1 - ppois(svData$ArmCountStart - 1, svData$ArmExpStart),4)
svData$ArmExpEnd = ifelse(svData$ArmExpEnd>0,svData$ArmExpEnd,0.1)
svData$StressedPPEnd = round(1 - ppois(svData$ArmCountEnd - 1, svData$ArmExpEnd),4)
svData$IsStressed = ifelse((svData$StressedPPStart<=0.001&svData$ArmCountStart>=10)|(svData$StressedPPEnd<=0.001&svData$ArmCountEnd>=10),1,0)
nrow(svData %>% filter(IsStressed==1))
View(svData %>% filter(ClusterCount==10))

# allocation of known / resolved types
svData$ResolvedType = 'NONE'

# RESOLVED TYPE: non-clustered SVs - simple types
svData$SampleClusterId = paste(svData$SampleId,svData$ClusterId,sep='_')
svData$ResolvedType = ifelse(svData$ClusterCount==1&(svData$Type=='DEL'|svData$Type=='DUP')&svData$ArmStart!=svData$ArmEnd,'NC_INVALID',svData$ResolvedType)
svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$ClusterCount==1&svData$Type=='DEL','NC_DEL',svData$ResolvedType)
svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$ClusterCount==1&svData$Type=='DUP','NC_DUP',svData$ResolvedType)
svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$ClusterCount==1&svData$Type=='INS','NC_INS',svData$ResolvedType)
svData$IsResolved = ifelse(svData$ResolvedType=='NONE',0,1)
View(svData %>% group_by(ResolvedType,IsStressed) %>% summarise(Count=n()) %>% arrange(ResolvedType,IsStressed)%>% spread(IsStressed,Count))
View(svData %>% group_by(IsResolved,IsStressed) %>% summarise(Count=n()) %>% arrange(IsResolved,IsStressed) )
View(svData %>% group_by(ResolvedType) %>% summarise(Count=n()) %>% arrange(ResolvedType))
View(svData %>% group_by(ResolvedType,IsStressed,ClusterSize) %>% summarise(Count=n()) %>% arrange(ResolvedType,IsStressed,ClusterSize) %>% spread(IsStressed,Count))


# RESOLVED TYPE: non-clustered SVs, dubious INVs and BNDs
dubiousNCSingleSVs = (svData %>% filter(ResolvedType=='NONE'&ClusterCount==1&(Type=='BND'|Type=='INV')&IsStressed==0&IsConsistent==0)
                      %>% filter(Ploidy<0.5&(AdjCNChgStart<0.5|AdjCNChgEnd<0.5)))

View(dubiousNCSingleSVs)
svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$Id %in% dubiousNCSingleSVs$Id,'NC_INVALID',svData$ResolvedType)

# remaining unresolved single SVs
View(svData %>% filter(ResolvedType=='NONE'&ClusterCount==1) %>% group_by(Type,IsStressed) %>% summarise(Count=n()))


# RESOLVED TYPE: span SVs
svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$IsSpan==1,'SPAN',svData$ResolvedType)
svData$IsResolved = ifelse(svData$ResolvedType=='NONE',0,1)
View(svData %>% group_by(ResolvedType,IsStressed) %>% summarise(Count=n()) %>% arrange(ResolvedType,IsStressed))


# svData$ArmExpStart = ifelse(svData$ArmExpStart>0,svData$ArmExpStart,0.1)
# svData$StressedPoissonProb = round(1 - ppois(svData$ArmCountStart - 1, svData$ArmExpStart),4)
# svData$IsStressed = ifelse(svData$ArmCountStart>=10 & svData$StressedPoissonProb <= 0.001,1,0)
# nrow(svData %>% filter(IsStressed==0))

# TI and DB analysis

# proportion of types of SVs in TIs, DBs or transitives

# nrow(svData %>% filter(IsSpan==1))
# View(head(svData))
# View(svData %>% filter(ShortTICount>1))
# nrow(svData %>% filter(IsDB==1))
# nrow(svData %>% filter(IsTI==1))
# View(svData %>% filter(ClusterCount==3&DoubleDupBE==1&IsSpan==0) %>% arrange(SampleId,ClusterId))
# nrow(svData %>% filter(ClusterCount==3&DoubleDupBE==1&IsSpan==1&MantaPrecise=='true'))
# nrow(svData %>% filter(ClusterCount>=3&ClusterCount<=50&DoubleDupBE==1&IsSpan==1&MantaPrecise=='true'))
# nrow(svData %>% filter(MantaPrecise=='true'))
# View(svData %>% group_by(TransType) %>% summarise(Count=n()))
# View(svData %>% group_by(LnkTypeStart) %>% summarise(Count=n()))

# Chain analysis and Transitive analysis
clusteredSvs = svData %>% filter(ClusterCount>1)
nrow(clusteredSvs)

allClusterData = (clusteredSvs %>% group_by(SampleId,ClusterId)
                   %>% summarise(SvCount=n(),
                                 SampleClusterId=first(ClusterId),
                                 ClusterCount=first(ClusterCount),
                                 ClusterDesc=first(ClusterDesc),
                                 Consistency=first(Consistency),
                                 TICount=sum(TICount),
                                 DBCount=sum(DBCount),
                                 TransCount=sum(IsTrans==1),
                                 ShortTICount=sum(ShortTICount),
                                 SpanCount=sum(IsSpan==1),
                                 DoubleDupBECount=sum(DoubleDupBE),
                                 SingleDupBECount=sum(SingleDupBE),
                                 DupBECount=first(DupBECount),
                                 DupBESiteCount=first(DupBESiteCount),
                                 IsConsistent=first(ifelse(Consistency==0,1,0)),
                                 CrossArmCount=sum(Type=='BND'|ArmStart!=ArmEnd),
                                 ArmCount=sum(first(ArmCount)),
                                 DelCount=sum(Type=='DEL'),
                                 DupCount=sum(Type=='DUP'),
                                 InsCount=sum(Type=='INS'),
                                 InvCount=sum(Type=='INV'),
                                 BndCount=sum(Type=='BND'),
                                 ChainCount=sum(ChainId>0),
                                 MaxChainCount=max(ChainCount),
                                 MaxChainTIs=max(ChainTICount),
                                 MaxChainDBs=max(ChainDBCount),
                                 StressedCount=sum(IsStressed),
                                 LineCount=sum(IsLINE))
                   %>% arrange(SampleId,ClusterId))

# RESOLVED TYPE: LINE and SVs in LINE clusters
lineClusters = allClusterData %>% filter(LineCount>0)
nrow(lineClusters)
sum(lineClusters$SvCount)
sum(lineClusters$LineCount)

# translate LINE info back to svData
svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$SampleClusterId %in% lineClusters$SampleClusterId,'LINE',svData$ResolvedType)
svData$ResolvedType = ifelse(svData$ResolvedType=='LINE'&svData$IsLINE==0,'LINE_CLUST',svData$ResolvedType)
svData$IsResolved = ifelse(svData$ResolvedType=='NONE',0,1)
View(svData %>% group_by(ResolvedType) %>% summarise(Count=n()) %>% arrange(ResolvedType))

# all other clusters
clusterData = allClusterData %>% filter(LineCount==0)
clusterData$IsStressed = ifelse(clusterData$StressedCount>0,1,0)
clusterData$ResolvedType = 'NONE'

# limit to cluster sizes and less since linking logic cuts out beyond that
clusterData = clusterData %>% filter(ClusterCount<=100)
nrow(clusterData)
sum(clusterData$ClusterCount)
sum(clusterData$SvCount)
sum(clusterData$StressedCount)

# RESOLVED TYPE: reciprocal inversions & translocations
svData$IsRecipInv = ifelse(svData$ClusterCount==2&svData$Type!='DEL'&svData$Type!='DUP'&svData$DoubleDB&svData$LnkSvStart==svData$LnkSvEnd,1,0)
View(svData %>% filter(IsRecipInv==1))
View(svData %>% filter(IsRecipInv==1&ClusterCount==2) %>% group_by(ClusterDesc) %>% summarise(Count=n()))
reciprocals = clusterData %>% filter(ClusterCount==2&(BndCount==2|InvCount==2)&DBCount==2&IsConsistent==1)
View(reciprocals)
svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$Type=='BND'&svData$SampleClusterId %in% reciprocals%SampleClusterId,'RECIP_TRANS',svData$ResolvedType)
svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$Type!='BND'&svData$SampleClusterId %in% reciprocals%SampleClusterId,'RECIP_INV',svData$ResolvedType)
clusterData$ResolvedType = ifelse(clusterData$ResolvedType=='NONE'&clusterData$BndCount==2&clusterData$SampleClusterId %in% reciprocals$SampleClusterId,'RECIP_TRANS',clusterData$ResolvedType)
clusterData$ResolvedType = ifelse(clusterData$ResolvedType=='NONE'&clusterData$InvCount==2&clusterData$SampleClusterId %in% reciprocals$SampleClusterId,'RECIP_INV',clusterData$ResolvedType)

View(svData %>% group_by(ResolvedType,IsStressed) %>% summarise(Count=n()) %>% arrange(ResolvedType,IsStressed))
View(svData %>% group_by(ResolvedType,IsStressed,IsConsistent) %>% summarise(Count=n()) %>% arrange(ResolvedType,IsStressed,IsConsistent))
View(svData %>% group_by(ResolvedType) %>% summarise(Count=n()) %>% arrange(ResolvedType))
View(clusterData %>% group_by(ResolvedType) %>% summarise(Count=n()) %>% arrange(ResolvedType))


# complex clusters
clusterData$ClusterBucket = ifelse(clusterData$ClusterCount<=10,clusterData$ClusterCount,ifelse(clusterData$ClusterCount<=50,round(clusterData$ClusterCount/5)*5,round(clusterData$ClusterCount/10)*10))
clusterData$ShortTIPerc = round(clusterData$ShortTICount/clusterData$SvCount,2)
clusterData$DBPerc = round(clusterData$DBCount/clusterData$SvCount,2)
clusterData$LinkPerc = round((clusterData$TICount+clusterData$DBCount+clusterData$SpanCount)/clusterData$SvCount,2)
View(clusterData)

# RESOLVED TYPE: all linked DBs
# dbChains = clusterData %>% filter(ClusterCount>=3&(DBCount+SpanCount)==ClusterCount&Consistency==0)
shortTIDBChains = clusterData %>% filter(ClusterCount>=2&(ShortTICount+DBCount+SpanCount)==ClusterCount&Consistency==0)
nrow(shortTIDBChains)
sum(shortTIDBChains$SvCount)
sum(shortTIDBChains$SpanCount)
svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$SampleClusterId %in% shortTIDBChains$SampleClusterId,'DB_CHAIN',svData$ResolvedType)
clusterData$ResolvedType = ifelse(clusterData$ResolvedType=='NONE'&clusterData$SampleClusterId %in% shortTIDBChains$SampleClusterId,'DB_CHAIN',clusterData$ResolvedType)
View(svData %>% group_by(ResolvedType) %>% summarise(Count=n()) %>% arrange(ResolvedType))
View(clusterData %>% group_by(ResolvedType) %>% summarise(Count=n()) %>% arrange(ResolvedType))


# RESOLVED TYPE: all TIs
tiChains = clusterData %>% filter(ClusterCount>=2&(TICount+SpanCount)>=(ClusterCount-1)&DBCount==0&Consistency==0)
View(tiChains)
View(tiChains %>% filter(SampleClusterId %in% shortTIDBChains$SampleClusterId))
svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$SampleClusterId %in% tiChains$SampleClusterId,'TI_CHAIN',svData$ResolvedType)
clusterData$ResolvedType = ifelse(clusterData$ResolvedType=='NONE'&clusterData$SampleClusterId %in% tiChains$SampleClusterId,'TI_CHAIN',clusterData$ResolvedType)
View(svData %>% group_by(ResolvedType) %>% summarise(Count=n()) %>% arrange(ResolvedType))
View(clusterData %>% group_by(ResolvedType) %>% summarise(Count=n(), SvCount=sum(SvCount)) %>% arrange(ResolvedType))
View(clusterData %>% group_by(ResolvedType,IsStressed) %>% summarise(Count=n(), SvCount=sum(SvCount)) %>% arrange(ResolvedType,IsStressed))

View(svData %>% filter(ClusterCount<20,IsSpan==0) %>% group_by(ClusterCount,ChainCount) %>%count() %>% spread(ChainCount,n))

View(svData %>% filter(ChainCount==0&svData$IsTI==1|svData$IsDB==1))

# RESOLVED TYPE: all chained SVs
allChainedClusters = clusterData %>% filter(ResolvedType=='NONE'&IsConsistent==1&DoubleDupBECount==SpanCount&(DBCount+TICount+SpanCount)>=(ClusterCount-1))
nrow(allChainedClusters)
sum(allChainedClusters$SvCount)
svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$SampleClusterId %in% allChainedClusters$SampleClusterId,'COMPLEX_CHAIN',svData$ResolvedType)
clusterData$ResolvedType = ifelse(clusterData$ResolvedType=='NONE'&clusterData$SampleClusterId %in% allChainedClusters$SampleClusterId,'COMPLEX_CHAIN',clusterData$ResolvedType)


# RESOLVED TYPE: other small clusters
cc2Unres = clusterData %>% filter(ClusterCount==2&ResolvedType=='NONE')
nrow(cc2Unres)
nrow(cc2Unres %>% filter(InsCount>=(ClusterCount-1))) # Inserts part of cluster
nrow(cc2Unres %>% filter(DelCount==2)) # 2 DELs clustered
nrow(cc2Unres %>% filter(IsConsistent==0&SingleDupBECount>0)) # Likely to be missing another SV
nrow(cc2Unres %>% filter(DBCount==2)) # Invalid from the types involved (ie not INVs or BNDs)
cc2Unres$Described = 0
cc2Unres$Described = ifelse(cc2Unres$Described==0&cc2Unres$InsCount>=(cc2Unres$ClusterCount-1),1,cc2Unres$Described)
cc2Unres$Described = ifelse(cc2Unres$Described==0&cc2Unres$DelCount==2,1,cc2Unres$Described)
cc2Unres$Described = ifelse(cc2Unres$Described==0&cc2Unres$IsConsistent==0&cc2Unres$SingleDupBECount>0,1,cc2Unres$Described)
cc2Unres$Described = ifelse(cc2Unres$Described==0&cc2Unres$DBCount==2,1,cc2Unres$Described)
cc2Unres$Described = ifelse(cc2Unres$Described==0&cc2Unres$DBCount==2,1,cc2Unres$Described)
sum(cc2Unres$Described)
sum(cc2Unres$Described==0)

# valid types
View(clusterData %>% filter(ResolvedType=='NONE'&IsConsistent==1&SingleDupBECount==0&DoubleDupBECount==0))
View(clusterData %>% filter(ResolvedType=='NONE'&IsConsistent==1&DoubleDupBECount==SpanCount))
View(clusterData %>% filter(ResolvedType=='NONE'&IsConsistent==1&DoubleDupBECount==SpanCount&(DBCount+TICount+SpanCount)==ClusterCount))

View(cc2Unres %>% filter(InvCount==2&TICount==1))

View(cc2Unres %>% filter(Described==0) %>% group_by(IsConsistent,IsStressed,ClusterDesc)
                 %>% summarise(Count=n())
                 %>% arrange(IsConsistent,IsStressed,ClusterDesc))

View(cc2Unres %>% filter(Described==0) %>% group_by(IsConsistent,ClusterDesc)
     %>% summarise(Count=n(),
                   TICount=sum(TICount),
                   DBCount=sum(DBCount),
                   DupBECount=sum(SingleDupBECount)/2)
     %>% arrange(IsConsistent,ClusterDesc))

View(cc2Unres %>% filter(TICount==2) %>% group_by(IsConsistent,ClusterDesc)
     %>% summarise(Count=n(),
                   TICount=sum(TICount),
                   DBCount=sum(DBCount),
                   DupBECount=sum(SingleDupBECount)/2)
     %>% arrange(IsConsistent,ClusterDesc))

# invalid due to an INS
nrow(clusterData %>% filter(ClusterCount>=2&(InsCount>=(ClusterCount-1))))
nrow(clusterData %>% filter(ClusterCount==2&InsCount>0))

cc2Unres = clusterData %>% filter(ClusterCount==2&Consistency==0)


# TAKING STOCK
View(svData %>% group_by(ResolvedType,IsStressed,ClusterSize)
     %>% summarise(Count=n()) %>% arrange(ResolvedType,IsStressed,ClusterSize))

View(svData %>% group_by(ResolvedType,IsStressed,ClusterSize)
     %>% summarise(Count=n()) %>% arrange(ResolvedType,IsStressed,ClusterSize) %>% spread(IsStressed,Count))


clusterData$LinkPercBucket = round(clusterData$LinkPerc/0.2)*0.2

View(clusterData %>% filter(IsConsistent==0&ResolvedType=='NONE') %>% group_by(LinkPercBucket,IsStressed)
     %>% summarise(BucketCount=n(),
                   SvCount=sum(SvCount),
                   TICount=sum(TICount),
                   DBCount=sum(DBCount),
                   ShortTICount=sum(ShortTICount))
     %>% arrange(IsStressed,LinkPercBucket))


# Characteristics of Stressed Large Clusters
scData = clusterData %>% filter(ResolvedType=='NONE'&StressedCount>0)
nrow(scData)
sum(scData$SvCount)
View(scData %>% filter(ClusterCount>=10&ClusterCount<=50))



# Re-evaluate Stressed Arms in light of Resolved SVs
urData = svData %>% filter(ResolvedType=='NONE')
nrow(urData)

svArmData = urData
svArmData$Chr = svArmData$ChrStart
svArmData$Arm = svArmData$ArmStart
svArmData$ArmCount = svArmData$ArmCountStart
svArmData$ArmExpected = svArmData$ArmExpStart

svArmBndData = urData %>% filter(Type=='BND')
svArmBndData$Chr = svArmBndData$ChrEnd
svArmBndData$Arm = svArmBndData$ArmEnd
svArmBndData$ArmCount = svArmBndData$ArmCountEnd
svArmBndData$ArmExpected = svArmBndData$ArmExpEnd

# merge rows prior to arm grouping
urArmData = rbind(svArmData, svArmBndData)
nrow(urArmData)

samArmData = urArmData %>% group_by(SampleId,Chr,Arm) %>% summarise(Count=n())
View(samArmData)

# determine IsStressed using poisson distribution using expected vs actual SV counts per arm
combinedArmData$ArmExpected = ifelse(combinedArmData$ArmExpected>0,combinedArmData$ArmExpected,0.1)





# group by cluster size, resolved and stressed
resolvedClusters = clusterData %>% filter(ResolvedType!='NONE')

resolvedClusters$ClusterBucket = ifelse(resolvedClusters$ClusterCount<=10,resolvedClusters$ClusterCount,
                                 ifelse(resolvedClusters$ClusterCount<=25,round(resolvedClusters$ClusterCount/5)*5,
                                 ifelse(resolvedClusters$ClusterCount<=100,round(resolvedClusters$ClusterCount/10)*10,
                                 round(resolvedClusters$ClusterCount/50)*50)))

resolvedClusters$IsStressed = ifelse(resolvedClusters$StressedCount>0,1,0)
resolvedClusters$IsResolved = ifelse(resolvedClusters$ResolvedCount>0,1,0)


resolvedClusterSummary = (resolvedClusters %>% group_by(ClusterBucket,IsStressed)
                  %>% summarise(BucketCount=n(),
                                SvCount=sum(SvCount),
                                TICount=sum(TICount),
                                DBCount=sum(DBCount),
                                ShortTICount=sum(ShortTICount))
                  %>% arrange(ClusterBucket))

View(resolvedClusterSummary)

View(resolvedClusters %>% filter(LineCount==0) %>% group_by(ClusterBucket,IsResolved,IsStressed)
                          %>% summarise(Count=sum(SvCount))
                          %>% spread(IsStressed,Count)
                          %>% arrange(ClusterBucket))


# clusters with high proportion of TIs or DBs
highLinkData = clusterData %>% filter(ClusterCount>=3&LinkPerc>=0.90)
View(highLinkData)

chainsDBShortTIs = clusterData %>% filter(ClusterCount>=3&(ShortTIPerc+DBPerc)>=0.8)
View(chainsDBShortTIs)

View(clusterData %>% filter(ShortTIPerc==1&ClusterCount>=3) %>% group_by(ClusterCount) %>% summarise(Count=n()))
View(clusterData %>% filter(DBPerc==1&ClusterCount>=3) %>% group_by(ClusterCount,StressedCount>0) %>% summarise(Count=n()))

View(svData %>% filter(ClusterCount==3&DBCount==3&ResolvedType=='NONE'))
View(svData %>% filter(ClusterCount==3&ResolvedType=='NONE'))


# clusters which are all/mostly chains of TIs or DBs



# summary stats across all cluster sizes
clusterSummary = (clusterData %>% group_by(ClusterBucket)
               %>% summarise(BucketCount=n(),
                             SvCount=sum(SvCount),
                             TICount=sum(TICount),
                             DBCount=sum(DBCount),
                             TransCount=sum(TransCount),
                             ShortTICount=sum(ShortTICount),
                             SpanCount=sum(SpanCount),
                             DoubleDupBECount=sum(DoubleDupBECount),
                             SingleDupBECount=sum(SingleDupBECount))
               %>% arrange(ClusterBucket))

View(clusterSummary)

write.csv(clusterSummary, "~/logs/r_output/clusterSummary.csv")


# Transitive analysis
clusterSummary = ((clusterData %>% group_by(ClusterBucket)
                   %>% summarise(SvCount=n(),
                                 ClusterCount=first(ClusterCount),
                                 TICount=sum(TICount),
                                 DBCount=sum(DBCount),
                                 ShortTICount=sum(ShortTICount),
                                 SpanCount=sum(SpanCount),
                                 LineCount=sum(LineCount))
                   %>% arrange(ClusterBucket)))

View(clusterSummary)

# number of TIs, DBs, double TIs and DBs, doubleDupBEs by cluster size
lnksByType = (svData %>% filter(ClusterSize!='None'&Type!='INS') %>% group_by(Type,ClusterSize)
                     %>% summarise(Count=n(),
                                   TICount=sum(IsTI==1),
                                   DBCount=sum(IsDB==1),
                                   TIPerc=round(sum(IsTI==1)/n(),2),
                                   DBPerc=round(sum(IsDB==1)/n(),2),
                                   SpanCount=sum(IsSpan==1),
                                   TransCount=sum(IsTrans==1),
                                   ShortTICount=sum(ShortTICount))
                     %>% arrange(Type,ClusterSize))

# proportion of SVs by type removed/converted by transitive logic
lnksByType$TransRemoved = round(lnksByType$TransCount/lnksByType$Count,3)
lnksByType$SpansPercent = round(lnksByType$SpanCount/lnksByType$Count,3)

# length of TIs and DBs
startLinks = svData %>% filter(LnkTypeStart=='TI'|LnkTypeStart=='DB')
startLinks$LinkType = startLinks$LnkTypeStart
startLinks$LinkLen = startLinks$LnkLenStart
endLinks = svData %>% filter(LnkTypeEnd=='TI'|LnkTypeEnd=='DB')
endLinks$LinkType = endLinks$LnkTypeEnd
endLinks$LinkLen = endLinks$LnkLenEnd
combinedLinks = rbind(startLinks, endLinks)

# switch short TIs back from negative DBs
combinedLinks$DubiousShortTI = ifelse(combinedLinks$LinkLen<0,1,0)
combinedLinks$LinkType = ifelse(combinedLinks$DubiousShortTI==1,"TI",ifelse(combinedLinks$LinkType=="DB","DB","TI"))
combinedLinks$LinkLen = ifelse(combinedLinks$DubiousShortTI==1,-combinedLinks$LinkLen,combinedLinks$LinkLen)

View(combinedLinks %>% filter(DubiousShortTI==1))
nrow(combinedLinks %>% filter(LinkType=="DB"))
View(combinedLinks %>% filter(LinkType=="DB"))

nrow(combinedLinks)

combinedLinks$LenBucketLog = ifelse(combinedLinks$LinkLen==0,0,
                              ifelse(combinedLinks$LinkLen>0,2**round(log(combinedLinks$LinkLen,2),0),
                                   -(2**round(log(-combinedLinks$LinkLen,2),0))))

nrow(combinedLinks %>% filter(LinkLen==0))

combinedLinks$LenBucket = ifelse(combinedLinks$LinkLen<28,combinedLinks$LinkLen,2**round(log(combinedLinks$LinkLen,2),0))

linkLengthStats = (combinedLinks %>% group_by(LenBucket)
                   %>% summarise(TICount=sum(LinkType=='TI'), DBCount=sum(LinkType=='DB'))
                   %>% arrange(LenBucket))


View(linkLengthStats)

linkLenPlot2 = (ggplot(data = linkLengthStats, aes(x = LenBucket))
               + geom_line(aes(y=TICount, colour='TI'))
               + geom_line(aes(y=DBCount, colour='DB'))
               + scale_x_log10()
               # + coord_cartesian(xlim = c(-1e4, 1e6))
               + ylab("SV Count") + labs(title = "Link Lengths - TIs and DBs"))

print(linkLenPlot2)

linkLengthActualStats = (combinedLinks %>% filter(LinkLen <= 500) %>% group_by(LinkLen)
                   %>% summarise(TICount=sum(LinkType=='TI'), DBCount=sum(LinkType=='DB'))
                   %>% arrange(LinkLen))

linkLenPlot3 = (ggplot(data = linkLengthActualStats, aes(x = LinkLen))
                + geom_line(aes(y=TICount, colour='TI'))
                + geom_line(aes(y=DBCount, colour='DB'))
                # + scale_x_log10()
                # + coord_cartesian(xlim = c(-1e4, 1e6))
                + scale_x_continuous()
                + ylab("SV Count") + labs(title = "Link Lengths - TIs and DBs"))

print(linkLenPlot3)


# linkLengthStats = (combinedLinks %>% group_by(LenBucketLog)
#                             %>% summarise(TICount=sum(LinkType=='TI'), DBCount=sum(LinkType=='DB'))
#                             %>% arrange(LenBucketLog))
View(linkLengthStats)

linkLenPlot = (ggplot(data = linkLengthStats %>% filter(LenBucketLog < 1e8), aes(x = LenBucketLog))
                      + geom_line(aes(y=TICount, colour='TI'))
                      + geom_line(aes(y=DBCount, colour='DB'))
                      # + scale_x_log10()
                      # + coord_cartesian(xlim = c(-1e4, 1e6))
                      + scale_x_continuous()
                      # + facet_wrap(as.formula(paste("~", facetWrap)))
                      + ylab("SV Count") + labs(title = "Link Lengths - TIs and DBs")
)

print(linkLenPlot)

combinedLinks$LenBucketFixed = round(combinedLinks$LinkLen/10)*10

linkLengthFixedStats = (combinedLinks %>% filter(LenBucketFixed>=-1000&LenBucketFixed<=1000) %>% group_by(LenBucketFixed)
                   %>% summarise(TICount=sum(LinkType=='TI'), DBCount=sum(LinkType=='DB'))
                   %>% arrange(LenBucketFixed))

View(linkLengthFixedStats)

linkFixedLenPlot = (ggplot(data = linkLengthFixedStats, aes(x = LenBucketFixed))
               + geom_line(aes(y=TICount, colour='TI'))
               + geom_line(aes(y=DBCount, colour='DB'))
               + scale_x_continuous()
               # + facet_wrap(as.formula(paste("~", facetWrap)))
               + ylab("SV Count") + labs(title = "Link Lengths - TIs and DBs"))
               
print(linkFixedLenPlot)

View(clusterData %>% filter(SampleId=='CPCT02060039T'))


# cluster-2 stats
cc2Data = clusterData %>% filter(ClusterCount==2)
nrow(cc2Data)
View(cc2Data)
nrow(cc2Data %>% filter(DoubleDupBECount==2))
nrow(cc2Data %>% filter(SingleDupBECount==2))
# cc2Data$InvalidCrossArm = ifelse(cc2Data$CrossArmCount==1|cc2Data$ArmCount==3,1,0)
cc2Data$IsComplete = ifelse(cc2Data$DoubleDupBECount==0&cc2Data$IsConsistent==1&cc2Data$SingleDupBECount==0,1,0)
sum(cc2Data$IsComplete==1)
sum(cc2Data$IsComplete==0)

cc2Incomplete = cc2Data %>% filter(IsComplete==0)
cc2Incomplete$HasDupBEs = ifelse(cc2Incomplete$DoubleDupBECount>0|cc2Incomplete$SingleDupBECount>0,1,0)
cc2Incomplete$HasLinks = ifelse(cc2Incomplete$TICount>0|cc2Incomplete$DBCount>0,1,0)
nrow(cc2Incomplete)

View(cc2Incomplete %>% group_by(ClusterDesc,HasDupBEs,IsConsistent,CrossArmCount,HasLinks)
     %>% summarise(Count=n())
     %>% arrange(ClusterDesc,HasDupBEs,IsConsistent,CrossArmCount,HasLinks))

cc2NoLinks = cc2Data %>% filter(DBCount==0&TICount==0)
nrow(cc2NoLinks)

cc2Complete = cc2Data %>% filter(DoubleDupBECount==0&IsConsistent==T&SingleDupBECount==0)
nrow(cc2Complete)

View(cc2Complete %>% group_by(ClusterDesc,DBCount,TICount)
     %>% summarise(Count=n())
     %>% arrange(ClusterDesc,DBCount,TICount))


# Chain analysis

max(svData$ChainCount)
max(svData$ChainTICount)
max(svData$ChainDBCount)

View(svData)

nrow(svData %>% filter(SingleDupBE==1&IsTI==0&IsDB==0))           
nrow(svData %>% filter(SingleDupBE==1&IsTI==0&IsDB==0&IsLINE==0))           
nrow(svData %>% filter(SingleDupBE==1&IsTI==0&IsDB==0&IsStressed==1))           

svData$ChainCount = ifelse(svData$ChainId==0,0,svData$ChainCount+1) # since the count is of links
View(svData %>% filter(ChainId==0&ChainCount>0))

lineClusters = clusterChaining %>% filter(LineCount>0)
sum(lineClusters$SvCount - lineClusters$LineCount)
sum(lineClusters$LineCount)


clusterChaining = (svData %>% filter(ClusterCount>1) %>% group_by(SampleId,ClusterId)
                   %>% summarise(SvCount=n(),
                                 ClusterCount=first(ClusterCount),
                                 ChainCount=n_distinct(ChainId)-1,
                                 ChainedCount=sum(ChainId>0), 
                                 LinkedCount=sum(ChainId>0|IsDB==1|IsTI==1), 
                                 # TICount=round(sum(ifelse(ChainId>0,ChainTICount/ChainCount,0)),2),
                                 # DBCount=round(sum(ifelse(ChainId>0,ChainDBCount/ChainCount,0)),2),
                                 DBCount=sum(IsDB==1),
                                 TICount=sum(IsTI==1),
                                 SpanDupBECount=sum(DoubleDupBE==1),
                                 SpanCount=sum(IsSpan==1),
                                 SingleDupBENoLinkCount=sum(SingleDupBE==1&IsTI==0&IsDB==0),
                                 DelCount=sum(Type=='DEL'),
                                 DupCount=sum(Type=='DUP'),
                                 InsCount=sum(Type=='INS'),
                                 InvCount=sum(Type=='INV'),
                                 BndCount=sum(Type=='BND'),
                                 LineCount=sum(IsLINE==1))
                   %>% arrange(SampleId,ClusterId))

View(clusterChaining)

regClusterChaining = clusterChaining %>% filter(LineCount==0&ClusterCount<=100)
View(regClusterChaining)

View(regClusterChaining %>% filter(SpanDupBECount==0))

regClusterChaining$RevClusterCount = regClusterChaining$ClusterCount - regClusterChaining$SpanDupBECount
clusteredChains = regClusterChaining %>% filter(RevClusterCount>1)
clusteredChains$ChainedPerc = round(clusteredChains$ChainedCount/clusteredChains$RevClusterCount,4)
clusteredChains$ChainedPercBucket = round(clusteredChains$ChainedPerc/0.1)*0.1

clusteredChains$LinkedPerc = round(clusteredChains$LinkedCount/clusteredChains$RevClusterCount,4)
clusteredChains$LinkedPercBucket = round(clusteredChains$LinkedPerc/0.1)*0.1

# clusteredChains$ClusterSize = ifelse(clusteredChains$RevClusterCount<=3,'Small',ifelse(clusteredChains$RevClusterCount<=10,'Med','Large'))
clusteredChains$ClusterSize = ifelse(clusteredChains$RevClusterCount<=5,clusteredChains$RevClusterCount,
                              ifelse(clusteredChains$RevClusterCount<=25,round(clusteredChains$RevClusterCount/5)*5,round(clusteredChains$RevClusterCount/10)*10))

View(clusteredChains)

clusteredChainStats = (clusteredChains %>% group_by(ChainedPercBucket,ClusterSize) 
                     %>% summarise(Count=n()) # AsPerc=round(n()/nrow(regClusterChaining),2) 
                     %>% arrange(ClusterSize,ChainedPercBucket))

View(clusteredChainStats)
View(clusteredChainStats %>% spread(ClusterSize,Count))

View(clusteredChains %>% group_by(ChainedPercBucket,ClusterSize) %>% summarise(SvCount=sum(SvCount)) %>% spread(ClusterSize,SvCount))
View(clusteredChains %>% group_by(LinkedPercBucket,ClusterSize) %>% summarise(SvCount=sum(SvCount)) %>% spread(ClusterSize,SvCount))

# remove unresolved single duplicate BEs and repeat
nrow(clusteredChains %>% filter(SingleDupBENoLinkCount==RevClusterCount))
cleanClusteredChains = clusteredChains %>% filter(SingleDupBENoLinkCount<RevClusterCount)

cleanClusteredChains$CleanLinkedPerc = round(cleanClusteredChains$LinkedCount/(cleanClusteredChains$RevClusterCount-cleanClusteredChains$SingleDupBENoLinkCount),4)
cleanClusteredChains$CleanLinkedPercBucket = round(cleanClusteredChains$CleanLinkedPerc/0.1)*0.1
View(cleanClusteredChains %>% group_by(CleanLinkedPercBucket,ClusterSize) %>% summarise(SvCount=sum(SvCount)) %>% spread(ClusterSize,SvCount))


View(clusteredChains %>% filter(RevClusterCount==4))

View(svData %>% filter(ClusterId==85&SampleId=='CPCT02010050T'))

# cosine similarity

vector1 = svData$AdjCNChgEnd
vector1 = svData$AdjCNChgEnd


cosine.similarity(.alpha, .beta,)




# bucket by percent chained
clusterChaining$ChainedBucket = round(clusterChaining$ChainedPerc/0.2)*0.2
clusterChaining$ClusterSize = ifelse(clusterChaining$ClusterCount<=4,'Small','Large')

chainedStats = (clusterChaining %>% group_by(ClusterSize,ChainedBucket)
                %>% summarise(ClusterCount=n(),
                              SvCount=sum(SvCount),
                              TIPerc=round(sum(TICount/SvCount),2),
                              DBPerc=round(sum(DBCount/SvCount),2),
                              ShortTIPerc=round(sum(ShortTICount/SvCount),2))
                %>% arrange(ClusterSize,ChainedBucket))

View(chainedStats)

chainData = (clusteredSvs %>% filter(ChainCount>0) %>% group_by(SampleId,ClusterId,ChainId)
             %>% summarise(SvCount=n(),
                           TICount=first(ChainTICount),
                           TIPerc=round(first(ChainTICount)/first(ChainCount),2),
                           DBCount=first(ChainDBCount),
                           DBPerc=round(first(ChainDBCount)/first(ChainCount),2),
                           ShortTICount=sum((LnkStartType=='TI'&LnkStartLen<=500)|(LnkEndType=='TI'&LnkEndLen<=500))/2
             )
             %>% arrange(SampleId,ClusterId,ChainId))


View(chainData)



# TI length
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




svData$LenBucketLog = ifelse(svData$Length>0,2**round(log(svData$Length,2),0),0)

resolvedLengthStats = (svData %>% filter(ResolvedType=='NC_DEL'|ResolvedType=='NC_DUP') %>% group_by(LenBucketLog,IsStressed)
                       %>% summarise(DelCount=sum(ResolvedType=='NC_DEL'), DupCount=sum(ResolvedType=='NC_DUP'))
                       %>% arrange(LenBucketLog))


svLenPlot = (ggplot(data = resolvedLengthStats, aes(x = LenBucketLog))
             + geom_line(aes(y=DelCount, colour='DELs'))
             + geom_line(aes(y=DupCount, colour='DUPs'))
             + scale_x_log10()
             + facet_wrap(~IsStressed)
             + ylab("SV Count") + labs(title = "Lengths - DELs and DUPs"))

print(svLenPlot)




# OLD CLUSTER-2 analysis




# templated insertions from cluster-count 2, back to same location
cc2Data = svData %>% filter(ClusterCount==2)
nrow(cc2Data)

# group by clusterId and then look at the BE pairs
cc2Clusters = (cc2Data %>% group_by(SampleId,ClusterId)
               %>% summarise(Id1=first(Id), Id2=last(Id), ClusterDesc=first(ClusterDesc),
                             ChrStart1=first(ChrStart),ChrEnd1=first(ChrEnd),
                             ChrStart2=last(ChrStart),ChrEnd2=last(ChrEnd),
                             CrossChr=ifelse(first(Type)=='BND'|last(Type)=='BND',1,0),
                             IsStressed=ifelse(first(IsStressed)==1,1,0),
                             IsLINE=ifelse(first(LEStart)!='false'|first(LEEnd)!='false'|last(LEStart)!='false'|last(LEEnd)!='false',1,0),
                             SS_IsTI=ifelse(first(ChrStart)==last(ChrStart),ifelse((first(PosStart)<last(PosStart) & first(OrientStart)==-1 & last(OrientStart)==1)|(first(PosStart)>last(PosStart) & first(OrientStart)==1 & last(OrientStart)==-1),1,0),0),
                             EE_IsTI=ifelse(first(ChrEnd)==last(ChrEnd),ifelse((first(PosEnd)<last(PosEnd) & first(OrientEnd)==-1 & last(OrientEnd)==1)|(first(PosEnd)>last(PosEnd) & first(OrientEnd)==1 & last(OrientEnd)==-1),1,0),0),
                             SE_IsTI=ifelse(first(ChrStart)==last(ChrEnd),ifelse((first(PosStart)<last(PosEnd) & first(OrientStart)==-1 & last(OrientEnd)==1)|(first(PosStart)>last(PosEnd) & first(OrientStart)==1 & last(OrientEnd)==-1),1,0),0),
                             ES_IsTI=ifelse(first(ChrEnd)==last(ChrStart),ifelse((first(PosEnd)<last(PosStart) & first(OrientEnd)==-1 & last(OrientStart)==1)|(first(PosEnd)>last(PosStart) & first(OrientEnd)==1 & last(OrientStart)==-1),1,0),0),
                             SS_IsDB=ifelse(first(ChrStart)==last(ChrStart),ifelse((first(PosStart)<=last(PosStart) & first(OrientStart)==1 & last(OrientStart)==-1)|(first(PosStart)>=last(PosStart) & first(OrientStart)==-1 & last(OrientStart)==1),1,0),0),
                             EE_IsDB=ifelse(first(ChrEnd)==last(ChrEnd),ifelse((first(PosEnd)<=last(PosEnd) & first(OrientEnd)==1 & last(OrientEnd)==-1)|(first(PosEnd)>=last(PosEnd) & first(OrientEnd)==-1 & last(OrientEnd)==1),1,0),0),
                             SE_IsDB=ifelse(first(ChrStart)==last(ChrEnd),ifelse((first(PosStart)<=last(PosEnd) & first(OrientStart)==1 & last(OrientEnd)==-1)|(first(PosStart)>=last(PosEnd) & first(OrientStart)==-1 & last(OrientEnd)==1),1,0),0),
                             ES_IsDB=ifelse(first(ChrEnd)==last(ChrStart),ifelse((first(PosEnd)<=last(PosStart) & first(OrientEnd)==1 & last(OrientStart)==-1)|(first(PosEnd)>=last(PosStart) & first(OrientEnd)==-1 & last(OrientStart)==1),1,0),0),
                             SS_Len=ifelse(first(ChrStart)==last(ChrStart),abs(first(PosStart)-last(PosStart)),1e7),
                             EE_Len=ifelse(first(ChrEnd)==last(ChrEnd),abs(first(PosEnd)-last(PosEnd)),1e7),
                             SE_Len=ifelse(first(ChrStart)==last(ChrEnd),abs(first(PosStart)-last(PosEnd)),1e7),
                             ES_Len=ifelse(first(ChrEnd)==last(ChrStart),abs(first(PosEnd)-last(PosStart)),1e7),
                             SS_SO=ifelse(first(ChrStart)==last(ChrStart)&first(OrientStart)==last(OrientStart),1,0),
                             SE_SO=ifelse(first(ChrStart)==last(ChrEnd)&first(OrientStart)==last(OrientEnd),1,0),
                             ES_SO=ifelse(first(ChrEnd)==last(ChrStart)&first(OrientEnd)==last(OrientStart),1,0),
                             EE_SO=ifelse(first(ChrEnd)==last(ChrEnd)&first(OrientEnd)==last(OrientEnd),1,0),
                             SS_CA=ifelse(first(ChrStart)!=last(ChrStart)|first(ArmStart)!=last(ArmStart),1,0),
                             SE_CA=ifelse(first(ChrStart)!=last(ChrEnd)|first(ArmStart)!=last(ArmEnd),1,0),
                             ES_CA=ifelse(first(ChrEnd)!=last(ChrStart)|first(ArmEnd)!=last(ArmStart),1,0),
                             EE_CA=ifelse(first(ChrEnd)!=last(ChrEnd)|first(ArmEnd)!=last(ArmEnd),1,0),
                             SS_DupBE=ifelse(first(ChrStart)==last(ChrStart)&first(OrientStart)==last(OrientStart)&abs(first(PosStart)-last(PosStart))<=50,1,0),
                             SE_DupBE=ifelse(first(ChrStart)==last(ChrEnd)&first(OrientStart)==last(OrientEnd)&abs(first(PosStart)-last(PosEnd))<=50,1,0),
                             ES_DupBE=ifelse(first(ChrEnd)==last(ChrStart)&first(OrientEnd)==last(OrientStart)&abs(first(PosEnd)-last(PosStart))<=50,1,0),
                             EE_DupBE=ifelse(first(ChrEnd)==last(ChrEnd)&first(OrientEnd)==last(OrientEnd)&abs(first(PosEnd)-last(PosEnd))<=50,1,0)
               )
               %>% arrange(SampleId,ClusterId))

View(cc2Clusters)
nrow(cc2Clusters)
cc2ClustersOrig = cc2Clusters
cc2Clusters = cc2ClustersOrig

#	Split NONE into ‘LONG TI, LONG DB, NONE, SAME OR'

cc2Clusters$ClusterDesc = stri_replace_all_fixed(cc2Clusters$ClusterDesc, "=1", "")
cc2Clusters$ClusterDesc = stri_replace_all_fixed(cc2Clusters$ClusterDesc, "=2", "_2")
cc2Desc = cc2Clusters %>% group_by(ClusterDesc) %>% summarise(Count=n())
View(cc2Desc)

cc2Clusters$DupBECount = cc2Clusters$SS_DupBE + cc2Clusters$SE_DupBE + cc2Clusters$ES_DupBE + cc2Clusters$EE_DupBE

# cc2Clusters$SS_Len = ifelse((cc2Clusters$SS_IsTI==1|cc2Clusters$SS_IsDB==1)&cc2Clusters$SS_Len<=1e4,cc2Clusters$SS_Len,1e7)
# cc2Clusters$SE_Len = ifelse((cc2Clusters$SE_IsTI==1|cc2Clusters$SE_IsDB==1)&cc2Clusters$SE_Len<=1e4,cc2Clusters$SE_Len,1e7)
# cc2Clusters$ES_Len = ifelse((cc2Clusters$ES_IsTI==1|cc2Clusters$ES_IsDB==1)&cc2Clusters$ES_Len<=1e4,cc2Clusters$ES_Len,1e7)
# cc2Clusters$EE_Len = ifelse((cc2Clusters$EE_IsTI==1|cc2Clusters$EE_IsDB==1)&cc2Clusters$EE_Len<=1e4,cc2Clusters$EE_Len,1e7)

# make note of long DBs and TIs which will fall out of the standard classifications below

# restrict DBs and TIs to within the same clustering base distance
cc2Clusters$SS_IsLongDB = ifelse(cc2Clusters$SS_IsDB==1&cc2Clusters$SS_Len>1e4,1,0)
cc2Clusters$SE_IsLongDB = ifelse(cc2Clusters$SE_IsDB==1&cc2Clusters$SE_Len>1e4,1,0)
cc2Clusters$ES_IsLongDB = ifelse(cc2Clusters$ES_IsDB==1&cc2Clusters$ES_Len>1e4,1,0)
cc2Clusters$EE_IsLongDB = ifelse(cc2Clusters$EE_IsDB==1&cc2Clusters$EE_Len>1e4,1,0)
cc2Clusters$SS_IsLongTI = ifelse(cc2Clusters$SS_IsTI==1&cc2Clusters$SS_Len>1e4,1,0)
cc2Clusters$SE_IsLongTI = ifelse(cc2Clusters$SE_IsTI==1&cc2Clusters$SE_Len>1e4,1,0)
cc2Clusters$ES_IsLongTI = ifelse(cc2Clusters$ES_IsTI==1&cc2Clusters$ES_Len>1e4,1,0)
cc2Clusters$EE_IsLongTI = ifelse(cc2Clusters$EE_IsTI==1&cc2Clusters$EE_Len>1e4,1,0)

cc2Clusters$SS_IsDB = ifelse(cc2Clusters$SS_IsDB==1&cc2Clusters$SS_Len<=1e4,1,0)
cc2Clusters$SE_IsDB = ifelse(cc2Clusters$SE_IsDB==1&cc2Clusters$SE_Len<=1e4,1,0)
cc2Clusters$ES_IsDB = ifelse(cc2Clusters$ES_IsDB==1&cc2Clusters$ES_Len<=1e4,1,0)
cc2Clusters$EE_IsDB = ifelse(cc2Clusters$EE_IsDB==1&cc2Clusters$EE_Len<=1e4,1,0)
cc2Clusters$SS_IsTI = ifelse(cc2Clusters$SS_IsTI==1&cc2Clusters$SS_Len<=1e4,1,0)
cc2Clusters$SE_IsTI = ifelse(cc2Clusters$SE_IsTI==1&cc2Clusters$SE_Len<=1e4,1,0)
cc2Clusters$ES_IsTI = ifelse(cc2Clusters$ES_IsTI==1&cc2Clusters$ES_Len<=1e4,1,0)
cc2Clusters$EE_IsTI = ifelse(cc2Clusters$EE_IsTI==1&cc2Clusters$EE_Len<=1e4,1,0)


# dubiously short TIs - convert these to DBs instead
dbVsTIs = cc2Clusters %>% filter((SS_IsTI==1&SS_Len<20)|(SE_IsTI==1&SE_Len<20)|(ES_IsTI==1&ES_Len<20)|(EE_IsTI==1&EE_Len<20))
View(dbVsTIs)

cc2Clusters$SS_IsDB = ifelse(cc2Clusters$SS_IsTI==1&cc2Clusters$SS_Len<20,1,ifelse(cc2Clusters$SS_IsDB==1,1,0))
cc2Clusters$SS_IsTI = ifelse(cc2Clusters$SS_IsDB==1,0,ifelse(cc2Clusters$SS_IsTI==1,1,0))
cc2Clusters$SE_IsDB = ifelse(cc2Clusters$SE_IsTI==1&cc2Clusters$SE_Len<20,1,ifelse(cc2Clusters$SE_IsDB==1,1,0))
cc2Clusters$SE_IsTI = ifelse(cc2Clusters$SE_IsDB==1,0,ifelse(cc2Clusters$SE_IsTI==1,1,0))
cc2Clusters$ES_IsDB = ifelse(cc2Clusters$ES_IsTI==1&cc2Clusters$ES_Len<20,1,ifelse(cc2Clusters$ES_IsDB==1,1,0))
cc2Clusters$ES_IsTI = ifelse(cc2Clusters$ES_IsDB==1,0,ifelse(cc2Clusters$ES_IsTI==1,1,0))
cc2Clusters$EE_IsDB = ifelse(cc2Clusters$EE_IsTI==1&cc2Clusters$EE_Len<20,1,ifelse(cc2Clusters$EE_IsDB==1,1,0))
cc2Clusters$EE_IsTI = ifelse(cc2Clusters$EE_IsDB==1,0,ifelse(cc2Clusters$EE_IsTI==1,1,0))

# note: count isn't mutually exclusive yet since haven't determined shortest lens and therefore SS_EE or SE_ES

# use the shortest length to infer what sort of link is made (ie a TI or DB)
# whether linked start-to-start or start-to-end positions

# possible combinations:
# all 4 have possible links
# 3 have possible links
# 2 have possible links but not the same pair
# 2 have possible links and the same pair
# only 1 possible link
# no possible links
cc2Clusters$SS_HasLink = ifelse(cc2Clusters$SS_IsTI==1|cc2Clusters$SS_IsDB==1,1,0)
cc2Clusters$SE_HasLink = ifelse(cc2Clusters$SE_IsTI==1|cc2Clusters$SE_IsDB==1,1,0)
cc2Clusters$ES_HasLink = ifelse(cc2Clusters$ES_IsTI==1|cc2Clusters$ES_IsDB==1,1,0)
cc2Clusters$EE_HasLink = ifelse(cc2Clusters$EE_IsTI==1|cc2Clusters$EE_IsDB==1,1,0)

cc2Clusters$SS_EE_MinLen = ifelse(cc2Clusters$SS_Len<cc2Clusters$EE_Len,cc2Clusters$SS_Len,cc2Clusters$EE_Len)
cc2Clusters$SS_EE_LenSum = cc2Clusters$SS_Len+cc2Clusters$EE_Len
cc2Clusters$SE_ES_MinLen = ifelse(cc2Clusters$SE_Len<cc2Clusters$ES_Len,cc2Clusters$SE_Len,cc2Clusters$ES_Len)
cc2Clusters$SE_ES_LenSum = cc2Clusters$SE_Len+cc2Clusters$ES_Len
View(cc2Clusters)

# first test if both pairs have links (comparing sum of lengths), then if either have a single link (comparing min of 2), then if any have a link
cc2Clusters$SS_EE = ifelse(cc2Clusters$SS_HasLink&cc2Clusters$EE_HasLink&cc2Clusters$SE_HasLink&cc2Clusters$ES_HasLink, ifelse(cc2Clusters$SS_EE_LenSum<cc2Clusters$SE_ES_LenSum,1,0),
                           ifelse((cc2Clusters$SS_HasLink|cc2Clusters$EE_HasLink)&(cc2Clusters$SE_HasLink|cc2Clusters$ES_HasLink), ifelse(cc2Clusters$SS_EE_MinLen<cc2Clusters$SE_ES_MinLen,1,0),
                                  ifelse(cc2Clusters$SS_HasLink|cc2Clusters$EE_HasLink,1,0)))

# now that side has been determined, re-evaluate the valid links
cc2Clusters$SS_IsDB = ifelse(cc2Clusters$SS_EE==1&cc2Clusters$SS_IsDB==1,1,0)
cc2Clusters$SS_IsTI = ifelse(cc2Clusters$SS_EE==1&cc2Clusters$SS_IsTI==1,1,0)
cc2Clusters$SE_IsDB = ifelse(cc2Clusters$SS_EE==0&cc2Clusters$SE_IsDB==1,1,0)
cc2Clusters$SE_IsTI = ifelse(cc2Clusters$SS_EE==0&cc2Clusters$SE_IsTI==1,1,0)
cc2Clusters$EE_IsDB = ifelse(cc2Clusters$SS_EE==1&cc2Clusters$EE_IsDB==1,1,0)
cc2Clusters$EE_IsTI = ifelse(cc2Clusters$SS_EE==1&cc2Clusters$EE_IsTI==1,1,0)
cc2Clusters$ES_IsDB = ifelse(cc2Clusters$SS_EE==0&cc2Clusters$ES_IsDB==1,1,0)
cc2Clusters$ES_IsTI = ifelse(cc2Clusters$SS_EE==0&cc2Clusters$ES_IsTI==1,1,0)

cc2Clusters$TICount = cc2Clusters$SS_IsTI+cc2Clusters$EE_IsTI+cc2Clusters$SE_IsTI+cc2Clusters$ES_IsTI
cc2Clusters$DBCount = cc2Clusters$SS_IsDB+cc2Clusters$EE_IsDB+cc2Clusters$SE_IsDB+cc2Clusters$ES_IsDB
View(cc2Clusters)

# validation
nrow(cc2Clusters %>% filter((SS_IsTI==1&SS_IsDB==1)|(SE_IsTI==1&SE_IsDB==1)|(ES_IsTI==1&ES_IsDB==1)|(EE_IsTI==1&EE_IsDB==1)))
nrow(cc2Clusters %>% filter(TICount>2|DBCount>2|(TICount+DBCount>2)))

# is this 2x templated insertions, or a TI and a DB
cc2Clusters$IsDoubleTI = ifelse(cc2Clusters$SS_EE==1,ifelse(cc2Clusters$SS_IsTI==1&cc2Clusters$EE_IsTI==1,1,0),ifelse(cc2Clusters$SE_IsTI==1&cc2Clusters$ES_IsTI==1,1,0))
cc2Clusters$IsDoubleDB = ifelse(cc2Clusters$SS_EE==1,ifelse(cc2Clusters$SS_IsDB==1&cc2Clusters$EE_IsDB==1,1,0),ifelse(cc2Clusters$SE_IsDB==1&cc2Clusters$ES_IsDB==1,1,0))
nrow(cc2Clusters %>% filter(IsDoubleDB==1))
nrow(cc2Clusters %>% filter(IsDoubleTI==1))
nrow(cc2Clusters %>% filter(DBCount==2))
nrow(cc2Clusters %>% filter(TICount==2))

cc2Clusters$S_Len = ifelse(cc2Clusters$SS_EE==1,cc2Clusters$SS_Len,cc2Clusters$SE_Len)
cc2Clusters$E_Len = ifelse(cc2Clusters$SS_EE==1,cc2Clusters$EE_Len,cc2Clusters$ES_Len)

cc2Clusters$S_IsShortTI = ifelse((cc2Clusters$SS_IsTI&cc2Clusters$SS_Len<=500)|(cc2Clusters$SE_IsTI&cc2Clusters$SE_Len<=500),1,0)
cc2Clusters$E_IsShortTI = ifelse((cc2Clusters$ES_IsTI&cc2Clusters$ES_Len<=500)|(cc2Clusters$EE_IsTI&cc2Clusters$EE_Len<=500),1,0)

# break the clusters into distinct groups for further analysis

#1 LINE
cc2g_LINE = cc2Clusters %>% filter(IsLINE==1)
cc2g_LINE$GroupType = "LINE"
nrow(cc2g_LINE)

#2 Duplicate breakend
cc2g_DupBE = cc2Clusters %>% filter(IsLINE==0&DupBECount>0)
cc2g_DupBE$GroupType = "DupBE"
nrow(cc2g_DupBE)

#3 None
cc2g_None = cc2Clusters %>% filter(IsLINE==0&DupBECount==0&TICount==0&DBCount==0)
cc2g_None$GroupType = "None"
View(cc2g_None)
nrow(cc2g_None)

cc2_Linked = cc2Clusters %>% filter(IsLINE==0&DupBECount==0&(TICount>0|DBCount>0))
nrow(cc2_Linked)

#4 2 x DSB
cc2g_DBDB = cc2_Linked %>% filter(IsDoubleDB==1)
cc2g_DBDB$GroupType = "DBx2"
nrow(cc2g_DBDB)

#5 2 x TI
cc2g_TITI = cc2_Linked %>% filter(IsDoubleTI==1)
cc2g_TITI$GroupType = "TIx2"
nrow(cc2g_TITI)

#6 DB and TI
cc2g_DBTI = cc2_Linked %>% filter(TICount==1&DBCount==1)
cc2g_DBTI$GroupType = "DBTI"
View(cc2g_DBTI)
nrow(cc2g_DBTI)
nrow(cc2_Linked %>% filter(TICount==1&DBCount==1))

#7 DB and None
cc2g_DBNone = (cc2_Linked %>% filter(DBCount==1&TICount==0))
cc2g_DBNone$GroupType = "DBNone"
View(cc2g_DBNone)
nrow(cc2g_DBNone)

#8 TI and None
cc2g_TINone = (cc2_Linked %>% filter(TICount==1&DBCount==0))
cc2g_TINone$GroupType = "TINone"
View(cc2g_TINone)
nrow(cc2g_TINone)

# check that all have been assigned
nrow(cc2Clusters)
print(nrow(cc2g_LINE)+nrow(cc2g_DupBE)+nrow(cc2g_TITI)+nrow(cc2g_DBDB)+nrow(cc2g_DBTI)+nrow(cc2g_TINone)+nrow(cc2g_DBNone)+nrow(cc2g_None))

# report on key stats for each group
cc2All = rbind(cc2g_LINE, cc2g_DupBE)
cc2All = rbind(cc2All, cc2g_None)
cc2All = rbind(cc2All, cc2g_DBDB)
cc2All = rbind(cc2All, cc2g_TITI)
cc2All = rbind(cc2All, cc2g_DBTI)
cc2All = rbind(cc2All, cc2g_DBNone)
cc2All = rbind(cc2All, cc2g_TINone)
nrow(cc2All)

# break down None BEs into further categories - Long TIs, Long DBs, same orientation and cross-arm

# for TINone and DBNone, we know which paired ends are not in the DB or TI, so can comment on what they are
# for the None (both) subset, at least one of the clustered ends must have same orientation
# so then report on what the other end has (ie longDB, longTI, crossArm or sameOrient)


# None can be SO_LongDB, SO_x2, SO_LongDB, SO_LongTI, SO_CA
cc2g_None$SS_Clustered = ifelse(cc2g_None$SS_Len<=1e4)


cc2g_None$HasLongDB = ifelse(cc2g_None$SS_IsLongDB|cc2g_None$SE_IsLongDB|cc2g_None$ES_IsLongDB|cc2g_None$EE_IsLongDB,1,0)
cc2g_None$HasLongTI = ifelse(cc2g_None$SS_IsLongTI|cc2g_None$SE_IsLongTI|cc2g_None$ES_IsLongTI|cc2g_None$EE_IsLongTI,1,0)

cc2g_None$LongLen = ifelse(cc2g_None$SS_EE==1,ifelse(cc2g_None$SE_IsLongDB==1|cc2g_None$SE_IsLongTI==1,cc2g_None$SE_Len,cc2g_None$ES_Len),
                           ifelse(cc2g_None$SS_IsLongDB==1|cc2g_None$SS_IsLongTI==1,cc2g_None$SS_Len,cc2g_None$EE_Len))

# 2xBNDs - will either both be clustered and same orientation or long DB or TI

nrow(cc2g_None)
nrow(cc2g_None %>% filter(SS_SO==1|SE_SO==1|ES_SO==1|EE_SO==1))
View(cc2g_None)

cc2g_None$SubGroupType = ifelse((cc2g_None$SS_SO==1&cc2g_None$EE_SO==1)|(cc2g_None$SE_SO==1&cc2g_None$ES_SO==1),"SOx2",
                                ifelse(cc2g_None$HasLongDB==1,'SO_LongDB',
                                       ifelse(cc2g_None$HasLongTI==1,'SO_LongTI','SO_CA')))

cc2g_None$SubGroupType = ifelse((cc2g_None$SS_CA==1&cc2g_None$EE_CA==1)|(cc2g_None$SE_CA==1&cc2g_None$ES_CA==1),"SO_CAx2",
                                ifelse((cc2g_None$SS_CA!=cc2g_None$EE_CA)|(cc2g_None$SE_CA!=cc2g_None$ES_CA),"SO_CA",
                                       ifelse(cc2g_None$HasLongTI==1,'SO_LongTI',
                                              ifelse(cc2g_None$HasLongDB==1,'SO_LongDB',
                                                     ifelse((cc2g_None$SS_SO==1&cc2g_None$EE_SO==1)|(cc2g_None$SE_SO==1&cc2g_None$ES_SO==1),"SOx2",'Unknown')))))


View(cc2g_None %>% group_by(SubGroupType,ClusterDesc) %>% summarise(Count=n()) %>% arrange(SubGroupType,ClusterDesc))

View(cc2g_None %>% group_by(SubGroupType) %>% summarise(Count=n(),AvgLongLen=round(sum(LongLen/n()),0)))

# repeat for TIs and DBs with None
cc2g_TIDBNones = rbind(cc2g_DBNone, cc2g_TINone)

cc2g_TIDBNones$HasLongDB = ifelse(cc2g_TIDBNones$SS_EE==1,ifelse(cc2g_TIDBNones$SE_IsLongDB==1|cc2g_TIDBNones$ES_IsLongDB==1,1,0),ifelse(cc2g_TIDBNones$SS_IsLongDB==1|cc2g_TIDBNones$EE_IsLongDB==1,1,0))
cc2g_TIDBNones$HasLongTI = ifelse(cc2g_TIDBNones$SS_EE==1,ifelse(cc2g_TIDBNones$SE_IsLongTI==1|cc2g_TIDBNones$ES_IsLongTI==1,1,0),ifelse(cc2g_TIDBNones$SS_IsLongTI==1|cc2g_TIDBNones$EE_IsLongTI==1,1,0))
cc2g_TIDBNones$HasOtherSO = ifelse(cc2g_TIDBNones$SS_EE==1,ifelse(cc2g_TIDBNones$SE_SO==1|cc2g_TIDBNones$ES_SO==1,1,0),ifelse(cc2g_TIDBNones$SS_SO==1|cc2g_TIDBNones$EE_SO==1,1,0))
cc2g_TIDBNones$HasOtherCA = ifelse(cc2g_TIDBNones$SS_EE==1,ifelse(cc2g_TIDBNones$SE_CA==1|cc2g_TIDBNones$ES_CA==1,1,0),ifelse(cc2g_TIDBNones$SS_CA==1|cc2g_TIDBNones$EE_CA==1,1,0))
cc2g_TIDBNones$LongLen = ifelse(cc2g_TIDBNones$SS_EE==1,ifelse(cc2g_TIDBNones$SE_IsLongDB==1|cc2g_TIDBNones$SE_IsLongTI==1,cc2g_TIDBNones$SE_Len,cc2g_TIDBNones$ES_Len),
                                ifelse(cc2g_TIDBNones$SS_IsLongDB==1|cc2g_TIDBNones$SS_IsLongTI==1,cc2g_TIDBNones$SS_Len,cc2g_TIDBNones$EE_Len))

cc2g_TIDBNones$SubGroupType = ifelse(cc2g_TIDBNones$GroupType=="TINone",
                                     ifelse(cc2g_TIDBNones$HasLongDB,"TI_LongDB",ifelse(cc2g_TIDBNones$HasLongTI,"TI_LongTI",ifelse(cc2g_TIDBNones$HasOtherCA,"TI_CA","TI_SO")),
                                            ifelse(cc2g_TIDBNones$HasLongDB,"DB_LongDB",ifelse(cc2g_TIDBNones$HasLongTI,"DB_LongTI",ifelse(cc2g_TIDBNones$HasOtherCA,"DB_CA","DB_SO")))))

View(cc2g_TIDBNones %>% group_by(SubGroupType) %>% summarise(Count=n(),AvgLongLen=round(sum(LongLen/n()),0)))


# merge these groups with a None BE pair
cc2g_Nones = rbind(cc2g_None, cc2g_TIDBNones)
nrow(cc2g_Nones)

cc2g_Nones$LongLenBucket = 2**round(log(cc2g_Nones$LongLen,2),0)

cc2g_NonesStats = (cc2g_Nones %>% group_by(GroupType,SubGroupType)
                   %>% summarise(Count=n())
                   %>% arrange(GroupType,SubGroupType))

View(cc2g_NonesStats)

cc2g_NonesDescStats = (cc2g_Nones %>% group_by(GroupType,SubGroupType,ClusterDesc)
                       %>% summarise(Count=n())
                       %>% arrange(GroupType,SubGroupType,ClusterDesc) %>% spread(ClusterDesc,Count))

View(cc2g_NonesDescStats)

cc2g_NonesLenStats = (cc2g_Nones %>% group_by(GroupType,SubGroupType,LongLenBucket)
                      %>% summarise(Count=n())
                      %>% arrange(GroupType,SubGroupType,LongLenBucket))

View(cc2g_NonesLenStats)


# overall stats
write.csv(cc2All, "~/logs/r_output/cluster2data.csv")

cc2Stats = (cc2All %>% group_by(GroupType)
            %>% summarise(Count=n(),
                          ShortTICount=sum(S_IsShortTI+E_IsShortTI),
                          StressedPercent=round(sum(IsStressed)/n(),2),
                          CrossArmPercent=round(sum(CrossChr)/n(),2)
            )
            %>% arrange(-Count))

View(cc2Stats)
write.csv(cc2Stats, "~/logs/r_output/cluster2stats.csv")

cc2StatsByDesc = (cc2All %>% group_by(GroupType,ClusterDesc)
                  %>% summarise(Count=n(),
                                ShortTICount=sum(S_IsShortTI+E_IsShortTI),
                                StressedPercent=round(sum(IsStressed)/n(),2),
                                CrossArmPercent=round(sum(CrossChr)/n(),2)
                  )
                  %>% arrange(GroupType,ClusterDesc))

View(cc2StatsByDesc)
write.csv(cc2StatsByDesc, "~/logs/r_output/cluster2statsDesc.csv")

# complete vs incomplete
cc2All$Status = ifelse(cc2All$GroupType=='LINE'|cc2All$GroupType=='DBx2'|cc2All$GroupType=='TIx2'|cc2All$GroupType=='DBTI','Complete',
                       ifelse(cc2All$GroupType=='DBNone'|cc2All$GroupType=='TINone'|cc2All$GroupType=='DupBE','Missing','UnLinked'))


cc2StatsByStatus = (cc2All %>% group_by(Status)
                    %>% summarise(Count=n(),
                                  ShortTICount=sum(S_IsShortTI+E_IsShortTI),
                                  StressedPercent=round(sum(IsStressed)/n(),2),
                                  CrossArmPercent=round(sum(CrossChr)/n(),2)
                    )
                    %>% arrange(Status))

View(cc2StatsByStatus)

cc2StatsByStatusAndDesc = (cc2All %>% group_by(Status,ClusterDesc)
                           %>% summarise(Count=n(),
                                         ShortTICount=sum(S_IsShortTI+E_IsShortTI),
                                         StressedPercent=round(sum(IsStressed)/n(),2),
                                         CrossArmPercent=round(sum(CrossChr)/n(),2)
                           )
                           %>% arrange(Status,ClusterDesc))

View(cc2StatsByStatusAndDesc)

# TI lengths
cc2gTILengths = (cc2All %>% filter(TICount>0) %>% group_by(TICount,S_IsShortTI,E_IsShortTI)
                 %>% summarise(Count=n(),
                               AvgLen1=round(sum(S_Len)/n(),0),
                               AvgLen2=round(sum(E_Len)/n(),0))
                 %>% arrange(TICount,S_IsShortTI,E_IsShortTI))
View(cc2gTILengths)
write.csv(cc2gTILengths, "~/logs/r_output/cluster2tiLengths.csv")

