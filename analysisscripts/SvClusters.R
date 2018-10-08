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

# Clustering and Resolved Logic
sv_set_common_fields<-function(svData)
{
  svData$IsLINE = ifelse(svData$LEStart!='false'|svData$LEEnd!='false',T,F)
  svData$IsFS = ifelse(svData$FSStart!='false'|svData$FSEnd!='false',T,F)
  svData$Length = ifelse(svData$Type=='BND'|svData$Type=='INS'|svData$ArmEnd!=svData$ArmStart, -1, svData$PosEnd-svData$PosStart)
  svData$DoubleDupBE = ifelse(svData$DupBEStart=='true'&svData$DupBEEnd=='true',T,F)
  svData$SingleDupBE = ifelse(svData$DoubleDupBE==0&(svData$DupBEStart=='true'|svData$DupBEEnd=='true'),T,F)
  svData$TICount = ifelse(svData$LnkTypeStart=='TI',0.5,0)+ifelse(svData$LnkTypeEnd=='TI',0.5,0)
  svData$DBCount = ifelse(svData$LnkTypeStart=='DB',0.5,0)+ifelse(svData$LnkTypeEnd=='DB',0.5,0)
  svData$IsSglTI = ifelse(svData$LnkTypeStart=='SGL',0.5,0)
  svData$IsTI = ifelse(svData$TICount>0,T,F)
  svData$IsDB = ifelse(svData$DBCount>0,T,F)
  svData$AssemblyLinked = ifelse(svData$AsmbMatchStart=='MATCH'|svData$AsmbMatchEnd=='MATCH',T,F)
  svData$InferLinked = ifelse(svData$AsmbMatchStart=='LINK_ONLY'|svData$AsmbMatchEnd=='LINK_ONLY',T,F)
  svData$ShortTICount=ifelse(svData$LnkTypeStart=='TI'&svData$LnkLenStart<=500,0.5,0)+ifelse(svData$LnkTypeEnd=='TI'&svData$LnkLenEnd<=500,0.5,0)
  svData$DoubleTI = ifelse(svData$TICount==1,T,F)
  svData$DoubleDB = ifelse(svData$DBCount==1,T,F)
  svData$IsSpan = ifelse(svData$TransType=='SPAN',T,F)
  svData$IsTrans = ifelse(svData$TransType=='TRANS',T,F)
  svData$ClusterSize = ifelse(svData$ClusterCount==1,'None',ifelse(svData$ClusterCount<=4,'Small','Large'))
  svData$IsConsistent = ifelse(svData$Consistency==0,T,F)
  svData$IsChained= (svData$ChainCount>0)
  svData$ChainCount = ifelse(svData$ChainCount>0,svData$ChainCount,ifelse(svData$IsTI==0&svData$IsDB==0,0,1)) # set ChainCount to 1 for single link
  svData$Filtered = F
  return (svData)
}

set_sv_stressed_state<-function(svData)
{
  # set stressed state
  svData$ArmExpStart = ifelse(svData$ArmExpStart>0,svData$ArmExpStart,0.1)
  svData$StressedPPStart = round(1 - ppois(svData$ArmCountStart - 1, svData$ArmExpStart),4)
  svData$ArmExpEnd = ifelse(svData$ArmExpEnd>0,svData$ArmExpEnd,0.1)
  svData$StressedPPEnd = round(1 - ppois(svData$ArmCountEnd - 1, svData$ArmExpEnd),4)
  svData$IsStressed = ifelse((svData$StressedPPStart<=0.001&svData$ArmCountStart>=10)|(svData$StressedPPEnd<=0.001&svData$ArmCountEnd>=10),T,F)

  return (svData)
}

set_sv_non_clustered_types<-function(svData)
{
  # RESOLVED TYPE: non-clustered SVs - simple types
  svData$SampleClusterId = paste(svData$SampleId,svData$ClusterId,sep='_')
  svData$ResolvedType = ifelse(svData$ClusterCount==1&(svData$Type=='DEL'|svData$Type=='DUP')&svData$ArmStart!=svData$ArmEnd,'NC_INVALID',svData$ResolvedType)
  svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$ClusterCount==1&svData$Type=='DEL','NC_DEL',svData$ResolvedType)
  svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$ClusterCount==1&svData$Type=='DUP','NC_DUP',svData$ResolvedType)
  svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$ClusterCount==1&svData$Type=='INS','NC_INS',svData$ResolvedType)
  svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$ClusterCount==1&svData$Type=='SGL','NC_SGL',svData$ResolvedType)

  # RESOLVED TYPE: non-clustered SVs, dubious INVs and BNDs
  dubiousNCSingleSVs = (svData %>% filter(ResolvedType=='NONE'&ClusterCount==1&(Type=='BND'|Type=='INV')&IsStressed==0&IsConsistent)
                        %>% filter(Ploidy<0.5&(AdjCNChgStart<0.5|AdjCNChgEnd<0.5)))

  svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$Id %in% dubiousNCSingleSVs$Id,'NC_INVALID',svData$ResolvedType)

  svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$ClusterCount==1,'NC_UNCLEAR',svData$ResolvedType)

  # RESOLVED TYPE: span SVs
  svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$IsSpan==1,'SPAN',svData$ResolvedType)

  return (svData)
}

get_sv_clustered_data<-function(svData)
{
  clusteredSvs = svData %>% filter(ClusterCount>1)

  allClusterData = (clusteredSvs %>% group_by(SampleId,ClusterId)
                    %>% summarise(SvCount=n(),
                                  SampleClusterId=first(SampleClusterId),
                                  ClusterCount=first(ClusterCount),
                                  ClusterDesc=first(ClusterDesc),
                                  Consistency=first(Consistency),
                                  TICount=sum(TICount),
                                  DBCount=sum(DBCount),
                                  TransCount=sum(IsTrans),
                                  ShortTICount=sum(ShortTICount),
                                  SpanCount=sum(IsSpan),
                                  DoubleDupBECount=sum(DoubleDupBE),
                                  SingleDupBECount=sum(SingleDupBE),
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
                                  StressedCount=sum(IsStressed),
                                  LineCount=sum(IsLINE),
                                  FilteredCount=sum(Filtered))
                    %>% arrange(SampleId,ClusterId))

  return (allClusterData)
}

set_sv_line_types<-function(svData, allClusterData)
{
  # RESOLVED TYPE: LINE and SVs in LINE clusters
  lineClusters = allClusterData %>% filter(LineCount>0)
  svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$SampleClusterId %in% lineClusters$SampleClusterId,'LINE',svData$ResolvedType)
  svData$ResolvedType = ifelse(svData$ResolvedType=='LINE'&svData$IsLINE,'LINE_CLUST',svData$ResolvedType)

  return (svData)
}

set_sv_clustered_types<-function(svData, allClusterData)
{
  # first exclude LINE, assuming it has already been resolved out
  clusterData = allClusterData %>% filter(LineCount==0)
  clusterData$IsStressed = ifelse(clusterData$StressedCount>0,T,F)
  clusterData$ResolvedType = 'NONE'

  # limit to cluster sizes and less since linking logic cuts out beyond that
  clusterData = clusterData %>% filter(ClusterCount<=100)

  # RESOLVED TYPE: reciprocal inversions & translocations
  svData$IsRecipInv = ifelse(svData$ClusterCount==2&svData$Type!='DEL'&svData$Type!='DUP'&svData$DoubleDB&svData$LnkSvStart==svData$LnkSvEnd,T,F)

  reciprocals = clusterData %>% filter(ClusterCount==2&(BndCount==2|InvCount==2)&DBCount==2&IsConsistent)

  svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$Type=='BND'&svData$SampleClusterId %in% reciprocals$SampleClusterId,'RECIP_TRANS',svData$ResolvedType)
  svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$Type!='BND'&svData$SampleClusterId %in% reciprocals$SampleClusterId,'RECIP_INV',svData$ResolvedType)
  clusterData$ResolvedType = ifelse(clusterData$ResolvedType=='NONE'&clusterData$BndCount==2&clusterData$SampleClusterId %in% reciprocals$SampleClusterId,'RECIP_TRANS',clusterData$ResolvedType)
  clusterData$ResolvedType = ifelse(clusterData$ResolvedType=='NONE'&clusterData$InvCount==2&clusterData$SampleClusterId %in% reciprocals$SampleClusterId,'RECIP_INV',clusterData$ResolvedType)

  # RESOLVED TYPE: simple TI from SV pairs
  simpleTIs = clusterData %>% filter(ResolvedType=='NONE'&ClusterCount==2&TICount>=1&Consistency==0)
  svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$SampleClusterId %in% simpleTIs$SampleClusterId,'TI_PAIR',svData$ResolvedType)
  clusterData$ResolvedType = ifelse(clusterData$ResolvedType=='NONE'&clusterData$SampleClusterId %in% simpleTIs$SampleClusterId,'TI_PAIR',clusterData$ResolvedType)

  # complex clusters
  clusterData$ClusterBucket = ifelse(clusterData$ClusterCount<=10,clusterData$ClusterCount,ifelse(clusterData$ClusterCount<=50,round(clusterData$ClusterCount/5)*5,round(clusterData$ClusterCount/10)*10))
  clusterData$ShortTIPerc = round(clusterData$ShortTICount/clusterData$SvCount,2)
  clusterData$DBPerc = round(clusterData$DBCount/clusterData$SvCount,2)
  clusterData$LinkPerc = round((clusterData$TICount+clusterData$DBCount+clusterData$SpanCount)/clusterData$SvCount,2)

  # RESOLVED TYPE: all linked DBs
  shortTIDBChains = clusterData %>% filter(ClusterCount>=2&(ShortTICount+DBCount+SpanCount)==ClusterCount&Consistency==0)
  svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$SampleClusterId %in% shortTIDBChains$SampleClusterId,'DB_CHAIN',svData$ResolvedType)
  clusterData$ResolvedType = ifelse(clusterData$ResolvedType=='NONE'&clusterData$SampleClusterId %in% shortTIDBChains$SampleClusterId,'DB_CHAIN',clusterData$ResolvedType)

  # RESOLVED TYPE: all TIs
  tiChains = clusterData %>% filter(ClusterCount>=2&(TICount+SpanCount)>=(ClusterCount-1)&DBCount==0&Consistency==0)
  svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$SampleClusterId %in% tiChains$SampleClusterId,'TI_CHAIN',svData$ResolvedType)
  clusterData$ResolvedType = ifelse(clusterData$ResolvedType=='NONE'&clusterData$SampleClusterId %in% tiChains$SampleClusterId,'TI_CHAIN',clusterData$ResolvedType)

  # RESOLVED TYPE: all chained SVs
  allChainedClusters = clusterData %>% filter(ResolvedType=='NONE'&IsConsistent&DoubleDupBECount==SpanCount&(DBCount+TICount+SpanCount)>=(ClusterCount-1))
  svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$SampleClusterId %in% allChainedClusters$SampleClusterId,'COMPLEX_CHAIN',svData$ResolvedType)
  clusterData$ResolvedType = ifelse(clusterData$ResolvedType=='NONE'&clusterData$SampleClusterId %in% allChainedClusters$SampleClusterId,'COMPLEX_CHAIN',clusterData$ResolvedType)

  # RESOLVED TYPE: partially chained SVs
  partiallyChainedClusters = clusterData %>% filter(ResolvedType=='NONE'&ChainCount>0&ChainCount<ClusterCount)
  svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$SampleClusterId %in% partiallyChainedClusters$SampleClusterId,'PARTIAL_CHAIN',svData$ResolvedType)
  clusterData$ResolvedType = ifelse(clusterData$ResolvedType=='NONE'&clusterData$SampleClusterId %in% partiallyChainedClusters$SampleClusterId,'PARTIAL_CHAIN',clusterData$ResolvedType)

  # RESOLVED TYPE: potentially chained SVs
  potentiallyChainedClusters = clusterData %>% filter(ResolvedType=='NONE'&TICount >= 0.6 * ClusterCount)
  svData$ResolvedType = ifelse(svData$ResolvedType=='NONE'&svData$SampleClusterId %in% potentiallyChainedClusters$SampleClusterId,'POTENTIAL_CHAIN',svData$ResolvedType)
  clusterData$ResolvedType = ifelse(clusterData$ResolvedType=='NONE'&clusterData$SampleClusterId %in% potentiallyChainedClusters$SampleClusterId,'POTENTIAL_CHAIN',clusterData$ResolvedType)

  return (svData)
}


rm(svData)

# SV data file
svData = read.csv('~/logs/CLUSTER_V25.csv')


nrow(svData)
View(head(svData,100))

svAllData = svData

# filter out multiple biopsy (approximately)
svData = svData %>% filter(!grepl("DRUP", SampleId)&!grepl("TIII", SampleId)&!grepl("TII", SampleId))
nrow(svData)

# FILTER FOR PONCount <2 for all subsequent analyses - no longer required since done already
View(svData %>% filter(PONCount>=2))
svData = svData %>% filter(PONCount<2)


svGridssMatched = svData %>% filter(SampleId %in% gridssSamples$SampleId)
nrow(svGridssMatched)
nrow(svGridssMatched %>% group_by(SampleId) %>% count())

# optimised version of verbose classification below:
svData = svGridssMatched

svData = sv_set_common_fields(svData)
svData = set_sv_stressed_state(svData)

# allocation of known / resolved types
svData$ResolvedType = 'NONE'
svData = set_sv_non_clustered_types(svData)
allClusterData = get_sv_clustered_data(svData)
svData = set_sv_line_types(svData, allClusterData)
svData = set_sv_clustered_types(svData, allClusterData)

# View(svData)

svResolvedSummary = (svData %>% group_by(ResolvedType,IsStressed,ClusterSize)
                     %>% summarise(Count=n()) %>% arrange(ResolvedType,IsStressed,ClusterSize))

View(svResolvedSummary)








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




