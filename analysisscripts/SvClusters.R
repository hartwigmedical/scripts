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

# Clustering and Resolved Logic
sv_set_common_fields<-function(svData)
{
  svData$IsLINE = ifelse(svData$LEStart!='None'|svData$LEEnd!='None',T,F)
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
  svData$InferLinked = ifelse(svData$AsmbMatchStart=='INFER_ONLY'|svData$AsmbMatchEnd=='INFER_ONLY',T,F)
  svData$ShortTICount=ifelse(svData$LnkTypeStart=='TI'&svData$LnkLenStart<=500,0.5,0)+ifelse(svData$LnkTypeEnd=='TI'&svData$LnkLenEnd<=500,0.5,0)
  svData$DoubleTI = ifelse(svData$TICount==1,T,F)
  svData$DoubleDB = ifelse(svData$DBCount==1,T,F)
  svData$IsSpan = F
  svData$IsTrans = F
  # svData$IsSpan = ifelse(svData$TransType=='SPAN',T,F)
  # svData$IsTrans = ifelse(svData$TransType=='TRANS',T,F)
  svData$ClusterSize = ifelse(svData$ClusterCount==1,'Single',ifelse(svData$ClusterCount<=4,'Small','Large'))
  svData$IsConsistent = ifelse(svData$Consistency==0,T,F)
  svData$IsChained= (svData$ChainCount>=1)
  svData$ChainCount = ifelse(svData$ChainCount>0,svData$ChainCount,ifelse(svData$IsTI==0&svData$IsDB==0,0,1)) # set ChainCount to 1 for single link
  svData$FoldbackCount = ifelse(svData$FoldbackLenStart>=0,0.5,0)+ifelse(svData$FoldbackLenEnd>=0,0.5,0)
  svData$Filtered = F
  svData$IsFoldback = (svData$ChainCount>0 & grepl(';',svData$ChainIndex))
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



