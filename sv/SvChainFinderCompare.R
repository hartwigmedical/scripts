

## ChainFinder comparison
cfData = read.csv('~/data/sv/LNX_CHAIN_FINDER_SVS.csv')
nrow(cfData)
View(cfData)
View(cfData %>% group_by(SampleId) %>% count)
nrow(cfData %>% group_by(SampleId) %>% count) # 1587

cfSamples = cfData %>% group_by(SampleId) %>% count
View(cfSamples)
nrow(cfSamples) # 3686 samples

## link to common Linx data
cfSvData = read.csv('~/data/sv/chainfinder/LNX_SVS.csv')
cfSvData = cfSvData %>% filter(SampleId %in% cfSamples$SampleId)
nrow(cfSvData)

cfSamplesWithChains = cfData %>% filter(CfChainId>=0) %>% group_by(SampleId) %>% count
nrow(cfSamplesWithChains) # 1587

cfSvData = cfSvData %>% filter(SampleId %in% cfSamplesWithChains$SampleId)
nrow(cfSvData)


sampleCounts = cfSvData %>% group_by(SampleId) %>% summarise(TotalSVs=n(),
                                                             SimpleSVs=sum(ClusterCount==1&(Type=='DUP'|Type=='INS'|Type=='DEL')&(PosEnd-PosStart<1e6)),
                                                             ClusteredSVs=sum(ClusterCount>1))
View(sampleCounts)

cfClusters = read.csv('~/data/sv/chainfinder/LNX_CLUSTERS.csv')
cfClusters = cfClusters %>% filter(SampleId %in% cfSamples$SampleId)


# 1a. Percentage of SVs in a chain/cluster 
cfSampleStats = cfData %>% filter(Type!='INF'&Type!='SGL') %>% group_by(SampleId) %>% 
  summarise(ClusteredSVs=n(),
            LinxClustered=sum(ClusterCount>1),
            CfChained=sum(ChainId>=1))

cfSampleStats = merge(cfSampleStats,sampleCounts %>% select(-ClusteredSVs),by='SampleId',all.x=T)

cfSampleStats = cfSampleStats %>% mutate(LinxClusterPerc=round(LinxClustered/TotalSVs,3),
                                         LinxNonSimpleClusterPerc=round(LinxClustered/(TotalSVs-SimpleSVs),3),
                                         CfChainedPerc=round(CfChained/TotalSVs,3),
                                         CfNonSimpleChainedPerc=round(CfChained/(TotalSVs-SimpleSVs),3))

View(cfSampleStats)

print(ggplot(cfSampleStats %>% filter(ClusteredSVs>50), aes(x=LinxClusterPerc, y=CfChainedPerc))
           + geom_point())

print(ggplot(cfSampleStats %>% filter(ClusteredSVs>50), aes(x=LinxNonSimpleClusterPerc, y=CfNonSimpleChainedPerc))
      + geom_point())


# 1b. Linx clusters more than CF regardless of total sample SV count
print(ggplot(cfSampleStats %>% mutate(NonSimpleSVs=TotalSVs-SimpleSVs,
                                      SampleMB=ifelse(NonSimpleSVs<50,'0: <50 SVs',ifelse(NonSimpleSVs<200,'1: 50-200 SVs','2: >200 SVs'))),
             aes(x=LinxNonSimpleClusterPerc, y=CfNonSimpleChainedPerc))
      + geom_point()
      + facet_wrap(~SampleMB)
      + labs(title = 'Percentage of non-simple SVs clustered by Linx and ChainFinder', x = 'Linx Clustered %', y = 'CF chained %'))


# 3. Linx clusters due to proximity much more than CF, followed by LOH
View(cfData %>% filter(CfChainId<0&ClusterCount>2&Type!='INF'&Type!='SGL'&Type!='INS') %>% 
       group_by(ClusterReason) %>% summarise(Total=n(),
                                             DELs=sum(Type=='DEL'),
                                             DUPs=sum(Type=='DUP'),
                                             INVs=sum(Type=='INV'),
                                             BNDs=sum(Type=='BND')))

# many of these unchained SVs are very close
View(cfData %>% filter(CfChainId<0&ClusterCount>2&Type!='INF'&Type!='SGL'&Type!='INS'&ClusterReason=='PROXIMITY') %>% 
      group_by(ProxDistance=2**round(log(ProxDistance,2))) %>% count)

# is it due to having only a few SVs on the arm?
armStats = rbind(cfSvData %>% filter(Type!='SGL'&Type!='INF') %>% select(SampleId,Chr=ChrStart,Arm=ArmStart,Pos=PosStart),
                 cfSvData %>% filter(Type!='SGL'&Type!='INF') %>% select(SampleId,Chr=ChrEnd,Arm=ArmEnd,Pos=PosEnd))
armStats = armStats %>% group_by(SampleId,Chr,Arm) %>% summarise(ArmBreakCount=n()) %>% ungroup()
View(armStats)
View(armStats %>% select(SampleId,ChrStart=Chr,ArmStart=Arm,ArmBreakCountStart=ArmBreakCount))

cfArmData = merge(cfData,armStats %>% select(SampleId,ChrStart=Chr,ArmStart=Arm,ArmBreakCountStart=ArmBreakCount),by=c('SampleId','ChrStart','ArmStart'),all.x=T)
cfArmData = merge(cfArmData,armStats %>% select(SampleId,ChrEnd=Chr,ArmEnd=Arm,ArmBreakCountEnd=ArmBreakCount),by=c('SampleId','ChrEnd','ArmEnd'),all.x=T)
View(cfArmData)

View(cfArmData %>% group_by(ArmSvCount=2**round(log(pmax(ArmBreakCountStart,ArmBreakCountEnd),2))) %>% count)

# doesn't appear to be
View(cfArmData %>% filter(CfChainId<0&ClusterCount>2&Type!='INF'&Type!='SGL'&Type!='INS'&ClusterReason=='PROXIMITY') %>% 
       group_by(ProximityDistance=2**round(log(ProximityDistance,2)),
                ArmSvCount=2**round(log(pmax(ArmBreakCountStart,ArmBreakCountEnd),2))) %>% count %>% spread(ArmSvCount,n,fill=0))

# Linx clusters SVs much further away
View(cfData %>% filter(CfProxDistance>1e5))

clusterDistances = cfData %>% filter(ClusterCount>3&Type!='INF'&Type!='SGL') %>% 
  group_by(ClusterReason,ClusterDistance=ifelse(ProxDistance==0,1,2**round(log(ProxDistance,2)))) %>% 
  summarise(LinxCount=n(),CfCount=sum(CfChainId>0))

View(clusterDistances)

ignoreCRs = c('LOH_CHAIN','OVERLAP_FOLDBACKS','CONSEC_BREAKS')

print(ggplot(clusterDistances %>% filter(ClusterReason!='PROXIMITY') %>% filter(!(ClusterReason %in% ignoreCRs)),aes(x=ClusterDistance))
      + geom_line(aes(y=LinxCount,color='LinxCount'))
      + geom_line(aes(y=CfCount,color='CfCount'))
      + scale_x_log10()
      + facet_wrap(~ClusterReason)
      + labs(title = "Clustering Reasons in Linx vs Chained by CF"))

print(ggplot(clusterDistances %>% filter(ClusterReason=='PROXIMITY'),aes(x=ClusterDistance))
      + geom_line(aes(y=LinxCount,color='LinxCount'))
      + geom_line(aes(y=CfCount,color='CfCount'))
      + labs(title = "Clustering by Proximity, Linx vs ChainFinder"))

combinedClusterDistances = rbind(cfData %>% filter(CfChainId>0) %>% mutate(Source='CF',ClusterReason='PROXIMITY') %>% 
                                   select(SampleId,Source,SvId,ClusterReason,ClusterDistance=CfProxDistance),
                                 cfData %>% filter(ClusterCount>1) %>% mutate(Source='LINX') %>% select(SampleId,Source,SvId,ClusterReason,ClusterDistance=ProxDistance))

View(combinedClusterDistances)

# much longer tail for Linx's clustering
print(ggplot(combinedClusterDistances %>% filter(ClusterDistance>0) %>% group_by(Source,DistanceBucket=2**round(log(ClusterDistance,2))) %>% count, 
             aes(x=DistanceBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~Source)
      + labs(title = "Clustering Distances Linx vs CF"))

# variety of reason for long-distance clustering in Linx
View(combinedClusterDistances %>% group_by(Source,ClusterReason) %>% count)
print(ggplot(combinedClusterDistances %>% 
               filter(ClusterDistance>0&Source=='LINX'&ClusterReason!='PROXIMITY'&ClusterReason!='LOH_CHAIN'&ClusterReason!='OVERLAP_FOLDBACKS'&ClusterReason!='CONSEC_BREAKS') %>% 
               group_by(ClusterReason,DistanceBucket=2**round(log(ClusterDistance,2))) %>% count, 
             aes(x=DistanceBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~ClusterReason)
      + labs(title = "Clustering Distances in Linx by Clustering Reason"))



# Size by type


# CF only clusters 3 or more SVs
View(cfData %>% filter(CfChainCount==1|CfChainCount==2)) # are these a mistake where it misses an SV?

# Cluster Summary stats
cfSampleSummary = cfData %>% filter(CfChainId>=0) %>% group_by(SampleId,ClusterId,ClusterCount,SglCount,ResolvedType,CfChainId,CfChainCount) %>% 
  summarise(SharedSvCount=n()) %>% mutate(LinxSvExcess=ClusterCount-SglCount-SharedSvCount,
                                          CfSvExcess=CfChainCount-SharedSvCount,
                                          MatchType=ifelse(LinxSvExcess==0&CfSvExcess==0,'EXACT',
                                                    ifelse(LinxSvExcess>0&CfSvExcess>0,'BOTH_EXCESS',
                                                    ifelse(LinxSvExcess>0,'LINX_EXCESS','CF_EXCESS'))))

View(cfSampleSummary)

## about 1/3 of clusters match exactly, and 1/3 have an excess in one or the other tool
View(cfSampleSummary %>% group_by(MatchType) %>% count)

## when clusters with SGLs are excluded, Linx has less than 1/2 the clusters with excess as CF
View(cfSampleSummary %>% filter(SglCount==0) %>% group_by(MatchType) %>% count)

## Linx typically forms larger clusters
View(cfSampleSummary %>% filter(MatchType=='LINX_EXCESS'|MatchType=='CF_EXCESS') %>%
       group_by(MatchType,ClusterExcess=ifelse(MatchType=='LINX_EXCESS',2**round(log(abs(LinxSvExcess),2)),2**round(log(abs(CfSvExcess),2)))) %>% 
       count %>% spread(MatchType,n,fill=0))


# Proximity differences
View(cfData %>% filter(CfChainId<0&ClusterCount>2&ClusterReason=='Prox'&Type!='INF'&Type!='SGL') %>% group_by(ClusterDistance=2**round(log(ClusterDistance,2))) %>% count)

View(cfData %>% filter(CfChainId<0&ClusterCount>2&ClusterReason=='Prox'))

# Large clusters with lots of Linx-only proxity clustering
View(cfData %>% filter(CfChainId<0&ClusterCount>2&ClusterReason=='Prox') %>% 
       filter(SampleId %in% cfSamples$SampleId) %>%
       group_by(SampleId,ClusterId,ClusterCount) %>%
       summarise(ProxSVs=n()) %>%
       mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_')) %>% mutate(SimpleChained=SampleClusterId %in% simpleChainedClusters$SampleClusterId))

simpleChainedClusters = cfClusters %>% filter(FullyChained=='true'&ChainCount==1&Replication=='false'&SglCount==0&InfCount==0) %>% 
  mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))
View(simpleChainedClusters)

## Characteristics of SVs for CF clustering when Linx doesn't
View(cfData %>% filter(ClusterCount==1&CfChainId>=0))

# about 1/3 of all Linx unclustered SVs were clustered, mostly due to proximity, but then dissolved
View(cfData %>% filter(ClusterCount==1&CfChainId>=0) %>% group_by(ClusterReason,Type) %>% count %>% spread(Type,n,fill=0))


#####
## Heatmap of cluster vs chain size
View(cfData %>% filter(CfChainId>0))

View(cfData %>% filter((ClusterCount>=3|CfChainCount>=3)&(CfChainCount==-1|CfChainCount>=3)) %>% 
       group_by(ClusterSize=ifelse(ClusterCount<=3,ClusterCount,2**round(log(ClusterCount,2))),
                ChainSize=ifelse(is.na(CfChainCount),0,ifelse(CfChainCount<=3,CfChainCount,2**round(log(CfChainCount,2))))) %>% count %>% spread(ChainSize,n,fill=0))

get_cf_group_size<-function(clusterSize,useLog=F)
{
  if(is.na(clusterSize))
    return (as.character('<3'))
  if(useLog)
    return (as.character(2**round(log(clusterSize,2))))
  
  if(clusterSize<=2)
    return ('<3')
  else if(clusterSize<=4)
    return (as.character(clusterSize))
  else if(clusterSize<=8)
    return (as.character('5-8'))
  else if(clusterSize<=16)
    return (as.character('9-16'))
  else if(clusterSize<=32)
    return (as.character('17-32'))
  else if(clusterSize<=64)
    return (as.character('33-64'))
  else
    return (as.character('>64'))
}

get_cf_group_index<-function(clusterSize)
{
  if(clusterSize=='<3')
    return (0)
  else if(clusterSize=='3')
    return (1)
  else if(clusterSize=='4')
    return (2)
  else if(clusterSize=='5-8')
    return (3)
  else if(clusterSize=='9-16')
    return (4)
  else if(clusterSize=='17-32')
    return (5)
  else if(clusterSize=='33-64')
    return (6)
  else
    return (7)
}

# testing
print(get_cf_group_size(-1))
print(get_cf_group_size(0))
print(get_cf_group_size(1))
print(get_cf_group_size(-1))

clusterChainCmp = cfData %>% filter(CfChainId==-1|CfChainCount>=3) %>% filter(ClusterCount>=3|CfChainId>0)

clusterChainCmp$LinxClusterSize=apply(clusterChainCmp[,c('ClusterCount'),drop=F], 1, function(x) get_cf_group_size(x[1]))
clusterChainCmp$CfChainSize=apply(clusterChainCmp[,c('CfChainCount'),drop=F], 1, function(x) get_cf_group_size(x[1]))
View(clusterChainCmp)

View(clusterChainCmp %>% group_by(LinxClusterSize,CfChainSize) %>% count %>% spread(CfChainSize, n,fill=0))

View(clusterChainCmp %>% group_by(LinxClusterSize=2**round(log(ClusterCount,2)),
                                  CfChainSize=ifelse(CfChainCount<=0,0,2**round(log(CfChainCount,2)))) %>% count %>% spread(CfChainSize, n,fill=0))

clusterChainSummary = clusterChainCmp %>% group_by(LinxClusterSize,CfChainSize) %>% summarise(Count=n(),CfChainCount=first(CfChainCount),LinxClusterCount=first(ClusterCount))
View(clusterChainSummary)

clusterChainSummary2 = clusterChainSummary %>% select(LinxClusterSize,CfChainSize,Count) %>%  spread(CfChainSize,Count,fill=0) %>% gather('CfChainSize','Count',2:9)
clusterChainSummary2$LinxSizeIndex=apply(clusterChainSummary2[,c('LinxClusterSize'),drop=F], 1, function(x) get_cf_group_index(x[1]))
clusterChainSummary2$CfSizeIndex=apply( clusterChainSummary2[,c('CfChainSize'),drop=F], 1, function(x) get_cf_group_index(x[1]))
View(clusterChainSummary2)

## Final Plot
print(ggplot(clusterChainSummary2, aes(x=reorder(LinxClusterSize,LinxSizeIndex),y=reorder(CfChainSize,CfSizeIndex))) 
      + geom_tile(aes(fill=Count),stat="identity") 
      + geom_text(aes(label=Count))
      + scale_fill_gradient(low="white",high="steelblue")
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), panel.background = element_blank())
      + theme(legend.position = "none")
      + labs(title='Linx vs ChainFinder Clustering', x='Linx Cluster Size',y='ChainFinder Chain Size'))


print(ggplot(clusterChainSummary, aes(x=reorder(LinxClusterSize,LinxClusterCount),y=reorder(CfChainSize,CfChainCount))) 
      + geom_hex(aes(fill=Count),stat="identity") 
      + geom_text(aes(label=Count))
      + scale_fill_gradient(low="white",high="steelblue")
      + labs(title='Linx vs ChainFinder Clustering', x='Linx Cluster Size',y='ChainFinder Chain Size'))

clusterChainCmp = merge(clusterChainCmp,clusterChainSummary %>% select(LinxClusterSize,CfChainSize,Count),by=c('LinxClusterSize','CfChainSize'),all.x=T)
View(clusterChainCmp %>% select(ClusterCount,CfChainCount,LinxClusterSize,CfChainSize,Count,everything()))

print(ggplot(clusterChainCmp, aes(x=pmin(ClusterCount,64),y=pmin(CfChainCount,64))) 
      + geom_hex() 
      #+ scale_x_log10()
      #+ scale_y_log10()
      + scale_fill_gradient(low="white",high="steelblue",trans="log10")
      + labs(title='Linx vs ChainFinder Clustering', x='Linx Cluster Size',y='ChainFinder Chain Size'))





#####
## Combined Group Data
#cfGroupData = read.csv('~/data/sv/chainfinder/LNX_CHAIN_FINDER_GROUPS.csv')
cfGroupData = read.csv('~/data/sv/LNX_CHAIN_FINDER_GROUPS.csv')

View(cfGroupData)
nrow(cfGroupData) # about 6100 in prod, from ~1540 samples
View(cfGroupData %>% group_by(SampleId) %>% count)

# 1. Most groups involve just 1 cluster and 1 chain, following by 1 extra of one or the other (ie 2 chains vs 1 cluster, or 1 cluster and 2 chains)
# to make the comparison fair, single and cluster-2s are ignored
cfGroupData = cfGroupData %>% mutate(Clusters3Plus=Clusters-SimpleSVs-Cluster2s)

View(cfGroupData %>% 
       # filter(SGLs==0) %>% 
       group_by(Clusters3Plus,Chains) %>% count %>% spread(Chains,n,fill=0))

# CF chaining lots of single-SVs and cluster-2s
View(cfGroupData %>% filter(Clusters3Plus==0))

# by SV count
View(cfGroupData %>% 
       filter(SGLs==0) %>% 
       group_by(Clusters3Plus,Chains) %>% summarise(SvCount=sum(TotalSVs)) %>% spread(Chains,SvCount,fill=0))

View(cfGroupData %>% 
       # filter(SGLs==0) %>% 
       group_by(Clusters3Plus,Chains) %>% summarise(SvCount=sum(TotalSVs)) %>% spread(Chains,SvCount,fill=0))

# Under-clustering by CF as Linx' clusters get bigger


View(cfGroupData %>% filter(Clusters==1&Chains==1&SGLs==0))
View(cfGroupData %>% filter(Clusters-SimpleSVs==1&Chains==1&SGLs==0))

View(cfGroupData %>% filter(Clusters-SimpleSVs==1&Chains==1&SGLs==0) %>% filter(ChainedSvCount>ClusteredSvCount))

singleCC = cfGroupData %>% filter(Clusters3Plus==1&Chains==1&SGLs==0) %>% 
  mutate(Excess=ifelse(ClusteredSvCount/TotalSVs>0.9&ChainedSvCount/TotalSVs>0.9,'MATCH',
                       ifelse(ClusteredSvCount>SharedSVs&ChainedSvCount>SharedSVs&ClusteredSvCount/TotalSVs>0.75&ChainedSvCount/TotalSVs>0.75,'BOTH',
                              ifelse(ClusteredSvCount>SharedSVs,'LINX','CF'))))

View(singleCC %>% group_by(Excess) %>% summarise(Count=n(),ExcessSVs=sum(abs(ClusteredSvCount-ChainedSvCount))))

## CF clustering where Linx keeps simple SVs separate - many of the larger ones are in fragile sites
View(cfGroupData %>% filter(Clusters-SimpleSVs==0))


View(cfData %>% filter(SampleId=='CPCT02110049T'&OverlapGroupId==4))

chainedSVs = cfData %>% filter(ClusterCount==1)
chainedSVs = merge(chainedSVs,cfSvData %>% select(SampleId,SvId=Id,FSStart,FSEnd),by=c('SampleId','SvId'),all.x=T)
View(chainedSVs %>% filter(Type!='INS'&!is.na(FSStart)) %>% group_by(Type,InFragileSite=(FSStart=='true'|FSEnd=='true')) %>% count %>% spread(InFragileSite,n,fill=0))

# 2 or more clusters in a single chain
cfOverClustering = cfGroupData %>% filter(SimpleSVs==0&Clusters>Chains)
View(cfOverClustering)
View(cfOverClustering)
View(cfOverClustering %>% summarise(Line=sum(grepl('LINE',ResolvedTypes)),
                                    Linex2=sum(grepl('LINE;LINE',ResolvedTypes)),
                                    RecipInv=sum(grepl('RECIP_INV',ResolvedTypes)),
                                    RecipTrans=sum(grepl('RECIP_TRANS',ResolvedTypes)),
                                    SimpleGroup=sum(grepl('SIMPLE_GRP',ResolvedTypes)),
                                    RecipTransDups=sum(grepl('RECIP_TRANS_DUPS',ResolvedTypes)),
                                    SynthDel=sum(grepl(';DEL',ResolvedTypes)|grepl('DEL;',ResolvedTypes)),
                                    SynthDup=sum(grepl(';DUP',ResolvedTypes)|grepl('DUP;',ResolvedTypes)),
                                    DelTI=sum(grepl('DEL_TI',ResolvedTypes))))

View(cfGroupData %>% filter(Clusters>Chains) %>% summarise(SimpleSVs=sum(SimpleSVs>0),
                                                           Cluster2s=sum(Cluster2s>0),
                                    Line=sum(grepl('LINE',ResolvedTypes)),
                                    Linex2=sum(grepl('LINE;LINE',ResolvedTypes)),
                                    RecipInv=sum(grepl('RECIP_INV',ResolvedTypes)),
                                    RecipTrans=sum(grepl('RECIP_TRANS',ResolvedTypes)),
                                    SimpleGroup=sum(grepl('SIMPLE_GRP',ResolvedTypes)),
                                    RecipTransDups=sum(grepl('RECIP_TRANS_DUPS',ResolvedTypes)),
                                    SynthDel=sum(grepl(';DEL',ResolvedTypes)|grepl('DEL;',ResolvedTypes)),
                                    SynthDup=sum(grepl(';DUP',ResolvedTypes)|grepl('DUP;',ResolvedTypes)),
                                    DelTI=sum(grepl('DEL_TI',ResolvedTypes))))







#######

## ChainFinder debug and experiments

cfCompareSample = read.csv('~/logs/LNX_CHAIN_FINDER_SVS.csv')
View(cfCompareSample)

View(cfCompareSample %>% filter(CfSvId==783898|CfSvId==783814))


cfCompare100 = read.csv('~/data/sv/chainfinder/LNX_CHAIN_FINDER_SVS.csv')
View(cfCompare100)

View(cfCompare %>% filter(ClusterCount>100&CfChainId>=0))


View(cfData %>% filter(SampleId=='CPCT02120063T'&ChainId>0))
View(cfData %>% filter(SampleId=='CPCT02120063T'))


# longest chained distance
View(cfCompare100 %>% filter(CfChainId>=0))

cfGroupDataSample = read.csv('~/logs/LNX_CHAIN_FINDER_GROUPS.csv')
View(cfGroupDataSample)

View(cfGroupData %>% filter(SampleId=='CPCT02120063T'))



# Local vs Remote
createClusterArmStats <- function(svData) 
{
  rbind(svData %>% mutate(Arm=ArmStart,Chr=ChrStart),svData %>% filter(ChrEnd!=0) %>% mutate(Arm=ArmEnd,Chr=ChrEnd)) %>%
    group_by(SampleId,ClusterId) %>% mutate(clusterSGLCount=sum(ChrEnd==0)) %>% ungroup() %>%
    group_by(SampleId,ClusterId,Arm,Chr,ResolvedType,ClusterBreakends=ClusterCount*2-clusterSGLCount) %>% 
    summarise(LocalCount=sum(ArmEnd==ArmStart&ChrEnd==ChrStart),CrossArmCount=sum(ArmEnd!=ArmStart&ChrEnd==ChrStart),SGLCount=sum(ChrEnd==0),RemoteCount=sum((ArmEnd!=ArmStart|ChrEnd!=ChrStart)&ChrEnd!=0)) %>%
    mutate(expectedLocalProportion = round((LocalCount+RemoteCount+SGLCount)/ClusterBreakends,4),actualLocalProportion=(round((LocalCount+SGLCount)/(LocalCount+RemoteCount+SGLCount),4)),diff=actualLocalProportion-expectedLocalProportion) %>%
    mutate(pValue=pmin(ppois(LocalCount+SGLCount,(LocalCount+RemoteCount+SGLCount)*expectedLocalProportion,FALSE,FALSE),ppois(LocalCount+SGLCount,(LocalCount+RemoteCount+SGLCount)*expectedLocalProportion,TRUE,FALSE)))
}

View(svData %>% filter(SampleId=='CPCT02010022T'))
View(svData %>% filter(SampleId=='CPCT02010022T') %>% group_by(ClusterId,ClusterCount,ResolvedType) %>% count)


clusters = read.csv('~/data/sv/chainfinder/LNX_CLUSTERS_PROD.csv')
clusters = clusters %>% filter(SampleId %in% cfSampleIds$SampleId)
View(clusters)

View(clusters %>% filter(SampleId %in% cfSamplesWithChains$SampleId) %>% filter(ClusterCount>100))

View(clusters %>% filter(SampleId %in% cfSamplesWithChains$SampleId) %>% group_by(SampleId) %>% summarise(SvCount=sum(ClusterCount)))





# CF Under-Clustering
View(cfData %>% filter(CfChainId<0) %>% separate(ClusterReason,c('ClusterReason','ClusterSvId',sep='_')) %>%
       group_by(Type,ClusterReason) %>% count %>% spread(Type,n,fill=0))

View(cfData %>% filter(CfChainId<0) %>% separate(ClusterReason,c('ClusterReason','ClusterSvId',sep='_')) %>%
       group_by(ClusterReason) %>% count)



specSvData = svData %>% filter(SampleId=='CPCT02010452T')
specSvData = svData %>% filter(SampleId=='CPCT02330070T')
cfSpecCompare = merge(cfData,specSvData %>% select(SampleId,SvId=Id,ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd,DBLenStart,DBLenEnd),by=c('SampleId','SvId'),all.x=T)
View(cfSpecCompare)

View(cfSpecCompare %>% filter(ClusterId %in% c(167,169,172,175,176,208,209)) %>% select(ClusterId,Type,ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd,Ploidy,everything()))


svSampleSummary = svData %>% group_by(SampleId) %>% summarise(SvCount=n(),
                                                              SglCount=sum(Type=='SGL'|Type=='INF'),
                                                              MaxPloidy=max(Ploidy))

View(svSampleSummary %>% filter(SvCount>100&SvCount<500&SglCount<15&MaxPloidy<20) %>% arrange(SglCount))

write.csv(head(svSampleSummary %>% filter(SvCount>100&SvCount<500&SglCount<15&MaxPloidy<20) %>% arrange(SglCount),100),
          '~/data/sv/chainfinder/cf_samples.txt',row.names = F,quote = F)

cfSampleIds = read.csv('~/data/sv/chainfinder/cf_samples.txt')
write.csv(cfSampleIds %>% select(SampleId),'~/data/sv/chainfinder/cf_sample_ids.csv',row.names = F,quote = F)




# CF proximity vs Linx clustering distances
clusterDistances = cfData %>% filter(CfChainId>0&CfProxDistance>=0) %>%
  mutate(LinxDistance=ifelse(ProxDistance==0,1,2**round(log(ProxDistance,2))),
         CfDistance=ifelse(CfProxDistance==0,1,2**round(log(CfProxDistance,2))))

View(cfData %>% filter(CfChainId>0))


print(ggplot(clusterDistances %>% 
               filter(ClusterReason!='LOH_CHAIN'&ClusterReason!='OVERLAP_FOLDBACKS'&ClusterReason!='CONSEC_BREAKS') %>% 
               group_by(ClusterReason,DistanceBucket=2**round(log(ClusterDistance,2))) %>% count, 
             aes(x=DistanceBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~ClusterReason)
      + labs(title = "Clustering Distances in Linx by Clustering Reason"))



combinedClusterDistances = rbind(cfData %>% filter(CfChainId>0) %>% mutate(Source='CF',ClusterReason='PROXIMITY') %>% 
                                   select(SampleId,Source,SvId,ClusterReason,ClusterDistance=CfProxDistance),
                                 cfData %>% filter(ClusterCount>1) %>% mutate(Source='LINX') %>% select(SampleId,Source,SvId,ClusterReason,ClusterDistance=ProxDistance))

View(combinedClusterDistances)



# much longer tail for Linx's clustering
print(ggplot(combinedClusterDistances %>% filter(ClusterDistance>0) %>% group_by(Source,DistanceBucket=2**round(log(ClusterDistance,2))) %>% count, 
             aes(x=DistanceBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~Source)
      + labs(title = "Clustering Distances Linx vs CF"))

# variety of reason for long-distance clustering in Linx
View(combinedClusterDistances %>% group_by(Source,ClusterReason) %>% count)
print(ggplot(combinedClusterDistances %>% 
               filter(ClusterDistance>0&Source=='LINX'&ClusterReason!='PROXIMITY'&ClusterReason!='LOH_CHAIN'&ClusterReason!='OVERLAP_FOLDBACKS'&ClusterReason!='CONSEC_BREAKS') %>% 
               group_by(ClusterReason,DistanceBucket=2**round(log(ClusterDistance,2))) %>% count, 
             aes(x=DistanceBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~ClusterReason)
      + labs(title = "Clustering Distances in Linx by Clustering Reason"))





View(cfData %>% filter((LnxDbLenStart>-500&LnxDbLenStart<0&CfLinkLenStart>0&abs(LnxDbLenStart+CfLinkLenStart)<=2)|
                         (LnxDbLenEnd>-500&LnxDbLenEnd<0&CfLinkLenEnd>0&abs(LnxDbLenEnd+CfLinkLenEnd)<=2)) %>%
       select(SampleId,SvId,LnxDbLenStart,CfLinkLenStart,LnxDbLenEnd,CfLinkLenEnd,everything()))


cfData = cfData %>% mutate(HasLinkStart=CfLinkLenStart>0,
                           HasLinkEnd=CfLinkLenEnd>0,
                           CfLinkDbOverhangMatchStart=HasLinkStart&LnxDbLenStart>-500&LnxDbLenStart<0&abs(LnxDbLenStart+CfLinkLenStart)<=2,
                           CfLinkDbOverhangMatchEnd=HasLinkEnd&LnxDbLenEnd>-500&LnxDbLenEnd<0&abs(LnxDbLenEnd+CfLinkLenEnd)<=2)

(sum(cfData$HasLinkStart)+sum(cfData$HasLinkEnd))/nrow(cfData %>% filter(CfChainId>0))
sum(cfData$HasLinkStart)

sum(cfData$HasLinkStart&cfData$CfLinkDbOverhangMatchStart)
sum(cfData$HasLinkEnd)

sum(cfData$CfLinkDbOverhangMatchStart)/sum(cfData$HasLinkStart)
sum(cfData$CfLinkDbOverhangMatchEnd)/sum(cfData$HasLinkEnd)

sum(cfData$HasLinkEnd&cfData$CfLinkDbOverhangMatchEnd)

View(cfData %>% filter((LnxDbLenStart>-500&LnxDbLenStart<0&CfLinkLenStart>0&abs(LnxDbLenStart+CfLinkLenStart)<=2)|
                         (LnxDbLenEnd>-500&LnxDbLenEnd<0&CfLinkLenEnd>0&abs(LnxDbLenEnd+CfLinkLenEnd)<=2)) %>%
       select(SampleId,SvId,LnxDbLenStart,CfLinkLenStart,LnxDbLenEnd,CfLinkLenEnd,everything()))


####
## Deletion Bridge stats

## Linx DBs
dbStats = rbind(cfCompare100 %>% filter(CfChainId>=0&DbMatchedStart!='NONE') %>% mutate(IsStart=T) %>% select(SampleId,SvId,DbMatched=DbMatchedStart,LnxDbLen=LnxDbLenStart,CfDbLen=CfDbLenStart),
                cfCompare100 %>% filter(CfChainId>=0&DbMatchedEnd!='NONE') %>% mutate(IsStart=F) %>% select(SampleId,SvId,DbMatched=DbMatchedEnd,LnxDbLen=LnxDbLenEnd,CfDbLen=CfDbLenEnd))

dbStats = rbind(cfCompare100 %>% filter(LnxDbLenStart>-1000) %>% mutate(Source='LINX') %>% select(SampleId,SvId,CfSvId,Source,DbMatched=DbMatchedStart,DbLen=LnxDbLenStart),
                cfCompare100 %>% filter(LnxDbLenEnd>-1000) %>% mutate(Source='LINX') %>% select(SampleId,SvId,CfSvId,Source,DbMatched=DbMatchedEnd,DbLen=LnxDbLenEnd),
                cfCompare100 %>% filter(CfDbLenStart>-1000) %>% mutate(Source='CF') %>% select(SampleId,SvId,CfSvId,Source,DbMatched=DbMatchedStart,DbLen=CfDbLenStart),
                cfCompare100 %>% filter(CfDbLenEnd>-1000) %>% mutate(Source='CF') %>% select(SampleId,SvId,CfSvId,Source,DbMatched=DbMatchedEnd,DbLen=CfDbLenEnd))

dbStats = dbStats %>% mutate(DbLengthBucket=ifelse(DbLen==0,0,ifelse(DbLen<0,-2**round(log(-DbLen,2)),2**round(log(DbLen,2)))))
View(dbStats)
View(dbStats %>% filter(is.na(DbLengthBucket)))

View(dbStats %>% group_by(Source,DbLengthBucket) %>% count %>% spread(Source,n,fill=0))
View(dbStats %>% group_by(Source) %>% summarise(Total=n(),ZeroDB=sum(DbLen==0),NegDBs=sum(DbLen<0)))

print(ggplot(dbStats %>% filter(DbLengthBucket>0&DbLengthBucket<1e4) %>% group_by(Source,DbLengthBucket) %>% count, aes(x=DbLengthBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~Source)
      + labs(title = "DB Lengths Linx vs CF"))

print(ggplot(dbStats %>% filter(DbLengthBucket<1e3) %>% group_by(Source,DbLengthBucket) %>% count, aes(x=DbLengthBucket, y=n))
      + geom_line()
      + facet_wrap(~Source)
      + labs(title = "DB Lengths Linx vs CF"))

View(dbStats %>% group_by(DbMatched) %>% summarise(Count=n(),
                                                   OverlapDBs=sum(LnxDbLen>-1000&LnxDbLen<0),
                                                   LnxAvg=round(mean(LnxDbLen)),LnxMin=min(LnxDbLen),LnxMax=max(LnxDbLen),
                                                   CfAvg=round(mean(CfDbLen)),CfMin=min(CfDbLen),CfMax=max(CfDbLen))) 




cfSamplesWithChains = cfCompare %>% filter(CfChainId>=0) %>% group_by(SampleId) %>% count
View(cfSamplesWithChains)
View(cfCompare %>% filter(SampleId %in% cfSamplesWithChains$SampleId) %>% filter(ClusterCount>100))
View(cfCompare %>% filter(SampleId %in% cfSamplesWithChains$SampleId) %>% filter(ClusterCount>100))





