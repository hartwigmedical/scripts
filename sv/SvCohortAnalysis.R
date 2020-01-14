# library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringi)
library(devtools)

svData = read.csv('~/data/sv/drivers/LNX_SVS.csv')
svData = read.csv('~/logs/LNX_SVS.csv')
clusters = read.csv('~/data/sv/drivers/LNX_CLUSTERS.csv')
clusters = read.csv('~/logs/LNX_CLUSTERS.csv')
nrow(svData)


#####
## Shattering characteristics
shClusters = clusters %>% filter(ClusterCount>=3&ClusterCount<=100&ResolvedType=='COMPLEX'&MinPloidy==1&MinPloidy==MaxPloidy&Foldbacks==0)
nrow(shClusters) # 14684

# fully chained and no single breakends
nrow(shClusters %>% filter(SglCount==0&InfCount==0&FullyChained=='true')) # 5357
ctClusters = shClusters %>% filter(SglCount==0&InfCount==0&FullyChained=='true')

View(ctClusters %>% group_by(ClusterCount=ifelse(ClusterCount<=10,ClusterCount,2**round(log(ClusterCount,2)))) %>% count)

print(ggplot(ctClusters %>% group_by(ClusterCount) %>% count,
             aes(x=ClusterCount, y=n))
      + geom_bar(stat='identity',colour='black')
      + labs(title = "Shattering - Cluster Count"))

# with gain, consecutive breakends or overlapping TIs - all very rare
nrow(ctClusters %>% filter(AcSameOrient>0)) # 446
nrow(ctClusters %>% filter(OverlapTIs>0)) # 377
nrow(ctClusters %>% filter(IntTIsCnGain>0)) # 244

# chain counts - 4K have only a single chain
View(ctClusters %>% group_by(ChainCount) %>% count)

print(ggplot(ctClusters %>% group_by(ChainCount) %>% count,
             aes(x=ChainCount, y=n))
      + geom_bar(stat='identity',colour='black')
      + labs(title = "Shattering - Chain Count"))

# deleted material
View(ctClusters %>% filter(ChainCount==1&TotalDeleted>0) %>% group_by(TotalDeleted=2**round(log(TotalDeleted,2))) %>% count)
View(ctClusters %>% filter(ChainCount==1&TotalDeleted==0))

print(ggplot(ctClusters %>% filter(ChainCount==1&TotalDeleted>0) %>% group_by(TotalDeleted=2**round(log(TotalDeleted,2))) %>% count,
             aes(x=TotalDeleted, y=n))
      + geom_line()
      + scale_x_log10()
      + labs(title = "Shattering - Deleted Material"))

print(ggplot(ctClusters %>% filter(ChainCount==1&TotalDeleted>0) %>% group_by(ClusterSize=2**round(log(pmin(ClusterCount,16),2)),
                                                                              TotalDeleted=2**round(log(TotalDeleted,2))) %>% count,
             aes(x=TotalDeleted, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~ClusterSize)
      + labs(title = "Shattering - Deleted Material by Cluster Size"))

# as a percent of total range
View(ctClusters %>% filter(ChainCount==1&TotalDeleted>0&TotalRange>0) %>% group_by(ClusterSize=2**round(log(pmin(ClusterCount,16),2)),
                                                                                   DeletedPerc=round(pmin(TotalDeleted/TotalRange,1),1)) %>% count)

print(ggplot(ctClusters %>% filter(ChainCount==1&TotalDeleted>0&TotalRange>0) %>% 
               group_by(ClusterSize=2**round(log(pmin(ClusterCount,16),2)),
                        DeletedPerc=round(pmin(TotalDeleted/TotalRange,1),1)) %>% count,
             aes(x=DeletedPerc, y=n))
      + geom_bar(stat='identity',colour='black')
      + facet_wrap(~ClusterSize)
      + labs(title = "Shattering - deleted material as percent of total"))

print(ggplot(ctClusters %>% filter(ChainCount==1&TotalDeleted>0&TotalRange>0&ClusterCount<=10) %>% 
               group_by(ClusterCount,DeletedPerc=round(pmin(TotalDeleted/TotalRange,1),1)) %>% count,
             aes(x=DeletedPerc, y=n))
      + geom_bar(stat='identity',colour='black')
      + facet_wrap(~ClusterCount)
      + labs(title = "Shattering - deleted material as percent of total"))

print(ggplot(ctClusters %>% filter(ChainCount==1&TotalDeleted>0&TotalRange>=1e3) %>% 
               group_by(MaterialLength=10**round(log(TotalRange,10)),
                        DeletedPerc=round(pmin(TotalDeleted/TotalRange,1),1)) %>% count,
             aes(x=DeletedPerc, y=n))
      + geom_bar(stat='identity',colour='black')
      + facet_wrap(~MaterialLength)
      + labs(title = "Shattering - Deleted Material as Percent of Cluster Range"))

print(ggplot(ctClusters %>% filter(ChainCount==1&TotalDeleted>0&TotalRange>=1e3) %>% 
               group_by(ClusterSize=2**round(log(pmin(ClusterCount,16),2)),
                        DeletedPerc=round(pmin(TotalDeleted/TotalRange,1),1)) %>% count,
             aes(x=DeletedPerc, y=n))
      + geom_bar(stat='identity',colour='black')
      + facet_wrap(~ClusterSize)
      + labs(title = "Shattering - Deleted Material Percent vs Cluster Size"))

nrow(clusters %>% filter(TotalDeleted!=TotalDBLength))
nrow(clusters %>% filter(AlleleValidPerc!=1))
nrow(clusters %>% filter(AlleleValidPerc==1))


# validation
View(ctClusters %>% filter(ChainCount==1&TotalDeleted>0&TotalRange>0&ClusterCount<=10&ChainEndsAway==1) %>% 
       mutate(DeletedPerc=round(pmin(TotalDeleted/TotalRange,1),3),
              DeletedPercAdj=round(pmin((TotalDeleted+TravDelLength)/TotalRange,1),3),
              DeletedPercDiff=DeletedPercAdj-DeletedPerc,
              RangeRatio=round(pmin(ChainedLength,TotalRange)/pmax(ChainedLength,TotalRange),3),
              DeleteRatio=round(pmin(TotalDeleted,TotalDBLength)/pmax(TotalDeleted,TotalDBLength),3),
              TravDelPerc=round(TravDelLength/TotalRange,3)) %>%
       filter(RangeRatio>0.5&DeleteRatio>0.5) %>% # DeletedPerc>0.75&
       select(SampleId,ClusterId,ClusterCount,ClusterDesc,DeletedPerc,DeletedPercAdj,DeletedPercDiff,TotalDeleted,TotalDBLength,TravDelCount,TravDelLength,TravDelPerc,TotalRange,ChainedLength,RangeRatio,OriginArms,FragmentArms,
              DBs,ShortDBs,TotalTIs,ExtTIs,IntTIs,everything()))

View(ctClusters %>% filter(SampleId=='CPCT02230001T'&ClusterId==74) %>% 
       mutate(DeletedPerc=round(pmin(TotalDeleted/TotalRange,1),3),
              RangeRatio=round(pmin(ChainedLength,TotalRange)/pmax(ChainedLength,TotalRange),3),
              DeleteRatio=round(pmin(TotalDeleted,TotalDBLength)/pmax(TotalDeleted,TotalDBLength),3)) %>%
       # filter(RangeRatio>0.5&DeleteRatio>0.5) %>% # DeletedPerc>0.75&
       select(SampleId,ClusterId,ClusterCount,ClusterDesc,DeletedPerc,TotalDeleted,DeleteRatio,TotalDBLength,TotalRange,RangeRatio,ChainedLength,OriginArms,FragmentArms,DBs,ShortDBs,TotalTIs,ExtTIs,IntTIs,everything()))

# traversed deletes and deleted lengths
View(ctClusters %>% group_by(ClusterCount=2**round(log(ClusterCount,2)),TravDelCount=2**round(log(TravDelCount,2))) %>% count %>% spread(TravDelCount,n,fill=0))

View(ctClusters %>% filter(TotalRange>0) %>% group_by(ClusterCount=2**round(log(ClusterCount,2)),TravDelPerc=round(TravDelLength/TotalRange,1)) %>% count %>% spread(TravDelPerc,n,fill=0))



View(ctClusters %>% filter(ChainCount==1&TotalDeleted>0&TotalRange>0&ClusterCount<=10) %>% 
  group_by(ClusterCount,DeletedPerc=round(pmin(TotalDeleted/TotalRange,1),1)) %>% count)



# clusters with almost all material retained
View(ctClusters %>% filter(ChainCount==1&ClusterCount>=10&TotalRange>0&ChainEndsAway==1) %>% mutate(DeletedPerc=round(TotalDeleted/TotalRange,2)) %>%
       select(SampleId,ClusterId,ClusterCount,DeletedPerc,TotalDeleted,TotalRange,OriginArms,FragmentArms,DBs,ShortDBs,TotalTIs,ExtTIs,IntTIs,everything()))

# restricted to those without BNDs or external TIs
View(ctClusters %>% filter(ChainCount==1&ClusterCount>=10&TotalRange>0&ChainEndsAway==1&BndCount==0&ExtTIs==0) %>% mutate(DeletedPerc=round(TotalDeleted/TotalRange,2)) %>%
       select(SampleId,ClusterId,ClusterCount,DeletedPerc,TotalDeleted,TotalRange,OriginArms,FragmentArms,DBs,ShortDBs,TotalTIs,ExtTIs,IntTIs,everything()))

View(ctClusters %>% filter(ClusterId==327) %>% mutate(DeletedPerc=round(TotalDeleted/TotalRange,2)) %>%
       select(SampleId,ClusterId,ClusterCount,DeletedPerc,TotalDeleted,TotalRange,OriginArms,FragmentArms,DBs,ShortDBs,TotalTIs,ExtTIs,IntTIs,everything()))

View(clusters %>% filter(SampleId=='CPCT02010687T'&ClusterId==327))

write.csv(ctClusters %>% filter(ChainCount==1&ClusterCount>=10&TotalRange>0&ChainEndsAway==1) %>% mutate(DeletedPerc=round(TotalDeleted/TotalRange,2)) %>%
       select(SampleId),'~/logs/shattering_samples.csv',row.names = F,quote = F)


## Comparison with shattering simulation
# simulations are suggesting that average/expected retained segments = 2/3
# cannot be sure of starting segments, but not an issue if perfect repair is rare, although consecutive segments can still be lost
# so fairer comparison is say a cluster with between 5-10 SVs

ctCleanClusters = ctClusters %>% filter(ChainCount==1&TotalDeleted>0&TotalRange>0&ClusterCount<=50&OriginArms==1)
nrow(ctCleanClusters) # 2814

ctSuperCleanClusters = ctClusters %>% filter(ChainCount==1&TotalDeleted>0&TotalRange>0&TotalRange>ChainedLength*1.1&ClusterCount<=50&OriginArms==1&ExtTIs-ExtShortTIs==0)
nrow(ctSuperCleanClusters)

View(ctCleanClusters %>% group_by(ClusterCount) %>% count)

ctCleanClusters = ctCleanClusters %>% mutate(DeletedPerc=round(pmin(TotalDeleted/TotalRange,1),3),
                                             AdjDeletedPerc=round(pmin((TotalDeleted+TravDelLength)/TotalRange,1),3),
                                             RangeRatio=round(pmin(ChainedLength,TotalRange)/pmax(ChainedLength,TotalRange),3),
                                             RetainedPerc=1-AdjDeletedPerc)

ctCleanClusters = ctCleanClusters %>% filter(RangeRatio>0.5)
nrow(ctCleanClusters) # 1405

ctCleanClusters = ctCleanClusters %>% mutate(ClusterSize=ifelse(ClusterCount<=6,ClusterCount,2**round(log(pmin(ClusterCount,16),2))))

print(ggplot(ctCleanClusters,aes(x=RetainedPerc)) 
      + stat_ecdf(geom='point')
      + facet_wrap(~ClusterSize))

# as a bar chart and in terms of retained percent
print(ggplot(ctCleanClusters %>% group_by(ClusterSize,RetainedPerc=round(RetainedPerc,1)) %>% count,aes(x=RetainedPerc, y=n))
      + geom_bar(stat='identity',colour='black')
      + facet_wrap(~ClusterSize)
      + labs(x='Retained %',y='Cluster Count',title='Shattering - Retained Material as % of Total by Cluster Size'))

View(ctCleanClusters %>% group_by(ClusterSize) %>% summarise(Count=n(),
                                                             AvgRetainedPerc=mean(RetainedPerc),
                                                             MedRetainedPerc=median(RetainedPerc)))

# by number of segments
ctCleanClusters = ctCleanClusters %>% mutate(DeletedSegs=DBs-ShortDBs+TravDelCount,
                                             TotalSegs=DeletedSegs+TotalTIs-ExtShortTIs+ImpliedTIs,
                                             DeletedSegPerc=round(DeletedSegs/TotalSegs,3),
                                             RetainedSegPerc=1-DeletedSegPerc,
                                             ImpliedClusterSize=ifelse(TotalSegs<=6,TotalSegs,2**round(log(pmin(TotalSegs,32),2))),
                                             RetainedRatio=round(pmin(RetainedPerc,RetainedSegPerc)/pmax(RetainedPerc,RetainedSegPerc),2),
                                             RetainedSegsVsMaterialDiff=abs(RetainedSegPerc-RetainedPerc))

View(ctCleanClusters %>% group_by(TotalSegs,DeletedSegs) %>% count() %>% spread(DeletedSegs,n,fill=0))

View(ctCleanClusters %>% select(SampleId,ClusterId,ClusterCount,DeletedSegs,DBs,ShortDBs,TotalSegs,TotalTIs,ShortTIs,ExtShortTIs,ImpliedTIs,TotalSegs,RetainedSegPerc,
                                TotalDeleted,TravDelCount,TravDelLength,TotalRange,RetainedPerc,everything()))

# comparing the 2 ratios - by material vs segments
View(ctCleanClusters %>% group_by(RetainedRatioDiff=round(RetainedPerc-RetainedSegPerc,1)) %>% count)

print(ggplot(ctCleanClusters %>% group_by(RetainedRatioDiff=round(RetainedPerc-RetainedSegPerc,1)) %>% count,aes(x=RetainedRatioDiff, y=n))
      + geom_bar(stat='identity',colour='black')
      + labs(x='Retained % Difference - +ve = more Material vs Segments',y='Cluster Count',title='Shattering - Retained Material vs Segments'))

nrow(ctCleanClusters)
write.csv(ctCleanClusters,'~/logs/shattering_clusters.csv',row.names = F,quote = F)

# most material lost
View(ctCleanClusters %>% 
       filter(TotalRange>TotalDeleted&OriginArms==1&ChainCount==1&(ExtTIs-ExtShortTIs)==0&ChainEndsAway==1&ChainedLength<1.1*TotalRange) %>% 
       filter(RetainedPerc<0.25) %>%
       select(SampleId,ClusterId,ClusterCount,ClusterDesc,RetainedPerc,TotalDeleted,TravDelCount,TravDelLength,TotalRange,
              DBs,ShortDBs,TotalTIs,ExtShortTIs,ExtTIs,everything()))


# discrepancies between total deleted and DeletedSegs
View(ctCleanClusters %>% filter(TotalRange>TotalDeleted&OriginArms==1&ChainCount==1&(ExtTIs-ExtShortTIs)==0&ChainEndsAway==1&ChainedLength<1.1*TotalRange) %>% filter(RetainedSegsVsMaterialDiff>0.5) %>%
       select(SampleId,ClusterId,ClusterCount,ClusterDesc,RetainedSegsVsMaterialDiff,DBs,ShortDBs,DeletedSegs,TotalSegs,TotalTIs,ExtShortTIs,ImpliedTIs,
              RetainedSegPerc,RetainedPerc,TotalDeleted,TravDelCount,TravDelLength,TotalRange,everything()))

# total range less than total deleted
View(ctCleanClusters %>% filter(TotalRange<0.8*TotalDeleted) %>%
       select(SampleId,ClusterId,ClusterCount,DBs,ShortDBs,TotalSegs,TotalTIs,ShortTIs,TotalDeleted,TotalDBLength,TravDelCount,TravDelLength,TotalRange,everything()))

# super clean clusters only

print(ggplot(ctCleanClusters %>% filter(ChainCount==1&TotalDeleted>0&TotalRange>0&TotalRange>ChainedLength*1.1&ClusterCount<=50&OriginArms==1&ExtTIs-ExtShortTIs==0) %>% 
               group_by(RetainedSegPerc=round(RetainedSegPerc,1)) %>% count,aes(x=RetainedSegPerc, y=n))
      + geom_bar(stat='identity',colour='black')
      + labs(x='% Retained Segments',y='Cluster Count', title='Shattering - Retained Inferred Segment %'))


View(ctCleanClusters %>% filter(RangeRatio<0.5) %>% 
       select(SampleId,ClusterId,ClusterCount,ClusterDesc,DeletedPerc,TotalDeleted,DeleteRatio,TotalDBLength,TotalRange,RangeRatio,ChainedLength,
              OriginArms,FragmentArms,DBs,ShortDBs,TotalTIs,ExtTIs,IntTIs,everything()))





## DEBUG
tmpClusters = read.csv('~/logs/LNX_CLUSTERS.csv')
View(tmpClusters %>% filter(IndelCount>0&MaxPloidy<2.5) %>% select(SampleId,ClusterId,ClusterCount,ClusterDesc,ResolvedType,TotalRange,TotalDeleted,
                                                     TotalTIs,DBs,ShortDBs,TotalDBLength,ChainedLength,IndelCount,FullyChained,ChainCount,Consistency,ChainEndsAway,
                                                     Foldbacks,MaxPloidy,everything()))

View(tmpClusters %>% filter(TotalRange!=ChainedLength) %>% select(SampleId,ClusterId,ClusterCount,ClusterDesc,ResolvedType,TotalRange,TotalDeleted,
                                                                  TotalTIs,DBs,ShortDBs,TotalDBLength,ChainedLength,IndelCount,FullyChained,ChainCount,Consistency,ChainEndsAway,
                                                                  Foldbacks,MaxPloidy,everything()))
                                                                   
View(tmpClusters %>% filter(TotalRange<TotalDeleted) %>% select(SampleId,ClusterId,ClusterCount,ClusterDesc,ResolvedType,TotalRange,TotalDeleted,
                                                                  TotalTIs,DBs,ShortDBs,TotalDBLength,ChainedLength,IndelCount,FullyChained,ChainCount,Consistency,ChainEndsAway,
                                                                  Foldbacks,MaxPloidy,everything()))

View(tmpClusters %>% filter((TotalDBLength>0|TotalDeleted>0)&abs(TotalDBLength-TotalDeleted)/pmax(TotalDeleted,TotalDBLength)>0.1) %>% 
       select(SampleId,ClusterId,ClusterCount,ClusterDesc,ResolvedType,TotalRange,TotalDeleted,TotalTIs,DBs,ShortDBs,TotalDBLength,ChainedLength,IndelCount,FullyChained,ChainCount,Consistency,
              ChainEndsAway,Foldbacks,MaxPloidy,everything()))

# deletion bridges by cluster type and size
ctClusters = ctClusters %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))
svData = svData %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))

View(ctClusters %>% group_by(SampleId) %>% count)

ctSvData = svData %>% filter(SampleClusterId %in% ctClusters$SampleClusterId)
ctSvData = ctSvData %>% mutate(IsFS=FSStart=='true'|FSEnd=='true')
nrow(ctSvData)
View(ctSvData)

svDbData = rbind(ctSvData %>% filter(DBLenStart>-1000) %>% select(SampleId,ClusterId,ClusterCount,DBLen=DBLenStart,IsFS),
                 ctSvData %>% filter(DBLenEnd>-1000) %>% select(SampleId,ClusterId,ClusterCount,DBLen=DBLenEnd,IsFS),
                 ctSvData %>% filter(Type=='DEL'&DBLenStart==-1000&DBLenEnd==-1000) %>% mutate(DBLen=PosEnd-PosStart) %>% 
                   select(SampleId,ClusterId,ClusterCount,DBLen,IsFS))

View(svDbData %>% group_by(SampleId,ClusterId,DBLen,IsFS) %>% count %>% group_by(n) %>% count)
View(svDbData %>% group_by(SampleId,ClusterId,DBLen) %>% count %>% group_by(n) %>% count)

# de-dup the pairs of DB lengths within a cluster - will only lose about 400 entries
svDbData = svDbData %>% group_by(SampleId,ClusterId,ClusterCount,DBLen) %>% summarise(IsFS=first(IsFS)|last(IsFS)) %>% ungroup()
View(svDbData)


View(svDbData %>% group_by(ClusterCount,DBLen) %>% count)
View(svDbData %>% group_by(ClusterCount) %>% count)
View(svDbData %>% group_by(ClusterCount,IsFS) %>% count %>% spread(IsFS,n,fill=0))

svDbData = svDbData %>% mutate(ClusterSize=as.character(ifelse(ClusterCount<=10,ClusterCount,2**round(log(ClusterCount,2)))))

print(ggplot(svDbData, aes(x=reorder(as.character(ClusterSize),ClusterCount),y=DBLen))
      + geom_violin(scale="count",fill="#6baed6")
      + scale_y_log10()
      # + facet_wrap(~IsFS)
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + labs(title = "Deletion Bridge length by Cluster Size",x='Cluster Size',y='DB Length'))


# deleted material taken from DELs and DBs - very similar picture to that from the cluster data
dbDeleted = svDbData %>% group_by(SampleId,ClusterId,ClusterCount) %>% summarise(TotalDeleted=sum(DBLen))

print(ggplot(dbDeleted %>% filter(TotalDeleted>0) %>% group_by(ClusterSize=2**round(log(pmin(ClusterCount,16),2)),
                                                               TotalDeleted=2**round(log(TotalDeleted,2))) %>% count,
             aes(x=TotalDeleted, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~ClusterSize)
      + labs(title = "Shattering - Deleted Material (from DBs) by Cluster Size"))


# clusters with very even DB lengths
clusterDbData = svDbData %>% group_by(SampleId,ClusterId,ClusterCount) %>% 
  summarise(DBCount=n(),TotalDeleted=sum(DBLen),ShortDBs=sum(DBLen<100),
            MinDB=min(DBLen),MaxDB=max(DBLen),MedDB=median(DBLen),StdDev=sd(DBLen)) %>% 
  mutate(MinMaxRatio=round(MinDB/MaxDB,2))

View(clusterDbData %>% filter(DBCount>=5&MinDB>1000))

View(ctSvData %>% filter(SampleId=='CPCT02220030T'&ClusterId==76))




#####
## Replication before Repair candidates
# - uniform ploidy & consistent
upClusters = clusters %>% filter(ClusterCount>=2&ClusterCount<=20&MinPloidy==1&MinPloidy==MaxPloidy&Consistency==0)

# - ignore single breakends for cleanness
upClusters = upClusters %>% filter(SglCount==0&InfCount==0)

# - keep to local types and skip shards for now
upClusters = upClusters %>% filter(BndCount==0)

# - restrict to local, simple cluster types
upClusters = upClusters %>% filter(ResolvedType %in% c('COMPLEX','DEL_TI','DEL','DUP','RECIP_INV','RECIP_INV_DEL_DUP','RESOLVED_FOLDBACK',
                                                       'PAIR_OTHER','FB_INV_PAIR','SIMPLE_GRP','RECIP_INV_DUPS'))

View(upClusters %>% group_by(ResolvedType) %>% count)
nrow(upClusters) # 13.8K

View(clusters %>% filter(grepl('REP_REP',Annotations)))

# - evidence of internal TIs
nrow(upClusters %>% filter(IntTIs>0))
nrow(upClusters %>% filter(IntTIsCnGain>0))
nrow(upClusters %>% filter(AcSameOrient>0))

# cluster 2 and 3s

# simple groups
View(upClusters %>% filter(ResolvedType=='SIMPLE_GRP'&ChainCount==1&ChainEndsAway==1))


# how many simple DELs and DUPs are in short DBs with other DELs and DUPs?
simpleDDs = svData %>% filter(Type %in% c('DEL','DUP')&ClusterCount==1)
# simpleDDs = simpleDDs %>% mutate(NearestType=ifelse(NearestLen<0,'',as.character(NearestType)))
nrow(simpleDDs)
View(simpleDDs %>% filter(NearestType!='') %>% group_by(NearestType,NearestLen=ifelse(NearestLen>0,2**round(log(NearestLen,2)),0)) %>% count %>% spread(NearestType,n,fill=0))
View(simpleDDs %>% filter(NearestType!=''&NearestLen>0&NearestLen<=50) %>% mutate(Length=PosEnd-PosStart) %>% select(SampleId,Id,Type,Length,NearestType,NearestLen,PosStart,PosEnd,
                                                                                                       DBLenStart,DBLenEnd,Ploidy,PloidyMin,PloidyMax,CNStart,CNEnd,everything()))


View(upClusters %>% filter(ClusterCount==3) %>% group_by(ClusterDesc) %>% count)
View(upClusters %>% filter(ClusterDesc=='DUP=1_INV=2') %>% 
       select(SampleId,ClusterId,ClusterCount,AcSameOrient,AcDsb,AcTIOnly,TotalTIs,OverlapTIs,IntTIs,IntTIsCnGain,Foldbacks,everything()))


upClusters = upClusters %>% filter(ClusterCount>=2&ClusterCount<=20&SglCount==0&InfCount==0&MinPloidy==1&MinPloidy==MaxPloidy&Consistency==0) # FullyChained=='true'


View(upClusters %>% filter(ResolvedType=='DEL'))
nrow(upClusters %>% filter(BndCount==0)) # 14174 with no shards



#####
## Non-BFB clusters with amplification
View(clusters %>% filter(ResolvedType=='COMPLEX') %>% group_by(Annotations) %>% count)
nonfbClusters = clusters %>% filter(ResolvedType=='COMPLEX'&Annotations==''&Replication=='true'&MaxPloidy>=4&Foldbacks==0)
nrow(nonfbClusters) # 2147
View(nonfbClusters)

# 35K complex clusters
nrow(clusters %>% filter(ResolvedType=='COMPLEX'))

# 1612 have a single breakend
nrow(nonfbClusters %>% filter(SglCount>0|InfCount>0))

# BFB clusters - 9700
nrow(clusters %>% filter(ResolvedType=='COMPLEX'&!grepl('DM',Annotations)&Foldbacks>0))

# 11K have amplification
nrow(clusters %>% filter(ResolvedType=='COMPLEX'&MinPloidy<MaxPloidy))
nrow(clusters %>% filter(ResolvedType=='COMPLEX'&MinPloidy<MaxPloidy&Foldbacks>0)) # 6K

# 535 have no foldbacks and no single breakends
nrow(nonfbClusters %>% filter(SglCount==0&InfCount==0))

View(nonfbClusters %>% filter(SglCount==0&InfCount==0))
View(nonfbClusters %>% filter(SglCount==0&InfCount==0&ChainCount==1&FullyChained=='true'))

View(nonfbClusters %>% filter(SglCount==0&InfCount==0&FullyChained=='true'&ChainCount<=4) %>% 
       mutate(OverlapTIPerc=round(OverlapTIs/TotalTIs,2)))


# shattering without any replication by way of comparison - 7.5K
nrow(clusters %>% filter(ResolvedType=='COMPLEX'&!grepl('DM',Annotations)&!grepl('BFB',Annotations)) %>%
       filter(SglCount==0&InfCount==0&MinPloidy==1&MinPloidy==MaxPloidy&Foldbacks==0))

# 5175 are fully chained
nrow(clusters %>% filter(ResolvedType=='COMPLEX'&!grepl('DM',Annotations)&!grepl('BFB',Annotations)) %>%
       filter(SglCount==0&InfCount==0&MinPloidy==1&MinPloidy==MaxPloidy&Foldbacks==0&FullyChained=='true'))

View(clusters %>% filter(ResolvedType=='COMPLEX'&!grepl('DM',Annotations)&!grepl('BFB',Annotations)) %>%
       filter(SglCount==0&InfCount==0&MinPloidy==1&MinPloidy==MaxPloidy&Replication=='false'
              &AcSameOrient>0&ClusterCount<50&ClusterCount>4) %>%
       select(SampleId,ClusterId,ClusterCount,AcSameOrient,AcDsb,AcTIOnly,TotalTIs,OverlapTIs,everything()))


# strict
View(clusters %>% filter(ResolvedType=='COMPLEX'&MinPloidy==1&MinPloidy==MaxPloidy&Foldbacks==0&ClusterCount<50&OverlapTIs>0&SglCount==0&InfCount==0
                         &ChainCount==1&FullyChained=='true'&ChainEndsAway==1&ClusterCount>3))

# looser
View(clusters %>% filter(ResolvedType=='COMPLEX'&MinPloidy==1&MinPloidy==MaxPloidy&ClusterCount<50&OverlapTIs>=2&SglCount<3&InfCount<3
                         &ClusterCount>3&AcSameOrient>2))


View(svData %>% filter(SampleId=='DRUP01080008T'&ClusterId==79))
View(svData %>% filter(SampleId=='CPCT02010419TII'&ClusterId==496))



#####
## DOUBLE MINUTES - see SvDoubleMinutes.R instead


#####
## TYFONAS
# Sum(FB CN) / Max Interval CN >= 0.61 AND >=27 variants with JCN >= ploidy

samplePloidy = read.csv('~/data/sample_ploidy_purity.csv')
nrow(samplePloidy)

svData = svData %>% mutate()

svData = svData %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'),Ploidy=(PloidyMin+PloidyMax)*0.5)
clusters = clusters %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))

tyfClusters = clusters %>% filter(Foldbacks>0&ClusterCount>=27)
nrow(tyfClusters) # 2094
nrow(tyfClusters %>% filter(Annotations=='BFB')) # 1651

tyfSVs = svData %>% filter(SampleClusterId %in% tyfClusters$SampleClusterId)
tyfSVs = tyfSVs %>% mutate(IsFoldback=FoldbackLenStart>0|FoldbackLenEnd>0)

tyfClusterMaxCN = tyfSVs %>% group_by(SampleId,ClusterId) %>% summarise(MaxCNStart=max(CNStart),MaxCNEnd=max(CNEnd))
tyfClusterMaxCN = tyfClusterMaxCN %>% mutate(MaxCN=pmax(MaxCNStart,MaxCNEnd))
View(tyfClusterMaxCN)

# tyfSVs = merge(tyfSVs,tyfClusterMaxCN %>% select(SampleId,ClusterId,MaxCN),by=c('SampleId','ClusterId'),all.x=T)
tyfSVs = merge(tyfSVs,samplePloidy %>% select(SampleId,SamplePloidy=Ploidy),by='SampleId',all.x=T)
# nrow(tyfSVs %>% filter(is.na(SamplePloidy)|is.na(MaxCN)))

tyfSvSummary = tyfSVs %>% group_by(SampleId,ClusterId,ClusterCount) %>% 
  summarise(FbPloidyTotal=sum(ifelse(IsFoldback,Ploidy,0)),
            MaxSvPloidy=max(Ploidy),
            HighPloidyCount=sum(Ploidy>=SamplePloidy),
            MaxCNStart=max(CNStart),MaxCNEnd=max(CNEnd)) %>% mutate(MaxCN=pmax(MaxCNStart,MaxCNEnd)) %>% select(-MaxCNStart,-MaxCNEnd)

tyfSvSummary = tyfSvSummary %>% mutate(FbPloidyRatio=pmin(round(FbPloidyTotal/MaxCN,3),1))
tyfSvSummary = tyfSvSummary %>% mutate(IsTyfonas=HighPloidyCount>=27&FbPloidyTotal/MaxCN>=0.61&MaxSvPloidy>=7)

tyfSvSummary = merge(tyfSvSummary,sampleCancerTypes,by='SampleId',all.x=T)
sampleCancerSubTypes = load_cancer_types('~/data/hpc_sample_cancer_types.csv',T,10)
# View(sampleCancerSubTypes)
tyfSvSummary = merge(tyfSvSummary,sampleCancerSubTypes %>% select(-CancerType),by='SampleId',all.x=T)

View(tyfSvSummary)
View(tyfSvSummary %>% filter(HighPloidyCount>=27&FbPloidyTotal/MaxCN>=0.61&MaxSvPloidy>=7))
View(tyfSvSummary %>% group_by(IsTyfonas,
                               FbPloidyRatio=round(FbPloidyRatio,1),
                               MaxCN=2**round(log(MaxCN,2))) %>% count %>% spread(MaxCN,n,fill=0))

View(tyfSvSummary %>% group_by(NHIGH=2**round(log(HighPloidyCount,2)),
                               FBIJCN=round(FbPloidyRatio,1)) %>% count %>% spread(NHIGH,n,fill=0))
tyfSvSummaryPlot = tyfSvSummary %>% filter(HighPloidyCount>0) %>% 
  group_by(HighPloidyCount=pmin(5*round(HighPloidyCount/5),100),
           FbPloidyRatio=round(FbPloidyRatio,1),
           MaxSvPloidy=ifelse(MaxSvPloidy<6.7,'Max SV Ploidy: 1-6','Max SV Ploidy: >7')) %>% 
  summarise(Clusters=n())

print(ggplot(tyfSvSummaryPlot,aes(x=HighPloidyCount,y=reorder(FbPloidyRatio,-FbPloidyRatio))) 
      + geom_tile(aes(fill=Clusters),colour="white",stat="identity",position="identity") 
      + geom_text(aes(label=Clusters))
      + scale_fill_gradient(low="white",high="steelblue")
      # + scale_x_log10()
      + facet_wrap(~MaxSvPloidy)
      + labs(x='# SVs with Ploidy > Sample Ploidy',y='Sum(Foldback Ploidy) / Max Cluster Copy Number')
      + ggtitle("TYFONAS criteria"))

View(tyfSvSummary %>% group_by(CancerType,IsTyfonas) %>% count %>% spread(IsTyfonas,n,fill=0))

View(tyfSvSummary %>% group_by(CancerType,Type=ifelse(IsTyfonas,'Tyfonas',ifelse(HighPloidyCount<27&FbPloidyTotal/MaxCN>=0.61&MaxSvPloidy>=7,'BFB','Other'))) %>% 
       count %>% spread(Type,n,fill=0))

View(tyfSvSummary %>% group_by(CancerType,NHIGH=ifelse(FbPloidyTotal/MaxCN>=0.61&MaxSvPloidy>=7,2**round(log(HighPloidyCount,2)),0)) %>% 
       count %>% spread(NHIGH,n,fill=0))

write.csv(tyfSvSummary %>% filter(HighPloidyCount>=27&FbPloidyTotal/MaxCN>=0.61&MaxSvPloidy>=7) %>% 
            filter(ClusterCount<=100) %>% select(SampleId,ClusterId),'~/logs/tyfonas_clusters.csv',row.names = F, quote = F)


#####
## PYRGO
# Enriched counts of DUP with length > 10kb and < 1Mbp within and JCN < ploidy within a 1Mbp sliding window

rm(pyrgoDups)
rm(pyrgoSampleData)

pyrgoDups = svData %>% filter(Type=='DUP'&ClusterCount==1) %>% mutate(Length=PosEnd-PosStart)
pyrgoDups = merge(pyrgoDups,samplePloidy %>% select(SampleId,SamplePloidy=Ploidy),by='SampleId',all.x=T)
pyrgoDups = pyrgoDups %>% filter(Ploidy<SamplePloidy&Length>1e4&Length<1e6)
nrow(pyrgoDups) # 70.3K

# bucket into 2MB windows and looks for counts of > 1
pyrgoDups = pyrgoDups %>% mutate(ChrPosBucket=paste(ChrStart,2e3*round(PosStart/2e6),'Kb',sep="_"))
pyrgoSampleData = pyrgoDups %>% group_by(SampleId,ChrPosBucket) %>% summarise(DupCount=n())
pyrgoSampleData = pyrgoSampleData %>% group_by(SampleId) %>% summarise(PyrgoDups=sum(ifelse(DupCount>1,DupCount,0)),Locations=sum(DupCount>1))
pyrgoSampleData = merge(pyrgoSampleData,pyrgoDups %>% group_by(SampleId) %>% summarise(DupCount=n()),by='SampleId',all.x=T)
pyrgoSampleData = merge(pyrgoSampleData,sampleCancerSubTypes,by='SampleId',all.x=T)
pyrgoSampleData = pyrgoSampleData %>% mutate(PyrgoRate=round(PyrgoDups/DupCount,2))
View(pyrgoSampleData)
View(pyrgoSampleData %>% group_by(PyrgoRate=round(PyrgoRate,1)) %>% count)
View(pyrgoSampleData %>% filter(DupCount>=5) %>% group_by(CancerType,PyrgosRate=0.05*round(PyrgosRate/0.05)) %>% count %>% spread(PyrgosRate,n,fill=0))

write.csv(pyrgoSampleData,'~/logs/pyrgo_dups.csv',row.names = F,quote = F)

print(ggplot(pyrgoSampleData,aes(x=DupCount,y=PyrgoDups))
      + geom_point()
      + labs(title = "PYRGO: DUPs (10kb-1MB) vs clustered DUPs within 2MB"))

print(ggplot(pyrgoSampleData %>% mutate(DupEnriched=(DupCount>=100)),aes(x=DupCount,y=PyrgoDups))
      + geom_point()
      + facet_wrap(~DupEnriched)
      + labs(title = "PYRGO: DUPs (10kb-1MB) vs clustered DUPs within 2MB - faceted by Enriched in DUPs (>100)"))

# enrichment by location - pull in oncogenes as well
oncoGenes = ampDrivers %>% group_by(Gene) %>% count %>% ungroup() %>% select(Gene)
View(oncoGenes)
View(ampDrivers)

geneData = read.csv('~/data/sv/ensembl_gene_data.csv')
View(geneData %>% filter(Chromosome==8&GeneStart>127e6&GeneEnd<129e6))
oncoGenes = merge(oncoGenes,geneData %>% select(Gene=GeneName,Chromosome,GeneStart,GeneEnd),by='Gene',all.x=T)
oncoGenes = oncoGenes %>% mutate(ChrPosBucket=paste(Chromosome,2e3*round((GeneStart+GeneEnd)*0.5/2e6),'Kb',sep="_"))

pyrgoDups = pyrgoDups %>% mutate(AmpGene=ChrPosBucket %in% oncoGenes$ChrPosBucket)
View(pyrgoDups)

pyrgoDriverGene = pyrgoDups %>% group_by(SampleId,ChrPosBucket) %>% summarise(DupCount=n(),AmpGene=first(AmpGene)) %>% filter(DupCount>=3&AmpGene)
pyrgoDriverGene = merge(pyrgoDriverGene,oncoGenes,by='ChrPosBucket',all.x=T)
View(pyrgoDriverGene)
View(pyrgoDriverGene %>% group_by(Gene,DupCount) %>% count %>% spread(DupCount,n,fill=0) %>% mutate(TotalSamples=`3`+`4`+`5`+`6`))

pyrgoLocations = pyrgoDups %>% group_by(ChrPosBucket,AmpGene) %>% summarise(Count=n()) %>% filter(Count>1)
View(pyrgoLocations)

print(ggplot(pyrgoLocations %>% filter(Count>=100),aes(x=reorder(ChrPosBucket,-Count),y=Count))
      + geom_bar(stat='identity',colour='black')
      + facet_wrap(~AmpGene)
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + labs(x='Location (Chr Pos (KB)',title = "PYRGO: Hotspot Locations - by overlaps an AMP Driver Oncogene"))


# validation
View(pyrgoDups %>% filter(SampleId=='CPCT02010240T'))

rm(pyrgoDriverGene)
rm(pyrgoDups)
rm(pyrgoSampleData)
rm(pyrgoLocations)

rm(geneData)





