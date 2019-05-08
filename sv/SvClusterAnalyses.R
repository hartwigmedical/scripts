

clusters = read.csv('~/logs/SVA_CLUSTERS.csv')

View(clusters %>% filter(ClusterCount>1&ClusterCount<100))
View(clusters %>% filter(ClusterCount>1&ClusterCount<100&DSBs>ClusterCount))
View(clusters %>% filter(TotalLinks<IntTIs+ExtTIs))
View(clusters %>% filter(SampleId=='COLO829T'))

ctClusters = clusters %>% filter(ClusterCount>1&ClusterCount<100) # approx 111K
ctClusters$SampleClusterId = paste(ctClusters$SampleId,ctClusters$ClusterId,sep='_')
nrow(ctClusters)
View(ctClusters)
rm(clusters)



# simple (no foldbacks), little to no replication, 1 or 2 chains = 5K
simpleClusters = ctClusters %>% filter(FullyChained=='true'&Foldbacks==0&ResolvedType=='SimpleChain'&MaxCopyNumber<=2.5&ChainCount<=2)
View(simpleClusters)

View(simpleClusters %>% select(SampleId,ClusterId,ClusterCount,ChainCount,IntTIs,ExtTIs,IntTIsWithGain,ExtTIsWithGain,OverlapTIs,DSBs,ShortDSBs,
                               ChainEndsFace,ChainEndsAway,ArmCount,OriginArms,FragmentArms,TotalLinks,AssemblyLinks,ShortTIRemotes,MinCopyNumber,MaxCopyNumber))

# 1. Single chain, no overlaps, chain ends facing away - count = 2300
# shows classic shattering, but also examples of failed replication, with no short DBs
simpleCT = ctClusters %>% filter(FullyChained=='true'&Foldbacks==0&ResolvedType=='SimpleChain'&MaxCopyNumber<=2.5&ChainCount==1&IntTIsWithGain==0&OverlapTIs==0&ChainEndsAway==1)
simpleCT = simpleCT %>% mutate(ShortDBPerc=round(ShortDSBs/DSBs,1),)

# further restrict to those having no external TIs
internalCT = simpleCT %>% filter(DSBs>0&ClusterCount>=4&ExtTIs==0)

nrow(simpleCT)
View(simpleCT %>% select(SampleId,ClusterId,ClusterCount,ChainCount,IntTIs,ExtTIs,IntTIsWithGain,ExtTIsWithGain,OverlapTIs,DSBs,ShortDSBs,
                               ChainEndsFace,ChainEndsAway,ArmCount,OriginArms,FragmentArms,TotalLinks,AssemblyLinks,ShortTIRemotes,MinCopyNumber,MaxCopyNumber))
View(simpleCT %>% group_by(ClusterCount) %>% count())
View(simpleCT %>% filter(DSBs>0&ClusterCount>4) %>% group_by(ShortDBPerc=round(ShortDSBs/DSBs,1)) %>% count())
View(simpleCT %>% filter(DSBs>0&ClusterCount>4&ExtTIs==0) %>% group_by(ShortDBPerc=round(ShortDSBs/DSBs,1)) %>% count())

# 2. Correlation with TI length
# showing a clear peak around 500K for TIs where no short DBs exist
tiLinks = read.csv('~/logs/SVA_LINKS.csv')
tiLinks$SampleClusterId = paste(tiLinks$SampleId,tiLinks$ClusterId,sep='_')
tiLinks$TILenBucket=2**round(log(tiLinks$TILength,2))
View(tiLinks)
# View(tiLinks %>% group_by(HasSGL=grepl('SGL',ClusterDesc),LocationType) %>% count())

View(tiLinks %>% group_by(LocationType) %>% count())
View(tiLinks %>% group_by(LocationType,HasOverlaps=OverlapCount>0) %>% count())
print(plot_length_facetted(tiLinks,"LocationType!='Unclear'&ClusterCount>3",'TILenBucket,LocationType','TILenBucket','LocationType','TI Length by LocationType'))

# factoring in whether TIs overlap other TIs or not
print(plot_length_facetted(tiLinks,"LocationType!='Unclear'&ClusterCount>3&OverlapCount==0",'TILenBucket,LocationType','TILenBucket','LocationType','TI Length by LocationType, no overlaps'))
print(plot_length_facetted(tiLinks,"LocationType!='Unclear'&ClusterCount>3&OverlapCount>0",'TILenBucket,LocationType','TILenBucket','LocationType','TI Length by LocationType, no overlaps'))
print(plot_length_facetted(tiLinks,"LocationType!='Unclear'&ClusterCount>3",'TILenBucket,HasOverlaps=OverlapCount>0','TILenBucket','HasOverlaps','TI Length by whether has overlaps or not'))


ctLinks = tiLinks %>% filter(SampleClusterId %in% internalCT$SampleClusterId)
ctLinks = merge(ctLinks,simpleCT %>% select(SampleClusterId,ShortDBPerc),by='SampleClusterId',all.x=T)
nrow(ctLinks)
View(ctLinks)

# plot_length_facetted<-function(data, filterStr, groupByStr, bucketStr, facetStr, titleStr, logScale=T)

print(plot_length_facetted(ctLinks,'ShortDBPerc>=0&ClusterCount>3','TILenBucket,ShortDBPerc','TILenBucket','ShortDBPerc','TI Length by Short DB Percent'))
print(plot_length_facetted(ctLinks,"LocationType!='Unclear'&ClusterCount>3",'TILenBucket,LocationType','TILenBucket','LocationType','TI Length by LocationType'))

# factoring in short DB prevalence
print(plot_length_facetted(ctLinks,"LocationType!='Unclear'&ClusterCount>3&ShortDBPerc==0",'TILenBucket,LocationType','TILenBucket','LocationType','TI Length by LocationType'))
print(plot_length_facetted(ctLinks,"LocationType!='Unclear'&ClusterCount>3&ShortDBPerc>0.9",'TILenBucket,LocationType','TILenBucket','LocationType','TI Length by LocationType'))


# analysis of DB lengths between TIs where predominantly not short
dbLengthsStart = ctLinks %>% filter(DBLenStart>-1000) %>% mutate(DBLen=DBLenStart)
dbLengthsEnd = ctLinks %>% filter(DBLenEnd>-1000) %>% mutate(DBLen=DBLenEnd)
dbLengthsCombined = rbind(dbLengthsStart,dbLengthsEnd) %>% mutate(DBLenBucket = ifelse(DBLen>0,2**round(log(DBLen,2)),0))
View(dbLengthsCombined)

print(plot_length_facetted(dbLengthsCombined,'ClusterCount>=3&(ShortDBPerc<=0.2|ShortDBPerc>=0.9)','DBLenBucket,MostlyShort=ShortDBPerc>=0.9','DBLenBucket','MostlyShort','DB Length'))
print(plot_length_facetted(dbLengthsCombined,'ClusterCount>=3&(ShortDBPerc<=0.1|ShortDBPerc>=0.9)','DBLenBucket,ShortDBPerc','DBLenBucket','ShortDBPerc','DB Length'))

dbLengths = ctLinks %>% group_by(SampleId,ClusterId) %>% summarise(LinksCount=n(),
                                                                   NegDB=sum(ifelse(DBLenStart>-1000&DBLenStart<0,1,0)+ifelse(DBLenEnd>-1000&DBLenEnd<0,1,0))/2,
                                                                   DB_0_100=sum(ifelse(DBLenStart>=0&DBLenStart<100,1,0)+ifelse(DBLenEnd>=0&DBLenEnd<100,1,0))/2,
                                                                   DB_100_1000=sum(ifelse(DBLenStart>=100&DBLenStart<1e3,1,0)+ifelse(DBLenEnd>=100&DBLenEnd<1e3,1,0))/2,
                                                                   DB_1_10K=sum(ifelse(DBLenStart>=1e3&DBLenStart<1e4,1,0)+ifelse(DBLenEnd>=1e3&DBLenEnd<1e4,1,0))/2,
                                                                   DB_10_100K=sum(ifelse(DBLenStart>=1e4&DBLenStart<1e5,1,0)+ifelse(DBLenEnd>=1e4&DBLenEnd<1e5,1,0))/2,
                                                                   DB_100K_Plus=sum(ifelse(DBLenStart>=1e5,1,0)+ifelse(DBLenEnd>=1e5,1,0))/2)

View(dbLengths)
View(dbLengths %>% filter(NegDB))
View(dbLengths %>% group_by())


# Shattering with Replication
repClusters = ctClusters %>% filter(FullyChained=='true'&Foldbacks==0&ResolvedType=='SimpleChain'&MaxCopyNumber<=6&ChainCount==1&OverlapTIs>0&ChainEndsAway==1)
nrow(repClusters)
View(repClusters)
View(repClusters %>% select(SampleId,ClusterId,ClusterCount,ChainCount,IntTIs,ExtTIs,IntTIsWithGain,ExtTIsWithGain,OverlapTIs,DSBs,ShortDSBs,
                         ChainEndsFace,ChainEndsAway,ArmCount,OriginArms,FragmentArms,TotalLinks,AssemblyLinks,ShortTIRemotes,MinCopyNumber,MaxCopyNumber))

# simpleCT = simpleCT %>% mutate(ShortDBPerc=round(ShortDSBs/DSBs,1),)


# allowing foldbacks but trying to exclude BFB
repComplexClusters = ctClusters %>% filter(FullyChained=='true'&ResolvedType=='ComplexChain'&MaxCopyNumber<=6&ChainCount==1&OverlapTIs>0&ChainEndsAway==1)
nrow(repComplexClusters)
View(repComplexClusters)



# clusters with lots of low-count overlaps
View(head(tiLinks,1000))
tiOverlaps = tiLinks %>% filter(ClusterCount>3&ClusterCount<=100) %>% group_by(SampleId,ClusterId,ResolvedType,FullyChained) %>% 
  summarise(LinkCount=n(),
            ChainCount=n_distinct(ChainId),
            OverlapTotal=sum(OverlapCount),
            OverlappingLinks=sum(OverlapCount>0),
            SampleClusterId=first(SampleClusterId))

tiOverlaps$OverlapRatio = round(tiOverlaps$OverlappingLinks/tiOverlaps$LinkCount,2)
View(tiOverlaps)
tiOverlaps = tiOverlaps %>% ungroup()

repComplexClusters2 = merge(ctClusters, tiOverlaps %>% filter(OverlapRatio>=0.5) %>% select(SampleClusterId,OverlapTotal,OverlappingLinks,OverlapRatio),
                            by='SampleClusterId',all.y=T)
View(repComplexClusters2)

View(tiLinks %>% filter(SampleId=='COLO829T'))




