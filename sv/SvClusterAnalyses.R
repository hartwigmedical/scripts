

clusters = read.csv('~/logs/SVA_CLUSTERS.csv')
# clusters = read.csv('~/data/sv/fusions/SVA_CLUSTERS.csv')
rm(clusters)


#####
## CLUSTER TYPES


cleanClusters = clusters %>% filter(Subclonal=='false'&SuperType!='ARTIFACT')
totalClusters = nrow(cleanClusters)
totaSVs = sum(cleanClusters$ClusterCount)

View(cleanClusters %>% group_by(SuperType,ResolvedType,Synthetic) %>% 
       summarise(Clusters=n(),
                 Percent=round(n()/totalClusters,3),
                 SvCount=sum(ClusterCount),
                 SvPercent=round(sum(ClusterCount)/totaSVs,3)))

resolvedTypes1 = cleanClusters %>% filter(!grepl('SGL_PAIR',ResolvedType)) %>%
  group_by(SuperType,ResolvedType,IsSynthetic=ifelse(Synthetic=='true','Synthetic','Simple')) %>% count() %>% spread(IsSynthetic,n)

resolvedTypes1[is.na(resolvedTypes1)] = 0
resolvedTypes1 = resolvedTypes1 %>% mutate(SyntheticPerc=round(Synthetic/(Simple+Synthetic),3))
View(resolvedTypes1)

View(resolvedTypes1 %>% filter(SuperType=='SIMPLE'|SuperType=='BREAK_PAIR'))


##### 
## Synthetic DEL and DUP clusters

simpleClusters = cleanClusters %>% filter(Synthetic=='true'&ResolvedType %in% c('DUP','DEL')) %>% 
  mutate(TILocations=ifelse(TotalTIs==IntTIs,'Internal',ifelse(TotalTIs==ExtTIs,'External','Mixed')))
nrow(simpleClusters)
View(simpleClusters)

View(simpleClusters %>% group_by(ResolvedType,ClusterCount) %>% count() %>% spread(ResolvedType,n))
View(simpleClusters %>% group_by(ResolvedType,TILocations,ClusterCount) %>% count() %>% spread(TILocations,n))
View(simpleClusters %>% group_by(ResolvedType,ArmCount,ClusterCount) %>% count() %>% spread(ArmCount,n))


####
## CHROMOTHRIPTIC-type clusters

ctClusters = cleanClusters %>% filter(ResolvedType %in% c('DUP','DEL','COMPLEX')) %>% 
  filter(SglCount==0&InfCount==0&ChainCount==1&(ChainEndsFace==1|ChainEndsAway==1)&HasReplicated=='false'&Foldbacks==0&MaxCopyNumber<=2)

View(ctClusters)
nrow(ctClusters)

# metrics - DSBs, short TIs, arm count, cluster count, external vs internal TIs
ctClusters = ctClusters %>% mutate(ClusterSize=2**round(log(ClusterCount,2)),
                                   TILocations=ifelse(TotalTIs==IntTIs,'Internal',ifelse(TotalTIs==ExtTIs,'External','Mixed')),
                                   DsbPerc=round(DSBs/ClusterCount,1),
                                   ShortDbPerc=round(ShortDSBs/ClusterCount,1),
                                   DsbBucket=2**round(log(DSBs,2)),
                                   ShortTiBucket=2**round(log(ShortTIs,2)),
                                   ExtTIPerc=round(ExtShortTIs/ShortTIs,1),
                                   ShortTIPerc=round(ShortTIs/TotalTIs,2))

View(ctClusters %>% group_by(ArmCount,ClusterSize) %>% count() %>% spread(ArmCount,n))
View(ctClusters %>% group_by(TILocations,ClusterSize) %>% count() %>% spread(TILocations,n))
View(ctClusters %>% group_by(ShortDbPerc,ShortTiPerc) %>% count() %>% spread(ShortDbPerc,n))
View(ctClusters %>% group_by(ArmCount,TILocations) %>% count() %>% spread(ArmCount,n))
View(ctClusters %>% filter(ShortTiPerc<=0.8) %>% group_by(ArmCount,TILocations) %>% count() %>% spread(ArmCount,n))

View(ctClusters %>% filter(ExtTIs==ExtShortTIs) %>% group_by(ClusterSize,TILocations,ChainEndsAway) %>% count() %>% spread(ClusterSize,n))
View(ctClusters %>% filter(ExtTIs==ExtShortTIs) %>% group_by(ClusterSize,TILocations,ChainEndsAway) %>% count() %>% spread(ClusterSize,n))

# 
nrow(ctClusters %>% filter(IntTIs+ExtTIs<TotalTIs))

# localise shattering, can have remote TIs
nrow(ctClusters %>% filter(ExtTIs==ExtShortTIs))

ctClusters = ctClusters %>% mutate(DsbBucket=2**round(log(DSBs,2)),
                                   ShortTiBucket=2**round(log(ShortTIs,2)),
                                   ExtTIPerc=round(ExtShortTIs/ShortTIs,1),
                                   ShortTIPerc=round(ShortTIs/TotalTIs,2))

View(ctClusters)

print(ggplot(ctClusters %>% filter(ExtTIs==ExtShortTIs), aes(x=ClusterCount, y=DSBs)) 
      + geom_point(aes(size=ShortTIs, color=ExtTIPerc)) 
      + theme_gray()
      + theme(panel.grid.major = element_line(), panel.grid.minor = element_line())
      + labs(subtitle="Shattering by #DSBs and short TIs", y="DSBs", x="Cluster Size"))


print(ggplot(ctClusters, aes(x=ClusterCount, y=ShortTIs)) 
      + geom_point() 
      + theme_gray()
      + theme(panel.grid.major = element_line(), panel.grid.minor = element_line())
      + labs(subtitle="Short TIs vs Cluster Size", y="Short TI Count", x="Cluster Size"))


# shattering involving more than 1 arm
print(ggplot(ctClusters, aes(x=ClusterCount, y=ArmCount)) 
      + geom_point(aes(size=DSBs, color=ShortTIPerc)) 
      + theme_gray()
      + theme(panel.grid.major = element_line(), panel.grid.minor = element_line())
      # + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0,size=10))
      # + scale_x_log10()
      + labs(subtitle="Shattering by #DSBs and short TIs", y="DSBs", x="Cluster Size"))



ctSummaryData2 = ctSummaryData %>% group_by(ClusterSize,ShortTiBucket,DsbBucket) %>% count()
View(ctSummaryData2)



#####
## Short TIs for different cluster type
print(ggplot(ctClusters, aes(x=ClusterCount, y=ShortTIs)) 
      + geom_point() 
      + theme_gray()
      + theme(panel.grid.major = element_line(), panel.grid.minor = element_line())
      + labs(subtitle="Short TIs vs Cluster Size", y="Short TI Count", x="Cluster Size"))


compClusters = cleanClusters %>% filter(ClusterCount<=500&ResolvedType %in% c('DUP','DEL','COMPLEX')) %>% 
  mutate(ClusterType=ifelse(Foldbacks==0&MaxCopyNumber<=2&HasReplicated=='false','Simple',ifelse(Foldbacks==0,'Replicated','BFB')),
         PloidyBucket=ifelse(MaxCopyNumber<=64,2**round(log(MaxCopyNumber,2)),64),
         FoldbackBucket=ifelse(Foldbacks<=32,2**round(log(Foldbacks,2)),32))

compClusterSummary = compClusters %>% group_by(SampleId,ClusterId,ClusterCount,ClusterType,PloidyBucket,FoldbackBucket) %>% 
  summarise(ShortTIs=sum(ShortTIs))

print(ggplot(compClusterSummary, aes(x=ClusterCount, y=ShortTIs)) 
      + geom_point() 
      + facet_wrap(~ClusterType)
      + labs(subtitle="Short TIs vs Cluster Size", y="Short TI Count", x="Cluster Size"))

print(ggplot(compClusterSummary, aes(x=ClusterCount, y=ShortTIs)) 
      + geom_point() 
      + facet_wrap(~FoldbackBucket)
      + labs(subtitle="Short TIs vs Cluster Size by Foldback count", y="Short TI Count", x="Cluster Size"))

print(ggplot(compClusterSummary %>% filter(PloidyBucket>=1), aes(x=ClusterCount, y=ShortTIs)) 
      + geom_point() 
      + facet_wrap(~PloidyBucket)
      + labs(subtitle="Short TIs vs Cluster Size by Ploidy Bucket", y="Short TI Count", x="Cluster Size"))

######
## Clusters with replication of any kind except foldbacks

repClusters = cleanClusters %>% filter(ResolvedType=='COMPLEX'&Foldbacks==0&(OverlapTIs>0|HasReplicated=='true')&SglCount==0&InfCount==0
                                       &ChainCount==1&(ChainEndsFace==1|ChainEndsAway==1)&MaxCopyNumber<=6)

View(repClusters)
nrow(repClusters) # only 731

#####
## Foldback cluster

fbClusters = cleanClusters %>% filter(ResolvedType=='COMPLEX'&Foldbacks>0&SglCount==0&InfCount==0&ChainCount==1)
nrow(fbClusters)

fbClusters = fbClusters %>% mutate(ClusterSize=2**round(log(ClusterCount,2)),
                                   ShortTIPerc=round(ShortTIs/TotalTIs,2))

print(ggplot(fbClusters, aes(x=ClusterCount, y=ArmCount)) 
      + geom_point(aes(size=Foldbacks, color=ShortTIPerc)) 
      + theme_gray()
      + theme(panel.grid.major = element_line(), panel.grid.minor = element_line())
      + labs(subtitle="Foldbacks and short TIs", y="Arm Count", x="Cluster Size"))

anyFbClusters = cleanClusters %>% filter(ResolvedType=='COMPLEX'&Foldbacks>0&SglCount==0&InfCount==0) %>% 
  mutate(ClusterSize=2**round(log(ClusterCount,2)),
         ShortTIPerc=round(ShortTIs/TotalTIs,2))

nrow(anyFbClusters)

print(ggplot(anyFbClusters %>% filter(ClusterCount<=500), aes(x=ClusterCount, y=ArmCount)) 
      + geom_point(aes(size=Foldbacks)) 
      + theme_gray()
      + scale_x_log10()
      + theme(panel.grid.major = element_line(), panel.grid.minor = element_line())
      + labs(subtitle="Foldbacks and short TIs"))

print(ggplot(anyFbClusters %>% filter(ClusterCount<=500), aes(x=ClusterCount, y=MaxCopyNumber)) 
      + geom_point(aes(size=Foldbacks)) 
      + theme_gray()
      + scale_x_log10()
      + scale_y_log10()
      + theme(panel.grid.major = element_line(), panel.grid.minor = element_line())
      + labs(subtitle="Foldbacks and short TIs"))

print(ggplot(anyFbClusters %>% filter(ClusterCount<=500&MaxCopyNumber<=50&Foldbacks<=20), aes(x=Foldbacks, y=MaxCopyNumber)) 
      + geom_point(aes(size=ClusterSize)) 
      + theme_gray()
      # + scale_x_log10()
      + scale_y_log10()
      + theme(panel.grid.major = element_line(), panel.grid.minor = element_line())
      + labs(subtitle="Foldbacks vs Max Ploidy", y="Max Ploidy"))




sampleCancerTypes = read.csv('~/data/sigs/sample_cancer_types.csv')
View(sampleCancerTypes)

View(svData %>% filter(SampleId=='CPCT02010762T'&ChrStart==12))

View(svData %>% filter(SampleId=='CPCT02010762T') %>% group_by(Type) %>% count())


# localised shattering allowing for some external or remote TI inserts, but all short


##### 
## CLUSTER STATS
compClusters = clusters %>% filter(ClusterCount>=3&ResolvedType=='COMPLEX'&ClusterCount<100) # approx 28.5K
nrow(compClusters)

View(compClusters %>% filter(DSBs>ClusterCount))
View(compClusters %>% filter(Foldbacks>ClusterCount))
View(compClusters %>% filter(TotalLinks>ClusterCount))

# prepare summary data
compClusters = compClusters %>% mutate(ClusterBucket=2**round(log(ClusterCount,2)),
                        SglPerc=round((SglCount+InfCount)/ClusterCount,1),
                        DbPerc=round(DSBs/ClusterCount,1),
                        FoldbackPerc=round(Foldbacks/ClusterCount,1),
                        ShortDbPerc=round(ShortDSBs/ClusterCount,1),
                        ShortTiPerc=round(AcTotalTIs/ClusterCount,1),
                        MaxPloidyBucket=2**round(log(MaxCopyNumber,2)),
                        OverlapTIs=ifelse(TotalLinks>0,round(OverlapTIs/TotalLinks,1),0))

View(compClusters %>% filter(SglPerc==0) %>% group_by(ShortDbPerc,ShortTiPerc) %>% count() %>% spread(ShortDbPerc,n))
View(compClusters %>% filter(SglPerc==0) %>% group_by(OverlapTIs,ShortTiPerc) %>% count() %>% spread(OverlapTIs,n))
View(compClusters %>% filter(SglPerc==0) %>% group_by(ClusterBucket,ShortTiPerc) %>% count() %>% spread(ClusterBucket,n))
View(compClusters %>% filter(SglPerc==0) %>% group_by(ClusterBucket,FoldbackPerc) %>% count() %>% spread(ClusterBucket,n))
View(compClusters %>% filter(SglPerc==0) %>% group_by(MaxPloidyBucket,FoldbackPerc) %>% count() %>% spread(MaxPloidyBucket,n))

# TIs with overlap (> 1kb) are fairly rare
View(compClusters %>% group_by(OverlapTIs) %>% count())

# about 2/3 have no SGLs or NONEs
View(compClusters %>% group_by(SglPerc) %>% count())


View(compClusters %>% filter(SglPerc==0&ClusterCount>50&ShortTiPerc>=0.9))

View(compClusters)

view_cluster_sv('CPCT02020398T',198)

# clusters with shattering and multiple short DSBS




ctClusters = clusters %>% filter(ClusterCount>1&ClusterCount<100) # approx 111K
ctClusters$SampleClusterId = paste(ctClusters$SampleId,ctClusters$ClusterId,sep='_')
nrow(ctClusters)
View(ctClusters)




# simple (no foldbacks), little to no replication, 1 or 2 chains = 5K
simpleClusters = ctClusters %>% filter(FullyChained=='true'&Foldbacks==0&MaxCopyNumber<=2.5&ChainCount<=2&HasReplicated=='false')
nrow(simpleClusters)
View(simpleClusters)

svData$SampleClusterId = paste(svData$SampleId,svData$ClusterId,sep='_')
scSvData = svData %>% filter(SampleClusterId %in% simpleClusters$SampleClusterId)
nrow(scSvData)
rm(scSvData)

# prepare info on DB lengths
dbDataStart = scSvData %>% filter(DBLenStart>-1000) %>% select(SampleId,ClusterId,Id,Type,ResolvedType,ClusterCount,DbLength=DBLenStart,SampleClusterId)
dbDataEnd = scSvData %>% filter(DBLenEnd>-1000) %>% select(SampleId,ClusterId,Id,Type,ResolvedType,ClusterCount,DbLength=DBLenEnd,SampleClusterId)
dbData = rbind(dbDataStart,dbDataEnd)
# dbData = dbData %>% mutate(DbLenBucket=ifelse(DbLength==0,0,2**round(log(DbLength,2))))
View(dbData)

print(plot_length_facetted(dbData,"ResolvedType!='SIMPLE'&DbLength>=-100&DbLength<=100",'DbLength,ResolvedType','DbLength','ResolvedType','DB Length by ResolvedType',F,T))

# prepare cluster level summaries
ctClusterSummary = dbData %>% group_by(SampleClusterId,SampleId,ClusterId,ResolvedType,ClusterCount) %>% 
  summarise(DbCount=n()/2,
            InvCount=sum(Type=='INV'),
            ShortDbCount=sum(DbLength<=100)/2) %>% 
  mutate(DbPerc=round(DbCount/ClusterCount,3),
         ShortDbPerc=round(ShortDbCount/ClusterCount,3),
         ShortDbBucket=round(ShortDbCount/ClusterCount/0.2)*0.2) %>% ungroup()

View(ctClusterSummary)

dbData = merge(dbData,ctClusterSummary %>% select(SampleClusterId,ShortDbPerc,ShortDbBucket),by='SampleClusterId',all.x=T)
dbData = dbData %>% mutate(ClusterSize=2**round(log(ClusterCount,2)))

print(plot_length_facetted(dbData,"ResolvedType=='COMPLEX'&DbLength>=-100&DbLength<=100",'DbLength,ShortDbBucket','DbLength','ShortDbBucket','DB Length by Short DB %',F,T))
print(plot_length_facetted(dbData,"ShortDbBucket>=0.95&DbLength>=-100&DbLength<=100",'DbLength,ClusterSize','DbLength','ClusterSize','DB Length by Short DB %',F,T))
print(plot_length_facetted(dbData,"ClusterCount>=3&ClusterCount<=10&ShortDbBucket>=1&DbLength>=-100&DbLength<=100",'DbLength,ClusterCount','DbLength','ClusterCount','DB Length by Short DB %',F,T))












View(simpleClusters %>% select(SampleId,ClusterId,ClusterCount,ChainCount,IntTIs,ExtTIs,IntTIsWithGain,ExtTIsWithGain,OverlapTIs,DSBs,ShortDSBs,
                               ChainEndsFace,ChainEndsAway,ArmCount,OriginArms,FragmentArms,TotalLinks,AssemblyLinks,ShortTIRemotes,MinCopyNumber,MaxCopyNumber))

# 1. Single chain, no overlaps, chain ends facing away - count = 2300
# shows classic shattering, but also examples of failed replication, with no short DBs
simpleCT = ctClusters %>% filter(FullyChained=='true'&Foldbacks==0&ResolvedType=='COMPLEX'&MaxCopyNumber<=2.5&ChainCount==1&IntTIsWithGain==0&OverlapTIs==0&ChainEndsAway==1)
simpleCT = simpleCT %>% mutate(ShortDBPerc=round(ShortDSBs/DSBs,1),)
nrow(simpleCT)

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

tiLinks = tiLinks %>% mutate(SampleClusterId = paste(SampleId,ClusterId,sep='_'),
                             TILenBucket=2**round(log(TILength,2)),
                             ClusterSize=ifelse(ClusterCount<=3,'Small',ifelse(ClusterCount<=10,'Medium','Large')),
                             NextSvDist=ifelse((NextSvDist==-1|NextSvDist<NextClusteredSvDist)&NextClusteredSvDist>-1,NextClusteredSvDist,NextSvDist),
                             NextDist=ifelse(NextSvDist<=0,NextSvDist,10**round(log(NextSvDist,10))),
                             NextClusterDist=ifelse(NextClusteredSvDist<=0,NextClusteredSvDist,10**round(log(NextClusteredSvDist,10))))

View(tiLinks %>% select(SampleId,ClusterId,NextSvDist,NextDist,NextClusteredSvDist,NextClusterDist))

View(tiNoSglLinks %>% group_by(NextDist,NextClusterDist) %>% count() %>% spread(NextClusterDist,n))
View(tiNoSglLinks %>% filter(LocationType=='External'|LocationType=='Remote') %>% group_by(NextDist,NextClusterDist) %>% count() %>% spread(NextClusterDist,n))
View(tiNoSglLinks %>% filter(ResolvedType=='DUP_EXT_TI'|ResolvedType=='DEL_EXT_TI') %>% group_by(NextDist,NextClusterDist) %>% count() %>% spread(NextClusterDist,n))

nrow(tiLinks %>% filter(NextSvDist==0&NextClusteredSvDist>0))

View(tiLinks %>% filter(ResolvedType!='LINE'&(LocationType=='External'|LocationType=='Remote')&NextDist==-1&NextClusterDist==-1))

       
View(tiLinks %>% filter(NextSvDist>-1&NextSvDist<NextClusteredSvDist) %>% 
       select(SampleId,ClusterId,NextSvDist,NextDist,NextClusteredSvDist,NextClusterDist,everything()))

View(svData %>% filter(SampleId=='CPCT02030387T'&ClusterId==30))
View(svData %>% filter(SampleId=='CPCT02030387T'&(ChrStart==3|ChrEnd==3)))

# exclude TIs involving a single for length analysis
sglSvData = svData %>% filter(Type=='NONE'|Type=='SGL')
nrow(tiLinks)
tiNoSglLinks = tiLinks %>% filter(!(Id1 %in% sglSvData$Id)&!(Id2 %in% sglSvData$Id))
nrow(tiNoSglLinks)

colnames(tiLinks)
View(tiLinks)
# View(tiLinks %>% group_by(HasSGL=grepl('SGL',ClusterDesc),LocationType) %>% count())

View(tiLinks %>% group_by(LocationType) %>% count())
View(tiLinks %>% group_by(LocationType,HasOverlaps=OverlapCount>0) %>% count())
print(plot_length_facetted(tiLinks,"LocationType!='Unclear'&ClusterCount>3",'TILenBucket,LocationType','TILenBucket','LocationType','TI Length by LocationType'))

# factoring in whether TIs overlap other TIs or not
print(plot_length_facetted(tiLinks,"LocationType!='Unclear'&ClusterCount>3&OverlapCount==0",'TILenBucket,LocationType','TILenBucket','LocationType','TI Length by LocationType, no overlaps'))
print(plot_length_facetted(tiLinks,"LocationType!='Unclear'&ClusterCount>3&OverlapCount>0",'TILenBucket,LocationType','TILenBucket','LocationType','TI Length by LocationType, no overlaps'))
print(plot_length_facetted(tiLinks,"LocationType!='Unclear'&ClusterCount>3",'TILenBucket,HasOverlaps=OverlapCount>0','TILenBucket','HasOverlaps','TI Length by whether has overlaps or not'))

# remote and external with no SVs near and no overlap
print(plot_length_facetted(tiLinks,"(LocationType=='Remote'|LocationType=='External')&ClusterCount>3&NextSVDistance>5e3&OverlapCount==0",
                           'TILenBucket,LocationType','TILenBucket','LocationType','Remote & isolated TI Length'))

# long TI peak grows with larger clusters
print(plot_length_facetted(tiLinks,"ResolvedType!='LINE'&(LocationType=='Remote'|LocationType=='External')&OverlapCount==0",
                           'TILenBucket,ClusterSize','TILenBucket','ClusterSize','Remote & isolated TI Length by ClusterSize'))

# no SVs on same arm is strongest identification of short TIs, followed by proximity (likely a DB)
print(plot_length_facetted(tiNoSglLinks,"ResolvedType!='LINE'&(LocationType=='Remote'|LocationType=='External')&OverlapCount==0&TraversedSVCount==0",
                           'TILenBucket,NearestDist','TILenBucket','NearestDist','Remote & isolated TI Length by Next SV Dist'))

print(plot_length_facetted(tiNoSglLinks,"ResolvedType!='LINE'&(LocationType=='Remote'|LocationType=='External')",
                           'TILenBucket,NextDist','TILenBucket','NextDist','Remote & isolated TI Length by Next Dist'))

print(plot_length_facetted(tiNoSglLinks,"ResolvedType!='LINE'&(LocationType=='Remote'|LocationType=='External')",
                           'TILenBucket,NextClusterDist','TILenBucket','NextClusterDist','Remote & isolated TI Length by Next clustered SV Dist'))

print(plot_length_facetted(tiNoSglLinks,"ResolvedType!='LINE'&LocationType=='Remote'",
                           'TILenBucket,NextDist','TILenBucket','NextDist','Remote & isolated TI Length by Next Dist'))

print(plot_length_facetted(tiNoSglLinks,"ResolvedType!='LINE'&LocationType=='Remote'",
                           'TILenBucket,NextClusterDist','TILenBucket','NextClusterDist','Remote & isolated TI Length by Next clustered SV Dist'))


View(tiNoSglLinks %>% filter(NextSVDistance>-1&NextSVDistance<=1))
View(tiNoSglLinks %>% filter(NextSVDistance>=1e6&NextSVDistance<=1e7&TILength>1e6))



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




