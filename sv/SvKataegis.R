## KATEGIS ANALYSIS                 

load('~/data/hmf_cohort_may_2019.RData')
sampleCancerTypes = highestPurityCohort %>% select(SampleId=sampleId,CancerType=cancerType)

snvKataegisData = read.csv('~/data/sv/prod_kataegis_data.csv')

View(snvKataegisData)
nrow(snvKataegisData)
nrow(snvKataegisData %>% group_by(SampleId,KataegisId) %>% count())
nrow(snvKataegisData %>% group_by(SampleId) %>% count())
View(snvKataegisData %>% group_by(LengthBucket,SnvCount) %>% count())

snvKataegisData = snvKataegisData %>% mutate(Length=MaxPosition-MinPosition,
                                             LengthBucket=2**round(log(Length,2))) 

View(snvKataegisData %>% group_by(LengthBucket) %>% count())


View(snvKataegisData %>% group_by(SampleId) %>% count())
View(snvKataegisData %>% group_by(SampleId) %>% summarise(Regions=n(),SNVs=sum(SnvCount)))

View(snvKataegisData %>% group_by(SampleId,Chromosome) %>% count())
View(snvKataegisData %>% group_by(SampleId,Chromosome) %>% count() %>% group_by(SampleId) %>% count())

# load Kataegis site SV annotations from LINX
svKatData = read.csv('~/data/sv/LNX_KATAEGIS.csv')

# temp to get the correct kat site positions until next Linx re-run
# svKatData = merge(svKatData,snvKataegisData %>% select(SampleId,KataegisId,KatLength=Length,KatLengthBucket=LengthBucket),by=c('SampleId','KataegisId'))

svKatData = svKatData %>% filter(SampleId %in% sampleCancerTypes$SampleId)
svKatData = merge(svKatData,sampleCancerTypes,by='SampleId',all.x=T)

# get SNV sample totals
sampleSnvTotals = read.csv('~/data/sample_snv_totals.csv')
View(sampleSnvTotals)
svKatData = merge(svKatData,sampleSnvTotals,by='SampleId',all.x=T)
# View(svKatData %>% filter(is.na(SampleSnvTotal)))

# join to SV data to look at local topology data and cluster info
svData = read.csv('~/data/sv/drivers/LNX_SVS.csv')
svDataKat = svData %>% filter(SampleId %in% sampleCancerTypes$SampleId)

katSvMatched = svKatData %>% filter(SvId>=0)
nrow(katSvMatched)

# merge in relevant SV fields
katSvMatched = merge(katSvMatched %>% mutate(Id=SvId), # ,SnvCount,SvIsStart,KatPosStart=PosStart,KatPosEnd=PosEnd,Distance,CancerType
                     svDataKat %>% select(SampleId,Id,Type,PosStart,PosEnd,LocTopTypeStart,LocTopTypeEnd,LocTopTIStart,LocTopTIEnd,
                                       ClusterId,ClusterCount,ResolvedType,HomologyStart,HomologyEnd),
                     by=c('SampleId','Id'),all.x=T)

katSvMatched = katSvMatched %>% mutate(DistBucket=2**round(log(Distance,2)),
                                       DistBucket2=round(Distance/50)*50,
                                       DistBucket3=round(Distance/1000)*10000,
                                       SvLength=ifelse(Type=='BND'|Type=='INF'|Type=='SGL',-1,PosEnd-PosStart),
                                       SvLengthBucket=ifelse(SvLength>0,2**round(log(SvLength,2)),2),
                                       SvLengthType = classify_del_dup_length(Type,SvLength),
                                       LocTopType=ifelse(SvIsStart=='true',as.character(LocTopTypeStart),as.character(LocTopTypeStart)),
                                       LocTopTiCount=ifelse(SvIsStart=='true',as.character(LocTopTIStart),as.character(LocTopTIEnd)),
                                       ClusterSize=2**round(log(ClusterCount,2)),
                                       HomologyLen=ifelse(SvIsStart=='true',stri_length(HomologyStart),stri_length(HomologyEnd)),
                                       HomologyLenBucket=ifelse(HomologyLen>10,10,HomologyLen),
                                       SnvCountBucket=2**round(log(SnvCount,2)))

View(katSvMatched)

# SV proximity to Kataegis site
View(katSvMatched %>% group_by(DistBucket) %>% count())

nrow(katSvMatched %>% filter(SampleSnvTotal<=highMutLoadThreshold&DistBucket<=1e3))/nrow(svKatData %>% filter(SampleSnvTotal<=highMutLoadThreshold))
nrow(katSvMatched)
nrow(katSvMatched %>% filter(DistBucket>1e3))
View(katSvMatched)
highMutLoadThreshold=5e4

print(ggplot(data = katSvMatched %>% group_by(DistBucket,MutationalLoad=ifelse(SampleSnvTotal>highMutLoadThreshold,'HIGH','LOW')) %>% count(), aes(x=DistBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~MutationalLoad)
      + labs(title = "SV Breakend Proximity to Kataegis Site - by Mutational Load (+/- 50K)",x='Distance',y='Count'))

# filter out high mutational load samples for further analysis
# svKatData = svKatData %>% filter(SampleSnvTotal<=highMutLoadThreshold)
katSvMatched = katSvMatched %>% filter(SampleSnvTotal<=highMutLoadThreshold)
nrow(katSvMatched)

## DEBUG

View(svKatData)
View(svKatData %>% group_by(Matched=SvId>=0) %>% count())
View(svKatData %>% filter(SnvCount>5) %>% group_by(CancerType,Matched=SvId>=0) %>% count() %>% spread(Matched,n) %>% mutate(MatchPerc=round(`TRUE`/(`TRUE`+`FALSE`),2)))
View(svKatData %>% filter(SnvCount>0&Distance<1e3) %>% group_by(CancerType,Matched=SvId>=0) %>% count() %>% spread(Matched,n) %>% mutate(MatchPerc=round(`TRUE`/(`TRUE`+`FALSE`),2)))
View(svKatData %>% filter(SnvCount>0&Distance<1e3) %>% group_by(pmin(SnvCount,12),Matched=SvId>=0) %>% count() %>% spread(Matched,n) %>% mutate(MatchPerc=round(`TRUE`/(`TRUE`+`FALSE`),2)))

View(svKatData %>% filter(CancerType!='Skin',SnvCount>4) %>% group_by(SnvTotalBucket=2**round(log(SampleTotal,2)),Matched=SvId>=0) 
     %>% count() %>% spread(Matched,n))



# by Resolved Type
View(svDataKat)
resolvedTypeCounts = svDataKat %>% group_by(ResolvedType) %>% summarise(ResolvedTypeTotal=n())
View(resolvedTypeCounts)

katByResolvedTypes = katSvMatched %>% filter(Distance<1e3) %>% group_by(ResolvedType) %>% summarise(Count=n()) %>% 
  mutate(Percent=round(Count/nrow(katSvMatched),2)) %>% arrange(-Count)

katByResolvedTypes = merge(katByResolvedTypes,resolvedTypeCounts,by='ResolvedType',all.x=T)
View(katByResolvedTypes)

katByResolvedTypes = katByResolvedTypes %>% mutate(Percent=round(Count/nrow(katSvMatched),2),
                                                   ResolvedTypePerc=round(Count/ResolvedTypeTotal,4)) %>% arrange(-Count)

print(ggplot(data = katByResolvedTypes %>% filter(ResolvedTypeTotal>5e3), 
             aes(x=reorder(ResolvedType,-ResolvedTypePerc),y=ResolvedTypePerc*100))
      + geom_bar(stat = "identity", colour = "black")
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + labs(title = "% of SVs by Event Type within 1K of Kataegis Site",x='Event Type',y='% of SVs'))


# percent of complex clusters with at least one breakend near a kat site
compClusterCount = nrow(svDataKat %>% filter(ResolvedType=='COMPLEX') %>% group_by(SampleId,ClusterId) %>% count())
print(compClusterCount)

compClusterWithKatCount = nrow(katSvMatched %>% filter(Distance<1e3&ResolvedType=='COMPLEX') %>% group_by(SampleId,ClusterId) %>% count())
print(compClusterWithKatCount)
print(compClusterWithKatCount/compClusterCount)


# DELs and DUPs by length
delDupTotals = svDataKat %>% filter(ClusterCount==1&(Type=='DEL'|Type=='DUP')) %>% 
  mutate(Length=PosEnd-PosStart,
         SvLengthType=classify_del_dup_length(Type,Length))

delDupTotals = delDupTotals %>% group_by(SvLengthType) %>% summarise(SvTotal=n())
View(delDupTotals)

katByDelDupTypes = katSvMatched %>% filter(Distance<1e3&(Type=='DEL'|Type=='DUP')&ClusterCount==1) %>% group_by(SvLengthType) %>% summarise(Count=n())
katByDelDupTypes = merge(katByDelDupTypes,delDupTotals,by='SvLengthType',all.x=T)
View(katByDelDupTypes)

katByDelDupTypes = katByDelDupTypes %>% mutate(SvLengthTypePerc=round(Count/SvTotal,4)) %>% arrange(-Count)
View(katByDelDupTypes)

print(ggplot(data = katByDelDupTypes %>% filter(grepl('DEL',SvLengthType)),aes(x=reorder(SvLengthType,-SvLengthTypePerc),y=SvLengthTypePerc*100))
      + geom_bar(stat = "identity", colour = "black")
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + labs(title = "% of DEL SVs by Length Category within 1K of Kataegis Site",x='Length Type',y='% of SVs'))


# pull in replication timing
svRepTiming = read.csv('~/logs/LNX_SVS_rep_timing.csv')
svRepTiming = svRepTiming %>% select(SampleId,Id,RepOriginStart,RepOriginEnd)
write.csv(svRepTiming,'~/data/sv/sv_rep_timing.csv', row.names = F, quote = F)

katSvMatched = merge(katSvMatched,svRepTiming,by=c('SampleId','Id'),all.x=T)

View(katSvMatched %>% filter(is.na(RepOriginEnd)))
katSvMatched = katSvMatched %>% mutate(RepTiming=ifelse(SvIsStart=='true',RepOriginStart,RepOriginEnd),
                                       RepTimingBucket=0.1*round(RepTiming/0.1))

repTimingCounts = rbind(svRepTiming %>% select(RepTiming=RepOriginStart),svRepTiming %>% filter(Type!='INF'&Type!='SGL') %>% select(RepTiming=RepOriginEnd))
repTimingCounts = repTimingCounts %>% group_by(RepTimingBucket=0.1*round(RepTiming/0.1)) %>% summarise(RepTimingTotal=n())
View(repTimingCounts)

katSvMatched = merge(katSvMatched,repTimingCounts,by='RepTimingBucket',all.x=T)

View(katSvMatched %>% filter(Distance<1e3) %>% group_by(RepTimingBucket,RepTimingTotal) %>% count() %>% 
       mutate(Percent=round(n/nrow(katSvMatched),2),
              RepTimingPerc=round(n/RepTimingTotal,3)) %>% arrange(-n))

# Homology correlation
homologyCounts = rbind(svData %>% select(Type,ClusterCount,Homology=HomologyStart),
                       svData %>% filter(Type!='INF'&Type!='SGL') %>% select(Type,ClusterCount,Homology=HomologyEnd))

homologyCounts = homologyCounts %>% mutate(HomologyLen=stri_length(Homology),
                                           HomologyLenBucket=ifelse(HomologyLen>10,10,HomologyLen))

homologyTotals = homologyCounts %>% group_by(HomologyLenBucket) %>% summarise(HomologyLenTotal=n())
View(homologyTotals)
katSvMatchedHom1 = merge(katSvMatched,homologyTotals,by='HomologyLenBucket',all.x=T)

katSvMatchedHomSummary = katSvMatchedHom1 %>% filter(Distance<1e3) %>%
       group_by(HomologyLenBucket,HomologyLenTotal) %>% count() %>% 
       mutate(Percent=round(n/nrow(katSvMatched),2),
              HomologyLenPerc=round(n/HomologyLenTotal,4)) %>% arrange(-n)

print(ggplot(data = katSvMatchedHomSummary,aes(x=reorder(as.character(HomologyLenBucket),HomologyLenBucket),y=HomologyLenPerc*100))
      + geom_bar(stat = "identity", colour = "black")
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + labs(title = "% of SVs by Homology Length within 1K of Kataegis Site",x='Homology Length',y='% of SVs'))


homologyTotalsByType = homologyCounts %>% group_by(ClusterType=ifelse(ClusterCount==1&Type=='DEL','DEL',ifelse(ClusterCount>=3,'COMPLEX','OTHER')),
                                                   HomologyLenBucket) %>% summarise(HomologyLenTotal=n())
View(homologyTotalsByType)
katSvMatchedHom2 = merge(katSvMatched %>% mutate(ClusterType=ifelse(ClusterCount==1&Type=='DEL','DEL',ifelse(ClusterCount>=3,'COMPLEX','OTHER'))),
                         homologyTotalsByType,by=c('ClusterType','HomologyLenBucket'),all.x=T)

View(katSvMatchedHom2 %>% filter(Distance<1e3) %>%
       group_by(ClusterType,HomologyLenBucket,HomologyLenTotal) %>% count() %>% 
       mutate(Percent=round(n/nrow(katSvMatched),2),
              HomologyLenPerc=round(n/HomologyLenTotal,4)) %>% arrange(ClusterType,-n))


# print(plot_length_facetted(katSvMatched,"ClusterCount==1&Distance<5e3",'DistBucket,Type','DistBucket','Type','Kataegis Distance vs Type, Cluster-1',F))
print(plot_length_facetted(katSvMatched,"ClusterCount==1&Distance<1e3",'DistBucket2,Type','DistBucket2','Type','Kataegis Distance vs Type, Cluster-1',F))
print(plot_length_facetted(katSvMatched,"Distance<5e3",'DistBucket,Type','DistBucket','Type','Kataegis Distance vs Type - all Cluster Counts',F))
print(plot_length_facetted(katSvMatched,"Distance<1e3",'DistBucket2,Type','DistBucket2','Type','Kataegis Distance vs Type - all Cluster Counts',F))
print(plot_length_facetted(katSvMatched,"Type!='INF'&ClusterCount>=3&Distance<1e3",'DistBucket2,Type','DistBucket2','Type','Kataegis Distance vs Type - Complex Clusters',F))

# by SV distance for DELs
print(plot_length_facetted(katSvMatched,"ClusterCount==1&Type=='DEL'&Distance<1e3",'DistBucket2,SvLengthType','DistBucket2','SvLengthType','Kataegis Distance vs DEL length',F))
# print(plot_length_facetted(katSvMatched,"ClusterCount==1&Type=='DEL'&Distance<1e3",'DistBucket2,SvLengthBucket','DistBucket2','SvLengthBucket','Kataegis Distance vs DEL length',F))
print(plot_length_facetted(katSvMatched,"ClusterCount==1&Type=='DUP'&Distance<1e3",'DistBucket2,SvLengthType','DistBucket2','SvLengthType','Kataegis Distance vs DUP length',F))
print(plot_length_facetted(katSvMatched,"ClusterCount>=3&Type=='INV'&Distance<1e3",'DistBucket2,SvLengthBucket','DistBucket2','SvLengthBucket','Kataegis Distance vs DEL length',F))
print(plot_length_facetted(katSvMatched,"SvLengthType=='DEL_SHORT'&Distance<1e3",'DistBucket2,ClusterSize','DistBucket2','ClusterSize','Kataegis Distance vs DEL length',F))

# by local topology
print(plot_length_facetted(katSvMatched,"Distance<5e3",'DistBucket2,LocTopType','DistBucket2','LocTopType','Kat-Distance vs LocTopType',F))
print(plot_length_facetted(katSvMatched,"Distance<1e3&LocTopType %in% c('TI_ONLY','DSB')",'DistBucket2,LocTopType','DistBucket2','LocTopType','Kat-Distance vs LocTopType',F))
print(plot_length_facetted(katSvMatched,"Distance<1e3&LocTopType!='ISOLATED_BE'",'DistBucket2,LocTopType','DistBucket2','LocTopType','Kat-Distance vs LocTopType',F))

# by local top TI count
print(plot_length_facetted(katSvMatched,"Distance<5e3",'DistBucket2,LocTopTiCount','DistBucket2','LocTopTiCount','Kat-Distance vs LocTop TI Count',F))

# by cluster size
print(plot_length_facetted(katSvMatched,"Distance<5e3&ClusterSize>=3&ClusterSize<1e3",'DistBucket2,ClusterSize','DistBucket2','ClusterSize','Kat-Distance ClusterSize',F))

# by cancer type
print(plot_length_facetted(katSvMatched,"Distance<1e3",'DistBucket2,CancerType','DistBucket2','CancerType','Kat-Distance vs CancerType',F))

# by kat-site length bucket
print(plot_length_facetted(katSvMatched,"Distance<1e3",'DistBucket2,KatLengthBucket','DistBucket2','KatLengthBucket','Kataegis Distance vs KatLengthBucket',F))

# by replication timing
print(plot_length_facetted(katSvMatched,"Distance<1e3",'DistBucket2,RepTimingBucket','DistBucket2',' RepTimingBucket','Kataegis Distance vs Rep-Timing',F))

# by homology
print(plot_length_facetted(katSvMatched,"Distance<1e3",'DistBucket2,HomologyLenBucket','DistBucket2',' HomologyLenBucket','Kataegis Distance vs Homology',F))

View(svData %>% filter(ClusterCount>=3&ResolvedType!='LINE') %>% group_by(LocTopTypeStart) %>% count())

# LINE only
View(katSvMatched)
print(plot_length_facetted(katSvMatched,"Distance<5e3&ResolvedType=='LINE'",'DistBucket2','DistBucket2','','Kat-Distance for LINE clusters',F))

print(ggplot(data = katSvMatched %>% filter(Distance<1e3&ResolvedType=='LINE'), aes(x=Distance,y=SnvCount))
      + geom_point())


print(ggplot(data = katSvMatched %>% filter(Distance<1e3&SvLengthType=='DEL_SHORT'), aes(x=Distance,y=SvLength))
      + geom_point())

View(katSvMatched %>% filter(ClusterCount==1&Distance<1e3&SvLengthType!='DEL_SHORT') %>% group_by(SampleId,Id) %>% count())


## Link to signatures to get samples with clear AID-APOBEC enrichment
cosmicSigs = as.matrix(read.csv("~/data/sigs/snv_cosmic_sigs.csv"))
cosmicSigsRnd = round(cosmicSigs,4)
View(cosmicSigsRnd)

snvMatrixData = as.matrix(read.csv('~/data/sigs/snv_prod_qc_pass_matrix_data_20190405.csv'))
nrow(snvMatrixData)
ncol(snvMatrixData)

cosmicSigFit = fit_to_signatures(snvMatrixData, cosmicSigs)

sampleSigContribs = as.data.frame(cosmicSigFit$contribution)
rownames(sampleSigContribs) = NULL
View(sampleSigContribs[,1:10])

rowIndex = data.frame(as.numeric(as.character(rownames(sampleSigContribs))))
colnames(rowIndex) <- c("Signature")
sampleSigContribs = cbind(rowIndex,sampleSigContribs)
View(sampleSigContribs[,1:10])

sampleSigContribs2 = gather(sampleSigContribs,"SampleId","Count",2:ncol(sampleSigContribs))
View(sampleSigContribs2)

sampleSnvTotals = sampleSigContribs2 %>% group_by(SampleId) %>% summarise(SampleTotal=round(sum(Count)))
View(sampleSnvTotals)


sampleSigContribs2 = merge(sampleSigContribs2,sampleSnvTotals,by='SampleId',all.x=T)
sampleSigContribs2 = sampleSigContribs2 %>% mutate(Count=round(Count),
                                                   SigPercent=round(Count/SampleTotal,2))

View(sampleSigContribs2 %>% filter((Signature==2|Signature==13)&SigPercent>0.2&SampleTotal>5000))


sig2And13Samples = sampleSigContribs2 %>% filter((Signature==2|Signature==13)&SigPercent>0.2&SampleTotal>5000)

svKatData = svKatData %>% mutate(AidApobec=SampleId %in% sig2And13Samples$SampleId)

svKatData = merge(svKatData,sampleSnvTotals,by='SampleId',all.x=T)

View(svKatData %>% filter(CancerType!='Skin',SnvCount>=3) %>% group_by(AidApobec,SnvTotalBucket=2**round(log(SampleTotal,2)),Matched=SvId>=0) 
     %>% count() %>% spread(Matched,n) %>% mutate(MatchPerc=round(`TRUE`/(`TRUE`+`FALSE`),2)))

# by whether has AID-APOBEC sigs
katSvMatched = katSvMatched %>% mutate(AidApobec=SampleId %in% sig2And13Samples$SampleId)
print(plot_length_facetted(katSvMatched,"ClusterCount==1&Type=='DEL'&Distance<1e3",'DistBucket2,AidApobec','DistBucket2','AidApobec','Kataegis Distance vs AID-APOBEC sigs',F))

View(svKatData %>% filter(is.na(SampleTotal)))

# matching by mutational load
View(svKatData %>% filter(SnvCount<=10&!is.na(SampleTotal)) %>% 
       group_by(MutLoad=ifelse(SampleTotal<10000,'Low',ifelse(SampleTotal<50000,'Medium','High')),SnvCount,Matched=SvId>=0) 
     %>% count() %>% spread(Matched,n) %>% mutate(MatchPerc=round(`TRUE`/(`TRUE`+`FALSE`),2)))


