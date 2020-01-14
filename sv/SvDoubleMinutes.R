# library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringi)
library(stringr)
library(devtools)




#####
## DOUBLE MINUTES

tmpDMs = read.csv('~/logs/LNX_DOUBLE_MINUTES.csv')
View(tmpDMs)
View(tmpDMs %>% filter(FullyChained=='true'&ClusterCount<=20))
doubleMinutes = tmpDMs

doubleMinutes = read.csv('~/data/sv/drivers/LNX_DOUBLE_MINUTES.csv')
View(doubleMinutes)
nrow(doubleMinutes)

#colnames(doubleMinutes)

# extract samples list to make subsequent LINX runs faster
dmSamples = doubleMinutes %>% group_by(SampleId) %>% count
nrow(dmSamples)
write.csv(dmSamples %>% select(SampleId),'~/logs/dm_samples.csv',row.names = F, quote = F)

sampleCancerTypes = load_cancer_types('~/data/hpc_sample_cancer_types.csv',F,10)
doubleMinutes = merge(doubleMinutes,sampleCancerTypes,by='SampleId',all.x=T)

doubleMinutes = doubleMinutes %>% mutate(HasSubclonalSVs=DMSvCount<ClusterCount,
                                         DmSVs=2**round(log(DMSvCount,2)),
                                         CompleteChain=(FullyChained=='true'),
                                         DMType=ifelse(ClusterDesc=='DUP'|DMSvTypes=='DUP=1','DUP',ifelse(CompleteChain,'COMPLEX_CHAINED','COMPLEX_PARTIAL')),
                                         HasAmpGenes=!is.na(AmpGenes)&AmpGenes!='',
                                         Ploidy=ifelse(MinPloidy>0,pmax(2**round(log(MinPloidy,2)),4),0),
                                         NormPloidy=ifelse(MinPloidy>0,pmax(2**round(log(MinPloidy/SamplePloidy,2)),4),0),
                                         DnaLength=ifelse(ChainLength>0,2**round(log(ChainLength,2)),0))

# summary stats

# types of DMs
nrow(doubleMinutes)
View(doubleMinutes %>% group_by(DMType,HasSubclonalSVs) %>% count)

# ploidy buckets
View(doubleMinutes %>% group_by(DMType,Ploidy=2**round(log(MinPloidy,2))) %>% count %>% spread(Ploidy,n,fill=0))
View(doubleMinutes %>% filter(HasDriver) %>% group_by(Type=ifelse(DMSvTypes=='DUP=1','DUP','COMPLEX'),
                                                      Ploidy=2**round(log(MinPloidy,2))) %>% count %>% spread(Ploidy,n,fill=0))

print(ggplot(doubleMinutes %>% filter(HasDriver&ClusterCount<=100&MinPloidy<=150) %>% 
               group_by(Type=ifelse(DMSvTypes=='DUP=1','DUP','COMPLEX'),Ploidy=5*round(MinPloidy/5)) %>% count,aes(x=Ploidy, y=n))
      + geom_bar(stat = "identity", colour = "black")
      + facet_wrap(~Type))

# occurrence
View(doubleMinutes %>% group_by(SampleId) %>% summarise(DMCount=n()))

# frequency per sample - possible evidence of underclustering
View(doubleMinutes %>% group_by(SampleId) %>% summarise(DMCount=n()) %>% group_by(DMCount) %>% summarise(Samples=n()))
View(doubleMinutes %>% group_by(SampleId) %>% summarise(DMCount=n()) %>% filter(DMCount>1) %>% arrange(-DMCount))

# frequency per cancer type
dmCancerTypes = doubleMinutes %>% group_by(SampleId,CancerType) %>% summarise(DmCount=n(),WithDriver=sum(HasDriver)) %>% 
  group_by(CancerType) %>% summarise(DmCount=n(),WithDriver=sum(WithDriver>0))

dmCancerTypes = merge(dmCancerTypes,sampleCancerTypes %>% group_by(CancerType) %>% summarise(SampleCount=n()),by='CancerType',all.x=T)
dmCancerTypes = dmCancerTypes %>% filter(!is.na(CancerType)) %>% mutate(Percent=round(DmCount/SampleCount,2),
                                                                        DriverPercent=round(WithDriver/DmCount,2))
# View(dmCancerTypes)

# DM SV count vs cluster size - most DMs have other low-ploidy material included
View(doubleMinutes %>% group_by(DmSVs,ClusterSize=2**round(log(ClusterCount,2))) %>% count %>% spread(ClusterSize,n,fill=0))
View(doubleMinutes %>% group_by(HighCN=(MaxCopyNumber>20),DmSVs,ClusterSize=2**round(log(ClusterCount,2))) %>% count %>% spread(ClusterSize,n,fill=0))

# occuring in the context of chromothripsis?? (St Judes)

# length of material
View(doubleMinutes %>% filter(DnaLength>0) %>% group_by(Type=ifelse(ClusterDesc=='DUP','DUP','COMPLEX'),DnaLength) %>% count %>% spread(Type,n,fill=0))

print(ggplot(doubleMinutes %>% filter(DnaLength>0) %>% group_by(DMType,DnaLength) %>% count,
             aes(x=DnaLength, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~DMType))

# as measured by length of copy number segments where CN is similar to the DM ploidy
doubleMinutes = doubleMinutes %>% mutate(MaterialLength=ifelse(FullyChained=='true',ChainLength,CNR_0))

View(doubleMinutes %>% mutate(DmRatio=round(MinPloidy/SamplePloidy,1),CnRatio=round(MaxCopyNumber/SamplePloidy,1)) %>% filter(MaterialLength==0) %>%
       select(SamplePloidy,MinPloidy,DmRatio,MaterialLength,CnSegmentLength,ChainLength,CnRatio,MaxCopyNumber,CNR_3,CNR_6,CNR_10,CNR_20,CNR_50,CNR_100,everything()))

View(doubleMinutes %>% mutate(MaterialLength=2**round(log(MaterialLength,2))) %>% group_by(Type=ifelse(ClusterDesc=='DUP','DUP','COMPLEX'),MaterialLength) %>% count %>% spread(Type,n,fill=0))

print(ggplot(doubleMinutes %>% mutate(MaterialLength=2**round(log(MaterialLength,2)),
                                      Type=paste(DMType,ifelse(HasDriver,'DRIVER','NO_DRIVER'),sep='_')) 
                                      %>% group_by(Type,MaterialLength) %>% count, aes(x=MaterialLength, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~Type))

print(ggplot(doubleMinutes %>% mutate(MaterialLength=2**round(log(MaterialLength,2))) %>% group_by(HasDriver,MaterialLength) %>% count,
             aes(x=MaterialLength, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~HasDriver))

# DM material mostly from a single chromosome - only ~10% having 2 chromosomes
View(doubleMinutes %>% mutate(ChrCount=stri_count_fixed(Chromosomes,';')+1) %>% group_by(ChrCount) %>% count)
View(doubleMinutes %>% mutate(ChrCount=stri_count_fixed(Chromosomes,';')+1) %>% group_by(ChrCount,FullyChained) %>% count %>% spread(FullyChained,n,fill=0))

# Adjacent MAP
View(doubleMinutes %>% group_by(MinAdjMAPRatio=2**round(log(MinAdjMAPRatio,2))) %>% count)

## Affected genes
View(doubleMinutes %>% group_by(HasAmpGenes,Ploidy=pmax(2**round(log(MinPloidy,2)),4)) %>% count %>% spread(Ploidy,n,fill=0))
View(doubleMinutes %>% group_by(HasAmpGenes,NormPloidy) %>% count %>% spread(NormPloidy,n,fill=0))
View(doubleMinutes %>% group_by(HasDriver,NormPloidy) %>% count %>% spread(NormPloidy,n,fill=0))
View(doubleMinutes %>% group_by(CancerType,HasDriver) %>% count %>% spread(HasDriver,n))
View(doubleMinutes %>% group_by(Chromosomes,HasDriver) %>% count %>% spread(HasDriver,n))
View(doubleMinutes %>% group_by(CancerType,Chromosomes,HasDriver) %>% count %>% spread(HasDriver,n))

# Link to oncogene drivers
drivers = read.csv('~/data/sv/drivers/LNX_DRIVERS.csv')
#View(drivers)
ampDrivers = drivers %>% filter(Category=='ONCO')
#View(ampDrivers)

doubleMinutes = within(doubleMinutes,rm(HasDriver))
doubleMinutes = within(doubleMinutes,rm(DriverAMPs))
colnames(doubleMinutes)
doubleMinutes$DriverAMPs=""
colIndex=ncol(doubleMinutes)

dmDriverGenes = data.frame(matrix(ncol = 4, nrow = 0))
colnames(dmDriverGenes) = c('SampleId','ClusterId','Gene','GeneMinCN')

for(i in 1:nrow(doubleMinutes))
{
  dmData = doubleMinutes[i,]
  sampleId = as.character(dmData$SampleId)
  clusterId = as.numeric(dmData$ClusterId)
  # dmGenes = ifelse(!is.na(dmData$AmpGenes),as.character(dmData$AmpGenes),'')
  driverData = ampDrivers %>% filter(SampleId==sampleId&ClusterId==clusterId) %>% group_by(Gene,GeneMinCN) %>% count %>% ungroup()
  ampStr = ""
  
  if(nrow(driverData) > 0) # dmGenes !=''
  {
    for(j in 1:nrow(driverData))
    {
      dd = driverData[j,]
      gene = as.character(dd$Gene)
      geneCN = dd$GeneMinCN
      
      # print(sprintf('sample(%s): searching for %s in %s',sampleId,gene,dmGenes))
      
      # if(!is.na(gene)) # grepl(gene,dmGenes)
      if(!is.na(gene))
      {
        ampStr = ifelse(ampStr=='',gene,paste(ampStr,gene,sep=';'))
        
        rowIndex = nrow(dmDriverGenes)+1
        dmDriverGenes[rowIndex,1] = sampleId
        dmDriverGenes[rowIndex,2] = clusterId
        dmDriverGenes[rowIndex,3] = gene
        dmDriverGenes[rowIndex,4] = geneCN
      }
    }
  
    doubleMinutes[i,colIndex] = ampStr
  }
}

doubleMinutes = doubleMinutes %>% mutate(HasDriver=DriverAMPs!='')
View(doubleMinutes)
View(dmDriverGenes)

View(dmDriverGenes %>% group_by(Gene) %>% count)


write.csv(doubleMinutes,'~/data/sv/DM_summary.csv',row.names = F, quote = F)


## MULTIPLE BIOPSY data

mbSvData = read.csv('~/data/sv/multi_biop/LNX_MB_SV_DATA.csv')
# mbMergeDataAll = read.csv('~/data/sv/multi_biop/LNX_MB_MERGE_DATA.csv')
mbClusterData = read.csv('~/data/sv/multi_biop/LNX_MB_CLUSTER_DATA.csv')
mbSVs = read.csv('~/data/sv/multi_biop/LNX_SVS.csv')
View(mbMergeDataAll)
mbSvCombined = merge(mbSVs,mbSvData %>% select(PatientId,SampleId,Id=SvId,MatchType,OtherSvId),by=c('SampleId','Id'),all.x=T)

View(mbSvData)

# not yet allowing multiple matches per SV
View(mbSvData %>% filter(MatchType=='Shared') %>% group_by(SampleId,SvId) %>% count %>% filter(n>1))
# View(mbSvData %>% group_by(MatchType) %>% count)
View(mbSvCombined %>% group_by(MatchType,DMSV) %>% count %>% spread(DMSV,n,fill=0))
View(mbSvCombined %>% filter(MatchType!='Partial') %>% group_by(MatchType,DMSV,Ploidy=ifelse(Ploidy>1,2**round(log(Ploidy,2)),1)) %>% 
       count %>% spread(Ploidy,n,fill=0))

nrow(mbSVs) # 27.6K
nrow(mbSvData) # 158K

# get all clusters marked as DM also with a driver AMP in either sample and 

mbDMs = doubleMinutes %>% filter(SampleId %in% doubleBiopsy$SampleId)
mbDMs = merge(mbDMs,doubleBiopsy,by='SampleId',all.x=T)
View(mbDMs)

mbDMs = mbDMs %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))
mbClusterData = mbClusterData %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))
mbSvCombined = mbSvCombined %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'),
                                       Ploidy=(PloidyMin+PloidyMax)*0.5)

write.csv(mbDMs,'~/logs/mb_dm_data.csv',row.names = F,quote = F)
write.csv(mbDMs %>% group_by(SampleId) %>% count %>% select(SampleId),'~/logs/mb_dm_sample_ids.csv',row.names = F,quote = F)

# 61 samples from 38 distinct patients and 58 distinct samples have a double minute, 3 of them have 2 distinct DMs
View(mbDMs %>% group_by(SampleId) %>% count)
View(mbDMs %>% group_by(PatientId) %>% count)

# 25 samples from 18 patients have a driver (and none have 2 per sample)
View(mbDMs %>% filter(HasDriver) %>% group_by(SampleId) %>% count)
View(mbDMs %>% filter(HasDriver) %>% group_by(PatientId) %>% count)

# 20 have drivers in both samples, 18 have a driver in only 1 sample (classified as a DM)
View(mbDMs %>% group_by(PatientId,SampleId) %>% summarise(DriverCount=n()) %>% group_by(PatientId) %>% summarise(DriverCount=sum(DriverCount>0)))

# from 18 patients, 11 have a driver in only one of the samples, 7 as DMs in both
View(mbDMs %>% filter(HasDriver) %>% group_by(PatientId,SampleId,DriverAMPs) %>% summarise(DriverCount=n()) %>% 
       group_by(PatientId,SampleId) %>% summarise(DriverCount=sum(DriverCount)) %>% group_by(PatientId) %>% count)

# with MB data
mbDMsMatchData = merge(mbDMs,mbClusterData %>% select(-ClusterCount,-ResolvedType),by=c('SampleId','ClusterId'),all.x=T)
View(mbDMsMatchData)

# 6 DMs have no matching clusters

# DM SVs matching
mbDMSVs = mbSvCombined %>% filter(DMSV=='true')
View(mbDMSVs)
View(mbDMSVs %>% filter(ClusterDesc=='DUP'))
View(mbDMSVs %>% group_by(MatchType) %>% count)
View(mbDMSVs %>% group_by(MatchType,ClusterCount=2**round(log(ClusterCount,2))) %>% count %>% spread(MatchType,n,fill=0))
View(mbDMSVs %>% group_by(MatchType,Ploidy=2**round(log(Ploidy,2))) %>% count %>% spread(MatchType,n,fill=0))

View(mbSvCombined %>% filter(SampleClusterId %in% mbDMs$SampleClusterId) %>% group_by(PatientId,SampleId,ClusterId,ClusterCount,ClusterDesc) %>%
       summarise(DMSvCount=sum(DMSV=='true'),
                 DMSvMatched=sum(DMSV=='true'&MatchType=='Shared')))

View(mbDMSVs %>% filter(MatchType=='Private') %>% 
       group_by(PatientId,MatchType,Type,ChrStart,ChrEnd,OrientStart,OrientStart,PosStartRnd=round(PosStart,-2),PosEndRnd=round(PosEnd,-2)) %>%
       summarise(Count=n(),Sample1=first(SampleId),Sample2=last(SampleId),Id1=first(Id),Id2=last(Id),
                 PosStart1=first(PosStart),PosStart2=last(PosStart),PosEnd1=first(PosEnd),PosEnd2=last(PosEnd),
                 PosDiffStart=abs(first(PosStart)-last(PosStart)),PosDiffEnd=abs(first(PosEnd)-last(PosEnd))) %>% filter(Count==2) %>%
       mutate(PosDiffTotal=PosDiffStart+PosDiffEnd))

View(mbDMSVs %>% filter(PatientId=='CPCT02070345'&Type=='BND'))

# identical DMs between the 2 samples - 8 patients / 16 samples/clusters
identicalDMs = mbDMs %>% group_by(PatientId,DMSvTypes) %>% 
  summarise(SampleCount=n(),Sample1=first(SampleId),Sample2=last(SampleId),
            Cluster1=first(ClusterId),Cluster2=last(ClusterId),
            ResolvedType1=first(ResolvedType),ResolvedType2=last(ResolvedType),
            ClusterCount1=first(ClusterCount),ClusterCount2=last(ClusterCount),
            DriverAMPs1=first(DriverAMPs),DriverAMPs2=last(DriverAMPs),
            MinPloidy1=first(MinPloidy),MinPloidy2=last(MinPloidy),
            MaxCopyNumber1=first(MaxCopyNumber),MaxCopyNumber2=last(MaxCopyNumber)) %>% filter(SampleCount==2)

View(identicalDMs)

identicalDMs = merge(identicalDMs,
                   mbClusterData %>% select(Sample1=SampleId,Cluster1=ClusterId,
                                            PrCount1=PrivateCount,ShCount1=SharedCount,MatchCC1=MatchingClustersCount,OtherCIds1=OtherClusterIds,OT1=OverlapType),
                   by=c('Sample1','Cluster1'),all.x=T)

identicalDMs = merge(identicalDMs,
                   mbClusterData %>% select(Sample2=SampleId,Cluster2=ClusterId,
                                            PrCount2=PrivateCount,ShCount2=SharedCount,MatchCC2=MatchingClustersCount,OtherCIds2=OtherClusterIds,OT2=OverlapType),
                   by=c('Sample2','Cluster2'),all.x=T)

View(identicalDMs)

View(mbDMs %>% filter(DMSvTypes=='DUP=1'&DMSvCount==1) %>% select(PatientId,SampleId,everything()))

# DRUP01330016T	245	CPCT02330074T	230 - exact matches but 1 DM SV has quite diff ploidy
View(mbSvCombined %>% filter(SampleClusterId=='DRUP01330016T_245'|SampleClusterId=='CPCT02330074T_230') %>% 
       select(DMSV,SampleId,Id,Type,MatchType,OtherSvId,ChrStart,ChrEnd,PosStart,PosEnd,everything()))

# DRUP01340001T	148	CPCT02340024T	175
View(mbSvCombined %>% filter(SampleClusterId=='DRUP01340001T_148'|SampleClusterId=='CPCT02340024T_175') %>% 
       select(DMSV,SampleId,Id,Type,MatchType,OtherSvId,ChrStart,ChrEnd,PosStart,PosEnd,everything()))

# similar DMs
similarDMs = mbDMs %>% group_by(PatientId) %>% summarise(SampleCount=n(),Sample1=first(SampleId),Sample2=last(SampleId),
                                                 Cluster1=first(ClusterId),Cluster2=last(ClusterId),
                                                 ResolvedType1=first(ResolvedType),ResolvedType2=last(ResolvedType),
                                                 ClusterCount1=first(ClusterCount),ClusterCount2=last(ClusterCount),
                                                 DriverAMPs1=first(DriverAMPs),DriverAMPs2=last(DriverAMPs),
                                                 MinPloidy1=first(MinPloidy),MinPloidy2=last(MinPloidy),
                                                 DMSvTypes1=first(DMSvTypes),DMSvTypes2=last(DMSvTypes)) %>% filter(SampleCount==2&DMSvTypes1!=DMSvTypes2)

View(similarDMs %>% select(PatientId,Sample1,Sample2,Cluster1,Cluster2,ClusterCount1,ClusterCount2,DMSvTypes1,DMSvTypes2,everything()))

# pull in MB match data
View(mbClusterData %>% filter(SampleClusterId %in% mbDMs$SampleClusterId))

similarDMs = merge(similarDMs,
                   mbClusterData %>% select(Sample1=SampleId,Cluster1=ClusterId,
                                            PrCount1=PrivateCount,ShCount1=SharedCount,MatchCC1=MatchingClustersCount,OtherCIds1=OtherClusterIds,OT1=OverlapType),
                   by=c('Sample1','Cluster1'),all.x=T)

similarDMs = merge(similarDMs,
                   mbClusterData %>% select(Sample2=SampleId,Cluster2=ClusterId,
                                            PrCount2=PrivateCount,ShCount2=SharedCount,MatchCC2=MatchingClustersCount,OtherCIds2=OtherClusterIds,OT2=OverlapType),
                   by=c('Sample2','Cluster2'),all.x=T)

View(similarDMs)

View(similarDMs %>% mutate(Linked=ifelse(OtherCIds2!='',stringr::str_detect(as.character(OtherCIds2),as.character(Cluster1)),F)
                           |ifelse(OtherCIds1!='',stringr::str_detect(as.character(OtherCIds1),as.character(Cluster2)),F)) %>%
       select(PatientId,Sample1,Sample2,Cluster1,Cluster2,Linked,ClusterCount1,ClusterCount2,DMSvTypes1,DMSvTypes2,OtherCIds1,OtherCIds2,PrCount1,PrCount2,ShCount1,ShCount2,everything()))

# investigations

# CPCT02240009T	DRUP01010072T	128	157
View(mbSvCombined %>% filter(SampleId=='CPCT02240009T'&ClusterId==128))

# exact matches
View(mbSvCombined %>% filter(SampleId %in% c('CPCT02010288T','CPCT02010288TII')&(ClusterId==47|ClusterId==25)) %>% 
       select(DMSV,SampleId,Id,Type,MatchType,OtherSvId,ChrStart,ChrEnd,PosStart,PosEnd,everything()))


View(mbSvCombined %>% filter(SampleClusterId %in% mbDMs$SampleClusterId) %>%
       #filter(SampleId=='CPCT02240009T'&ClusterId==128) %>%
       filter(PatientId=='CPCT02240009') %>% 
       filter(DMSV=='true'|Ploidy>8) %>%
       select(PatientId,SampleId,ClusterId,Id,Type,DMSV,Ploidy,MatchType,OtherSvId,OverlapType))


       

# plotData = rbind(similarDMs %>% select(SampleId=Sample1,ClusterId=Cluster1),similarDMs %>% select(SampleId=Sample2,ClusterId=Cluster2))
# write.csv(plotData,'~/logs/mb_dm_similar_sample_ids.csv',row.names = F,quote = F)


# unmatched DMs with drivers
loneDMDrivers = mbDMs %>% filter(HasDriver) %>% group_by(PatientId) %>% summarise(SampleCount=n(),
                                                                                  SampleId=first(SampleId),
                                                                                  ClusterId=first(ClusterId),
                                                                                  ClusterCount=first(ClusterCount),
                                                                                  ClusterDesc=first(ClusterDesc),
                                                                                  Gene=first(DriverAMPs),
                                                                                  DMSvTypes=first(DMSvTypes),
                                                                                  MinPloidy=first(MinPloidy),
                                                                                  MaxCopyNumber=first(MaxCopyNumber)) %>% 
  filter(SampleCount==1) %>% select(-SampleCount)

View(loneDMDrivers)

loneDMDrivers = merge(loneDMDrivers,doubleBiopsy %>% filter(!(SampleId %in% loneDMDrivers$SampleId)) %>% select(PatientId,OtherSampleId=SampleId),
                      by='PatientId',all.x=T)

# look for matched clusters in the alternate sample
loneDMDrivers = merge(loneDMDrivers,mbClusterData %>% select(-ClusterCount,-ResolvedType),by=c('SampleId','ClusterId'),all.x=T)
View(loneDMDrivers)

# look for a driver gene in the other sample
loneDMDrivers = merge(loneDMDrivers,
                      ampDrivers %>% select(OtherSampleId=SampleId,Gene,OtherClusterId=ClusterId,GeneMinCN,Chromosome,
                                            OtherClusterCount=ClusterCount,OtherResolvedType=ResolvedType),by=c('OtherSampleId','Gene'),all.x=T)

View(loneDMDrivers)
write.csv(loneDMDrivers,'~/logs/mb_dm_driver_unmatched.csv',quote = F, row.names = F)

plotData = loneDMDrivers %>% select(SampleId,ClusterId)
plotData = rbind(plotData,loneDMDrivers %>% filter(!is.na(OtherClusterId)) %>% select(SampleId=OtherSampleId,ClusterId=OtherClusterId))
View(plotData)
write.csv(plotData,'~/logs/mb_dm_driver_sample_ids.csv',row.names = F,quote = F)


View(mbDMsByPatient)
missingSamples = mbDMsByPatient %>% filter(SampleCount==1)
missingSamples = merge(missingSamples,
                       doubleBiopsy %>% filter(!(SampleId %in% missingSamples$Sample1)) %>% select(PatientId,OtherSampleId=SampleId),
                                               by='PatientId',all.x=T)

missingDrivers = missingSamples %>% filter(DriverAMPs1!='') %>% mutate(Gene=DriverAMPs1)

missingDrivers = merge(missingDrivers,ampDrivers,by.x=c('OtherSampleId','Gene'),by.y=c('SampleId','Gene'),all.x=T)
View(missingDrivers)

View(missingSamples)
View(doubleBiopsy)

ampDrivers = merge(ampDrivers,doubleBiopsy,by='SampleId',all.x=T)
View(ampDrivers %>% filter(PatientId=='CPCT02040071'))
View(ampDrivers)



## DM FILTERING

## Suspect DMs
# low copy-number points within DM chain
# high incidence of foldbacks and foldback ploidy
# high number of SGLs suggesting further unknown SVs and possible foldbacks
# not fully chained
# low-end copy-number and ploidy
# no gene or driver gene
# incidence of non-DM SVs with high ploidy


# min CN of non-DM SV segments
View(doubleMinutes %>% filter(ChainCount>0) %>% group_by(MinCnBucket=round(ChainMinCnPercent,1)) %>% count)
View(doubleMinutes %>% filter(ChainCount>0) %>% arrange(ChainMinCnPercent))
View(doubleMinutes %>% filter(ClusterCount>1&DMSvTypes=='DUP=1'&ChainMinCnPercent<0.5&ClusterCount<100))


# number and ploidy of foldbacks
View(doubleMinutes %>% group_by(FbSumBucket=round(FbSumPloidy/MaxPloidy,1)) %>% count)
View(doubleMinutes %>% group_by(FbCount=ifelse(FbCount>1,2**round(log(FbCount,2)),FbCount),FbSumBucket=round(FbSumPloidy/MaxPloidy,1)) 
     %>% count %>% spread(FbCount,n,fill=0))

# number and ploidy of SGLs
View(doubleMinutes %>% group_by(SglSumBucket=round(SglSumPloidy/MaxPloidy,1)) %>% count)
View(doubleMinutes %>% group_by(SglCount=ifelse(SglCount>1,2**round(log(SglCount,2)),SglCount),SglSumBucket=round(SglSumPloidy/MaxPloidy,1)) 
     %>% count %>% spread(SglCount,n,fill=0))

# non-DM SVs whose ploidy is > the DM min ploidy
View(doubleMinutes %>% group_by(NonDMSvPerc=ifelse(ClusterCount-DMSvCount>0,round(NonDmSvsGtPloidy/(ClusterCount-DMSvCount),1),0)) %>% count)

View(doubleMinutes %>% group_by(NonDMSvFullPerc=ifelse(ClusterCount-DMSvCount>0,round(NonDmSvsGtPloidy/(ClusterCount-DMSvCount),1),0),
                                NonDMSvHalfPerc=ifelse(ClusterCount-DMSvCount>0,round(NonDmSvsGtHalfPloidy/(ClusterCount-DMSvCount),1),0)) 
     %>% count %>% spread(NonDMSvFullPerc,n,fill=0))


View(doubleMinutes %>% filter(DMSvTypes=='INV=2'&FbCount>=2))
View(doubleMinutes %>% filter(DMSvTypes=='INF=2'))
View(doubleMinutes %>% filter(DMSvTypes=='SGL=2'))

# possibly missing a variant
View(doubleMinutes %>% filter(DMSvCount<=5&NonDmSvsGtPloidy>0))

View(dmSvs %>% filter(SampleId=='CPCT02030269T'&ClusterId==94))

doubleMinutes = doubleMinutes %>% mutate(LowCNScore=round(ifelse(MinChainCnPercent>0.5,10,round(MinChainCnPercent*10))),
                                         FbScore=round(((1-pmin(FbSumPloidy/MaxPloidy,1))^FbCount)*10),
                                         SglScore=round(((1-pmin(SglSumPloidy/MaxPloidy,1))^SglCount)*10),
                                         ChainedScore=ifelse(FullyChained=='true',10,pmin(ClusterCount,10)),
                                         NonDMSvScore=ifelse(ClusterCount-DMSvCount==0,10,round((1-(NonDmSvsGtPloidy/(ClusterCount-DMSvCount)))*10)),
                                         DMScore=LowCNScore+FbScore+SglScore+ChainedScore+NonDMSvScore)

View(doubleMinutes %>% select(SampleId,ClusterId,ClusterCount,ClusterDesc,DMSvTypes,DMScore,LowCNScore,FbScore,SglScore,ChainedScore,NonDMSvScore,
                              everything()))

View(doubleMinutes %>% group_by(DMScore=5*round(DMScore/5),ClusterCount=2**round(log(ClusterCount,2))) %>% count %>% spread(DMScore,n,fill=0))

View(doubleMinutes %>% mutate(FbRatio=pmin(FbSumPloidy/MaxPloidy,1)) %>% select(FbCount,FbSumPloidy,MaxPloidy,FbRatio,FoldbackScore,SampleId,ClusterId))
View(doubleMinutes %>% mutate(SglRatio=pmin(SglSumPloidy/MaxPloidy,1)) %>% select(SglCount,SglSumPloidy,MaxPloidy,SglRatio,SglScore,
                                                                                  SampleId,ClusterId,DMSvTypes,ClusterCount,ClusterDesc))
View(doubleMinutes %>% select(MinChainCnPercent,LowCNScore,MaxCopyNumber,MaxPloidy,ClusterCount,DMSvTypes,SampleId,ClusterId))
View(doubleMinutes %>% select(ClusterCount,DMSvCount,NonDmSvsGtPloidy,NonDMSvScore,SampleId,ClusterId,DMSvTypes,ClusterCount,ClusterDesc))



# ploidy buckets and ratios
View(doubleMinutes %>% mutate(PloidyToCN=round(MaxPloidy/MaxCopyNumber,1),
                              CNBucket=2**round(log(MaxCopyNumber,2))) %>% 
       group_by(PloidyToCN,CNBucket) %>% count %>% spread(PloidyToCN,n,fill=0))
       

View(doubleMinutes %>% mutate(PloidyToCN=round(MaxPloidy/MaxCopyNumber,1)) %>% 
       group_by(PloidyToCN,DMType) %>% count %>% spread(PloidyToCN,n,fill=0))


View(doubleMinutes)
View(doubleMinutes %>% filter(HasDriver))
View(doubleMinutes %>% filter(HasDriver&grepl(';',DriverAMPs)))

# plot clusters with between 2 and 100 SVs
write.csv(doubleMinutes %>% filter(ClusterCount>1&ClusterCount<100) %>% select(SampleId,ClusterId),'~/logs/dm_clusters.csv',row.names = F, quote = F)


View(doubleMinutes)
View(doubleMinutes %>% filter(!HasAmpGenes|NearDMSvCount>0|FullyChained=='false'|MaxPloidy<0.6*MaxCopyNumber) %>%
       mutate(PloidyToCN=round(MaxPloidy/MaxCopyNumber,2)) %>%
       select(SampleId,DMSvCount,DMSvTypes,HasAmpGenes,NearDMSvCount,FullyChained,PloidyToCN,MaxPloidy,MaxCopyNumber,ClusterCount,ClusterDesc,everything()))

View(doubleMinutes %>% filter(HasAmpGenes&NearDMSvCount==0&FullyChained=='true') %>%
       mutate(PloidyToCN=round(MaxPloidy/MaxCopyNumber,2)) %>%
       select(SampleId,DMSvCount,DMSvTypes,DriverAMPs,AmpGenes,PloidyToCN,MaxPloidy,MaxCopyNumber,ClusterCount,ClusterDesc,everything()))


#####
## SV properties of DMs

dmSvs = read.csv('~/data/sv/DM_LNX_SVS.csv')
doubleMinutes = doubleMinutes %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))
dmSvs = dmSvs %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))
dmSvs = dmSvs %>% filter(SampleClusterId %in% doubleMinutes$SampleClusterId)
nrow(dmSvs)
View(dmSvs)
sum(doubleMinutes$ClusterCount)
nrow(dmSvs %>% filter(DMSV=='true'))
sum(doubleMinutes$DMSvCount) # 4433

dmSvs = dmSvs %>% mutate(PloidyMid=(PloidyMax+PloidyMin)*0.5)

# evidence of shattering prior to formation of DM
dmSvSummary = dmSvs %>% group_by(SampleId,ClusterId,ClusterCount,ClusterDesc) %>%
  summarise(DMSvCount=sum(DMSV=='true'),
            DMMinPloidy=max(ifelse(DMSV=='true',PloidyMin,0)),
            CDBDmSvCount=sum(ifelse(DMSV=='true',ifelse(LocTopTypeStart=='DSB',0.5,0)+ifelse(LocTopTypeEnd=='DSB',0.5,0),0)),
            CDBNonDmSvCount=sum(ifelse(DMSV=='false'&PloidyMax<0.5*DMMinPloidy,ifelse(LocTopTypeStart=='DSB',0.5,0)+ifelse(LocTopTypeEnd=='DSB',0.5,0),0)),
            DBDmSvCount=sum(ifelse(DMSV=='true',ifelse(LocTopTypeStart=='DSB',0.5,0)+ifelse(LocTopTypeEnd=='DSB',0.5,0),0)),
            DBNonDmSvCount=sum(ifelse(DMSV=='false'&PloidyMax<0.5*DMMinPloidy,ifelse(LocTopTypeStart=='DSB',0.5,0)+ifelse(LocTopTypeEnd=='DSB',0.5,0),0)))

View(dmSvSummary)            
View(dmSvSummary %>% group_by(HasDMSvDBs=CDBDmSvCount>0,HasNonDMSvDBs=CDBNonDmSvCount>0) %>% count)
dmSvSummary = dmSvSummary %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))

# DMs formed from multiple breaks and no foldbacks
cleanMultiSvDMs = doubleMinutes %>% filter(DMSvCount>1&ClusterCount<100&FbCount==0&SglCount==0&FullyChained=='true'&ChainMinCnPercent>0.6)
View(cleanMultiSvDMs)
View(dmSvSummary %>% filter(SampleClusterId %in% cleanMultiSvDMs$SampleClusterId & CDBDmSvCount>=1))

cleanFbSvDMs = doubleMinutes %>% filter(DMSvCount>1&ClusterCount<100&FbCount>0&SglCount==0&FullyChained=='true'&ChainMinCnPercent>0.6)
View(cleanFbSvDMs)

View(dmSvs %>% filter(SampleId=='CPCT02010440T'&ClusterId==136) %>% select(SampleId,ClusterId,Id,Type,DMSV,PloidyMid,PloidyMin,DBLenStart,DBLenEnd,
                                                                          LocTopTypeStart,LocTopTypeEnd,LocTopIdStart,LocTopIdEnd))

View(dmSvs %>% filter(SampleId=='CPCT02120029T'&ClusterId==81) %>% select(SampleId,ClusterId,Id,Type,DMSV,PloidyMid,PloidyMin,DBLenStart,DBLenEnd,
                                                                          LocTopTypeStart,LocTopTypeEnd,LocTopIdStart,LocTopIdEnd))

nrow(dmDriverGenes)
View(doubleMinutes)

# shattering prior:
# CPCT02120029T'&ClusterId==81

dmSvs = merge(dmSvs,doubleMinutes %>% select(SampleClusterId,DMSvCount,DMType,DMSvTypes,MaxPloidy,MaxCopyNumber,DriverAMPs),by='SampleClusterId',all.x=T)
dmSvs = merge(dmSvs,doubleMinutes %>% select(SampleClusterId,DriverAMPs),by='SampleClusterId',all.x=T)
View(head(dmSvs,100))
View(dmSvs %>% filter(SampleId=='CPCT02050186T'&ChrStart==8))
View(dmSvs %>% filter(is.na(MaxCopyNumber)))
dmSvs = dmSvs %>% mutate(IsDMSV=DMSV=='true',
                         PloidyMid=(PloidyMax+PloidyMin)*0.5,
                         PloidyRatio=pmin(round(PloidyMid/MaxPloidy,2),1),
                         PloidyRatioBucket=paste('PRB',0.25*round(PloidyRatio/0.25),ifelse(IsDMSV,'DMSv','Sv'),sep='_'))

dmClusters = dmSvs %>% group_by(SampleId,ClusterId,DriverAMPs,DMSvCount,ClusterCount,PloidyRatioBucket) %>% count %>% spread(PloidyRatioBucket,n,fill=0)
View(dmClusters)
View(dmClusters %>% filter(!IsDMSV&(PRB_0.75>0|PRB_1>0)))

View(doubleMinutes)

View(dmSvs %>% group_by(SampleId,ClusterId,DMSvCount,ClusterCount,PloidyRatioBucket) %>% count %>% spread(PloidyRatioBucket,n,fill=0) %>%
       mutate(HighPBR=PRB_0.75+PRB_1) %>% filter(HighPBR>DMSvCount) %>% mutate(ExcessPerc=round((HighPBR-DMSvCount)/ClusterCount,2)))


# length by copy number / sample ploidy ratio



## DEBUG

sampleDMs = read.csv('~/logs/LNX_DOUBLE_MINUTES.csv')
View(sampleDMs)

View(doubleMinutes %>% filter(DriverAMPs!='') %>% group_by(DriverAMPs) %>% count)
View(doubleMinutes %>% filter(AmpGenes=='AR'|grepl('AR;',AmpGenes)|grepl(';AR',AmpGenes)) %>% filter(CancerType=='Prostate'))
View(doubleMinutes %>% filter(grepl('EGFR',AmpGenes)&CancerType=='CNS'))

prevDoubleMinutes = doubleMinutes
nrow(prevDoubleMinutes)
newDoubleMinutes = read.csv('~/logs/LNX_DOUBLE_MINUTES.csv')
nrow(newDoubleMinutes)
View(newDoubleMinutes)

View(newDoubleMinutes2 %>% filter(MaxTeloCentroCn>2.3*MaxPloidy))

colnames(newDoubleMinutes)





View(doubleMinutes %>% group_by(DMSvs=pmin(DMSvCount,10)) %>% count)
View(doubleMinutes %>% group_by(ResolvedType=ifelse(SingleDUP,'DUP',as.character(ResolvedType)),DMSvs=pmin(DMSvCount,10)) %>% count)


View(doubleMinutes %>% group_by(ResolvedType=ifelse(ClusterDesc=='DUP','DUP',as.character(ResolvedType)),DMSvs=pmin(DMSvCount,10)) %>% count)




