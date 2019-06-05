library(tidyr)
library(dplyr)
library(RMySQL)
library(purple)
#library(GenomicRanges)
#library(BSgenome.Hsapiens.UCSC.hg19)


## DATA PREP
doubleBiopsy = multipleBiopsyCohort %>% group_by(patientId) %>% filter(n() == 2) %>% mutate(n = row_number())
nrow(doubleBiopsy)
View(doubleBiopsy)
dbSvs = mbcSvs %>% filter(SampleId %in% doubleBiopsy$sampleId)
nrow(dbSvs)

View(doubleBiopsy %>% group_by(patientId) %>% count()) #  %>% group_by(n) %>% count()
write.csv(doubleBiopsy %>% select(patientId,sampleId), '~/data/sv/multiple_biopsy_samples.csv', row.names = F, quote = F)
doubleBiopsy = read.csv('~/data/sv/multiple_biopsy_samples.csv')

svData = read.csv('~/logs/SVA_SVS.csv')
mbSvDataActuals = svData %>% filter(SampleId %in% doubleBiopsy$sampleId)
write.csv(mbSvDataActuals, '~/data/sv/SVA_SVS_MB.csv', row.names = F, quote = F)
nrow(mbSvDataActuals)
rm(svData)
nrow(mbSvData)
View(mbSvData)
View(mbSvData %>% group_by(SampleId) %>% count())

# write out a subset of SV data for the MB samples for use in generating the MB data
write.csv(mbSvDataActuals %>% select(SampleId,Id,Type,ChrStart,PosStart,OrientStart,ChrEnd,PosEnd,OrientEnd,
                                     ClusterId,ClusterCount,ClusterReason,ResolvedType,Ploidy,CNStart,CNEnd) %>% arrange(SampleId,ClusterId),
          '~/data/sv/multiple_biopsy_sv_data.csv', row.names = F, quote = F)

## ANALYSIS
mbSvData = read.csv('~/data/sv/SVA_MB_SV_DATA.csv')
mbMergeDataAll = read.csv('~/data/sv/SVA_MB_MERGE_DATA.csv')
mbClusterData = read.csv('~/data/sv/SVA_MB_CLUSTER_DATA.csv')
View(mbSvData)
nrow(mbSvData)
nrow(mbMergeDataAll)
View(mbMergeData)
View(mbClusterData)
View(mbMergeData %>% group_by(ClusterReason) %>% count())
View(mbClusterData %>% group_by(OverlapType) %>% count())
View(mbClusterData %>% group_by(OverlapType,ResolvedType) %>% count() %>% spread(OverlapType,n))

mbSvData = merge(mbSvData,mbClusterData %>% select(SampleId,ClusterId,OverlapType),by=c('SampleId','ClusterId'),all.x=T)

defaultPlotColours = c("#fb8072","#bc80bd","#bebada", "#fdb462","#80b1d3","#8dd3c7","#b3de69","#fccde5","#ffffb3","#d9d9d9",
                       'saddlebrown','slateblue1','indianred3','limegreen')


# 0. Overall results by Match Type
# consider private and exactly shared as fine
mbClusteringSummary = mbClusterData %>% 
  mutate(ResolvedType=ifelse(ResolvedType=='None'|ResolvedType=='ComplexChain','Complex',as.character(ResolvedType))) %>%
  group_by(ResolvedType) %>% 
  summarise(ResolvedTypeCount=n(),
            SvCount=sum(ClusterCount),
            CorrectClusters=sum(ifelse(OverlapType=='Exact'|OverlapType=='Private',1,0)),
            IncorrectClusters=sum(ifelse(OverlapType!='Exact'&OverlapType!='Private',1,0))) %>%
  mutate(CorrectPerc=round(CorrectClusters/ResolvedTypeCount,3)*100,
         IncorrectPerc=round(IncorrectClusters/ResolvedTypeCount,3)*100)

View(mbClusteringSummary)         

csPlot = (ggplot(mbClusteringSummary, 
                   aes(x=ResolvedType, y=IncorrectPerc, fill=ResolvedType))
            + geom_bar(stat = "identity", colour = "black")
            + labs(x = "Resolved Type", y="Incorrect %")
            + scale_fill_manual(values = defaultPlotColours)
            + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
            + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
            + theme(axis.text.x = element_text(angle=90, hjust=1,size=10))
            + theme(legend.position="none")
            + coord_flip())

plot(csPlot)

# shared-private merges by resolved type, split by proportion of cluster with over-clustering
View(mbClusterData)

mbClusterSummary = mbClusterData %>% 
  mutate(MixedPerc=1-abs(SharedCount-ClusterCount/2)/(ClusterCount/2),
         MixedPercBucket=round(MixedPerc/0.1)*0.1,
         ClusterSize=ifelse(ClusterCount<=3,'Small',ifelse(ClusterCount<=10,'Medium','Large')))

View(mbClusterSummary %>% filter(ResolvedType=='SimpleChain') %>% group_by(MixedPercBucket,ClusterSize) %>% count() %>% spread(ClusterSize,n))
View(mbClusterSummary %>% filter(ResolvedType=='ComplexChain'|ResolvedType=='None') %>% group_by(MixedPercBucket,ClusterSize) %>% count() %>% spread(ClusterSize,n))
View(mbClusterSummary %>% filter(ResolvedType=='Line') %>% group_by(MixedPercBucket,ClusterSize) %>% count() %>% spread(ClusterSize,n))


# 0b. Comparion for Shared variants of resultant cluster type and size
mbSvData = merge(mbSvData,mbClusterData %>% select(SampleId,ClusterId,OverlapType),by=c('SampleId','ClusterId'),all.x=T)
sharedSvData = mbSvData %>% filter(MatchType=='Shared')
View(sharedSvData)
sharedCombinedData = merge(sharedSvData,sharedSvData %>% select(OtherSvId=SvId,OtherClusterId=ClusterId,OtherSampleId=SampleId,
                                                                OtherResolvedType=ResolvedType,OtherClusterCount=ClusterCount),by=c('OtherSvId'),all.x=T)
View(sharedCombinedData)

# distribution by clusterSize
sharedCombinedData = sharedCombinedData %>% mutate(ClusterSize=2**round(log(ClusterCount,2)),OtherClusterSize=2**round(log(OtherClusterCount,2)))

View(sharedCombinedData %>% group_by(ClusterSize,OtherClusterSize) %>% count() %>% spread(OtherClusterSize,n))
View(sharedCombinedData %>% filter(ResolvedType!='DUP_BE'&OtherResolvedType!='DUP_BE') %>% group_by(ResolvedType,OtherResolvedType) %>% count() %>% spread(OtherResolvedType,n))



# 1. Shared - Private cluster merges in super sets and mixed clusters
mbMergeData = mbMergeDataAll %>% filter(!((MatchType1=='Shared'&MatchType2=='Partial')|(MatchType1=='Partial'&MatchType2=='Shared')))
nrow(mbMergeData %>% filter((MatchType1=='Shared'&MatchType2=='Partial')|(MatchType1=='Partial'&MatchType2=='Shared')))
View(mbMergeData %>% group_by(ClusterReason) %>% count())

mbMergeData = merge(mbMergeData,mbClusterData,by=c('SampleId','ClusterId'),all.x=T)
View(mbMergeData)
View(mbMergeData %>% group_by(ClusterReason,ResolvedType) %>% count() %>% spread(ClusterReason,n))

# compared with counts of merge types by cluster
allClusterReasons = mbSvDataActuals %>% filter(ResolvedType!='LowQual') %>%
  mutate(ResolvedType=ifelse(ResolvedType=='None'|ResolvedType=='ComplexChain','Complex',as.character(ResolvedType)),
         ClusterSize=ifelse(ClusterCount<=3,'Small',ifelse(ClusterCount<=10,'Medium','Large'))) %>%
  group_by(ResolvedType,ClusterSize) %>% 
  summarise(#SvCount=n(),
            ArmEndPloidy=sum(grepl('ArmEndPloidy',ClusterReason)),
            BEPloidy=sum(grepl('BEPloidy',ClusterReason)),
            ComArm=sum(grepl('ComArm',ClusterReason)),
            DelDupInv=sum(grepl('DelDupInv',ClusterReason)),
            Foldback=sum(grepl('Foldback',ClusterReason)),
            LOH=sum(grepl('LOH',ClusterReason)),
            LOHChain=sum(grepl('LOHChain',ClusterReason)),
            LooseOverlap=sum(grepl('LooseOverlap',ClusterReason)),
            Prox=sum(grepl('Prox',ClusterReason)),
            Single=sum(grepl('Single',ClusterReason)))

View(allClusterReasons)

allClusterReasonsData = allClusterReasons %>% gather('ClusterReason','Count', 3:ncol(allClusterReasons))
View(allClusterReasonsData)

allClusterReasonsSummary = allClusterReasonsData %>% group_by(ResolvedType,ClusterReason) %>% summarise(Count=sum(Count))

acrPlot = (ggplot(allClusterReasonsSummary, aes(x=ResolvedType, y=Count, fill=ClusterReason))
                       + geom_bar(stat = "identity", colour = "black")
                       + labs(x = "Cluster Type", y="SV Count")
                       + scale_fill_manual(values = defaultPlotColours)
                       + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                       + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                       + theme(axis.text.x = element_text(angle=90, hjust=1,size=10))
                       + coord_flip())

plot(acrPlot)

acrPlot2 = (ggplot(allClusterReasonsSummary %>% filter(ResolvedType=='Complex'|ResolvedType=='SimpleChain'), aes(x=ResolvedType, y=Count, fill=ClusterReason))
           + geom_bar(stat = "identity", colour = "black")
           + labs(x = "Cluster Type", y="SV Count")
           + scale_fill_manual(values = defaultPlotColours)
           + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
           + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
           + theme(axis.text.x = element_text(angle=90, hjust=1,size=10))
           + coord_flip())

plot(acrPlot2)


# express as a percentage of all clustering reasons
# exclude any cluster with a line type

svClusterReasons = allClusterReasonsData %>% filter(ResolvedType!='Line'&ClusterReason!='Single') %>% 
  group_by(ClusterReason,ClusterSize) %>% summarise(TotalCount=sum(Count))
View(svClusterReasons)


View(mbMergeData)
tmp = merge(mbMergeData,mbClusterData %>% select(SampleId,ClusterId,ClusterCount),by=c('SampleId','ClusterId'),all.x=T)
View(tmp)
privSharedClusterReasons = merge(mbMergeData,mbClusterData %>% select(SampleId,ClusterId,ClusterCount),by=c('SampleId','ClusterId'),all.x=T) %>% 
  mutate(ClusterSize=ifelse(ClusterCount<=3,'Small',ifelse(ClusterCount<=10,'Medium','Large'))) %>%
  group_by(ClusterReason,ClusterSize) %>% summarise(PrivateSharedCount=n())

# privSharedClusterReasons = mbMergeData %>% group_by(ClusterReason) %>% summarise(PrivateSharedCount=n())
View(privSharedClusterReasons)

# 2. The background false-clustering rate is around 5% if proximity is a fair gauge for that
clusterReasonCompare = merge(svClusterReasons,privSharedClusterReasons,by=c('ClusterReason','ClusterSize'),all.x=T)
clusterReasonCompare[is.na(clusterReasonCompare)] = 0
clusterReasonCompare$OverclusteringPerc = round(clusterReasonCompare$PrivateSharedCount/clusterReasonCompare$TotalCount,3) * 100
View(clusterReasonCompare)

ocrPlot2 = (ggplot(clusterReasonCompare %>% filter(ClusterReason!='BEPloidy'&ClusterReason!='LOHChain'), 
                   aes(x=reorder(ClusterReason,OverclusteringPerc), y=OverclusteringPerc, fill=ClusterReason))
            + geom_bar(stat = "identity", colour = "black")
            + labs(x = "Cluster Reason", y="Overclustering %")
            + scale_fill_manual(values = defaultPlotColours)
            + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
            + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
            + theme(axis.text.x = element_text(angle=90, hjust=1,size=10))
            + theme(legend.position="none")
            + facet_wrap(~ClusterSize)
            + coord_flip())

plot(ocrPlot2)



# 3. Many subset clusters are single LowQual SVs, about 5% (484 of 9800)
View(mbSvData %>% filter(ResolvedType=='LowQual') %>% group_by(OverlapType) %>% count())
View(mbSvData %>% filter(ResolvedType=='LowQual') %>% group_by(MatchType) %>% count())
View(mbSvData %>% filter(ResolvedType=='LowQual'))

lowQual = mbSvData %>% filter(ResolvedType=='LowQual'&MatchType!='Private')
lowQual = merge(lowQual,mbSvData %>% select(OtherSvId=SvId,OtherClusterId=ClusterId,OtherSampleId=SampleId,OtherResolvedType=ResolvedType),by=c('OtherSvId'),all.x=T)
View(lowQual %>% group_by(OtherResolvedType) %>% count())
nrow(lowQual %>% filter(OtherResolvedType!='LowQual'))

# make-up of other cluster where low-qual is a subset
lowQualSubsets = mbSvData %>% filter(ResolvedType=='LowQual'&OverlapType=='Subset')
nrow(lowQualSubsets)
lowQualSubsets = merge(lowQualSubsets,mbSvData %>% select(OtherSvId=SvId,OtherClusterId=ClusterId,OtherSampleId=SampleId),by=c('OtherSvId'),all.x=T)
lowQualSubsets = lowQualSubsets %>% distinct()
lowQualSubsets = merge(lowQualSubsets,mbClusterData %>% 
                         select(OtherSampleId=SampleId,OtherClusterId=ClusterId,OtherClusterCount=ClusterCount,OtherResolvedType=ResolvedType,
                                OtherOverlapType=OverlapType,PrivateCount,SharedCount,MatchingClustersCount,OtherClusterIds),
                       by=c('OtherSampleId','OtherClusterId'),all.x=T)
View(lowQualSubsets)
View(lowQualSubsets %>% group_by(OtherResolvedType,
                                 ClusterSize=ifelse(OtherClusterCount<=3,'Small',ifelse(OtherClusterCount<=10,'Medium','Large'))) %>% count())

View(lowQualSubsets %>% filter(PrivateCount==0))

View(mbSvData %>% filter(MatchType=='Partial') %>% group_by(Type,ResolvedType=ifelse(ResolvedType=='LowQual','LowQual','NotLowQual'),Type) %>% 
       count() %>% spread(Type,n))


View(mbSvData %>% filter((SampleId=='CPCT02010347T'&ClusterId==43)|(SampleId=='CPCT02010347TII'&ClusterId %in% c(54,55))))

View(mbMergeDataAll %>% filter(Type))

# 4. Rate for shared SVs being subsets of a larger cluster by reason
complexSuperSets = mbClusterData %>% filter(OverlapType=='ComplexSuperset'&PrivateCount==0&MatchingClustersCount==2)
View(complexSuperSets)
View(complexSuperSets %>% group_by(MatchingClustersCount) %>% count())
complexSuperSets = complexSuperSets %>% separate(OtherClusterIds,c('OtherCluster1','OtherCluster2'),sep = ';')

mixedSuperSets = mbClusterData %>% filter(OverlapType=='Mixed'&PrivateCount==0)
View(mixedSuperSets)
View(mbSvData %>% filter((SampleId=='CPCT02120047T'&ClusterId==121)|(SampleId=='CPCT02120047TII'&(ClusterId==128|ClusterId==129))))

cssSvData = merge(complexSuperSets,mbSvDataActuals %>% select(SampleId,ClusterId,Id,ClusterReason),by=c('SampleId','ClusterId'),all.x=T)
View(cssSvData)
cssSvData = merge(cssSvData,mbSvData %>% filter(MatchType!='Private') %>% select(Id=SvId,OtherSvId),by='Id',all.x=T)
cssSvData = merge(cssSvData,mbSvDataActuals %>% select(OtherSvId=Id,OtherClusterReason=ClusterReason,OtherResolvedType=ResolvedType),by='OtherSvId',all.x=T)
cssSvData = merge(cssSvData,mbSvDataActuals %>% select(OtherSvId=Id,OtherResolvedType=ResolvedType),by='OtherSvId',all.x=T)
View(cssSvData)

sharedDiffClustering = cssSvData %>% group_by(SampleId,ClusterId) %>% summarise(ClusterCount=first(ClusterCount),
                                                              LowQualCount=sum(OtherResolvedType=='LowQual'),
                                                              ArmEndPloidy=sum(grepl('ArmEndPloidy',ClusterReason)),
                                                              OtherArmEndPloidy=sum(grepl('ArmEndPloidy',OtherClusterReason)),
                                                              BEPloidy=sum(grepl('BEPloidy',ClusterReason)),
                                                              OtherBEPloidy=sum(grepl('BEPloidy',OtherClusterReason)),
                                                              ComArm=sum(grepl('ComArm',ClusterReason)),
                                                              OtherComArm=sum(grepl('ComArm',OtherClusterReason)),
                                                              DelDupInv=sum(grepl('DelDupInv',ClusterReason)),
                                                              OtherDelDupInv=sum(grepl('DelDupInv',OtherClusterReason)),
                                                              Foldback=sum(grepl('Foldback',ClusterReason)),
                                                              OtherFoldback=sum(grepl('Foldback',OtherClusterReason)),
                                                              LOH=sum(grepl('LOH',ClusterReason)),
                                                              OtherLOH=sum(grepl('LOH',OtherClusterReason)),
                                                              LOHChain=sum(grepl('LOHChain',ClusterReason)),
                                                              OtherLOHChain=sum(grepl('LOHChain',OtherClusterReason)),
                                                              LooseOverlap=sum(grepl('LooseOverlap',ClusterReason)),
                                                              OtherLooseOverlap=sum(grepl('LooseOverlap',OtherClusterReason)),
                                                              Prox=sum(grepl('Prox',ClusterReason)),
                                                              OtherProx=sum(grepl('Prox',OtherClusterReason)))


sharedDiffClustering = sharedDiffClustering %>% mutate(ArmEndPloidyDiff=ArmEndPloidy!=OtherArmEndPloidy,
                                                       BEPloidyDiff=BEPloidy!=OtherBEPloidy,
                                                       ComArmDiff=ComArm!=OtherComArm,
                                                       DelDupInvDiff=DelDupInv!=OtherDelDupInv,
                                                       FoldbackDiff=Foldback!=OtherFoldback,
                                                       LOHDiff=LOH!=OtherLOH,
                                                       LooseOverlapDiff=LooseOverlap!=OtherLooseOverlap,
                                                       ProxDiff=Prox!=OtherProx)

View(sharedDiffClustering)
nrow(sharedDiffClustering %>% filter(LowQualCount>0)) # 30/112
nrow(sharedDiffClustering %>% filter(LowQualCount==0&ArmEndPloidyDiff)) # 15
nrow(sharedDiffClustering %>% filter(LowQualCount==0&BEPloidyDiff)) # 1
nrow(sharedDiffClustering %>% filter(LowQualCount==0&ComArmDiff)) # 1
nrow(sharedDiffClustering %>% filter(LowQualCount==0&DelDupInvDiff)) # 18
nrow(sharedDiffClustering %>% filter(LowQualCount==0&FoldbackDiff)) # 5
nrow(sharedDiffClustering %>% filter(LowQualCount==0&LOHDiff)) # 31
nrow(sharedDiffClustering %>% filter(LowQualCount==0&LooseOverlapDiff)) # 2
nrow(sharedDiffClustering %>% filter(LowQualCount==0&ProxDiff)) # 3


# 5. Exploration of over-clustering for certain reasons
aepMerges = read.csv('~/logs/cr_arm_end_ploidy.csv')
nrow(aepMerges)
aepMerges = aepMerges %>% filter(SampleId %in% doubleBiopsy$sampleId) 
aepMerges = aepMerges %>% select(-Time)
View(aepMerges)

aepSvData = merge(aepMerges,mbSvData %>% select(SvId,MatchType1=MatchType,OverlapType,ClusterCount,ResolvedType),by.x='SvId1',by.y='SvId',all.x=T)
aepSvData = merge(aepSvData,mbSvData %>% select(SvId,MatchType2=MatchType),by.x='SvId2',by.y='SvId',all.x=T)

aepSvData = aepSvData %>% filter(!(MatchType1=='Private'&MatchType2=='Private')) %>%
    mutate(SharedPrivate=(MatchType1=='Private'&MatchType2!='Private')|(MatchType2=='Private'&MatchType1!='Private'),
           ArmEndMAPBucket=ifelse(ArmEndMAP>0,2**round(log(ArmEndMAP,2)),0),
           ClusterBMPBucket=ifelse(ClusterBoundaryMinPloidy>0,2**round(log(ClusterBoundaryMinPloidy,2)),0),
           ClusterVsArmDiff=ClusterBoundaryMinPloidy-ArmEndMAP,
           ClusterVsArmDiffBucket=ifelse(ClusterVsArmDiff>0,2**round(log(ClusterVsArmDiff,2)),0))

View(aepSvData)
View(aepSvData %>% group_by(SharedPrivate,FacingCentromere,ClusterVsArmDiffBucket) %>% count() %>% spread(SharedPrivate,n))
View(aepSvData %>% group_by(SharedPrivate,ArmEndMAPBucket) %>% count() %>% spread(SharedPrivate,n))
View(aepSvData %>% group_by(SharedPrivate,ClusterBMPBucket) %>% count() %>% spread(SharedPrivate,n))

beMerges = read.csv('~/logs/cr_be_ploidy.csv')
nrow(beMerges)
beMerges = beMerges %>% filter(SampleId %in% doubleBiopsy$sampleId) 
beMerges = beMerges %>% select(-Time)
View(beMerges)

beSvData = merge(beMerges,mbSvData %>% select(SvId,MatchType1=MatchType,OverlapType,ClusterCount,ResolvedType),by.x='SvId1',by.y='SvId',all.x=T)
beSvData = merge(beSvData,mbSvData %>% select(SvId,MatchType2=MatchType),by.x='SvId2',by.y='SvId',all.x=T)

beSvData = beSvData %>% filter(!(MatchType1=='Private'&MatchType2=='Private')) %>%
  mutate(SharedPrivate=(MatchType1=='Private'&MatchType2!='Private')|(MatchType2=='Private'&MatchType1!='Private'),
         BPBucket=ifelse(BreakendPloidy>0,2**round(log(BreakendPloidy,2)),0),
         MAPBucket=ifelse(MajorAP>0,2**round(log(MajorAP,2)),0),
         Diff=BreakendPloidy-MajorAP,
         DiffBucket=ifelse(Diff>0,2**round(log(Diff,2)),0))

View(beSvData)
View(beSvData %>% group_by(SharedPrivate,DiffBucket) %>% count() %>% spread(SharedPrivate,n))
View(beSvData %>% group_by(SharedPrivate,BPBucket) %>% count() %>% spread(SharedPrivate,n))
View(beSvData %>% group_by(SharedPrivate,MAPBucket) %>% count() %>% spread(SharedPrivate,n))


# View(mbSvDataActuals %>% filter(Id %in% c(14976201,14976203,15907324,15907327)))


# complexSuperSets = complexSuperSets %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))

mbSvDataActuals = mbSvDataActuals %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))

cssSVs = mbSvData %>% filter(OverlapType=='ComplexSuperset') %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))
subsSVs = mbSvData %>% filter(OverlapType=='Subset') %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))
supersetSVs = mbSvDataActuals %>% filter(Id %in% complexSuperSets$ClusterId)
View(mbSvData)

View() %>% filter(ResolvedType!='Line') %>% group_by(ClusterReason) %>% summarise(TotalCount=sum(Count))







# INVESTIGATIONS
View(mbSvData)
View(mbSvData %>% filter(MatchType=='Partial') %>% arrange(PatientId,SvId))
View(mbSvData %>% filter(MatchType=='Partial') %>% group_by(PatientId) %>% summarise(Count=n(),SglCount=sum(Type=='SGL')))

View(mbSvData %>% filter(PatientId=='CPCT02050120'))

# example1
View(mbSvDataActuals %>% filter(SampleId=='DRUP01050018T'&ClusterId==84))
View(mbSvDataActuals %>% filter(SampleId=='CPCT02050120T'&ClusterId==292))
View(mbClusterData %>% filter(SampleId=='DRUP01050018T'&ClusterId==84))
View(mbClusterData %>% filter(SampleId=='CPCT02050120T'&(ClusterId==82|ClusterId==547|ClusterId==292)))
View(mbSvData %>% filter(SampleId=='DRUP01050018T'&ClusterId==84))
View(mbSvData %>% filter(SampleId=='CPCT02050120T'&ClusterId==547))
View(mbSvData %>% filter(SampleId=='CPCT02050120T'&ClusterId==82))

# example 2
View(mbSvData %>% filter(SampleId=='CPCT02010894T'&ClusterId==1))
View(mbSvData %>% filter(SvId %in% c(16422010,16422068)))
View(mbSvDataActuals %>% filter(Id %in% c(16320712,16320711,16320702)))

View(mbSvDataActuals %>% filter(Id %in% c(16422068,16422010)))

# example 3
View(mbSvData %>% filter(SampleId=='CPCT02110042T'&ClusterId==181))
View(mbSvData %>% filter(SvId %in% c(16422010,16422068)))
View(mbSvDataActuals %>% filter(Id %in% c(15516993,15516995,15516996,15517514)))
View(mbSvData %>% filter(SampleId=='DRUP01110007T'&ClusterId==197))
View(mbSvData %>% filter(SampleId=='DRUP01110007T'&ClusterId==142))
View(mbClusterData %>% filter(SampleId=='DRUP01110007T'&ClusterId==197))
View(mbClusterData %>% filter(SampleId=='DRUP01110007T'&ClusterId==142))

View(mbSvDataActuals %>% filter(SampleId=='CPCT02110042T'&ClusterId==181))


View(mbClusterData %>% filter(ClusterType=='SimpleSuperset'))
supersets = mbClusterData %>% filter(ClusterType=='SimpleSuperset')
nrow(supersets)
View(supersets)
subsets = mbClusterData %>% filter(ClusterType=='Subset')
nrow(subsets)
View(subsets %>% filter(!(ClusterId %in% supersets$OtherClusterIds)))





####### OLD R DATA PREP ROUTINES
load(file = "~/data/sv/cohort.RData")
View(multipleBiopsyCohort)

allSvs = read.csv("~/data/sv/SVA_SVS.csv")
mbcSvs = allSvs %>% filter(SampleId %in% multipleBiopsyCohort$sampleId)
nrow(mbcSvs)
save(mbcSvs, file = "~/data/sv/mbcSvs.RData")



reformattedSVS = dbSvs %>% select(sampleId = SampleId, everything()) %>% left_join(doubleBiopsy %>% select(sampleId, n), by = "sampleId")

query = reformattedSVS %>% filter(n == 1)
subject = reformattedSVS %>% filter(n == 2)
nrow(query)
nrow(subject)
View(query)

# tidy up preparation objects
rm(reformattedSVS)
rm(dbSvs)
rm(allSvs)
rm(mbcSvs)


save(query, subject, file = "~/data/sv/multipleBiopsy.RData")

# find matching SVs by chromosome, position and orientation, and store the results in a table mapping an index in the query data to an index in the subject data
sv_overlaps <- function(query, subject, maxgap = -1) 
{
  require(tidyr)
  require(dplyr)
  require(GenomicRanges)
  
  queryStartRange <- GRanges(paste0(query$patientId, query$ChrStart, query$OrientStart), IRanges(query$PosStart, query$PosStart))
  subjectStartRange <- GRanges(paste0(subject$patientId, subject$ChrStart, subject$OrientStart), IRanges(subject$PosStart, subject$PosStart))
  startOverlaps = data.frame(findOverlaps(queryStartRange, subjectStartRange, type="any", select="all", maxgap = maxgap))
  
  queryEndRange <- GRanges(paste0(query$patientId, query$ChrEnd, query$OrientEnd), IRanges(query$PosEnd, query$PosEnd))
  subjectEndRange <- GRanges(paste0(subject$patientId, subject$ChrEnd, subject$OrientEnd), IRanges(subject$PosEnd, subject$PosEnd))
  endOverlaps = data.frame(findOverlaps(queryEndRange, subjectEndRange, type="any", select="all", maxgap = maxgap))
  
  overlaps = inner_join(startOverlaps, endOverlaps, by = c("queryHits", "subjectHits"))
  
  overlapQueryData = query[overlaps$queryHits, ] %>%
    mutate(queryHits = overlaps$queryHits) %>%
    select(queryHits, sampleId, ChrStart, ChrEnd, PosStart, PosEnd, OrientStart, OrientEnd)
  
  overlapSubjectData = subject[overlaps$subjectHits, ] %>%
    mutate(subjectHits = overlaps$subjectHits) %>%
    select(subjectHits, subjectPosStart = PosStart, subjectPosEnd = PosEnd, subjectOrientStart = OrientStart, subjectOrientEnd = OrientEnd)
  
  overlapsData = bind_cols(overlapQueryData, overlapSubjectData) %>%
    filter(OrientStart == subjectOrientStart, OrientEnd == subjectOrientEnd) %>%
    select(-subjectOrientStart, -subjectOrientEnd) %>%
    mutate(PosStartDiff = abs(PosStart - subjectPosStart), PosEndDiff = abs(PosEnd - subjectPosEnd), positionDiff = PosStartDiff + PosEndDiff) %>%
    group_by(ChrStart, ChrEnd, PosStart, PosEnd,OrientStart,OrientEnd) %>%
    top_n(1, -positionDiff) %>%
    group_by(queryHits) %>%
    top_n(1, -subjectHits)
  
  return (overlapsData %>% select(queryHits, subjectHits))
}

####### ANALYSIS


load(file = "~/data/sv/multipleBiopsy.RData")
nrow(query)
nrow(subject)
overlaps = sv_overlaps(query, subject, 10)
nrow(overlaps)
View(overlaps)

names(subject)[names(subject) == 'sampleId'] <- 'SampleId'
names(query)[names(query) == 'sampleId'] <- 'SampleId'
names(subject)[names(subject) == 'scope'] <- 'Scope'
names(query)[names(query) == 'scope'] <- 'Scope'

# any matched SV will have the scope set to 'Shared' and will then show the matching cluster-related data from the other SV
# unmatched SVs will haev scope set to 'Private' and have 'NA' for all cluster-related fields 
subject$Scope <- "Private"
query$Scope <- "Private"
query$ClusterCountMatch <- NA
query$ClusterCountDesc <- NA

subject[overlaps$subjectHits, "Scope"] <- "Shared"
query[overlaps$queryHits, "Scope"] <- "Shared"

query[overlaps$queryHits, "ClusterCountMatch"] <- query[overlaps$queryHits, "ClusterCount"] == subject[overlaps$subjectHits, "ClusterCount"]
query[overlaps$queryHits, "ClusterDescMatch"] <- query[overlaps$queryHits, "ClusterDesc"] == subject[overlaps$subjectHits, "ClusterDesc"]
query[overlaps$queryHits, "ResolvedTypeMatch"] <- query[overlaps$queryHits, "ResolvedType"] == subject[overlaps$subjectHits, "ResolvedType"]
query[overlaps$queryHits, "ResolvedTypeOther"] <- subject[overlaps$subjectHits, "ResolvedType"]
query[overlaps$queryHits, "SubjectClusterReason "] <- subject[overlaps$subjectHits, "ClusterReason"]

query[overlaps$queryHits, "SubjectClusterCount"] <- subject[overlaps$subjectHits, "ClusterCount"]
query[overlaps$queryHits, "SubjectClusterId"] <- subject[overlaps$subjectHits, "ClusterId"]
subject[overlaps$subjectHits, "QueryClusterCount"] <- query[overlaps$queryHits, "ClusterCount"]
subject[overlaps$subjectHits, "QueryClusterId"] <- query[overlaps$queryHits, "ClusterId"]
subject[overlaps$subjectHits, "QueryClusterReason"] <- query[overlaps$queryHits, "ClusterReason"]

subjectClusterCount = subject %>% select(patientId, SubjectClusterId = ClusterId, SubjectClusterCount = ClusterCount) %>% distinct()
subjectClusterCount %>% group_by(patientId, SubjectClusterId) %>% count() %>% filter(n > 1)

# classify each cluster as: Private - all SVs are only in one sample, Exact - all shared SVs match,
# SimpleSuperset - one sample has all the shared of a single other cluster, plus some private SVs
# Subset - one sample has no private SVs, and some but not all the shared of a single other cluster
# ComplexSuperset = one sample has no or some private, shared overlapping with more than 1 other cluster
# otherwise mixed
queryClusterMap = query %>% 
  filter(Scope == 'Shared') %>% 
  select(SampleId, patientId, ClusterId, ClusterCount, SubjectClusterId) %>%
  distinct() %>%
  left_join(subjectClusterCount, by = c("patientId","SubjectClusterId")) %>%
  group_by(SampleId, ClusterId, ClusterCount) %>% 
  summarise(SubjectClusterIdCount = n(), SubjectClusterTotalCount = sum(SubjectClusterCount))

queryClusterStatus = query %>% 
  group_by(SampleId, ClusterId, ClusterCount, ResolvedType, Scope) %>% 
  summarise(n = n()) %>% ungroup() %>%
  spread(Scope, n, fill = 0) %>%
  left_join(queryClusterMap, by = c("SampleId", "ClusterId", "ClusterCount")) %>%
  mutate(
    SubjectClusterIdCount = ifelse(is.na(SubjectClusterIdCount), 0, SubjectClusterIdCount),
    SubjectClusterTotalCount = ifelse(is.na(SubjectClusterTotalCount), 0, SubjectClusterTotalCount)
  ) %>%
  mutate(
    Status = "Mixed",
    Status = ifelse(Shared == 0, "Private", Status),
    Status = ifelse(Private == 0 & SubjectClusterIdCount == 1 & ClusterCount == SubjectClusterTotalCount, "Exact", Status),
    Status = ifelse(Private > 0 & SubjectClusterIdCount == 1 & Shared == SubjectClusterTotalCount, "SimpleSuperset", Status),
    Status = ifelse(Private >= 0 & SubjectClusterIdCount > 1 & Shared == SubjectClusterTotalCount, "ComplexSuperset", Status),
    Status = ifelse(Private == 0 & SubjectClusterIdCount == 1 & ClusterCount < SubjectClusterTotalCount, "Subset", Status)
  ) 

View(queryClusterStatus)
View(queryClusterStatus %>% group_by(Status) %>% count())
queryClusterStatusSummary = queryClusterStatus %>% group_by(ResolvedType,Status) %>% count() %>% spread(Status,n)
View(queryClusterStatusSummary)

# same for the subject data
queryClusterCount = query %>% select(patientId, QueryClusterId = ClusterId, QueryClusterCount = ClusterCount) %>% distinct()
queryClusterCount %>% group_by(patientId, QueryClusterId) %>% count() %>% filter(n > 1)


subjectClusterMap = subject %>% 
  filter(Scope == 'Shared') %>% 
  select(SampleId, patientId, ClusterId, ClusterCount, QueryClusterId) %>%
  distinct() %>%
  left_join(queryClusterCount, by = c("patientId","QueryClusterId")) %>%
  group_by(SampleId, ClusterId, ClusterCount) %>% 
  summarise(QueryClusterIdCount = n(), QueryClusterTotalCount = sum(QueryClusterCount))

subjectClusterStatus = subject %>% 
  group_by(SampleId, ClusterId, ClusterCount, ResolvedType, Scope) %>% 
  summarise(n = n()) %>% ungroup() %>%
  spread(Scope, n, fill = 0) %>%
  left_join(subjectClusterMap, by = c("SampleId", "ClusterId", "ClusterCount")) %>%
  mutate(
    QueryClusterIdCount = ifelse(is.na(QueryClusterIdCount), 0, QueryClusterIdCount),
    QueryClusterTotalCount = ifelse(is.na(QueryClusterTotalCount), 0, QueryClusterTotalCount)
  ) %>%
  mutate(
    Status = "Mixed",
    Status = ifelse(Shared == 0, "Private", Status),
    Status = ifelse(Private == 0 & QueryClusterIdCount == 1 & ClusterCount == QueryClusterTotalCount, "Exact", Status),
    Status = ifelse(Private > 0 & QueryClusterIdCount == 1 & Shared == QueryClusterTotalCount, "SimpleSuperset", Status),
    Status = ifelse(Private >= 0 & QueryClusterIdCount > 1 & Shared == QueryClusterTotalCount, "ComplexSuperset", Status),
    Status = ifelse(Private == 0 & QueryClusterIdCount == 1 & ClusterCount < QueryClusterTotalCount, "Subset", Status)
  ) 


# look into the clusters which are either a Subset or Superset 
# these will be reverse of each other between the 2 datasets (ie queryClusterStatus and subjectClusterStatus)

queryClusterStatusSummary = queryClusterStatus %>% group_by(ResolvedType,Status) %>% count() %>% spread(Status,n)
View(queryClusterStatusSummary)
View(queryClusterStatus %>% filter(ResolvedType!='Line') %>% group_by(ResolvedType,Status) %>% count() %>% spread(Status,n))
View(subjectClusterStatus %>% filter(ResolvedType!='Line') %>% group_by(ResolvedType,Status) %>% count() %>% spread(Status,n))

View(queryClusterStatus %>% filter(ResolvedType!='Line') %>% group_by(Status) %>% count())
View(subjectClusterStatus %>% filter(ResolvedType!='Line')%>% group_by(Status) %>% count())


View(queryClusterStatus)
View(subjectClusterStatus)

View(query %>% filter(SampleId=='CPCT02390005T'&ClusterId %in% c(48,49,50)))
View(subject %>% filter(patientId=='CPCT02390005'&ClusterId %in% c(53)))

View(query %>% filter(patientId=='CPCT02010503'&ClusterId %in% c(248,249)))
View(subject %>% filter(SampleId=='CPCT02010256TII'&ClusterId %in% c(283)))
View(subject %>% filter(SampleId=='CPCT02010256TII'))
View(query %>% filter(patientId=='CPCT02010256'))
View(query %>% filter(patientId=='CPCT02010256'&(ChrStart==1|ChrEnd==1|ChrStart==3|ChrEnd==3)&(Type=='BND'|Type=='SGL')))


View(query %>% filter(Id==15671393))
View(subject %>% filter(Id==16047403|Id==16047402))





# linking Subset and Simple Superset clusters where no SVs are private

# set Status onto each SV
query2 = merge(query,queryClusterStatus %>% select(SampleId,ClusterId,Status,Private,Shared),by=c('SampleId','ClusterId'),all.x=T)
subject2 = merge(subject,subjectClusterStatus %>% select(SampleId,ClusterId,Status,Private,Shared),by=c('SampleId','ClusterId'),all.x=T)
# nrow(query2 %>% filter(is.na(Status))) all were matched
# nrow(subject2 %>% filter(is.na(Status)))

# determine the likely cause of over or additional clustering by comparing the set of clustering reasons
# if set of cluster reasons differs in shared SVs, this is likely a difference in the sample's make-up - eg copy number profile, LOH or long DEL-DUP threshold
# if not, the cause is likely in a private SV
subjectMixedSVs = subject2 %>% filter(Status=='ComplexSuperset'|Status=='Mixed') %>% filter(ResolvedType!='Line') %>%
  select(SampleId,ClusterId,Id,ClusterCount,ClusterReason,ResolvedType,QueryClusterReason,QueryClusterId,QueryClusterCount,Shared,Private,Scope)

subjectMixedSVs$ClusterReason = stringi::stri_replace_all_fixed(subjectMixedSVs$ClusterReason, 'InvOverlap', 'InvDelDup')
subjectMixedSVs$ClusterReason = stringi::stri_replace_all_fixed(subjectMixedSVs$ClusterReason, 'LongDelDup', 'InvDelDup')

View(subjectMixedSVs %>% filter(Shared>0&Private>0))
nrow(subjectMixedSVs %>% filter(Shared>0&Private>0))
nrow(subjectMixedSVs %>% filter(Shared>0&Private>0) %>% group_by(SampleId,ClusterId) %>% count())

# find examples of shared and private variants being linked
# need to merge 2 SV records by their IDs
clustersWithMixed = subjectMixedSVs %>% filter(Shared>0&Private>0)

sharedPrivateClusterReasons = data.frame(matrix(nrow=0,ncol=5))
colnames(sharedPrivateClusterReasons) = c("SampleId", "ClusterId", "SvId1", "SvId2", "ClusterReason")

find_shared_private_clustering<-function(svData)
{
  rowCount = nrow(svData)
  for(i in 1:rowCount)
  {
    rowData = svData[i,]
    sampleId = rowData$SampleId
    clusterId = rowData$ClusterId
    currentId = as.character(rowData$Id)
    clusterReason = rowData$ClusterReason

    # print(sprintf("i=%d svId(%s from %d)", i, currentId, rowData$Id))
    
    for(j in i+1:rowCount)
    {
      if(j > rowCount)
        break
      
      nextRowData = svData[j,]
      nextSampleId = nextRowData$SampleId
      nextClusterId = nextRowData$ClusterId
      
      if(nextSampleId!=sampleId|nextClusterId!=clusterId)
      {
        break
      }
      
      # check for a link between the SVs by clustering reason
      if(rowData$Scope != nextRowData$Scope)
      {
        nextClusterReason = nextRowData$ClusterReason
        nextSvId = as.character(nextRowData$Id)
        
        if(grepl(currentId,nextClusterReason))
        {
          #print(sprintf("sample(%s) cluster(%s) id(%s) reason(%s) vs id(%s) reason(%s)", 
          #              sampleId, clusterId, currentId, clusterReason, nextSvId, nextClusterReason))
          
          reasons = strsplit(nextClusterReason,';')
          
          for(i in 1:length(reasons[[1]]))
          {
            cr = reasons[[1]][i]
            if(grepl(currentId,cr))
            {
              reasonId = strsplit(cr,'_')
              reason = reasonId[[1]][1]
              # print(paste("i=",i,", reasonId=",cr,", reason=",reason,sep=''))

              sharedPrivateClusterReasons = sharedPrivateClusterReasons %>% 
                add_row(SampleId=sampleId, ClusterId=clusterId, SvId1=currentId, SvId2=nextSvId, ClusterReason=reason)

              break;
            }
          }
        }
      }
    }
  }
  
  return (sharedPrivateClusterReasons)
}

# sharedPrivateClusterReasons = find_shared_private_clustering(head(clustersWithMixed,20))
sharedPrivateClusterReasons = find_shared_private_clustering(clustersWithMixed)
View(sharedPrivateClusterReasons)
View(sharedPrivateClusterReasons %>% group_by(ClusterReason) %>% count())
write.csv(sharedPrivateClusterReasons, '~/data/sv/shared_private_cluster_reasons.csv', row.names = F, quote = F)

View(sharedPrivateClusterReasons %>% group_by(SampleId,ClusterId,ClusterReason) %>% count() %>% group_by(ClusterReason) %>% count())

spClusterReasons = merge(sharedPrivateClusterReasons,subject %>% select(SampleId,ClusterId,ClusterCount,ClusterDesc),by=c('SampleId','ClusterId'),all.x=T)
View(spClusterReasons)

View(spClusterReasons %>% filter(ClusterReason=='Prox') %>% group_by(SampleId,ClusterId) 
     %>% summarise(ClusterSize=first(2**round(log(ClusterCount,2),0)))
     %>% group_by(ClusterSize) %>% count())


View(subjectMixedSVs %>% group_by(ClusterReason) %>% count())
View(subjectMixedSVs %>% filter(Id==15230026|Id==15230025))



nrow(subjectMixedSVs)

subjectMixedClusterReasons = subjectMixedSVs %>% group_by(SampleId,ClusterId) %>%
  summarise(SharedCount=first(Shared),
            PrivateCount=first(Private),
            ResolvedType=first(ResolvedType),
            AcrPrivateLOH=sum(Scope=='Private'&grepl('LOH',ClusterReason)&!grepl('LOH',QueryClusterReason)),
            AcrSharedLOH=sum(Scope=='Shared'&grepl('LOH',ClusterReason)&!grepl('LOH',QueryClusterReason)),
            AcrPrivateComArm=sum(Scope=='Private'&grepl('ComArm',ClusterReason)&!grepl('ComArm',QueryClusterReason)),
            AcrSharedComArm=sum(Scope=='Shared'&grepl('ComArm',ClusterReason)&!grepl('ComArm',QueryClusterReason)),
            AcrPrivateArmEndPloidy=sum(Scope=='Private'&grepl('ArmEndPloidy',ClusterReason)&!grepl('ArmEndPloidy',QueryClusterReason)),
            AcrSharedArmEndPloidy=sum(Scope=='Shared'&grepl('ArmEndPloidy',ClusterReason)&!grepl('ArmEndPloidy',QueryClusterReason)),
            AcrPrivateBEPloidy=sum(Scope=='Private'&grepl('BEPloidy',ClusterReason)&!grepl('BEPloidy',QueryClusterReason)),
            AcrSharedBEPloidy=sum(Scope=='Shared'&grepl('BEPloidy',ClusterReason)&!grepl('BEPloidy',QueryClusterReason)),
            AcrPrivateInvDelDup=sum(Scope=='Private'&grepl('InvDelDup',ClusterReason)&!grepl('InvDelDup',QueryClusterReason)),
            AcrSharedInvDelDup=sum(Scope=='Shared'&grepl('InvDelDup',ClusterReason)&!grepl('InvDelDup',QueryClusterReason)),
            AcrPrivateLohChain=sum(Scope=='Private'&grepl('LohChain',ClusterReason)&!grepl('LohChain',QueryClusterReason)),
            AcrSharedLohChain=sum(Scope=='Shared'&grepl('LohChain',ClusterReason)&!grepl('LohChain',QueryClusterReason)),
            AcrPrivateSingle=sum(Scope=='Private'&grepl('Single',ClusterReason)&!grepl('Single',QueryClusterReason)),
            AcrSharedSingle=sum(Scope=='Shared'&grepl('Single',ClusterReason)&!grepl('Single',QueryClusterReason)),
            AcrPrivateFoldback=sum(Scope=='Private'&grepl('Foldback',ClusterReason)&!grepl('Foldback',QueryClusterReason)),
            AcrSharedFoldback=sum(Scope=='Shared'&grepl('Foldback',ClusterReason)&!grepl('Foldback',QueryClusterReason)),
            AcrPrivateLooseOverlap=sum(Scope=='Private'&grepl('LooseOverlap',ClusterReason)&!grepl('LooseOverlap',QueryClusterReason)),
            AcrSharedLooseOverlap=sum(Scope=='Shared'&grepl('LooseOverlap',ClusterReason)&!grepl('LooseOverlap',QueryClusterReason)))

View(subjectMixedClusterReasons)
write.csv(subjectMixedClusterReasons, '~/data/sv/multiple_biopsy_cluster_reasons.csv', row.names = F, quote = F)

subjectMixedClusterReasonsSummary = subjectMixedClusterReasons %>% group_by(ResolvedType) %>%
  summarise(LOH=sum(AcrPrivateLOH>0|AcrSharedLOH>0),
            LOHSO=sum(AcrPrivateLOH==0&AcrSharedLOH>0),
            ComArm=sum(AcrPrivateComArm>0|AcrSharedComArm>0),
            ComArmSO=sum(AcrPrivateComArm==0&AcrSharedComArm>0),
            ArmEndPloidy=sum(AcrPrivateArmEndPloidy>0|AcrSharedArmEndPloidy>0),
            ArmEndPloidySO=sum(AcrPrivateArmEndPloidy==0&AcrSharedArmEndPloidy>0),
            BEPloidy=sum(AcrPrivateBEPloidy>0|AcrSharedBEPloidy>0),
            BEPloidySO=sum(AcrPrivateBEPloidy==0&AcrSharedBEPloidy>0),
            InvDelDup=sum(AcrPrivateInvDelDup>0|AcrSharedInvDelDup>0),
            InvDelDupSO=sum(AcrPrivateInvDelDup==0&AcrSharedInvDelDup>0),
            LohChain=sum(AcrPrivateLohChain>0|AcrSharedLohChain>0),
            LohChainSO=sum(AcrPrivateLohChain==0&AcrSharedLohChain>0),
            Single=sum(AcrPrivateSingle>0|AcrSharedSingle>0),
            SingleSO=sum(AcrPrivateSingle==0&AcrSharedSingle>0),
            Foldback=sum(AcrPrivateFoldback>0|AcrSharedFoldback>0),
            FoldbackSO=sum(AcrPrivateFoldback==0&AcrSharedFoldback>0),
            LooseOverlap=sum(AcrPrivateLooseOverlap>0|AcrSharedLooseOverlap>0),
            LooseOverlapSO=sum(AcrPrivateLooseOverlap==0&AcrSharedLooseOverlap>0))
            
View(subjectMixedClusterReasonsSummary)

View(query2 %>% filter(Status=='Subset'))

# look at supersets with no privates to see the reason for clustering despite have all the same SVs
View(subject2 %>% filter(Status=='ComplexSuperset'|Status=='Mixed') %>% select(SampleId,ClusterId,ClusterCount,ClusterReason,QueryClusterReason,QueryClusterId,QueryClusterCount,Shared,everything())) 


View(subjectClusterStatus %>% filter(Status=='ComplexSuperset'))





View(query %>% filter(Scope=='Shared') %>% group_by(ResolvedType,Scope,ResolvedTypeOther) %>% count())
View(query %>% filter(Scope=='Shared') %>% group_by(ClusterCount,Scope,ClusterCountOther) %>% count())

View(query %>% filter(Scope=='Shared',ResolvedTypeMatch==F) %>% group_by(ResolvedType,Scope,ResolvedTypeOther) %>% count() %>% spread(ResolvedTypeOther,nn))
View(query %>% filter(Scope=='Shared') %>% group_by(ResolvedType,Scope, ResolvedTypeOther) %>% count() %>% spread(ResolvedTypeOther,nn))
View(query %>% filter(Scope=='Shared') %>% group_by(ClusterCount,Scope,SubjectClusterCount) %>% count() %>% spread(SubjectClusterCount,nn))
View(query %>% group_by(SampleId,ClusterId,ClusterCount,Scope) %>% tally %>% spread(Scope,n) %>% filter(!is.na(Private),!is.na(Shared)))



## incorrect copy/paste?
subjectClusterMap = subject %>% 
  filter(Scope == 'Shared') %>% 
  select(SampleId, ClusterId, ClusterCount, QueryClusterId) %>% 
  distinct() %>% 
  arrange(SampleId, ClusterId, ClusterCount, QueryClusterId) %>%
  group_by(SampleId, ClusterId, ClusterCount) %>% summarise(OtherClusterIdCount = n(), QueryClusterId = paste0(QueryClusterId, collapse = ",")) %>% 
  select(SampleId,  ClusterId,  OtherClusterIdCount)

queryClusterStatus = query %>% 
  group_by(SampleId,ClusterId,ClusterCount,SubjectClusterCount, ResolvedType,Scope) %>% 
  count() %>% ungroup() %>%
  spread(Scope, n, fill = 0) %>% mutate(SubjectClusterCount = ifelse(is.na(SubjectClusterCount), 0, SubjectClusterCount)) %>%
  left_join(queryClusterMap, by = c("SampleId", "ClusterId")) %>%
  mutate(OtherClusterIdCount = ifelse(is.na(OtherClusterIdCount), 0, OtherClusterIdCount)) %>%
  mutate(
    Status = "Mixed",
    Status = ifelse(Shared == 0, "Private", Status),
    Status = ifelse(Private == 0 & OtherClusterIdCount == 1 & ClusterCount == SubjectClusterCount, "Exact", Status),
    Status = ifelse(Private > 0 & OtherClusterIdCount == 1 & Shared == SubjectClusterCount, "Superset", Status),
    Status = ifelse(Private == 0 & OtherClusterIdCount == 1 & ClusterCount < SubjectClusterCount, "Subset", Status))
## 




####### DB retrieval

query_patient_id_lookup<-function(dbConnect) {
  query = paste(
    "SELECT CPCTCNT.patientId AS sampleId, concat('CPCT02', CPCTCNT.itemValue, LPAD(RIGHT(CPCTPN.itemValue,4), 4, '0')) as patientId",
    "  FROM drupEcrf CTCT2YN,  drupEcrf CPCTCNT, drupEcrf CPCTPN ",
    " WHERE CPCTCNT.patientId = CPCTPN.patientId AND CTCT2YN.patientId = CPCTCNT.patientId ",
    "   AND CPCTCNT.item = 'FLD.CPCTCNT' AND CPCTCNT.itemValue != ''",
    "   AND CPCTPN.item = 'FLD.CPCTPN' AND CPCTPN.itemValue != ''",
    "   AND CTCT2YN.item = 'FLD.CTCT2YN' AND CTCT2YN.itemValue = 'Yes'",
    sep = "")
  
  result = dbGetQuery(dbConnect, query)
  return (result)
}


sample_to_patient_id<-function(sampleId, lookup) {
  colnames(lookup) <- c("truncatedSampleIds", "patientIds")
  
  lookup = rbind(manual_patient_id(), lookup)
  
  index = match(substr(sampleId, 1, 12) , substr(lookup[[1]], 1, 12))
  if (is.na(index)) {
    substr(sampleId, 1, 12)
  } else {
    lookup[index, c(2)]
  }
}

manual_patient_id<-function() {
  truncatedSampleIds  = c("CPCT02020192", "CPCT02030224", "DRUP01010007", "DRUP01070024", "DRUP01050008",
                          "DRUP01010065", "DRUP01330002", "DRUP01340004", "DRUP01340003", "DRUP01340002", "DRUP01070008")
  patientIds = c("CPCT02020438", "CPCT02030292", "DRUP01010044", "CPCT02070110", "CPCT02050116",
                 "CPCT02010639", "CPCT02330049", "CPCT02340029", "CPCT02340014", "CPCT02340026", "CPCT02070023")
  return (data.frame(truncatedSampleIds, patientIds, stringsAsFactors = FALSE))
}

