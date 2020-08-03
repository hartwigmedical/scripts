library(ggplot2)
library(cowplot)
library(gridGraphics)


## Load Sample Set

sampleData = read.csv('~/data/sample_cancer_types_4500.csv') # taken Jun 15th 2020
nrow(sampleData) # 4513
View(sampleData)
View(sampleData %>% group_by(CancerType) %>% count)

write.csv(sampleData %>% select(SampleId),'~/data/sv/cohort/sample_ids_4512.csv',row.names = F,quote = F)
write.csv(sampleData,'~/data/sv/cohort/sample_data_4512.csv',row.names = F,quote = F)

# Merge Minor Cancer Types
sampleData = convert_to_major_cancer_types(sampleData,25)

# Sample de-dup
View(sampleData %>% group_by(HmfPatientId) %>% summarise(Samples=n(),FirstSample=first(SampleId),SecondSample=last(SampleId)) %>% filter(Samples>=2))
sampleDedupSummary = sampleData %>% group_by(HmfPatientId) %>% summarise(Samples=n(),SampleId=first(SampleId))
nrow(sampleDedupSummary) # 4169
View(sampleDedupSummary)

samplesDD = sampleData %>% filter(SampleId %in% sampleDedupSummary$SampleId)
View(samplesDD)

write.csv(samplesDD,'~/data/sv/cohort/sample_data_dedup.csv',row.names = F,quote = F)
samplesDD = read.csv('~/data/sv/cohort/sample_data_dedup.csv')

cancerTypeTotals = samplesDD %>% group_by(CancerType) %>% summarise(SampleCount=n())
View(cancerTypeTotals)



## Load SV Data
svData = read.csv('~/data/sv/cohort/LNX_SVS.csv')
# svData = read.csv('~/data/sv/cohort/logs/LNX_SVS.csv')
svData = svData %>% filter(SampleId %in% samplesDD$SampleId)
svData = svData %>% filter(!(ResolvedType %in% c('LOW_VAF','DUP_BE'))) # remove artefect types
nrow(svData)
nrow(svData %>% group_by(SampleId) %>% count)


get_cluster_size<-function(clusterSize)
{
  if(clusterSize<=4)
    return (as.character(clusterSize))
  else if(clusterSize<=8)
    return (as.character('5-8'))
  else if(clusterSize<=16)
    return (as.character('9-16'))
  else if(clusterSize<=32)
    return (as.character('17-32'))
  else if(clusterSize<=64)
    return (as.character('33-64'))
  else if(clusterSize<=128)
    return (as.character('65-128'))
  else
    return (as.character('>128'))
}




## Load Cluster Data
clusters = read.csv('~/data/sv/cohort/LNX_CLUSTERS.csv')
clusters = clusters %>% filter(SampleId %in% samplesDD$SampleId)
nrow(clusters %>% group_by(SampleId) %>% count)

# remove artefect types
artefects = c('LOW_VAF','DUP_BE')
excludedTypes = c('INF','PAIR_INF')
clusters = clusters %>% filter(!(ResolvedType %in% c('LOW_VAF','DUP_BE')))
View(clusters %>% group_by(ResolvedType) %>% count)


#####
## LINE and PseduoGene frequency

## Load Links Data
links = read.csv('~/data/sv/cohort/LNX_LINKS.csv')
links = links %>% filter(SampleId %in% samplesDD$SampleId)
nrow(links)
nrow(links %>% group_by(SampleId) %>% count)

View(links)
View(links %>% filter(ExonMatch!=''))

# 268 samples with pseudogenes
nrow(links %>% filter(ExonMatch!='') %>% group_by(SampleId) %>% count) 
psdGeneData = links %>% filter(ExonMatch!='') %>% select(SampleId,ClusterId,ExonMatch)
psdGeneData = psdGeneData %>% separate(ExonMatch,c('GeneName','TransId'),sep=';')
psdGeneSamples = psdGeneData %>% group_by(SampleId,GeneName) %>% count %>% group_by(SampleId) %>% summarise(Pseudogenes=n())
View(psdGeneSamples)
View(links %>% filter(ExonMatch!=''))

View(psdGeneData %>% group_by(SampleId,GeneName) %>% count %>% filter(SampleId=='WIDE01010564T'|SampleId=='WIDE01010575T'))

lineSamples = svData %>% group_by(SampleId) %>% summarise(LineBEs=sum(ResolvedType=='LINE'&(LEStart!='None'|LEEnd!='None')))
nrow(lineSamples)
View(lineSamples)

linePsdgSampleData = merge(lineSamples,psdGeneSamples,by='SampleId',all=T)
linePsdgSampleData[is.na(linePsdgSampleData)] = 0

linePsdgSampleData = linePsdgSampleData %>% mutate(LineGroup=ifelse(LineBEs==0,'0',ifelse(LineBEs<=10,'1-10',ifelse(LineBEs<=100,'11-100','>100'))),
                                                   LineBucket=ifelse(LineBEs==0,0,ifelse(LineBEs<=10,10,ifelse(LineBEs<=100,100,1000))),
                                                   PseudogeneGroup=ifelse(Pseudogenes==0,'0',ifelse(Pseudogenes==1,'1',ifelse(Pseudogenes<=3,'2-3',
                                                                    ifelse(Pseudogenes<=8,'4-8',ifelse(Pseudogenes<=16,'9-16','>16'))))),
                                                   PseudogeneBucket=ifelse(Pseudogenes==0,0,ifelse(Pseudogenes==1,1,ifelse(Pseudogenes<=3,2,
                                                                    ifelse(Pseudogenes<=8,4,ifelse(Pseudogenes<=16,8,16))))))

lineGroupTotals = linePsdgSampleData %>% group_by(LineGroup) %>% summarise(LineSampleCount=n())
View(lineGroupTotals)

linePsdgSampleData = merge(linePsdgSampleData,lineGroupTotals,by='LineGroup',all.x=T)
View(linePsdgSampleData)
View(linePsdgSampleData %>% #filter(PseudogeneBucket!=0) %>%
       mutate(LineLabel=sprintf('%s (n=%d)',LineGroup,LineSampleCount)) %>%
       group_by(LineLabel,LineBucket,PseudogeneGroup,PseudogeneBucket) %>% summarise(LineSampleCount=first(LineSampleCount),
                                                                                     Samples=n(),
                                                                                     PseudogenePercent=n()/first(LineSampleCount)))
print(ggplot(linePsdgSampleData %>% filter(PseudogeneBucket!=0) %>%
               mutate(LineLabel=sprintf('%s (n=%d)',LineGroup,LineSampleCount)) %>%
               group_by(LineLabel,LineBucket,PseudogeneGroup,PseudogeneBucket) %>% summarise(PseudogenePercent=n()/first(LineSampleCount)), 
             aes(x=reorder(LineLabel,LineBucket), y=PseudogenePercent, fill=reorder(PseudogeneGroup,-PseudogeneBucket)))
      + geom_bar(stat='identity',colour='black') 
      + scale_fill_manual(values = c('skyblue4','skyblue3','skyblue2','skyblue1','skyblue','white'))
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      # + theme(legend.title = element_text(text='Count of Pseudogenes per Sample' ))
      + labs(title = 'LINE Break Junctions vs Pseudogenes per Sample', x = 'LINE Break Junctions', y = 'Percent of Samples',fill='# of Pseudogene Insertions'))


## Homozygous Disruptions
drivers = read.csv('~/data/sv/cohort/LNX_DRIVERS.csv')
nrow(drivers)
nrow(drivers %>% filter(DriverType %in% c('DEL','HOM_DISRUPTION','AMP')))
nrow(drivers %>% group_by(SampleId) %>% count)

drivers = drivers %>% filter(SampleId %in% sampleData$SampleId)

delHomDis = drivers %>% filter(DriverType %in% c('DEL','HOM_DISRUPTION'))
View(delHomDis)

# 2607 samples from 4500 have a DEL driver
nrow(delHomDis %>% filter(DriverType=='DEL') %>% group_by(SampleId) %>% count) 

# 854 samples from 4500 have a HOM_DISRUPTION driver
nrow(delHomDis %>% filter(DriverType=='HOM_DISRUPTION') %>% group_by(SampleId) %>% count) 

topDelHomDisGenes = delHomDis %>% group_by(SampleId,Gene) %>% summarise(DelCount=sum(DriverType=='DEL'),
                                                                        HomDisCount=sum(DriverType=='HOM_DISRUPTION')) %>%
  group_by(Gene) %>% summarise(Samples=n(),
                               DelCount=sum(DelCount>0),
                               HomDisCount=sum(HomDisCount>0)) %>% arrange(-Samples)

totalSampleCount=nrow(samplesDD)
topDelHomDisGenes = topDelHomDisGenes %>% mutate(DelCohortPercent=round(DelCount/totalSampleCount,3),
                                                 HomDisCohortPercent=round(HomDisCount/totalSampleCount,3))

View(topDelHomDisGenes)

homDisTopGenes = head(topDelHomDisGenes %>% arrange(-HomDisCount),4)
delTopGenes = head(topDelHomDisGenes %>% arrange(-DelCount),10)
View(homDisTopGenes)

print(ggplot(homDisTopGenes %>% mutate(Gene=sprintf('%s n=%d',Gene,HomDisCount)),aes(x=reorder(Gene,-HomDisCount),y=paste(HomDisCohortPercent*100,'%',sep='')))
      + geom_bar(stat = "identity", colour = "black")
      + labs(title='Percent of Samples With Homozygous Disruption', x='', y='')
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))


# types of clusters causing HDs
homDisClusters = delHomDis %>% filter(DriverType=='HOM_DISRUPTION') %>% 
  mutate(ClusterType=ifelse(ClusterCount>=10|ResolvedType=='COMPLEX','COMPLEX',
                     ifelse(ResolvedType=='RECIP_INV'|ResolvedType=='RECIP_TRANS',as.character(ResolvedType),
                     ifelse(ResolvedType=='DEL'|ResolvedType=='DUP',as.character(ResolvedType),
                     ifelse(ResolvedType %in% twoBreakTypes,'TWO_BREAK','OTHER')))))


print(ggplot(homDisClusters %>% group_by(ClusterType) %>% count,aes(x=reorder(ClusterType,-n),y=n))
      + geom_bar(stat = "identity", colour = "black")
      + labs(title='Cluster Types Causing Homozygous Disruption', x='Gene', y='')
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))


delHomDisExp = read.csv('~/data/sv/cohort/isofox_del_hom_dis.sample_gene_perc_data.csv')
View(delHomDisExp)
View(delHomDisExp %>% group_by(GeneName) %>% count)
View(delHomDisExp %>% group_by(SampleId) %>% count)

View(homDisExp)
View(topDelHomDisGenes)

tsgNonDelHomDisDrivers = read.csv('~/data/sv/cohort/tsg_drivers.csv')
tsgNonDelHomDisDrivers = tsgNonDelHomDisDrivers %>% mutate(DriverType='OTHER')
View(delHomDis)
delHomDisExp2 = merge(delHomDisExp %>% select(SampleId,Gene=GeneName,TPM),rbind(delHomDis %>% select(SampleId,Gene,DriverType),tsgNonDelHomDisDrivers),
                      by=c('SampleId','Gene'),all.x=T)
delHomDisExp2 = delHomDisExp2 %>% mutate(DriverType=ifelse(is.na(DriverType),'WILD_TYPE',as.character(DriverType)))

driverSampleCounts = delHomDisExp2 %>% group_by(DriverType,Gene) %>% summarise(SampleCount=n())
delHomDisExp2 = merge(delHomDisExp2,driverSampleCounts,by=c('DriverType','Gene'),all.x=T)

View(delHomDisExp2)
View(delHomDisExp2 %>% group_by(DriverType) %>% count)


print(ggplot(delHomDisExp2 %>% filter(Gene %in% homDisTopGenes$Gene) %>% mutate(GeneLabel=sprintf('%s (%d)',Gene,SampleCount)), aes(Gene, TPM))
      + scale_y_log10()
      + theme(axis.text.x = element_text(angle = 90, hjust=1,size=7))
      + geom_boxplot(varwidth=F, fill="plum") 
      + labs(title="Expression by Driver Type", x="Cancer Type", y="TPM")
      + facet_wrap(~DriverType))

print(ggplot(delHomDisExp2 %>% filter(Gene %in% homDisTopGenes$Gene) %>% mutate(GeneLabel=sprintf('%s (%d)',Gene,SampleCount)), aes(DriverType, TPM))
      + scale_y_log10()
      + theme(axis.text.x = element_text(angle = 90, hjust=1,size=7))
      + geom_boxplot(varwidth=F, fill="plum") 
      # + geom_text(label="Look at this!", x=1,y=1,label.size = 0.35,color="black",fill="#69b3a2")
      + labs(title="Expression by Driver Type", x="", y="TPM")
      + facet_wrap(~Gene))

# Mann-Whitney / Wilcox test on distributions

View(homDisTopGenes)
for(gene in homDisTopGenes$Gene)
{
  geneExp = delHomDisExp2 %>% filter(Gene==gene) %>% filter(DriverType %in% c('WILD_TYPE','HOM_DISRUPTION')) %>% select(SampleId,Gene,DriverType,TPM)
  mwwGt = wilcox.test(TPM ~ DriverType, data=geneExp, alternative="greater") 
  mwwLt = wilcox.test(TPM ~ DriverType, data=geneExp, alternative="less")
  pValMin = pmin(mwwGt$p.value, mwwLt$p.value)
  print(sprintf('gene(%s) GT=%f LT=%f P-Val=%g', gene, mwwGt$p.value, mwwLt$p.value, pValMin))
}



## Driver AMP Gene

## TO-DO
# make CN Gain Bucket clear - eg 2-4, 5-8 etc
# EGFR change y-axis to count rather than percent, and drop Lung


View(drivers)
ampDrivers = drivers %>% filter(DriverCategory=='AMP') %>%
  mutate(AmpType=ifelse(DriverCategory!='AMP','NONE',
                        ifelse(EventType=='GAIN_ARM','GAIN_ARM',ifelse(EventType=='GAIN_CHR','GAIN_CHR',
                        ifelse(grepl('BFB',Annotations),'BFB',ifelse(grepl('DM',Annotations),'DM',ifelse(ResolvedType=='COMPLEX','COMPLEX',
                        ifelse(ResolvedType=='DUP','DUP','COMPLEX'))))))),
         CnGainBucket=ifelse(GeneMinCN>128,128,2**round(log(GeneMinCN/SamplePloidy,2))),
         GeneLoc=paste(Gene,':',Chromosome,Arm,sep=''))

ampDrivers = ampDrivers %>% filter(!(!is.na(ClusterId)&ClusterId>=0&is.na(Foldbacks)))

#View(ampDrivers %>% group_by(CnGainBucket=2**round(log(CNGain,2))) %>% count())
#View(ampDrivers %>% group_by(CnGainBucket=2**round(log(GeneMinCN/SamplePloidy,2))) %>% count())

ampSummary = ampDrivers %>% group_by(CancerType,SampleId,Gene,GeneLoc,DriverCategory,CnGainBucket) %>% 
  summarise(ClusterCount=n(),
            AmpDM=sum(AmpType=='DM'),
            AmpArm=sum(AmpType=='GAIN_ARM'),
            AmpChr=sum(AmpType=='GAIN_CHR'),
            AmpBFB=sum(AmpType=='BFB'),
            AmpComplex=sum(AmpType=='COMPLEX'),
            AmpDup=sum(AmpType=='DUP')) %>% ungroup()

ampSummary = ampSummary %>% mutate(AmpMainType=ifelse(AmpDM==ClusterCount,'DM',ifelse(AmpBFB==ClusterCount,'BFB',
                                               ifelse(AmpComplex==ClusterCount,'COMPLEX',ifelse(AmpDup==ClusterCount,'DUP',
                                               ifelse(AmpArm==ClusterCount,'ARM',ifelse(AmpChr==ClusterCount,'CHR',
                                               ifelse(AmpDM>0,'DM',ifelse(AmpBFB>0,'BFB','COMPLEX')))))))))

ampTypes = ampSummary %>% group_by(AmpMainType) %>% count() %>% arrange(AmpMainType) %>% select(AmpMainType)
ampTypes = unique(ampTypes$AmpMainType)
print(ampTypes)
ampGenes = c('MYC','ERBB2','FGFR1','AR','MDM2','CCNE1','EGFR','ZNF217','KRAS','CDK4','MET','TERT')
ampColours = getFactorColours(ampTypes,ampTypes)

ampDataSummary = ampSummary %>% group_by(Gene,ClusterType=AmpMainType) %>% summarise(Count=n()) %>% arrange(ClusterType)
ampDataSummary = merge(ampDataSummary,ampSummary %>% group_by(Gene) %>% summarise(GeneSampleCount=n()),by='Gene',all.x=T)
ampDataSummary = ampDataSummary %>% mutate(Percent=round(Count/GeneSampleCount,3))
View(ampDataSummary)

# create an ALL-gene category
ampAllGenes = ampSummary %>% group_by(Gene='ALL',ClusterType=AmpMainType) %>% summarise(Count=n()) %>% arrange(ClusterType)
ampAllGenes = merge(ampAllGenes,ampSummary %>% group_by(Gene='ALL') %>% summarise(GeneSampleCount=n()),by='Gene',all.x=T)
ampAllGenes = ampAllGenes %>% mutate(Percent=round(Count/GeneSampleCount,3))
View(ampAllGenes)
ampDataSummary = rbind(ampDataSummary,ampAllGenes) %>% filter(Gene %in% ampGenes | Gene =='ALL')
ampDataSummary = ampDataSummary %>% filter(Gene %in% ampGenes | Gene =='ALL')
ampDataSummary = ampDataSummary %>% mutate(GeneLabel=sprintf('%s (n=%d)',Gene,GeneSampleCount))

View(ampDataSummary)

print(ggplot(ampDataSummary, aes(x=reorder(GeneLabel,-GeneSampleCount), y=Percent, fill=ClusterType))
      + geom_bar(stat = "identity", colour = "black")
      + scale_fill_manual(values = ampColours)
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9))
      + labs(title='Frequency of driver amplification mechanisms overall and for selected driver genes', y='Cluster Type',x=''))

## EGFR by Cancer Type to show the difference
specificGene = 'EGFR'
specGeneData = ampSummary %>% filter(Gene==specificGene)
specGeneSummary = specGeneData %>% group_by(ClusterType=AmpMainType,
                                            CancerType=ifelse(CancerType %in% c('Lung','Nervous system'),as.character(CancerType),'Other')) %>% 
             summarise(Count=n()) %>% arrange(ClusterType)

specGeneSummary = merge(specGeneSummary,specGeneSummary %>% group_by(CancerType) %>% summarise(CancerDriverCount=sum(Count)),by='CancerType',all.x=T)
specGeneSummary = specGeneSummary %>% mutate(Percent=round(Count/nrow(specGeneData),3),
                                             CancerTypeLabel=sprintf('%s (n=%d)',CancerType,CancerDriverCount)) 
View(specGeneSummary)

geneAmpColours = getFactorColours(unique(specGeneSummary$ClusterType),ampTypes)
print(geneAmpColours)

print(ggplot(specGeneSummary, aes(x=CancerTypeLabel, y=Percent, fill=ClusterType))
      + geom_bar(stat = "identity", colour = "black")
      + scale_fill_manual(values = geneAmpColours)
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9))
      + labs(title='EGFR Amplification by Cluster Type', y='Cluster Type',x=''))


## AMPs by copy number gain
ampCnGainSummary = ampSummary %>% group_by(CnGainBucket,ClusterType=AmpMainType) %>% summarise(Count=n()) %>% arrange(ClusterType)
ampCnGainSummary = merge(ampCnGainSummary,ampSummary %>% group_by(CnGainBucket) %>% summarise(CnGainSampleCount=n()),by='CnGainBucket',all.x=T)
ampCnGainSummary = ampCnGainSummary %>% mutate(Percent=round(Count/CnGainSampleCount,3),
                                               CnGainLabel=sprintf('%s (n=%d)',CnGainBucket,CnGainSampleCount))

View(ampCnGainSummary)

print(ggplot(ampCnGainSummary %>% filter(CnGainBucket>2), aes(x=reorder(CnGainLabel,CnGainBucket), y=Percent, fill=ClusterType))
      + geom_bar(stat = "identity", colour = "black")
      + scale_fill_manual(values = ampColours)
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9))
      + labs(title='Amplification by Copy Number Gain and Cluster Type', y='Cluster Type',x=''))


cnGainPlot = driverFactorByGroupPlot(summaryData %>% select(Group=CnGainBucket,SampleId,Factor=AmpMainType,Count),ampTypes,
                                     'Gene','AMP Type','AMP Cluster Type Type by CN Gain',0)



#####
## COMPLEX Clusters

# 1. Frequency of complex events per sample
complexClusters = clusters %>% filter(ResolvedType=='COMPLEX'&ClusterCount>=3)
complexClusters = complexClusters %>% mutate(ClusterSize=ifelse(ClusterCount<=4,'3-4',ifelse(ClusterCount<=10,'5-10',ifelse(ClusterCount<=20,'11-20','>20'))))

ccSummary = complexClusters %>% group_by(SampleId,ClusterSize) %>% summarise(CountOfClusters=n(),ClusterSizeN=min(ClusterCount)) %>% 
  group_by(ClusterSize,CountOfClusters=pmin(CountOfClusters,20)) %>% summarise(Count=n(),ClusterSizeN=first(ClusterSizeN))

View(ccSummary)

print(ggplot(ccSummary) +
        geom_bar(aes(x=CountOfClusters,y=Count,fill=reorder(ClusterSize,-ClusterSizeN)),stat='identity',alpha=1) +
        labs(title='Frequency of Complex Clusters by Size',x='Clusters per Sample',y='# Clusters',fill='SVs per Cluster'))


ccSummary2 = complexClusters %>% group_by(SampleId) %>% summarise(CountOfClusters=n(),
                                                                  Clusters3_4=sum(ClusterSize=='3-4'),
                                                                  Clusters5_10=sum(ClusterSize=='5-10'),
                                                                  Clusters11_20=sum(ClusterSize=='11-20'),
                                                                  Clusters20p=sum(ClusterSize=='>20'))
View(ccSummary2)
ccSummary3 = ccSummary2 %>% gather('ClusterSize','ClusterSizeFreq',3:6)
View(ccSummary3)
ccSummary4 = ccSummary3 %>% group_by(ClusterSize,CountOfClusters=pmin(CountOfClusters,50)) %>% summarise(ClusterSizeFreq=sum(ClusterSizeFreq))
ccSummary4 = ccSummary4 %>% mutate(ClusterSizeN=ifelse(ClusterSize=='Clusters3_4',3,ifelse(ClusterSize=='Clusters5_10',5,ifelse(ClusterSize=='Clusters11_20',11,20))),
                                   ClusterSizeLabel=ifelse(ClusterSize=='Clusters3_4','3-4',ifelse(ClusterSize=='Clusters5_10','5-10',ifelse(ClusterSize=='Clusters11_20','11-20','>20'))))
View(ccSummary4)

print(ggplot(ccSummary4) +
        geom_bar(aes(x=CountOfClusters,y=ClusterSizeFreq,fill=reorder(ClusterSizeLabel,-ClusterSizeN)),stat='identity',alpha=1) +
        labs(title='Frequency of Complex Clusters by Size',x='Clusters per Sample',y='# Clusters',fill='SVs per Cluster'))


# 2. Complexity of complex events
complexSVs = svData %>% filter(ResolvedType=='COMPLEX')
sampleChrData = rbind(complexSVs %>% select(SampleId,ClusterId,ClusterCount,Chromosome=ChrStart), 
                      complexSVs %>% filter(Type=='BND') %>% select(SampleId,ClusterId,ClusterCount,Chromosome=ChrEnd))

sampleChrSummary = sampleChrData %>% group_by(SampleId,ClusterId,ClusterCount,Chromosome) %>% count() %>% group_by(SampleId,ClusterId,ClusterCount) %>% summarise(ChrCount=n())

sampleChrSummary = sampleChrSummary %>% mutate(ClusterSizeLog=2**round(log(ClusterCount,2)),
                                               ClusterSize=ifelse(ClusterCount<=2,as.character(ClusterCount),ifelse(ClusterCount<=5,'3-5',
                                                           ifelse(ClusterCount<=10,'6-10',ifelse(ClusterCount<=25,'11-25','>25')))))
View(sampleChrSummary)
View(sampleChrSummary %>% group_by(ClusterSize,ChrCount) %>% count)

print(ggplot(sampleChrSummary %>% group_by(ClusterSize,ChrCount) %>% summarise(Count=n(),ClusterCount=first(ClusterCount)), 
             aes(x=reorder(ClusterSize,ClusterCount),y=Count,fill=ChrCount))
      + geom_bar(stat = "identity", colour = "black")
      # + scale_x_log10()
      + scale_fill_gradient(low="white",high="darkblue")
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9))
      + labs(title='Complexity of Complex Events', y='Sample Count',x=''))


# 3. Amplification in complex events 
#X: Bucketed max JCN relative to sample ploidy 
#Y: Count of events
#Stack: Bucketed # of variants per complex event

samplePloidies = read.csv('~/data/sv/cohort/sample_purity_ploidy.csv')
View(samplePloidies)
clusterJcn = complexSVs %>% group_by(SampleId,ClusterId,ClusterCount) %>% summarise(MaxSvJcn=max(Jcn)) 
clusterJcn = merge(clusterJcn,samplePloidies %>% select(SampleId,Ploidy),by='SampleId',all.x=T)

clusterJcn = clusterJcn %>% mutate(RelativeJcn=round(2*MaxSvJcn/Ploidy,3),
                                   MaxJcn=ifelse(RelativeJcn<0.5,'0.5',ifelse(RelativeJcn<1.5,'1.0',ifelse(RelativeJcn<=3,'2-3',
                                           ifelse(RelativeJcn<=8,'4-8',ifelse(RelativeJcn<-16,'9-16','>16'))))),
                                   ClusterSize=ifelse(ClusterCount<=2,as.character(ClusterCount),ifelse(ClusterCount<=5,'3-5',
                                               ifelse(ClusterCount<=10,'6-10',ifelse(ClusterCount<=25,'11-25','>25')))))


clusterJcnSummary = clusterJcn %>% group_by(ClusterSize,MaxJcn) %>% summarise(Count=n(),RelativeJcn=first(RelativeJcn),ClusterCount=first(ClusterCount))

View(clusterJcnSummary)

print(ggplot(clusterJcnSummary, aes(reorder(MaxJcn,RelativeJcn),y=Count,fill=ClusterSize))
      + geom_bar(stat = "identity", colour = "black")
      # + scale_x_log10()
      # + scale_fill_gradient(low="lightblue",high="darkblue")
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9))
      + labs(title='Amplification in complex events ', y='Cluster Count',x=''))


## Relative contribution of resolved types by cancer type
# X: cancer type
# Y: % of break junctions resolved as resolved type
# Stack: Resolved type

svResolvedTypes = svData %>% filter(!(ResolvedType %in% excludedTypes)&ResolvedType!='INS') %>% group_by(SampleId,ResolvedType) %>% summarise(ClusterCount=n())
nrow(svResolvedTypes)

twoBreakTypes = c('DUP_TI','RECIP_INV_DEL_DUP','RECIP_TRANS','RECIP_INV','DEL_TI','RESOLVED_FOLDBACK','FB_INV_PAIR','RECIP_TRANS_DEL_DUP','RECIP_TRANS_DUPS','RECIP_INV_DUPS','PAIR_OTHER')
recips = c('RECIP_INV_DEL_DUP','RECIP_TRANS','RECIP_INV','RECIP_TRANS_DEL_DUP','RECIP_TRANS_DUPS','RECIP_INV_DUPS')
otherPairs = c('DUP_TI','DEL_TI','RESOLVED_FOLDBACK','FB_INV_PAIR','PAIR_OTHER')
simpleSVs = c('DEL','DUP','INS')
unbalTrans = c('UNBAL_TRANS','UNBAL_TRANS_TI')
simpleSVTypes = c('DEL','DUP','INS','SGL_PAIR_DEL','SGL_PAIR_DUP','SIMPLE_GRP','SGL_PAIR_INS')
unresolvedSVs = c('INV','SGL','INF','PAIR_INF')
complexTypes = c('DOUBLE_MINUTE','COMPLEX')

get_cluster_type<-function(resolvedType,useSimple=F)
{
  if(resolvedType=='LINE')
  {
    return (resolvedType)
  }
  else if(resolvedType %in% complexTypes)
  {
    return ('COMPLEX')
  }
  else if(resolvedType %in% simpleSVTypes)
  {
    if(useSimple)
      return ('SIMPLE_SV')
    else if(resolvedType %in% simpleSVs)
      return (resolvedType)
    else
      return ('OTHER_PAIR')
  }
  else if(resolvedType %in% recips)
  {
    return ('RECIPROCAL')
  }
  else if(resolvedType %in% otherPairs)
  {
    return ('OTHER_PAIRS')
  }
  else if(resolvedType %in% unbalTrans)
  {
    return ('UNBAL_TRANS')
  }
  else if(resolvedType %in% unresolvedSVs)
  {
    return ('UNRESOLVED')
  }
  else
    return (resolvedType)
}

# testing
print(get_cluster_type('COMPLEX'))
print(get_cluster_type('LINE'))
print(get_cluster_type('DEL'))
print(get_cluster_type('DEL',T))
print(get_cluster_type('DUP',T))
print(get_cluster_type('UNBAL_TRANS',T))
print(get_cluster_type('SGL_PAIR_DEL',T))
print(get_cluster_type('SGL_PAIR_DEL'))
print(get_cluster_type('RECIP_INV'))
print(get_cluster_type('DOUBLE_MINUTE'))
print(get_cluster_type('INF'))

svResolvedTypes$ClusterType = apply(svResolvedTypes[,c('ResolvedType'),drop=F], 1, function(x) get_cluster_type(x[1]))
View(svResolvedTypes)
View(svResolvedTypes %>% group_by(ClusterType,ResolvedType) %>% count)

svResolvedTypes = merge(svResolvedTypes,samplesDD %>% select(SampleId,CancerType),by='SampleId',all.x=T)
svResolvedTypes = merge(svResolvedTypes,cancerTypeTotals,by='CancerType',all.x=T)
svResolvedTypes = merge(svResolvedTypes,svResolvedTypes %>% group_by(CancerType) %>% summarise(CancerSvTotal=sum(ClusterCount)),by='CancerType',all.x=T)

svResolvedSummary = svResolvedTypes %>% group_by(CancerType,SampleCount,ClusterType,CancerSvTotal) %>% summarise(SvCount=sum(ClusterCount)) %>% 
  mutate(ClusterTypePercent=round(SvCount/CancerSvTotal,3),
         CancerLabel=sprintf('%s (n=%d)',CancerType,SampleCount))

svResolvedSummary = merge(svResolvedSummary,ctSampleTotals %>% select(CancerType,MedianSampleSvCount),by='CancerType',all.x=T)

View(svResolvedSummary)

print(ggplot(svResolvedSummary, aes(reorder(CancerLabel,SampleCount),y=ClusterTypePercent,fill=ClusterType))
      + geom_bar(stat = "identity", colour = "black")
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9))
      + theme(legend.title = element_blank())
      + labs(title='Relative Contribution of Resolved Types by Cancer Type ', y='',x='')
      + coord_flip())

ctSampleTotals = svResolvedTypes %>% group_by(CancerType,SampleId) %>% summarise(SampleSvTotal=sum(ClusterCount))

ctSampleTotals = merge(ctSampleTotals,ctSampleTotals %>% group_by(CancerType) %>% summarise(SampleCount=n(),MedianSampleSvCount=median(SampleSvTotal)),
                       by='CancerType',all.x=T)
View(ctSampleTotals)

ctSampleTotals = ctSampleTotals %>% mutate(CancerLabel=sprintf('%s (n=%d)',CancerType,SampleCount))

print(ggplot(ctSampleTotals,aes(x=reorder(CancerLabel,SampleCount),y=SampleSvTotal))
      + geom_violin(scale="area",fill="#6baed6")
      + stat_summary(fun.y="median", geom="point")
      + scale_y_log10()
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + labs(title = "",y='',x='')
      + coord_flip())


## Relative contribution of resolved types by sample
# X: sample type (arranged by cancer type and mutational load
# Y: % of break junctions resolved as resolved type
# Y2: Count of break junctions
# Stack: Resolved type

sampleResolvedTypes = svResolvedTypes %>% group_by(SampleId,ClusterType,CancerType,SampleCount) %>% summarise(ClusterTypeSvTotal=sum(ClusterCount))
sampleResolvedTypes = merge(sampleResolvedTypes,sampleResolvedTypes %>% group_by(SampleId) %>% summarise(SampleSvTotal=sum(ClusterTypeSvTotal)),by='SampleId',all.x=T)
View(sampleResolvedTypes)
View(sampleResolvedTypes %>% group_by(CancerType) %>% count)

sampleResolvedTypes = sampleResolvedTypes %>% mutate(ClusterTypePercent=ClusterTypeSvTotal/SampleSvTotal,
                                                     SampleLabel=sprintf('%s (SVs=%d)',SampleId,SampleSvTotal))

sampleSvTotals = svResolvedTypes %>% group_by(SampleId,CancerType,SampleCount) %>% summarise(SampleSvTotal=sum(ClusterCount)) %>% ungroup()

print(ggplot(sampleSvTotals
             # %>% filter(CancerType %in% c('Esophagus','Breast','Prostate','Uterus'))
             # %>% filter(CancerType %in% c('Uterus','Prostate','Breast'))
             # %>% filter(CancerType=='Uterus')
             %>% arrange(-SampleCount) %>% mutate_at(vars(CancerType), funs(factor(., levels=unique(.)))))
      + geom_bar(aes(x=reorder(SampleId,-SampleSvTotal),y=SampleSvTotal),stat="identity",width=1)
      + scale_y_log10()
      + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
              axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks = element_blank(),
              strip.text.x = element_text(angle = 90, hjust = 1,size=7),
              panel.spacing.x=unit(0.1, "lines")#, panel.spacing.y=unit(1, "lines")
              )
      + labs(title='Relative Contribution of Resolved Types by Sample', x='', y='')
      + facet_grid(~CancerType, scales="free", space="free")
)

print(ggplot(sampleResolvedTypes
             # %>% filter(CancerType %in% c('Esophagus','Breast','Prostate','Uterus'))
             # %>% filter(CancerType %in% c('Uterus','Prostate','Breast'))
             # %>% filter(CancerType=='Uterus')
             %>% arrange(-SampleCount) %>% mutate_at(vars(CancerType), funs(factor(., levels=unique(.)))))
      + geom_bar(aes(x=reorder(SampleId,-SampleSvTotal),y=ClusterTypePercent,fill=ClusterType),stat="identity",width=1)
      + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
              axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks = element_blank(),
              # strip.text.x = element_text(angle = 90, hjust = 1,size=7),
              strip.text.x = element_blank(),
              panel.spacing.x=unit(0.1, "lines")#, panel.spacing.y=unit(1, "lines")
      )
      + labs(title='', x='', y='SV % by Cluster Type')
      # + labs(title='Relative Contribution of Resolved Types by Sample', x='', y='SV % by Cluster Type')
      + facet_grid(~CancerType, scales="free", space="free")
)


## Frequency of drivers vs complexity of event
# X: Bucketed # of variants per complex event
# Y: Count of events
# Stack: # of drivers per event

View(complexClusters)

driversInComplex = clusters %>% filter(ResolvedType=='COMPLEX') %>% select(SampleId,ClusterId,ClusterCount,Foldbacks,ArmCount,MaxJcn) %>% 
  mutate(ClusterSize=ifelse(ClusterCount<=2,as.character(ClusterCount),ifelse(ClusterCount<=4,'3-4',ifelse(ClusterCount<=8,'5-8',
                     ifelse(ClusterCount<=16,'9-16',ifelse(ClusterCount<=32,'17-32',ifelse(ClusterCount<=64,'33-64',ifelse(ClusterCount<=128,'65-128','>128'))))))),
         HighJcnBucked=ifelse(MaxJcn<=2,'2',ifelse(MaxJcn<=4,'2-4',ifelse(MaxJcn<=8,'4-8','>8'))),
         HighJcn=MaxJcn>=8)

driversInComplex = merge(driversInComplex,driversInComplex %>% group_by(ClusterSize,HighJcn) %>% summarise(ClusterSizeCounts=n()),by=c('ClusterSize','HighJcn'),all.x=T)

driversByCluster = drivers %>% filter(ClusterId>=0) %>% group_by(SampleId,ClusterId) %>% summarise(DriverCount=n())
driversInComplex = merge(driversInComplex,driversByCluster,by=c('SampleId','ClusterId'),all.x=T)

driversInComplex = driversInComplex %>% mutate(DriverCountBucket=ifelse(is.na(DriverCount),0,ifelse(DriverCount<=2,as.character(DriverCount),
                                                                 ifelse(DriverCount<=5,'3-5','>5'))))
View(driversInComplex)
View(driversInComplex %>% filter(DriverCount>0) %>% group_by(DriverCountBucket,ClusterSize) %>% count)
View(driversInComplex %>% filter(DriverCount>0) %>% group_by(DriverCountBucket,ClusterSize) %>% 
       summarise(Percent=n()/first(ClusterSizeCounts),ClusterCount=first(ClusterCount),DriverCount=first(DriverCount)))

print(ggplot(driversInComplex %>% filter(DriverCount>0) %>% group_by(DriverCountBucket,ClusterSize,HighJcn) %>% 
               summarise(Percent=n()/first(ClusterSizeCounts),ClusterCount=first(ClusterCount),DriverCount=first(DriverCount)), 
             aes(reorder(ClusterSize,ClusterCount),y=Percent,fill=reorder(DriverCountBucket,-DriverCount)))
      + geom_bar(stat="identity", colour="black")
      + facet_wrap(~HighJcn)
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9))
      + labs(title='Frequency of Drivers vs Complexity of Event', y='% of Clusters with Driver Count',x='Cluster Size'))




#####
## BFB Characteristics


ampClusters = clusters %>% filter(ResolvedType=='COMPLEX'&MaxJcn>=2)
ampClusters = merge(ampClusters,samplePloidies %>% select(SampleId,Ploidy),by='SampleId',all.x=T)
ampClusters = ampClusters %>% mutate(RelativeJcn=round(2*MaxJcn/Ploidy,3))
nrow(ampClusters)

View(ampClusters %>% select(SampleId,ClusterId,ClusterCount,Foldbacks,MaxJcn,MinJcn,ArmCount,OriginArms,FragmentArms,ShortTIs,everything()))

View(ampClusters %>% filter(ArmCount==1&!grepl('DM',Annotations)&Ploidy<=3&Foldbacks>=0) %>%
       select(SampleId,ClusterId,ClusterCount,Foldbacks,RelativeJcn,MaxJcn,MinJcn,ArmCount,OriginArms,FragmentArms,ShortTIs,everything()))


View(ampClusters %>% filter(ArmCount==1&OriginArms==1&!grepl('DM',Annotations)&Ploidy<=3&Foldbacks==0) %>%
       select(SampleId,ClusterId,ClusterCount,Foldbacks,RelativeJcn,MaxJcn,MinJcn,ArmCount,OriginArms,FragmentArms,ShortTIs,everything()))

topAmpClusters = head(ampClusters %>% filter(ArmCount==1&OriginArms==1&!grepl('DM',Annotations)&Ploidy<=3&Foldbacks==0) %>% arrange(-MaxJcn),10) %>% 
  mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))

ampSVs = svData %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_')) %>% filter(SampleClusterId %in% topAmpClusters$SampleClusterId)
nrow(ampSVs)
View(ampSVs)
write.csv(ampSVs %>% group_by(SampleId,ClusterId) %>% summarise(Chromosome=first(ChrStart)) %>% select(SampleId,Chromosome),'~/logs/amp_clusters.csv',row.names = F,quote = F)

View(svData %>% filter(SampleId=='CPCT02020665T'))





#####
## AMP Driver characteristics

## number of 
View(ampDrivers)

ampDrivers = merge(ampDrivers,clusters %>% select(SampleId,ClusterId,ArmCount),by=c('SampleId','ClusterId'),all.x=T)

View(ampDrivers %>% filter(ClusterId>=0) %>% group_by(Gene,SampleId) %>% summarise(Clusters=n(),
                                                                                   ArmCount=max(ArmCount),
                                                                                   Foldbacks=max(Foldbacks)) %>%
  group_by(Gene) %>% summarise(Samples=n(),
                               AvgClusters=round(mean(Clusters),1),
                               AvgArmCount=round(mean(ArmCount),1),
                               MaxArmCount=max(ArmCount),
                               AvgFoldbacks=round(mean(Foldbacks),1),
                               MaxFoldbacks=max(Foldbacks)))

ampSummary = ampDrivers %>% group_by(CancerType,SampleId,Gene,GeneLoc,DriverCategory,CnGainBucket) %>% 
  summarise(ClusterCount=n(),
            AmpDM=sum(AmpType=='DM'),
            AmpArm=sum(AmpType=='GAIN_ARM'),
            AmpChr=sum(AmpType=='GAIN_CHR'),
            AmpBFB=sum(AmpType=='BFB'),
            AmpComplex=sum(AmpType=='COMPLEX'),
            AmpDup=sum(AmpType=='DUP')) %>% ungroup()



#####
## Shard Frequencies

# merge resolved types
shardExcluded = c('LINE','SGL','INF','PAIR_INF','INV','SGL_PAIR_DUP','SGL_PAIR_DEL','SGL_PAIR_INS','INS')
recipTypes = c('RECIP_INV_DUPS','RECIP_INV','RECIP_INV_DEL_DUP','RECIP_TRANS_DUPS','RECIP_TRANS_DEL_DUP')
miscTypes = c('DOUBLE_MINUTE','FB_INV_PAIR','SIMPLE_GRP','RESOLVED_FOLDBACK','PAIR_OTHER')

shardSvSummary = svData %>% 
  filter(!(ResolvedType %in% shardExcluded)) %>%
  mutate(ClusterType=ifelse(ResolvedType %in% c('DEL','DEL_TI'),'SYNTHETIC_DEL',
                            ifelse(ResolvedType %in% c('DUP','DUP_TI'),'SYNTHETIC_DUP',
                                   ifelse(ResolvedType %in% c('UNBAL_TRANS','UNBAL_TRANS_TI'),'UNBAL_TRANS',
                                          ifelse(ResolvedType %in% recipTypes,'RECIPROCAL',
                                                 ifelse(ResolvedType %in% miscTypes,'OTHER',as.character(ResolvedType))))))) %>%
  mutate(InShard=(ClusterCount>1&(LnkLenStart>0&LnkLenStart<=1e3)|(LnkLenEnd>0&LnkLenEnd<=1e3))) %>%
  group_by(ClusterType) %>% summarise(TotalSVs=n(),ShardSVs=sum(InShard)) %>% mutate(ShardPercent=round(ShardSVs/TotalSVs,3))

View(shardSvSummary)

plot1 = (ggplot(shardSvSummary,aes(x=reorder(ClusterType,TotalSVs),y=ShardSVs))
         + geom_bar(stat = "identity")
         + scale_y_log10()
         + labs(title = "Shards Counts",x='',y='')
         + coord_flip())

plot2 = (ggplot(shardSvSummary,aes(x=reorder(ClusterType,TotalSVs),y=ShardPercent))
      + geom_bar(stat = "identity")
      # + scale_y_log10()
      + labs(title = "% of SVs in Shards",x='',y='')
      + theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
      + coord_flip())

# Shard lengths
shardTIs = links %>% filter(TILength<=1e3)
shardTIs = shardTIs %>% filter(TILength>0) %>%  # has since been solved
  filter(!(ResolvedType %in% shardExcluded)) %>%
  mutate(ClusterType=ifelse(ResolvedType %in% c('DEL','DEL_TI'),'SYNTHETIC_DEL',
                                                  ifelse(ResolvedType %in% c('DUP','DUP_TI'),'SYNTHETIC_DUP',
                                                         ifelse(ResolvedType %in% c('UNBAL_TRANS','UNBAL_TRANS_TI'),'UNBAL_TRANS',
                                                                ifelse(ResolvedType %in% recipTypes,'RECIPROCAL',
                                                                       ifelse(ResolvedType %in% miscTypes,'OTHER',as.character(ResolvedType)))))))
shardTIs = merge(shardTIs,shardSvSummary %>% select(ClusterType,TotalSVs),by='ClusterType',all.x=T)

plot3 = ggplot(shardTIs,aes(y=TILength,x=reorder(ClusterType,TotalSVs))) +
  geom_violin(scale="area",fill="#6baed6") +
  stat_summary(fun.y="median", geom="point") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10)) +
  labs(title='Shard Lengths',x='',y='') +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
  coord_flip()

# print(plot3)
plot_grid(plot1,plot2,plot3,nrow=1,ncol=3)


##### 
## ChainFinder Comparison

cfData = read.csv('~/data/sv/LNX_CHAIN_FINDER_SVS.csv')
nrow(cfData) # 113K
nrow(cfData %>% group_by(SampleId) %>% count) # 1587

# 1. Heatmap of Linx cluster vs CF chain size

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

clusterChainCmp = cfData %>% filter(CfChainId==-1|CfChainCount>=3) %>% filter(ClusterCount>=3|CfChainId>0)
clusterChainCmp$LinxClusterSize=apply(clusterChainCmp[,c('ClusterCount'),drop=F], 1, function(x) get_cf_group_size(x[1]))
clusterChainCmp$CfChainSize=apply(clusterChainCmp[,c('CfChainCount'),drop=F], 1, function(x) get_cf_group_size(x[1]))

clusterChainSummary = clusterChainCmp %>% group_by(LinxClusterSize,CfChainSize) %>% summarise(Count=n(),CfChainCount=first(CfChainCount),LinxClusterCount=first(ClusterCount))
# View(clusterChainSummary)

clusterChainSummary2 = clusterChainSummary %>% select(LinxClusterSize,CfChainSize,Count) %>%  spread(CfChainSize,Count,fill=0) %>% gather('CfChainSize','Count',2:9)
clusterChainSummary2$LinxSizeIndex=apply(clusterChainSummary2[,c('LinxClusterSize'),drop=F], 1, function(x) get_cf_group_index(x[1]))
clusterChainSummary2$CfSizeIndex=apply( clusterChainSummary2[,c('CfChainSize'),drop=F], 1, function(x) get_cf_group_index(x[1]))
#View(clusterChainSummary2)

# plot 1
print(ggplot(clusterChainSummary2, aes(x=reorder(LinxClusterSize,LinxSizeIndex),y=reorder(CfChainSize,CfSizeIndex))) 
      + geom_tile(aes(fill=Count),stat="identity") 
      + geom_text(aes(label=Count))
      + scale_fill_gradient(low="white",high="steelblue")
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), panel.background = element_blank())
      + theme(legend.position = "none")
      + labs(title='Linx vs ChainFinder Clustering', x='Linx Cluster Size',y='ChainFinder Chain Size'))


# Clustering Reasons and Proximity

# Plot 1
combinedDistances = rbind(clusterDistances %>% mutate(Source='LINX') %>% select(ProxDistance,Source),
                          clusterDistances %>% filter(CfChainId>0) %>% mutate(Source='BOTH') %>% select(ProxDistance,Source))

print(ggplot(combinedDistances,aes(y=ProxDistance,x=Source))
      + geom_violin(scale="count",fill="#6baed6")
      + stat_summary(fun.y="median", geom="point")
      + scale_y_log10()
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + labs(x='',y='Clustering Distance'))

# Plot 2
print(ggplot(clusterDistances %>% 
               mutate(ClusterReason=ifelse(ClusterReason %in% c('LOH_CHAIN','OVERLAP_FOLDBACKS','CONSEC_BREAKS'),'OTHER',as.character(ClusterReason))) %>%
               filter(ClusterReason!='PROXIMITY'),aes(x=ClusterReason,y=ProxDistance))
      + geom_violin(scale="count",fill="#6baed6")
      + stat_summary(fun.y="median", geom="point")
      + scale_y_log10()
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + labs(x='',y='Clustering Distance by Reason'))

# Plot 3 Counts by type
View(cfData %>% group_by(Type) %>% count)
View(cfData %>% filter(Type=='INV'|Type=='BND') %>% filter(ClusterCount!=2))
View(cfData %>% filter(Type=='INV'|Type=='BND') %>% filter(OverlapType=='CHAIN-FINDER'))


cfData = cfData %>% mutate(OverlapType=ifelse(ClusterCount==1&CfChainCount>=3,'CHAIN-FINDER',
                                              ifelse(ClusterCount>=3&CfChainId==-1,'LINX',ifelse(ClusterCount==2,'CLUSTER-2','BOTH'))))

print(ggplot(cfData %>% filter(Type!='INS') %>% group_by(Type,OverlapType) %>% count,aes(x=Type,y=n,fill=OverlapType))
      + geom_bar(stat = "identity")
      # + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + labs(title='Overlap in Clustering by SV Type',x='',y='SV Count'))

print(ggplot(cfData %>% filter(Type!='INS'&OverlapType!='CLUSTER-2') %>% group_by(Type,OverlapType) %>% count,aes(x=OverlapType,y=n,fill=Type))
      + geom_bar(stat = "identity")
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + labs(title='Overlap in Clustering by SV Type',x='',y='SV Count'))

# Types and Lengths when Linx doesn't cluster
cfOnly = cfData %>% filter(ClusterCount==1&CfChainCount>=3&Type!='INS')
cfOnly = cfOnly %>% mutate(Length=ifelse(Type=='BND',0,PosEnd-PosStart))
cfOnly = merge(cfOnly,svData %>% select(SampleId,SvId=Id,FSStart,FSEnd),by=c('SampleId','SvId'),all.x=T)
View(cfOnly)
nrow(cfOnly)
View(cfOnly %>% group_by(Type) %>% count)
View(cfOnly %>% filter(is.na(FSStart)))
cfOnly = cfOnly %>% filter(!is.na(FSStart)) %>% mutate(Location=ifelse(FSStart=='true'|FSEnd=='true','FragileSite','Other'),
                                                       IsFoldback=FoldbackLnkStart>=0|FoldbackLnkEnd>=0)


print(ggplot(cfOnly %>% filter(Type!='BND'),
             aes(x=Type,y=Length))
      + geom_violin(scale="count",fill="#6baed6")
      + stat_summary(fun.y="median", geom="point")
      + scale_y_log10()
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + facet_wrap(~Location)
      + labs(title='Unclusterted SVs chained by ChainFinder',x='',y='SV Lengths'))

linxOnly = cfData %>% filter(ClusterCount>=3&CfChainId==-1&Type %in% c('DEL','DUP','INV'))
linxOnly = merge(linxOnly,svData %>% select(SampleId,SvId=Id,FoldbackLnkStart,FoldbackLnkEnd),by=c('SampleId','SvId'),all.x=T)
linxOnly = linxOnly %>% mutate(Length=PosEnd-PosStart,
                               IsFoldback=(!is.na(FoldbackLnkStart)&FoldbackLnkStart>=0)|(!is.na(FoldbackLnkEnd)&FoldbackLnkEnd>=0)) %>% filter(Length>0)
View(linxOnly)
nrow(linxOnly)

print(ggplot(linxOnly %>% mutate(Type=ifelse(Type=='INV'&IsFoldback,'FOLDBACK',as.character(Type))),aes(x=Type,y=Length))
      + geom_violin(scale="count",fill="#6baed6")
      + stat_summary(fun.y="median", geom="point")
      + scale_y_log10()
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + labs(title='Clusterted SVs, unchained by ChainFinder',x='',y='SV Lengths'))


## CF only for DELs and DUPs follows sample counts of simple DELs and DUPs
cfSamples = cfData %>% group_by(SampleId) %>% count

cfSvData = svData %>% filter(SampleId %in% cfSamples$SampleId)
cfSvSamples = cfSvData %>% group_by(SampleId) %>% count
View(cfSvData %>% group_by(SampleId) %>% count)

cfOnlyDelsDups = cfOnly %>% filter(SampleId %in% cfSvSamples$SampleId) %>% group_by(SampleId) %>% summarise(CfDelCount=sum(Type=='DEL'),CfDupCount=sum(Type=='DUP'))
View(cfOnlyDelsDups)
svDelDups = cfSvData %>% filter(Type %in% c('DUP','DEL') & ClusterCount==1) %>% group_by(SampleId) %>% summarise(SimpleDelCount=sum(Type=='DEL'),SimpleDupCount=sum(Type=='DUP'))

delDupsComp = merge(svDelDups,cfOnlyDelsDups,by='SampleId',all=T)
delDupsComp[is.na(delDupsComp)] = 0
View(delDupsComp)

print(ggplot(delDupsComp, aes(x=SimpleDelCount, y=CfDelCount))
      + geom_point()
      + labs(title = "Simple DELs vs ChainFinder Chained DELs"))

print(ggplot(delDupsComp, aes(x=SimpleDupCount, y=CfDupCount))
      + geom_point()
      + labs(title = "Simple DUPs vs ChainFinder Chained DUPs"))


get_prox_bucket<-function(distance,asLabel=T)
{
  if(distance<=100)
  {
    if(asLabel)
      return ('<0.1K')
    else
      return (100)
  }
  else if(distance<=5e3)
  {
    if(asLabel)
      return ('0.1-5K')
    else
      return (5e3)
  }
  else if(distance<=25e3)
  {
    if(asLabel)
      return ('5-25K')
    else
      return (25e3)
  }
  else if(distance<=1e5)
  {
    if(asLabel)
      return ('25-100K')
    else
      return (1e5)
  }
  else if(distance<=1e6)
  {
    if(asLabel)
      return ('100K-1M')
    else
      return (1e6)
  }
  else
  {
    if(asLabel)
      return ('>1M')
    else
      return (1e7)
  }
}

print(get_prox_bucket(250))
print(get_prox_bucket(1e7))
print(get_prox_bucket(0,T))

clusterDistances = cfData %>% filter(ClusterCount>=3&Type!='INF'&Type!='SGL')

clusterDistances$ProxDistLabel = apply(clusterDistances[,c('ProxDistance'),drop=F], 1, function(x) get_prox_bucket(x[1]))
clusterDistances$ProxDistIndex = apply(clusterDistances[,c('ProxDistance'),drop=F], 1, function(x) get_prox_bucket(x[1],F))
View(clusterDistances)
clusterDistanceSummary = clusterDistances %>% group_by(ClusterReason,ProxDistLabel,ProxDistIndex) %>% summarise(Linx=n(),ChainFinder=sum(CfChainId>0))

View(clusterDistanceSummary)

combinedDistances2 = rbind(clusterDistances %>% 
                             mutate(ClusterReason=ifelse(ClusterReason %in% c('LOH_CHAIN','OVERLAP_FOLDBACKS','CONSEC_BREAKS'),'OTHER',as.character(ClusterReason))) %>%
                             mutate(Source=ClusterReason,ProxDistance=pmax(ProxDistance,1)) %>% select(ProxDistance,Source),
                          clusterDistances %>% filter(CfChainId>0) %>% mutate(Source='BOTH') %>% select(ProxDistance,Source))

print(ggplot(combinedDistances2,aes(y=ProxDistance,x=Source))
      + geom_violin(scale="count",fill="#6baed6")
      + stat_summary(fun.y="median", geom="point")
      + scale_y_log10()
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + labs(x='',y='Clustering Distance'))


combinedDistances2$ProxDistLabel = apply(combinedDistances2[,c('ProxDistance'),drop=F], 1, function(x) get_prox_bucket(x[1]))
combinedDistances2$ProxDistIndex = apply(combinedDistances2[,c('ProxDistance'),drop=F], 1, function(x) get_prox_bucket(x[1],F))
clusterDistanceSummary3 = combinedDistances2 %>% group_by(Source,ProxDistLabel,ProxDistIndex) %>% summarise(Count=n())
View(clusterDistanceSummary3)

ignoreCRs = c('LOH_CHAIN','OVERLAP_FOLDBACKS','CONSEC_BREAKS')

print(ggplot(clusterDistanceSummary3 %>% filter(ProxDistIndex>5e3),aes(x=reorder(ProxDistLabel,ProxDistIndex),y=Count,fill=Source))
      + geom_bar(stat = "identity",position='dodge')
      # + scale_y_log10()
      + labs(title = "Clustering Reasons in Linx vs Chained by CF"))

ignoreCRs = c('LOH_CHAIN','OVERLAP_FOLDBACKS','CONSEC_BREAKS')

print(ggplot(clusterDistanceSummary %>% filter(ClusterReason!='PROXIMITY') %>% filter(!(ClusterReason %in% ignoreCRs)),
             aes(x=reorder(ProxDistLabel,ProxDistIndex)))
      + geom_bar(aes(y=LinxCount,color='LinxCount'),stat = "identity",position='dodge')
      + geom_bar(aes(y=CfCount,color='CfCount'),stat = "identity",position='dodge')
      # + scale_x_log10()
      + facet_wrap(~ClusterReason)
      + labs(title = "Clustering Reasons in Linx vs Chained by CF"))

clusterDistanceSummary2 = clusterDistanceSummary %>% gather('Tool','Count',4:5)
View(clusterDistanceSummary2)

print(ggplot(clusterDistanceSummary2  %>% filter(ClusterReason=='PROXIMITY') %>% filter(!(ClusterReason %in% ignoreCRs)),
             aes(x=reorder(ProxDistLabel,ProxDistIndex),y=Count,fill=Tool))
      + geom_bar(stat='identity',position='dodge')
      # + scale_y_log10()
      + facet_wrap(~ClusterReason)
      + labs(title = "Clustering Reasons in Linx vs Chained by CF"))





