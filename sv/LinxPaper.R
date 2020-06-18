library(ggplot2)
library(cowplot)
library(gridGraphics)


## Load Sample Set

sampleData = read.csv('~/data/sample_cancer_types_4500.csv') # taken Jun 15th 2020
nrow(sampleData) # 4513
View(sampleData)
View(sampleData %>% group_by(CancerType) %>% count)

# Merge Minor Cancer Types
sampleData = convert_to_major_cancer_types(sampleData,25)

# Sample de-dup
View(sampleData %>% group_by(HmfPatientId) %>% summarise(Samples=n(),FirstSample=first(SampleId),SecondSample=last(SampleId)) %>% filter(Samples>=2))
sampleDedupSummary = sampleData %>% group_by(HmfPatientId) %>% summarise(Samples=n(),SampleId=first(SampleId))
nrow(sampleDedupSummary) # 4169
View(sampleDedupSummary)

samplesDD = sampleData %>% filter(SampleId %in% sampleDedupSummary$SampleId)
View(samplesDD)

cancerTypeTotals = samplesDD %>% group_by(CancerType) %>% summarise(SampleCount=n())
View(cancerTypeTotals)

write.csv(sampleData %>% select(SampleId),'~/data/sv/cohort/sample_ids_4512.csv',row.names = F,quote = F)


## Load SV Data
svData = read.csv('~/data/sv/cohort/LNX_SVS.csv')
svData = svData %>% filter(SampleId %in% samplesDD$SampleId)
svData = svData %>% filter(!(ResolvedType %in% c('LOW_VAF','DUP_BE'))) # remove artefect types
nrow(svData)
nrow(svData %>% group_by(SampleId) %>% count)





## Load Cluster Data
clusters = read.csv('~/data/sv/cohort/LNX_CLUSTERS.csv')
clusters = clusters %>% filter(SampleId %in% samplesDD$SampleId)
nrow(clusters %>% group_by(SampleId) %>% count)

# remove artefect types
clusters = clusters %>% filter(!(ResolvedType %in% c('LOW_VAF','DUP_BE')))


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

lineSamples = svData %>% group_by(SampleId) %>% summarise(LineBEs=sum(ResolvedType=='LINE'&LEStart!='None')+sum(ResolvedType=='LINE'&LEEnd!='None'))
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


View(linePsdgSampleData)


write.csv(linePsdgSampleSummary,'~/logs/linePsdgSampleSummary.csv',row.names = F,quote = F)

# 'Count of LINE Break Junctions'
# add n = sample count

linePsdgSampleSummary = linePsdgSampleSummary %>% group_by(CancerType) %>% mutate(MedianLineBE=median(LineBEs),
                                                                                  PseudogeneTotal=sum(Pseudogenes)) %>% ungroup()
linePsdgSampleSummary = merge(linePsdgSampleSummary,cancerTypeTotals,by='CancerType',all.x=T)
View(linePsdgSampleSummary)

print(ggplot(linePsdgSampleSummary %>% mutate(LineBEs=LineBEs+0.5,CancerTypeLabel=sprintf('%s (n=%d)',CancerType,SampleCount)), 
             aes(x=reorder(CancerTypeLabel,-MedianLineBE),y=LineBEs))
      + geom_violin(scale="area",fill="#6baed6")
      + stat_summary(fun.y="median", geom="point")
      + scale_y_log10()
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + labs(title = "LINE Break Junctions per Sample",y='LINE Break Junctions',x=''))

print(ggplot(linePsdgSampleSummary %>% mutate(Pseudogenes=Pseudogenes+0.5,CancerTypeLabel=sprintf('%s (n=%d)',CancerType,SampleCount)), 
             aes(x=reorder(CancerTypeLabel,-PseudogeneTotal),y=Pseudogenes))
      + geom_violin(scale="area",fill="#6baed6")
      + scale_y_log10()
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + labs(title = "Pseudogenes per Sample",y='Pseudogenes',x=''))


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

sgls = svData %>% filter(Type=='SGL')
View(sgls %>% filter(DBLenStart>-1000))
View(sgls %>% filter(LnkLenStart>0))


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

complexClusters = clusters %>% filter(ResolvedType=='COMPLEX'&ClusterCount>=10)
complexSVs = svData %>% filter(ResolvedType=='COMPLEX')

# 1. Frequency of complex events per sample
ccPerSample = complexClusters %>% group_by(SampleId) %>% summarise(SampleCC=n())
View(ccPerSample)

ccZeroSamples = samplesDD %>% filter(!(SampleId %in% ccPerSample$SampleId)) %>% mutate(SampleCC=0) %>% select(SampleId,SampleCC)
nrow(ccZeroSamples)
ccPerSample = rbind(ccPerSample,ccZeroSamples)

ccPerSample = ccPerSample %>% mutate(CCBucket=ifelse(SampleCC<=2,as.character(SampleCC),ifelse(SampleCC<=5,'3-5',
                                              ifelse(SampleCC<=10,'6-10',ifelse(SampleCC<=25,'11-25','>25')))))
View(ccPerSample %>% group_by(CCBucket) %>% summarise(SampleCC=first(SampleCC),Count=n()))

print(ggplot(ccPerSample %>% group_by(CCBucket) %>% summarise(SampleCC=first(SampleCC),Count=n()), aes(x=reorder(CCBucket,SampleCC), y=Count))
      + geom_bar(stat = "identity", colour = "black")
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9))
      + labs(title='Number of Complex Clusters per Sample', y='Sample Count',x=''))

# 2. Complexity of complex events
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


Relative contribution of resolved types by cancer type
X: cancer type
Y: % of break junctions resolved as resolved type
Stack: Resolved type

Relative contribution of resolved types  by sample
X: sample type (arranged by cancer type and mutational load(
  Y: % of break junctions resolved as resolved type
  Y2: Count of break junctions
  Stack: Resolved type
  
  Frequency of drivers vs complexity of event
  X: Bucketed # of variants per complex event
  Y: Count of events
  Stack: # of drivers per event



#####
## Cohort Summary
sampleTotals = clusters %>% group_by(SampleId) %>% summarise(SampleSvCount=sum(ClusterCount))
View(sampleTotals)

clusterSummary = clusters %>% group_by(SampleId,ResolvedType,HasFoldbacks=Foldbacks>0,Replication) %>% 
  summarise(Clusters=n(),
            SvCount=sum(ClusterCount),
            MaxSvCount=max(ClusterCount),
            MedianSvCount=median(ClusterCount),
            MaxCN=max(MaxJcn))

View(clusterSummary)

twoBreakTypes = c('DUP_TI','RECIP_INV_DEL_DUP','RECIP_TRANS','RECIP_INV','DEL_TI','RESOLVED_FOLDBACK','FB_INV_PAIR','RECIP_TRANS_DEL_DUP','RECIP_TRANS_DUPS','RECIP_INV_DUPS','PAIR_OTHER')
simpleSVs = c('DEL','DUP','INS','SGL_PAIR_DEL','SGL_PAIR_DUP','SIMPLE_GRP','SGL_PAIR_INS')
unresolvedSVs = c('INV','UNBAL_TRANS','UNBAL_TRANS_TI','SGL','INF','PAIR_INF')
complexTypes = c('DOUBLE_MINUTE','COMPLEX')

clusters = clusters %>% mutate(ClusterType=ifelse(ResolvedType=='LINE','LINE',
                                                  ifelse(ResolvedType %in% complexTypes | ClusterCount>=10,'COMPLEX',
                                                         ifelse(ResolvedType %in% twoBreakTypes,'TWO_BREAK',
                                                                ifelse(ResolvedType %in% unresolvedSVs,'UNRESOLVED',
                                                                       ifelse(ResolvedType %in% simpleSVs,'SIMPLE_SV','UNKONWN'))))))

sampleClusterTypes = clusters %>% group_by(SampleId,ClusterType) %>% summarise(ClusterTypeCount=sum(ClusterCount))
sampleClusterTypes = rbind(sampleClusterTypes,clusters %>% group_by(SampleId,ClusterType='SHORT_TI') %>% summarise(ClusterTypeCount=-sum(ShortTIs)))
sampleClusterTypes = rbind(sampleClusterTypes,clusters %>% group_by(SampleId,ClusterType='SHORT_DB') %>% summarise(ClusterTypeCount=-sum(ShortDBs)))
sampleClusterTypes = merge(sampleClusterTypes,sampleTotals,by='SampleId',all.x=T)
sampleClusterTypes = merge(sampleClusterTypes,sampleData %>% select(SampleId,CancerType),by='SampleId',all.x=T)

View(sampleClusterTypes)

cancerTypes = unique(sampleData$CancerType)
# length(cancerTypes)
plotIndex = 1
cancerPlots = list()
topSampleCount=50

for(cancerType in cancerTypes)
{
  cancerData = sampleClusterTypes %>% filter(CancerType==cancerType) %>% arrange(-SampleSvCount)
  cancerData = cancerData %>% mutate(SampleLabel=sprintf('%s (%d)',SampleId,SampleSvCount))
  
  if(nrow(cancerData) > 0)
  {
    print(sprintf('CancerType(%s) has %d records', cancerType, nrow(cancerData)))
    topSamples = head(cancerData %>% group_by(SampleId) %>% summarise(SampleSvCount=first(SampleSvCount)) %>% arrange(-SampleSvCount),topSampleCount)
    
    plot = (ggplot(cancerData %>% filter(SampleId %in% topSamples$SampleId), 
                   aes(x=reorder(SampleLabel,-SampleSvCount),y=ClusterTypeCount,fill=ClusterType))
            + geom_bar(stat = "identity", colour = "black")
            + labs(title=paste('SV Counts by Type for ',cancerType,sep=''), x='', y='SV by Cluster Type')
            + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
            + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
            + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))
    
    cancerPlots[[plotIndex]] = plot
    plotIndex = plotIndex + 1
  }
}

plot_grid(cancerPlots[[1]],cancerPlots[[2]],nrow=2,ncol=1)
plot_grid(cancerPlots,ncol=1)




