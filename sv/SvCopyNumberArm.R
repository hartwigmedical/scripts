########################################
## Chromosome & Arm Copy Number Analysis

# LINX provides the following data per chromosome and arm:
# - min, max, avg and median copy number
# - CN at telomere and centromere
# - whether has an arm-level LOH

cnChrData = read.csv('~/data/sv/CN_CHR_ARM_DATA.csv')
nrow(cnChrData)

View(cnChrData %>% group_by(SampleId) %>% count())

# restrict to HPC deduped cohort
nrow(highestPurityCohort)
cnChrData = cnChrData %>% filter(SampleId %in% highestPurityCohort$sampleId)


#sampleCancerTypes = read.csv('~/data/hpc_sample_cancer_types.csv')
sampleCancerTypes = load_cancer_types('~/data/hpc_sample_cancer_types.csv',F,10)
sampleCancerTypes = sampleCancerTypes %>% select(SampleId,CancerType)
View(sampleCancerTypes)

cnChrData = cnChrData %>% filter(SampleId %in% sampleCancerTypes$SampleId)

# cnChrData = within(cnChrData, rm(CancerType))
cnChrData = merge(cnChrData,sampleCancerTypes,by='SampleId',all.x=T)
View(cnChrData %>% group_by(CancerType) %>% count())

sampleSummary = cnChrData %>% group_by(CancerType,SampleId) %>% summarise(ChrCount=n()) %>% 
  group_by(CancerType) %>% summarise(SampleCount=n(),SampleChrCount=sum(ChrCount))

View(sampleSummary)

chromosomes = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y')
chrNumbers = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
#View(chromosomes)
chrIndex = cbind(chromosomes, chrNumbers)
colnames(chrIndex) = c("Chromosome","ChrIndex")
chrIndex = as.data.frame(chrIndex)
chrIndex = chrIndex %>% mutate(ChrIndex=as.numeric(as.character(ChrIndex)))
View(chrIndex)

cnChrData = merge(cnChrData,sampleSummary,by='CancerType',all.x=T)
cnChrData = merge(cnChrData,chrIndex,by='Chromosome',all.x=T)

View(cnChrData %>% group_by(WholeChrLoh,LohP,LohQ) %>% count())

write.csv(cnChrData,'~/data/sv/chromosome_cn_data.csv', row.names = F, quote = F)
cnChrData = read.csv('~/data/sv/chromosome_cn_data.csv')
View(head(cnChrData,100))


gainThreshold = 1.4
lossThreshold = 0.7
cnCompPerc = 0.15
cnCompDiff = 0.5

cnChrData = cnChrData %>% mutate(AdjPloidy=ifelse(IsMale=='true'&(Chromosome=='X'|Chromosome=='Y'),Ploidy/2,Ploidy),
                                 ShortArm=(Chromosome==13|Chromosome==14|Chromosome==15|Chromosome==21|Chromosome==22|Chromosome=='Y'),
                                 PArmTCGain=TeloCnP>gainThreshold*AdjPloidy&CentroCnP>gainThreshold*AdjPloidy,
                                 QArmTCGain=TeloCnQ>gainThreshold*AdjPloidy&CentroCnQ>gainThreshold*AdjPloidy,
                                 PArmGain=AvgCnP>gainThreshold*AdjPloidy,
                                 QArmGain=AvgCnQ>gainThreshold*AdjPloidy,
                                 WholeChrTCGain=PArmTCGain&QArmTCGain,
                                 WholeChrAvgGain=PArmGain&QArmGain,
                                 WholeChrGain=WholeChrTCGain&WholeChrAvgGain,
                                 PvsQGain=AvgCnP>AvgCnQ*gainThreshold&(AvgCnP+AvgCnQ)/2>AdjPloidy,
                                 QvsPGain=AvgCnQ>AvgCnP*gainThreshold&(AvgCnP+AvgCnQ)/2>AdjPloidy,
                                 PArmLoss=AvgCnP<lossThreshold*AdjPloidy,
                                 QArmLoss=AvgCnQ<lossThreshold*AdjPloidy,
                                 WholeChrLoss=PArmLoss&QArmLoss,
                                 PvsQLoss=AvgCnP<AvgCnQ*lossThreshold&(AvgCnP+AvgCnQ)/2<AdjPloidy,
                                 QvsPLoss=AvgCnQ<AvgCnP*lossThreshold&(AvgCnP+AvgCnQ)/2<AdjPloidy,
                                 PCentroGain=(CentroCnP>CentroCnQ*(1+cnCompPerc)&CentroCnP>CentroCnQ+cnCompDiff),
                                 QCentroGain=(CentroCnQ>CentroCnP*(1+cnCompPerc)&CentroCnQ>CentroCnP+cnCompDiff))

cnChrData = cnChrData %>% mutate(HasCentroCNChg=PCentroGain|QCentroGain,
                                 CentroCNChg=ifelse(HasCentroCNChg,ifelse(PCentroGain,CentroCnP-pmax(CentroCnQ,0),CentroCnQ-pmax(CentroCnP,0)),0))

# 1. % with 1.4x whole chr gain
# 2. % with 1.4x one arm gain only
# 3. % with 1.4x one arm gain only AND pCentroCN =~ qCentroCN
# 4. % with <0.6x whole chr loss
# 5. % with <0.6x one arm loss only
# 6. % with <0.6x one arm loss only AND pCentroCN =~ qCentroCN
# 7. % with centromere gain or loss
# 8. % with whole chromosome LOH
# 9. % with pLOH/qLOH only

View(cnChrData)


#data = cnChrData %>% mutate(SumField=PvsQGain,FacetField=PCentroGain)
#sampleCount = nrow(data %>% group_by(SampleId) %>% count)
#panCancerData = data %>% group_by(Chromosome,ChrIndex,FacetField) %>% summarise(Count=n(),
#                                                                                SumFieldCount=sum(SumField),
#                                                                                Percent=sum(SumField)/n()) %>% 
#  mutate(CancerType='PanCancer',SampleCount=sampleCount)

chrCancerTypePlot<-function(data,title='',facetStr='',countsByFacet=F,checkFiltered=F)
{
  # if no faceting is requiring, the denominator is the sample count per cancer type
  sampleCount = nrow(data %>% group_by(SampleId) %>% count)
  
  cancerSampleCounts = data %>% group_by(CancerType,SampleId) %>% count %>% group_by(CancerType) %>% summarise(CancerSampleCount=n())
  data = merge(data,cancerSampleCounts,by='CancerType',all.x=T)
  
  if(checkFiltered)
    data = data %>% filter(!Filtered)

  if(facetStr=='')
  {
    panCancerData = data %>% group_by(Chromosome,ChrIndex) %>% summarise(Percent=sum(SumField)/n()) %>% mutate(CancerType='PanCancer',
                                                                                                               CancerSampleCount=sampleCount)
    
    perCancerData = data %>% group_by(Chromosome,ChrIndex,CancerType,CancerSampleCount) %>% summarise(Percent=sum(SumField)/n())
  }
  else if(countsByFacet)
  {
    cancerFacetSampleCounts = data %>% group_by(CancerType,SampleId,FacetField) %>% count %>% group_by(CancerType,FacetField) %>% summarise(CancerFacetSampleCount=n())
    data = merge(data,cancerFacetSampleCounts,by=c('CancerType','FacetField'),all.x=T)
    
    facetSampleCounts = data %>% group_by(FacetField,SampleId) %>% count %>% group_by(FacetField) %>% summarise(FacetSampleCount=n())
    
    panCancerData = data %>% group_by(Chromosome,ChrIndex,FacetField) %>% summarise(Percent=sum(SumField))
    panCancerData = merge(panCancerData,facetSampleCounts,by='FacetField',all.x=T)

    panCancerData = panCancerData %>% mutate(Percent=Percent/FacetSampleCount,
                                             CancerType='PanCancer',
                                             CancerSampleCount=sampleCount)
    
    panCancerData = panCancerData %>% select(-FacetSampleCount)
    panCancerData = panCancerData %>% ungroup()
    
    perCancerData = data %>% group_by(Chromosome,ChrIndex,CancerType,CancerSampleCount,CancerFacetSampleCount,FacetField) %>% summarise(Percent=sum(SumField)/first(CancerFacetSampleCount))
    perCancerData = perCancerData %>% ungroup()
    perCancerData = perCancerData %>% select(-CancerFacetSampleCount)
  }
  else
  {
    panCancerData = data %>% group_by(Chromosome,ChrIndex,FacetField) %>% summarise(Percent=sum(SumField)/sampleCount) %>% mutate(CancerType='PanCancer',
                                                                                                                                  CancerSampleCount=sampleCount)
    
    perCancerData = data %>% group_by(Chromosome,ChrIndex,CancerType,CancerSampleCount,FacetField) %>% summarise(Percent=sum(SumField)/first(CancerSampleCount))
  }

  allData = rbind(panCancerData,perCancerData)
  
  allData = allData %>% mutate(CancerTypeLabel=paste(CancerType,' (',CancerSampleCount,')',sep=''))
  
  plot = (ggplot(allData, aes(x=reorder(Chromosome,as.numeric(ChrIndex)),y=reorder(CancerTypeLabel,CancerSampleCount))) 
          + geom_tile(aes(fill = Percent),colour="white",stat = "identity",position="identity") 
          + geom_text(aes(label=ifelse(Percent==0,'',round(Percent*100, 0))))
          + scale_fill_gradient(low="white",high="steelblue")
          + labs(x='Chromosome',y='Cancer Type')
          + ggtitle(title))
  
  if(facetStr != "")
  {
    plot = plot + facet_wrap(~FacetField)
  }

  return (plot)
}

# see below for plot debug and validation



#####
## Key plots for chromosomal and arm level gain and loss

print(chrCancerTypePlot(cnChrData %>% mutate(SumField=WholeChrGain),'Whole Chromosome Gain by Chromosome'))

print(chrCancerTypePlot(cnChrData %>% mutate(SumField=WholeChrLoss),'Whole Chromosome Loss By Chromosome'))

print(chrCancerTypePlot(cnChrData %>% filter(!ShortArm) %>% mutate(SumField=PvsQGain,FacetField=PCentroGain),'P Arm Gain by Gain in Centromere','PCentroGain'))
print(chrCancerTypePlot(cnChrData %>% filter(!ShortArm) %>% mutate(SumField=QvsPGain,FacetField=QCentroGain),'Q Arm Gain by Gain in Centromere','QCentroGain'))

print(chrCancerTypePlot(cnChrData %>% filter(!ShortArm) %>% mutate(SumField=PvsQLoss,FacetField=QCentroGain),'P Arm Loss by Loss in Centromere','QCentroGain'))
print(chrCancerTypePlot(cnChrData %>% filter(!ShortArm) %>% mutate(SumField=QvsPLoss,FacetField=PCentroGain),'Q Arm Loss by Loss in Centromere','PCentroGain'))

print(chrCancerTypePlot(cnChrData %>% filter(!ShortArm) %>% mutate(SumField=QvsPLoss|PvsQGain),'P Arm Gain / Q Arm Loss'))
print(chrCancerTypePlot(cnChrData %>% filter(!ShortArm) %>% mutate(SumField=PvsQLoss|QvsPGain),'Q Arm Gain / P Arm Loss'))

print(chrCancerTypePlot(cnChrData %>% filter(!ShortArm) %>%
                          mutate(SumField=HasCentroCNChg,FacetField=CentroCNChgBucket,Filtered=CentroCNChgBucket<1),
                        'Centromere Copy Number Change by Scale of Change','CentroCNChgBucket',F,T))

print(chrCancerTypePlot(cnChrData %>% filter(!ShortArm) %>% mutate(SumField=HasCentroCNChg,FacetField=TP53Driver),
                        'Centromere Copy Number Change by Samples with TP53 Driver','TP53Driver',T))


#####
## NET pancreatic and low-ploidy samples

sampleCancerAndSubTypes = read.csv('~/data/hpc_sample_cancer_types.csv')
netPancreasSamples = sampleCancerAndSubTypes %>% filter(CancerType=='NET'&CancerSubtype=='Pancreatic')

View(sampleCancerAndSubTypes %>% filter(is.na(CancerSubtype)))

highChrLohSamples = cnChrData %>% group_by(SampleId) %>% summarise(ChrLohLoss=sum(WholeChrLoh=='true')) %>% filter(ChrLohLoss>=10)
View(highChrLohSamples)

altTertGenes = read.csv('~/logs/sample_alt_tert_genes.csv')
men1Drivers = altTertGenes %>% filter(Gene=='MEN1')
daxxDrivers = altTertGenes %>% filter(Gene=='DAXX')
atrxDrivers = altTertGenes %>% filter(Gene=='ATRX')
View(men1Drivers)

specChrData = cnChrData %>% filter(SampleId %in% netPancreasSamples$SampleId | SampleId %in% highChrLohSamples$SampleId) %>% 
  mutate(LowPloidy=Ploidy<=1.5,
         HighChrLoh=SampleId %in% highChrLohSamples$SampleId,
         NetPancreatic=SampleId %in% netPancreasSamples$SampleId,
         HasMen1Driver=SampleId %in% men1Drivers$SampleId,
         HasDaxxDriver=SampleId %in% daxxDrivers$SampleId,
         HasAtrxDriver=SampleId %in% atrxDrivers$SampleId)

specChrData = specChrData %>% mutate(CopyNumber=pmax((AvgCnP+AvgCnQ)/2,0),
                                     CnRatio=CopyNumber,
                                     SampleLabel=sprintf("%s %s Ploidy=%.1f",SampleId,CancerType,Ploidy))

View(specChrData)

View(specChrData %>% filter(NetPancreatic) %>% group_by(SampleId,HasMen1Driver,HasDaxxDriver,HasAtrxDriver) %>% 
       summarise(ChrLohLoss=sum(WholeChrLoh=='true')))


print(ggplot(specChrData %>% filter(LowPloidy&!NetPancreatic), aes(x=reorder(Chromosome,as.numeric(ChrIndex)),y=reorder(SampleLabel,Ploidy))) 
        + geom_tile(aes(fill=CnRatio),colour="white",stat = "identity",position="identity") 
        + geom_text(aes(label=sprintf("%.1f",CopyNumber)))
        + scale_fill_gradient(low="white",high="steelblue")
        + labs(x='Chromosome',y='Sample')
        + ggtitle("Low Ploidy Samples (NET-Pancreatic excluded) Avg Copy Number per Chromosome"))

print(ggplot(specChrData %>% filter(NetPancreatic&HighChrLoh), aes(x=reorder(Chromosome,as.numeric(ChrIndex)),y=reorder(SampleLabel,Ploidy))) 
        + geom_tile(aes(fill=CopyNumber/Ploidy),colour="white",stat = "identity",position="identity") 
        + geom_text(aes(label=sprintf("%.1f",CopyNumber/Ploidy)))
        + scale_fill_gradient(low="white",high="steelblue")
        + labs(x='Chromosome',y='Sample')
        + ggtitle("NET-Pancreatic Samples with 10+ LOH Chromsome Events, Avg Copy Number per Chromosome"))


#####
## Chr 7 and 10 investigations

# CHECK: is chr7 extra copy highly associated with DM in EGFR for CNS
chr7GainSamples = cnChrData %>% filter(CancerType=='CNS'&Chromosome==7&WholeChrGain)
View(chr7GainSamples)

egfrDrivers = ampSummary %>% filter(Gene=='EGFR'&CancerType=='CNS')
View(egfrDrivers)

ptenSamples = read.csv('~/logs/pten_samples.csv')
nrow(ptenSamples)

View(cnChrData %>% filter(CancerType=='CNS'&(Chromosome==10|Chromosome==7)) %>% group_by(SampleId) %>%
       summarise(HasLoss=sum(WholeChrLoss&Chromosome==10),
                 HasGain=sum(WholeChrGain&Chromosome==7)) %>% group_by(HasLoss,HasGain) %>% count)

cnsSamples = sampleCancerTypes %>% filter(CancerType=='CNS') %>% select(SampleId)

calc_fisher_et(cnsSamples,
               chr7GainSamples %>% select(SampleId),
               cnChrData %>% filter(CancerType=='CNS'&Chromosome==10&WholeChrLoss) %>% select(SampleId),
               'Chr7Gain','Chr10Loss')

calc_fisher_et(cnsSamples,
               chr7GainSamples %>% select(SampleId),
               egfrDrivers %>% select(SampleId),
               'Chr7Gain','EGFR')

calc_fisher_et(cnsSamples,
               cnChrData %>% filter(CancerType=='CNS'&Chromosome==10&WholeChrLoss) %>% select(SampleId),
               ptenSamples %>% filter(SampleId %in% cnsSamples$SampleId),
               'Chr10Loss','PTEN')

#####
## NET subtypes for chr 18 loss
chr18NET = cnChrData %>% filter(Chromosome==18&CancerType=='NET')
chr18NET = merge(chr18NET,sampleCancerAndSubTypes %>% select(SampleId,CancerSubtype),by='SampleId',all.x=T)
View(chr18NET)
print(chrCancerTypePlot(chr18NET %>% mutate(CancerType=CancerSubtype,SumField=WholeChrLoss),'NET Sub-types: Whole Chromosome 18 Loss'))

View(chr18NET %>% select(CancerSubtype,SampleId,WholeChrLoss,everything()))

smadSamples = read.csv('~/logs/smad4_samples.csv')
View(smadSamples)

netSamples = cnChrData %>% filter(CancerType=='NET')
netSamples = merge(netSamples,sampleCancerAndSubTypes %>% select(SampleId,CancerSubtype),by='SampleId',all.x=T)
netSmallIntData = netSamples %>% filter(CancerSubtype=='Small Intestinal')
netSmallIntSamples = netSmallIntData %>% group_by(SampleId) %>% count
View(netSmallIntData)
View(netSmallIntSamples)

calc_fisher_et(netSmallIntSamples,
               netSmallIntData %>% filter(Chromosome==18&WholeChrLoss) %>% select(SampleId),
               smadSamples %>% filter(SampleId %in% netSmallIntSamples$SampleId),
               'Chr18Loss','SMAD4')

crData = cnChrData %>% filter(CancerType=='Colon/Rectum')
crSamples = crData %>% group_by(SampleId) %>% count
View(crSamples)

calc_fisher_et(crSamples,
               crData %>% filter(Chromosome==18&WholeChrLoss) %>% select(SampleId),
               smadSamples %>% filter(SampleId %in% crSamples$SampleId),
               'Chr18Loss','SMAD4')



wgdSamples = read.csv('~/logs/wgd_samples.csv')
nrow(chr7GainSamples %>% filter(SampleId %in% wgdSamples$SampleId))
nrow(sampleCancerTypes %>% filter(CancerType=='CNS') %>% filter(!(SampleId %in% chr7GainSamples$SampleId)) 
       %>% filter(SampleId %in% wgdSamples$SampleId))


#####
## TP53 vs centromeric CN change
tp53Samples = driverGenes %>% filter(Gene=='TP53')
View(tp53Samples)
View(cnChrData)
highChrLohSamples = cnChrData %>% group_by(SampleId) %>% summarise(ChrLohLoss=sum(WholeChrLoh=='true')) %>% filter(ChrLohLoss>=10)

cnChrData = cnChrData %>% mutate(TP53Driver=SampleId %in% tp53Samples$SampleId)

# set SampleCount to be the number of samples with a TP53 driver by cancer type
tp53CancerCounts = cnChrData %>% group_by(CancerType,SampleId) %>% summarise(TP53Driver=first(TP53Driver)) %>% 
  group_by(CancerType) %>% summarise(Tp53SampleCount=sum(TP53Driver))

View(tp53CancerCounts)
View(cnChrData %>% group_by(CancerType,SampleId) %>% summarise(TP53Driver=first(TP53Driver)))
tp53CnChrData = merge(cnChrData,tp53CancerCounts,by='CancerType',all.x=T)

View(tp53CnChrData %>% filter(CancerType=='Skin'&Chromosome==1) %>% group_by(TP53Driver,HasCentroCNChg) %>% count)

print(chrCancerTypePlot(tp53CnChrData %>% filter(!ShortArm) %>% mutate(SampleCount=Tp53SampleCount,SumField=HasCentroCNChg,FacetField=TP53Driver),
                        'Centromere Copy Number Change by Samples with TP53 Driver','TP53Driver'))


View(cnChrData %>% filter(HasCentroCNChg) %>% select(CentroCnP,CentroCnQ,CentroCNChg,CentroCNChgBucket,everything()))
View(cnChrData %>% filter(HasCentroCNChg) %>% group_by(pmin(CentroCNChgBucket,3)) %>% count)


#####
## Co-occurrence for centromeric gain/loss between chromosomes by cancer-type

ctStr = 'Breast'
tmpCtData = cnChrData %>% filter(CancerType==ctStr)
sampleCount = nrow(sampleCancerTypes %>% filter(CancerType==ctStr))
print(sampleCount)
chr1=5
chr1ArmGain = tmpCtData %>% filter(Chromosome==chr1&QCentroGain) %>% select(SampleId)
nrow(chr1ArmGain)
chr2=7
chr2ArmGain = tmpCtData %>% filter(Chromosome==chr2&PCentroGain) %>% select(SampleId)
nrow(chr2ArmGain)

# Breast,5PGain,9PGain
calc_fisher_et(tmpCtData %>% group_by(SampleId) %>% count,chr1ArmGain,chr2ArmGain,'Chr1ArmGain','Chr2ArmGain',log=T)

chromosomes = c('1','2','3','4','5','6','7','8','9','10','11','12','16','17','18','19','20','X')
cancerTypes = c('Skin') 
cancerTypes = unique(sampleCancerTypes$CancerType)
for(cancerType in cancerTypes)
{
  print(sprintf('processing %s',cancerType))
  cnCtData = cnChrData %>% filter(CancerType==cancerType)
  cnCtSamples = cnCtData %>% group_by(SampleId) %>% count %>% ungroup()
  
  for(i in 1:length(chromosomes))
  {
    chr1 = chromosomes[i]
    # print(sprintf('processing %s',chr1))
    
    for(j in i+1:length(chromosomes))
    {
      if(j > length(chromosomes))
        break
      
      chr2 = chromosomes[j]
      
      # test P and Q arm gains independently
      
      # print(sprintf('%s and %s: prob=%f', label1, label2, prob))
      
      # gain on both
      label1 = sprintf('%sPGain',chr1)
      label2 = sprintf('%sPGain',chr2)
      prob = calc_fisher_et(cnCtSamples,cnCtData %>% filter(Chromosome==chr1&PCentroGain),cnCtData %>% filter(Chromosome==chr2&PCentroGain),label1,label2,F)

      if(prob < 0.001)
      {
        print(sprintf('%s,%s,%s,%f', cancerType, label1, label2, prob))
      }
      
      # gain on P, gain on Q
      label1 = sprintf('%sPGain',chr1)
      label2 = sprintf('%sQGain',chr2)
      prob = calc_fisher_et(cnCtSamples,cnCtData %>% filter(Chromosome==chr1&PCentroGain),cnCtData %>% filter(Chromosome==chr2&QCentroGain),label1,label2,F)
      
      if(prob < 0.001)
      {
        print(sprintf('%s,%s,%s,%f', cancerType, label1, label2, prob))
      }

      label1 = sprintf('%sQGain',chr1)
      label2 = sprintf('%sPGain',chr2)
      prob = calc_fisher_et(cnCtSamples,cnCtData %>% filter(Chromosome==chr1&QCentroGain),cnCtData %>% filter(Chromosome==chr2&PCentroGain),label1,label2,F)
      
      if(prob < 0.001)
      {
        print(sprintf('%s,%s,%s,%f', cancerType, label1, label2, prob))
      }

      label1 = sprintf('%sQGain',chr1)
      label2 = sprintf('%sQGain',chr2)
      prob = calc_fisher_et(cnCtSamples,cnCtData %>% filter(Chromosome==chr1&QCentroGain),cnCtData %>% filter(Chromosome==chr2&QCentroGain),label1,label2,F)
      
      if(prob < 0.001)
      {
        print(sprintf('%s,%s,%s,%f', cancerType, label1, label2, prob))
      }
      
    }
  }

}
    


# write out data for LINX stats to process using its 3-variable co-occurrence routine
View(tsgSampleData %>% select(SampleId,Gene,LohType,CancerType) %>% arrange(Gene,SampleId))

View(tsgSampleData)
write.csv(tsgSampleData %>% select(SampleId,Gene,CancerType,LohType) %>% arrange(Gene,SampleId), 
          '~/data/sv/coc_3var_gene_loh_data.csv', row.names = F, quote = F)

cocGeneLohCancerResults = read.csv('~/data/sv/SVA_STATS_DATA.csv')
View(cocGeneLohCancerResults)





## Elevated CN segments

# at each successive 3M base junction, record the absolute value of CN change vs the previous 3M base junction
# then plot / analyse frequency of each (bucket of) CN change
# plot rates for centromeric vs non-centromeric segments

cnChgSegments = read.csv('~/logs/CN_CHANGE_SEGMENTS.csv')
View(cnChgSegments)
View(cnChgSegments %>% group_by(Chromosome) %>% summarise(Segments=sum(Frequency)))
cnChgSegments = cnChgSegments %>% filter(SampleId %in% sampleCancerTypes$SampleId)
cnChgSegments = cnChgSegments %>% mutate(CnChgBucket=ifelse(CnChange<=1,CnChange,ifelse(CnChange>=4,'4+','2-3')))

View(cnChgSegments %>% group_by(CnChgBucket) %>% summarise(Count=n(),FreqCount=sum(Frequency)))

View(cnChgSegments %>% filter(CnChgBucket!='0') %>% group_by(SampleId) %>% summarise(Count=n(),FreqCount=sum(Frequency)))


sampleCount=nrow(cnChrData %>% group_by(SampleId) %>% count)
print(sampleCount)
cohortSegmentCount=sampleCount*1e3
print(cohortSegmentCount)
#View(cnChgSegments %>% group_by(CnChgBucket) %>% summarise(Count=n(),Rate=round(Count/cohortSegmentCount,6)))

cnChgSegmentSummary = cnChgSegments %>% group_by(CnChgBucket) %>% summarise(Count=sum(Frequency)) %>% ungroup() 
View(cnChgSegmentSummary)

# tally up counts and a rate for centromere CN changes
centroCnChgSegments = cnChrData %>% filter(!ShortArm) %>% select(SampleId,Chromosome,CentroCNChg) %>% 
  mutate(CnChgBucket=round(CentroCNChg),
         CnChgBucket=ifelse(CnChgBucket<=1,CnChgBucket,ifelse(CnChgBucket>=4,'4+','2-3')))

# View(centroCnChgSegments %>% group_by(CnChgBucket) %>% summarise(Count=n(),Rate=round(Count/centromereCount,6)))

centroCnChgSegmentSummary = centroCnChgSegments %>% group_by(CnChgBucket) %>% summarise(Count=n()) %>% 
  ungroup() %>% mutate(Type='Centromere')

centromereCount=sum(centroCnChgSegmentSummary$Count)
centroCnChgSegmentSummary = centroCnChgSegmentSummary %>% mutate(Rate=round(Count/centromereCount,6))
View(centroCnChgSegmentSummary)

# remove CN changes in the centromere
cnChgSegmentSummary = merge(cnChgSegmentSummary,centroCnChgSegmentSummary %>% select(CnChgBucket,CentroCount=Count),all.x=T)
# View(cnChgSegmentSummary)
cnChgSegmentSummary = cnChgSegmentSummary %>% mutate(Count=pmax(Count-CentroCount,0))
cnChgSegmentSummary = within(cnChgSegmentSummary,rm(CentroCount))
cnChgSegmentSummary = cnChgSegmentSummary %>% mutate(Type='All')
View(cnChgSegmentSummary)

# now calculate a rate per CN change bucket
cnChgSegmentSummary = cnChgSegmentSummary %>% mutate(Rate=round(Count/sum(cnChgSegmentSummary$Count),4))

allCnChgSegmentSummary = rbind(cnChgSegmentSummary,centroCnChgSegmentSummary)
View(allCnChgSegmentSummary)

print(ggplot(allCnChgSegmentSummary %>% filter(CnChgBucket!='0'), aes(x=Type, y=Rate, fill=CnChgBucket))
      + geom_bar(stat = "identity", colour = "black")
      + scale_fill_manual(values = tcColours)
      + labs(x='', y='Rate of Copy Number Change')
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10)))

View(allCnChgSegmentSummary %>% filter(CnChgBucket!='0') %>% select(CnChgBucket,Type,Rate) %>% spread(Type,Rate))

print(ggplot(allCnChgSegmentSummary %>% filter(CnChgBucket!='0') %>% select(CnChgBucket,Type,Rate) %>% spread(Type,Rate))
      + geom_bar(aes(x='All',y=All,fill=CnChgBucket), stat = "identity", colour = "black")
      + geom_bar(aes(x='Centromere',y=Centromere,fill=CnChgBucket), stat = "identity", colour = "black")
      + scale_fill_manual(values = tcColours)
      + labs(x='', y='Rate of Copy Number Change')
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10)))

#####
## repeated by chromosome
cnChgSegments = cnChgSegments %>% filter(SampleId %in% sampleCancerTypes$SampleId)
# cnChgSegments = cnChgSegments %>% mutate(CnChgBucket=pmax(pmin(round(CnChange),5),1))
cnChgSegments = cnChgSegments %>% mutate(CnChgBucket=ifelse(CnChange<=1,CnChange,ifelse(CnChange>=4,'4+','2-3')))

cnChgSegmentSummary = cnChgSegments %>% group_by(Chromosome,CnChgBucket) %>% summarise(Count=sum(Frequency)) %>% ungroup() 
View(cnChgSegmentSummary)

# tally up counts and a rate for centromere CN changes
centroCnChgSegments = cnChrData %>% filter(!ShortArm) %>% select(SampleId,Chromosome,CentroCNChg) %>% 
  mutate(CnChgBucket=round(CentroCNChg),
         CnChgBucket=ifelse(CnChgBucket<=1,CnChgBucket,ifelse(CnChgBucket>=4,'4+','2-3')))

centroCnChgSegmentSummary = centroCnChgSegments %>% group_by(Chromosome,CnChgBucket) %>% summarise(Count=n()) %>% 
  ungroup() %>% mutate(Type='Centromere')

centromereCounts=centroCnChgSegmentSummary %>% group_by(Chromosome) %>% summarise(ChrCount=sum(Count))
View(centromereCounts)
print(centromereCount)
print(sampleCount*18)
centroCnChgSegmentSummary = centroCnChgSegmentSummary %>% mutate(Rate=round(Count/sampleCount,6))
View(centroCnChgSegmentSummary)
View(centroCnChgSegmentSummary %>% filter(CnChgBucket!='0'))

# remove CN changes in the centromere
cnChgSegmentSummary = merge(cnChgSegmentSummary,centroCnChgSegmentSummary %>% select(CnChgBucket,Chromosome,CentroCount=Count),by=c('Chromosome','CnChgBucket'),all.x=T)
# View(cnChgSegmentSummary)
cnChgSegmentSummary = cnChgSegmentSummary %>% mutate(Count=pmax(Count-CentroCount,0))
cnChgSegmentSummary = within(cnChgSegmentSummary,rm(CentroCount))
cnChgSegmentSummary = cnChgSegmentSummary %>% mutate(Type='All')
View(cnChgSegmentSummary)

# now calculate a rate per CN change bucket
chromosomeSegCounts = chrLengths %>% mutate(SegCount=round(Length/3e6))
View(chromosomeSegCounts)

cnChgSegmentSummary = merge(cnChgSegmentSummary,chromosomeSegCounts %>% select(Chromosome,SegCount),by='Chromosome',all.x=T)
cnChgSegmentSummary = cnChgSegmentSummary %>% mutate(Rate=round(Count/(SegCount*sampleCount),4))
View(cnChgSegmentSummary)
cnChgSegmentSummary[is.na(cnChgSegmentSummary)] = 0
View(cnChgSegmentSummary %>% filter(CnChgBucket!='0'))
cnChgSegmentSummary = within(cnChgSegmentSummary,rm(SegCount))

allCnChgSegmentSummary = rbind(cnChgSegmentSummary,centroCnChgSegmentSummary)
View(allCnChgSegmentSummary)

View(allCnChgSegmentSummary %>% filter(CnChgBucket!='0'&!(Chromosome %in% c(13,14,15,21,22,'Y'))) %>% select(CnChgBucket,Chromosome,Type,Rate) %>% spread(Type,Rate))

print(ggplot(allCnChgSegmentSummary %>% filter(CnChgBucket!='0'&!(Chromosome %in% c(13,14,15,21,22,'Y'))) %>% 
               select(CnChgBucket,Chromosome,Type,Rate) %>% spread(Type,Rate))
      + geom_bar(aes(x='All',y=All,fill=CnChgBucket), stat = "identity", colour = "black")
      + geom_bar(aes(x='Centromere',y=Centromere,fill=CnChgBucket), stat = "identity", colour = "black")
      + scale_fill_manual(values = tcColours)
      + facet_wrap(~Chromosome)
      + labs(x='', y='Rate of Copy Number Change')
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + ggtitle('Rate of Copy Number Change by Chromomsome in Centromeric vs non-Centromeric regions'))



## DEBUG

print(chrCancerTypePlot(cnChrData %>% filter(!ShortArm) %>% mutate(SumField=PvsQGain,FacetField=PCentroGain),'P Arm Gain by Gain in Centromere','PCentroGain'))

View(cnChrData %>% filter(CancerType=='Skin'&Chromosome==6) %>% group_by(PvsQGain,PCentroGain) %>% count)


print(chrCancerTypePlot(cnChrData %>% mutate(SumField=WholeChrTCGain&WholeChrAvgGain),'Whole Chromosome Gain by Chromosome'))


print(chrCancerTypePlot(cnChrData %>% filter(!ShortArm) %>% mutate(SumField=PvsQGain,FacetField=PCentroGain),
                        'P Arm Gain by Gain in Centromere','PCentroGain'))

print(chrCancerTypePlot(cnChrData %>% filter(!ShortArm) %>% mutate(SumField=PvsQGain,FacetField=PCentroGain),
                        'P Arm Gain by Gain in Centromere','PCentroGain',T))


nrow(cnChrData %>% filter(CancerType=='Lung'&Chromosome==5&HasCentroCNChg&TP53Driver)) # 180
nrow(cnChrData %>% filter(CancerType=='Lung'&Chromosome==5&TP53Driver)) # 286
nrow(cnChrData %>% filter(CancerType=='Lung'&Chromosome==5&HasCentroCNChg&!TP53Driver)) # 35
nrow(cnChrData %>% filter(CancerType=='Lung'&Chromosome==5&!TP53Driver)) # 93
print(180/286)
print(35/93)


sampleCount = 3782
cancerSampleCounts = cnChrData %>% group_by(CancerType,SampleId) %>% count %>% group_by(CancerType) %>% summarise(CancerSampleCount=n())
View(cancerSampleCounts)
tmp = merge(cnChrData,cancerSampleCounts,by='CancerType',all.x=T)
tmp = tmp %>% mutate(SumField=WholeChrTCGain&WholeChrAvgGain)

# non-faceted version
panCancerData = tmp %>% group_by(Chromosome,ChrIndex) %>% summarise(Percent=sum(SumField)/n()) %>% mutate(CancerType='PanCancer',
                                                                                                          CancerSampleCount=sampleCount)

perCancerData = tmp %>% group_by(Chromosome,ChrIndex,CancerType,SampleCount) %>% summarise(Percent=sum(SumField)/n())

# faceted version 
tmp = tmp %>% mutate(TP53Driver=SampleId %in% tp53Samples$SampleId,
                     FacetField=TP53Driver,
                     SumField=HasCentroCNChg)

cancerSampleCounts = tmp %>% group_by(CancerType,SampleId,FacetField) %>% count %>% group_by(CancerType,FacetField) %>% summarise(CancerSampleCount=n())
View(cancerSampleCounts)
View(tmp)
tmp = within(tmp,rm(CancerSampleCount))
tmp = merge(tmp,cancerSampleCounts,by=c('CancerType','FacetField'),all.x=T)
View(tmp %>% filter(CancerType=='Skin'&Chromosome==1))

facetSampleCounts = tmp %>% group_by(FacetField,SampleId) %>% count %>% group_by(FacetField) %>% summarise(FacetSampleCount=n())

panCancerData = tmp %>% group_by(Chromosome,ChrIndex,FacetField) %>% summarise(Percent=sum(SumField))
panCancerData = merge(panCancerData,facetSampleCounts,by='FacetField',all.x=T)
View(panCancerData)

panCancerData = panCancerData %>% mutate(Percent=Percent/FacetSampleCount,
                                         CancerType='PanCancer',
                                         CancerSampleCount=sampleCount)

View(panCancerData)

perCancerData = tmp %>% group_by(Chromosome,ChrIndex,CancerType,CancerSampleCount,FacetField) %>% summarise(Percent=sum(SumField)/first(CancerSampleCount))


allData = rbind(panCancerData %>% select(-FacetSampleCount) %>% ungroup(),perCancerData %>% ungroup())
View(allData)

allData = allData %>% mutate(CancerTypeLabel=paste(CancerType,' (',CancerSampleCount,')',sep=''))
View(allData)
View(allData %>% filter(CancerType=='Skin'))

print(chrCancerTypePlot(cnChrData %>% mutate(SumField=WholeChrTCGain&WholeChrAvgGain),'Whole Chromosome Gain by Chromosome'))
# print(chrCancerTypePlot(cnChrData %>% mutate(SumField=PvsQGain,FacetField=PCentroGain),'P-vs-Q Arm Gain by Gain in Centromere','PCentroGain'))







# chrCtPArmGain = cnChrData %>% group_by(Chromosome,ChrIndex,CancerType,SampleCount) %>% summarise(Percent=sum(PvsQGain)/n())
print(chrCancerTypePlot(cnChrData %>% mutate(SumField=PvsQGain),'P-Arm Gain by Chromosome'))

# chrCtPArmGain2 = cnChrData %>% group_by(Chromosome,ChrIndex,CancerType,SampleCount,PCentroGain) %>% summarise(Percent=sum(PvsQGain)/first(SampleCount))

print(chrCancerTypePlot(cnChrData %>% filter(!ShortArm) %>% mutate(SumField=PArmLoss),'P-Arm Loss by Chromosome'))


# chrCtQArmGain = cnChrData %>% group_by(Chromosome,ChrIndex,CancerType,SampleCount) %>% summarise(Percent=sum(QvsPGain)/n())
print(chrCancerTypePlot(cnChrData %>% mutate(SumField=),'Q vs P Arm Gain'))

print(chrCancerTypePlot(cnChrData %>% filter(!ShortArm) %>% mutate(SumField=QvsPGain,FacetField=QCentroGain),'Q vs P Arm Gain by Gain in Centromere','QCentroGain'))

chrCtQArmLoss = cnChrData %>% group_by(Chromosome,ChrIndex,CancerType,SampleCount) %>% summarise(Percent=sum(QArmLoss)/n())
print(chrCancerTypePlot(cnChrData %>% filter(!ShortArm) %>% mutate(SumField=),'Q-Arm Loss'))

chrCtPArmLoh = cnChrData %>% group_by(Chromosome,ChrIndex,CancerType,SampleCount) %>% summarise(Percent=sum(LohP=='true')/n())
print(chrCancerTypePlot(cnChrData %>% mutate(SumField=),'P-Arm LOH'))

chrCtWholeChrLoh = cnChrData %>% group_by(Chromosome,ChrIndex,CancerType,SampleCount,WholeChrLoss) %>% summarise(Percent=sum(WholeChrLoh=='true')/first(SampleCount))
print(chrCancerTypePlot(cnChrData %>% mutate(SumField=),'Whole Chr LOH by Total Loss','WholeChrLoss'))
View(chrCtPWholeChrLoh)

chrCtPArmLoh = cnChrData %>% group_by(Chromosome,ChrIndex,CancerType,SampleCount,QCentroGain) %>% summarise(Percent=sum(LohP=='true')/first(SampleCount))
print(chrCancerTypePlot(cnChrData %>% mutate(SumField=),'P-Arm LOH','QCentroGain'))


View(drivers %>% filter(CancerType=='Skin'&Chromosome==6) %>% group_by(Gene,DriverType) %>% count)


# LOH data to check copy number neutrality
lohEvents = read.csv('~/data/sv/CN_LOH_EVENTS.csv')
View(lohEvents)
centroLohs = lohEvents %>% filter(SegStart=='CENTROMERE'|SegEnd=='CENTROMERE')
View(centroLohs)
rm(lohEvents)

chrCentromeres = read.csv('~/hmf/repos/hmftools/hmf-common/src/main/resources/refgenome/hg19_centromere.tsv',sep='\t',header=F)
View(chrCentromeres)
colnames(chrCentromeres) = c('Chromosome','Centromere')
chrLengths = read.csv('~/hmf/repos/hmftools/hmf-common/src/main/resources/refgenome/hg19_len.tsv',sep='\t',header=F)
colnames(chrLengths) = c('Chromosome','Length')
chrArmLengths = merge(chrCentromeres,chrLengths,by='Chromosome',all=T)
chrArmLengths = chrArmLengths %>% mutate(PArmLength=Centromere,QArmLength=Length-Centromere)
View(chrLengths)
# convert to arm
chrArms = chrArmLengths %>% select(Chromosome,PArmLength,QArmLength) %>% gather('Arm','ArmLength',2:3)
chrArms = chrArms %>% mutate(Arm=ifelse(Arm=='PArmLength','P','Q'))
chrArms = chrArms %>% mutate(ChrArm=paste(Chromosome,Arm,sep='_'))
View(chrArms)


View(chrArmLengths)
write.csv(chrArmLengths,'~/data/sv/chr_cento_lengths.csv',row.names = F, quote = F)
write.csv(chrArms,'~/data/sv/chr_arm_lengths.csv',row.names = F, quote = F)

centroLohs = merge(centroLohs,chrArms %>% select(ChrArm,ArmLength),by='ChrArm',all.x=T)

centroLohs = centroLohs %>% mutate(LengthBucket=2**round(log(Length,2)),
                                   LengthPerc=pmin(Length/ArmLength,1),
                                   LengthPercBucket=round(LengthPerc,1),
                                   Arm=ifelse(SegStart=='CENTROMERE','Q','P'),
                                   ChrArm=paste(Chromosome,Arm,sep='_'),
                                   WholeArm=(SegStart=='TELOMERE'|SegEnd=='TELOMERE'))


View(centroLohs)
View(centroLohs %>% group_by(Chromosome,Arm,WholeArm) %>% count() %>% spread(WholeArm,n))
View(centroLohs %>% group_by(Chromosome,Arm,LengthBucket) %>% count())
plot_length_facetted(centroLohs,'!WholeArm','LengthBucket,ChrArm','LengthBucket','ChrArm','LOH Length by ChrArm',F)
plot_length_facetted(centroLohs,'!WholeArm','LengthPercBucket,ChrArm','LengthPercBucket','ChrArm','LOH Length % by ChrArm',F)

View(lohEvents %>% group_by(SegStart,SegEnd) %>% count())
View(lohEvents %>% group_by(LohType,TypeStart,TypeEnd,SegStart,SegEnd) %>% count())

lohEvents = lohEvents %>% mutate(HighCn=CnStart>5|CnEnd>5,
                                 LengthBucket=2**round(log(Length,2)),
                                 TypeStart=ifelse(SegStart=='CENTROMERE','CENTRO',ifelse(SegStart=='TELOMERE','TELO','SV')),
                                 TypeEnd=ifelse(SegEnd=='CENTROMERE','CENTRO',ifelse(SegEnd=='TELOMERE','TELO','SV')),
                                 LohType=ifelse(TypeStart=='TELO'&TypeEnd=='TELO','LOH_CHR',
                                         ifelse((TypeStart=='TELO'&TypeEnd=='CENTRO')|(TypeEnd=='TELO'&TypeStart=='CENTRO'),'LOH_ARM',
                                         ifelse(TypeStart=='TELO'|TypeEnd=='TELO','LOH_SV_TELO',ifelse(TypeStart=='CENTRO'|TypeEnd=='CENTRO','LOH_SV_CENTRO','LOH_FOCAL')))))

View(lohEvents %>% group_by(LohType,CnNeutral,LengthBucket) %>% count() %>% spread(LengthBucket,n))
View(lohEvents %>% group_by(CnNeutral) %>%  count())


# incoporate SV data

# percentrimeric = HSATII

centroTeloLohs = lohEvents %>% filter((SegStart=='CENTROMERE'&SegEnd=='SGL')|(SegEnd=='CENTROMERE'&SegStart=='SGL')
                                      |(SegStart=='TELOMERE'&SegEnd=='SGL')|(SegEnd=='TELOMERE'&SegStart=='SGL'))

centroTeloLohs = centroTeloLohs %>% mutate(SvId=ifelse(StartSV!=-1,StartSV,EndSV))
centroTeloLohs = merge(centroTeloLohs,svData %>% filter(Type=='SGL') %>% select(SampleId,SvId=Id,RepeatClass,RepeatType,ResolvedType),by=c('SampleId','SvId'),all.x=T)
View(centroTeloLohs)
View(centroTeloLohs %>% group_by(SegStart,SegEnd,RepeatClass,RepeatType,ResolvedType) %>% count)



# elevCNS = read.csv('~/logs/CN_ELEV_SEGMENTS.csv')
View(cnChgSegments)

elevCNS = merge(elevCNS,centroLengths %>% select(Chromosome,Length),by='Chromosome',all.x=T)
elevCNS = merge(elevCNS,sampleCancerAndSubTypes,by='SampleId',all.x=T)
elevCNS = elevCNS %>% mutate(AdjElevSegCount=ElevSegCount/(Length/3e6))

View(elevCNS %>% group_by(Chromosome) %>% summarise(Avg=round(mean(AdjElevSegCount),2),
                                                    Median=round(median(AdjElevSegCount),2),
                                                    Max=round(max(AdjElevSegCount),2),
                                                    ElevationRate=ElevSegCount/1e3))


