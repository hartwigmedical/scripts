library(grid)
library(gridExtra)
library(ggplot2)
library(cowplot)

################
## Driver Analysis

getAllColours<-function()
{
  colours = c("cornsilk3","seagreen3", "tomato3","hotpink", "darkorange", "thistle2", "steelblue2", "darkgreen", "indianred", "honeydew2",
              "turquoise3", "lightpink2", "goldenrod2", "darkslateblue", "yellowgreen", "wheat2", "violetred2", "ivory3", "coral1", "springgreen2")
  
  return (colours)
}

getSubsetColours<-function(requiredCount)
{
  allColours = getAllColours()
  colours = c()
  count = 1
  
  for(i in 1:requiredCount)
  {
    colours[count] = allColours[count]
    count = count + 1
  }
  
  return (colours)
}

# print(getSubsetColours(3))

getFactorColours<-function(subsetFactors,allFactors)
{
  factorColours = getSubsetColours(length(allFactors))

  # return required colours
  colours = c()
  count = 1
  
  for(factor in subsetFactors)
  {
    for(i in 1:length(allFactors))
    {
      if(allFactors[i] == factor)
      {
        colours[count] = factorColours[i]
        count = count + 1
        break
      }
    }
  }
  
  return (colours)
}

driverFactorByGroupPlot<-function(data,factorAllTypes,groupLabel,splitLabel,title='',countLimit=0,asPercent=T)
{
  # calculate sample count
  sampleCounts = data %>% group_by(Group,SampleId) %>% count() %>% group_by(Group) %>% summarise(SampleCount=n())

  data = data %>% arrange(Factor)
  subsetFactors = unique(data$Factor)
  # subsetFactors = data %>% group_by(Factor) %>% count() %>% arrange(Factor) %>% select(Factor)
  
  # print(subsetFactors)
  # print(factorAllTypes)
  
  colours = getAllColours()

  if(length(factorAllTypes) > 0)
    colours = getFactorColours(subsetFactors,factorAllTypes)
  
  # print(colours)

  # group data and convert to percentage
  
  dataSummary = data %>% group_by(Group,Factor) %>% summarise(Count=sum(Count)) %>% arrange(Factor)
  dataTotals = data %>% group_by(Group) %>% summarise(TotalCount=sum(Count))
  dataSummary = merge(dataSummary,dataTotals,by='Group',all.x=T)
  
  dataSummary = merge(dataSummary,sampleCounts,by='Group',all.x=T)
  
  if(countLimit > 0)
    dataSummary = dataSummary %>% filter(TotalCount>=countLimit)
  
  dataSummary = dataSummary %>% mutate(GroupLabel=paste(Group,' (',SampleCount,')',sep=''),
                                       Percent=round(Count/TotalCount,6))
  
  if(asPercent)
  {
    plot = (ggplot(dataSummary, aes(x=reorder(GroupLabel,-SampleCount), y=Percent, fill=Factor))
             + geom_bar(stat = "identity", colour = "black")
             + labs(x="", y=paste(splitLabel,'Percent',sep=' ')))
  }
  else
  {
    plot = (ggplot(dataSummary, aes(x=reorder(GroupLabel,-TotalCount), y=Count, fill=Factor))
             + geom_bar(stat = "identity", colour = "black")
             + labs(x="", y=paste(splitLabel,'Count',sep=' ')))
  }
  
  plot = (plot + scale_fill_manual(values = colours)
           + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
           + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
           + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9))
           + ggtitle(title))
  
  return (plot)  
}

## Common Data prep
drivers = read.csv('~/data/sv/drivers/LNX_DRIVERS.csv')
clusters = read.csv('~/data/sv/drivers/LNX_CLUSTERS.csv')

#load('~/data/hmf_cohort_may_2019.RData')
#sampleCancerTypes = highestPurityCohort %>% select(SampleId=sampleId,CancerType=cancerType)

sampleCancerTypes = load_cancer_types('~/data/hpc_sample_cancer_types.csv',F,10)
nrow(sampleCancerTypes)

# LINX HPC cohort of 3784 samples, corresponds to those in cancer types file, not de-duped
hpcSamples = read.csv('~/data/sv/hpc_sample_ids.csv')
nrow(hpcSamples)

View(hpcSamples %>% filter(SampleId %in% sampleCancerTypes$SampleId))

drivers = drivers %>% filter(SampleId %in% sampleCancerTypes$SampleId)
drivers = merge(drivers,sampleCancerTypes,by='SampleId',all.x=T)

# get data from clusters
drivers = merge(drivers,clusters %>% select(SampleId,ClusterId,Foldbacks,Annotations,SuperType),by=c('SampleId','ClusterId'),all.x=T)

# filter incompletes
drivers = drivers %>% filter(!(EventType=='GAIN'&ClusterId==-1))

drivers = drivers %>% mutate(DriverCategory=ifelse(DriverType=='AMP','AMP',ifelse(grepl('LOH',EventType),'LOH','DEL')))
View(drivers)
View(drivers %>% group_by(DriverCategory,Category,DriverType,EventType) %>% count())

# Resolved Types
driverResolvedTypes = drivers %>% filter(ClusterId>=0) %>% group_by(CancerType,SampleId,Gene,DriverCategory,Chromosome,Arm) %>% 
  summarise(Count=n(),
            DlCount=sum(Likelihood),
            SuperType=first(SuperType),
            ResolvedType=first(ResolvedType)) %>% ungroup()

driverResolvedTypes = driverResolvedTypes %>% mutate(ResolvedType=ifelse(Count==1,as.character(ResolvedType),'MULTIPLE'),
                                                     GeneLoc=paste(Gene,':',Chromosome,Arm,sep=''))

View(driverResolvedTypes)
View(driverResolvedTypes %>% group_by(SuperType,ResolvedType) %>% count())

## LOH EVENTS

driverLOHReport<-function(lohEvents,resolvedTypes,fileName)
{
  outputFile = paste("~/logs/r_output/pdfs/", fileName, ".pdf", sep = "")
  print(paste("writing output to file: ", outputFile, sep=''))
  
  # DATA OUTPUT TO PDF
  pdf(file=outputFile, height = 14, width = 20)
  
  par(mar=c(1,1,1,1))
  
  print(sprintf("plotting LOH analysis"))

  # LOH reports
  # lohTypes = c('LOH_ARM','LOH_SV_CENTRO','LOH_CHR','LOH_FOCAL','LOH_SV_TELO')
  lohEvents = lohEvents %>% arrange(LohType)
  lohTypes = unique(lohEvents$LohType)

  # events by cancer type
  panCancer = lohEvents %>% mutate(CancerType='PanCancer')
  panCancer = rbind(panCancer,lohEvents)
  
  ctPlot = driverFactorByGroupPlot(panCancer %>% select(Group=CancerType,SampleId,Factor=LohType,Count=DlCount),lohTypes,
                                'Cancer Type','LOH Type','LOH Type by Cancer Type')
  
  # top X genes
  genePlot = driverFactorByGroupPlot(lohEvents %>% select(Group=GeneLoc,SampleId,Factor=LohType,Count=DlCount),lohTypes,
                                'Gene','LOH Type','LOH Type by Gene',50)

  delTypePlot = driverFactorByGroupPlot(lohEvents %>% select(Group=TsgType,SampleId,Factor=LohType,Count=DlCount),lohTypes,
                                     'TSG Type','LOH Type','LOH Type by TSG Type')

  grid.arrange(ctPlot, genePlot, delTypePlot,ncol=1, nrow=3, newpage = TRUE)
  
  # top X genes by cancer type  
  topGeneCount = 10
  geneEventCounts = lohEvents %>% group_by(Gene,SampleId) %>% count() %>% group_by(Gene) %>% summarise(SampleCount=n()) %>% arrange(-SampleCount)
  topGenes = head(geneEventCounts,topGeneCount)
  
  print(sprintf("plotting top %d genes", topGeneCount))
  
  topGenesList = unique(topGenes$Gene)
  topGenesCount = length(topGenesList)
  
  plotsPerPage = 2
  geneIndex = 1

  for(i in 1:ceiling(topGenesCount/2))
  {
    gene1 = topGenesList[geneIndex]
    plot1 = driverFactorByGroupPlot(lohEvents %>% filter(Gene==gene1) %>% 
                                    select(Group=CancerType,SampleId,Factor=LohType,Count=DlCount),lohTypes,
                                  'Cancer Type','LOH Type',paste(gene1,': LOH Type by Cancer Type',sep=''),0,F)

    if(geneIndex+1 > topGenesCount)
    {
      grid.arrange(plot1, ncol=1, nrow=1, newpage = TRUE)
      break
    }

    gene2 = topGenesList[geneIndex+1]
    plot2 = driverFactorByGroupPlot(lohEvents %>% filter(Gene==gene2) %>% 
                                      select(Group=CancerType,SampleId,Factor=LohType,Count=DlCount),lohTypes,
                                    'Cancer Type','LOH Type',paste(gene2,': LOH Type by Cancer Type',sep=''),0,F)
      
    grid.arrange(plot1, plot2, ncol=1, nrow=2, newpage = TRUE)
    geneIndex = geneIndex + 2
  }
  
  print("plotting resolved types")

  # resolved type plots
  # keep types for top X resolved types and merge the rest into 'Other'
  resolvedTypeCounts = resolvedTypes %>% group_by(ResolvedType) %>% count() %>% arrange(-n)
  topResolvedTypes = head(resolvedTypeCounts,8)

  resolvedTypes = resolvedTypes %>% mutate(ClusterType=ifelse(ResolvedType %in% topResolvedTypes$ResolvedType,ResolvedType,'OTHER'))
  
  plotRtCt = driverFactorByGroupPlot(resolvedTypes %>% select(Group=CancerType,SampleId,Factor=ClusterType,Count=DlCount),c(),
                                'Cancer Type','Resolved Type','LOH Resolved Type by Cancer Type')
  
  plotRtGene = driverFactorByGroupPlot(resolvedTypes %>% select(Group=GeneLoc,SampleId,Factor=ClusterType,Count=DlCount),c(),
                                'Gene','LOH Type','LOH Type by Gene',50)

  grid.arrange(plotRtCt, plotRtGene, ncol=1, nrow=2, newpage = TRUE)
  
  print(sprintf("PDF complete"))

  dev.off()
}

View(drivers %>% filter(DriverCategory=='LOH'))

lohSummary = drivers %>% filter(DriverCategory=='LOH') %>% group_by(CancerType,SampleId,Gene,DriverCategory,Chromosome,Arm) %>% 
  summarise(Count=n(),
            DlCount=sum(Likelihood),
            LohChr=sum(EventType=='LOH_CHR'),
            LohArm=sum(EventType=='LOH_ARM'),
            LohSvCentro=sum(EventType=='LOH_SV_CENTRO'),
            LohSvTelo=sum(EventType=='LOH_SV_TELO'),
            LohFocal=sum(EventType=='LOH')) %>% ungroup()

lohSummary = lohSummary %>% mutate(LohType=ifelse(LohChr>0,'LOH_CHR',ifelse(LohArm>0,'LOH_ARM',
                                           ifelse(LohSvCentro>0,'LOH_SV_CENTRO',ifelse(LohSvTelo>0,'LOH_SV_TELO','LOH_FOCAL')))),
                                   GeneLoc=paste(Gene,':',Chromosome,Arm,sep=''))

View(lohSummary %>% filter(DriverCategory=='LOH') %>% group_by(LohChr,LohArm,LohFocal,LohSvCentro,LohSvTelo) %>% count())
View(lohSummary)

lohResolvedTypes = driverResolvedTypes %>% filter(DriverCategory=='LOH')

fragileSiteGenes = read.csv('~/data/sv/drivers/fragile_site_genes.csv')
topTsgGenes = read.csv('~/data/sv/drivers/KnownTSGs.csv')
# View(topTsgGenes %>% filter(Gene %in% fragileSiteGenes$Gene))

lohSummary = lohSummary %>% mutate(TsgType=ifelse(Gene %in% fragileSiteGenes$GeneName,'FRAGILE_STE',ifelse(Gene %in% topTsgGenes$GeneName,'TOP_TSG','OTHER')))

View(lohSummary %>% group_by(TsgType,LohType) %>% count() %>% spread(TsgType,n))

driverLOHReport(lohSummary,lohResolvedTypes,'driver_loh')

lohTypes = c('LOH_ARM','LOH_CHR','LOH_FOCAL','LOH_SV_CENTRO','LOH_SV_TELO')



## DELETIONS
delResolvedTypes = driverResolvedTypes %>% filter(DriverCategory=='DEL'&!is.na(ResolvedType))
# View(delResolvedTypes)

# keep types for top X resolved types and merge the rest into 'Other'
delResolvedTypeCounts = delResolvedTypes %>% group_by(ResolvedType) %>% count() %>% arrange(-n)
topDelResolvedTypes = head(delResolvedTypeCounts,8)
View(topDelResolvedTypes)

delResolvedTypes = delResolvedTypes %>% mutate(ClusterType=ifelse(ResolvedType %in% topDelResolvedTypes$ResolvedType,ResolvedType,'OTHER'),
                                               TsgType=ifelse(Gene %in% fragileSiteGenes$Gene,'FRAGILE_SITE',ifelse(Gene %in% topTsgGenes$Gene,'TOP_TSG','OTHER')))
View(delResolvedTypes)
# View(delResolvedTypes %>% group_by(SampleId,Gene) %>% count())

driverDELReport<-function(resolvedTypes,fileName)
{
  outputFile = paste("~/logs/r_output/pdfs/", fileName, ".pdf", sep = "")
  print(paste("writing output to file: ", outputFile, sep=''))
  
  # DATA OUTPUT TO PDF
  pdf(file=outputFile, height = 14, width = 20)
  
  par(mar=c(1,1,1,1))
  
  print(sprintf("producing DEL driver report"))
  
  resolvedTypes = resolvedTypes %>% arrange(ClusterType)
  clusterTypes = unique(resolvedTypes$ClusterType)
  
  panCancer = resolvedTypes %>% mutate(CancerType='PanCancer')
  panCancer = rbind(panCancer,resolvedTypes)
  
  ctPlot = driverFactorByGroupPlot(panCancer %>% select(Group=CancerType,SampleId,Factor=ClusterType,Count=DlCount),clusterTypes,
                                'Cancer Type','Resolved Type','DEL Resolved Type by Cancer Type')
  
  genePlot = driverFactorByGroupPlot(delResolvedTypes %>% select(Group=GeneLoc,SampleId,Factor=ClusterType,Count=DlCount),clusterTypes,
                                'Gene','DEL Type','DEL Type by Gene',50)
  
  delTypePlot = driverFactorByGroupPlot(delResolvedTypes %>% select(Group=TsgType,SampleId,Factor=ClusterType,Count=DlCount),clusterTypes,
                                        'TSG Type','DEL Type','DEL Type by TSG Type')
  
  grid.arrange(ctPlot, genePlot, delTypePlot, ncol=1, nrow=3, newpage = TRUE)
  
  # top X genes by cancer type  
  topGeneCount = 10
  geneEventCounts = delResolvedTypes %>% group_by(Gene,SampleId) %>% count() %>% group_by(Gene) %>% summarise(SampleCount=n()) %>% arrange(-SampleCount)
  topGenes = head(geneEventCounts,topGeneCount)
  
  print(sprintf("plotting top %d genes", topGeneCount))
  
  topGenesList = unique(topGenes$Gene)
  topGenesCount = length(topGenesList)
  
  plotsPerPage = 2
  geneIndex = 1
  
  for(i in 1:ceiling(topGenesCount/2))
  {
    gene1 = topGenesList[geneIndex]
    
    plot1 = driverFactorByGroupPlot(delResolvedTypes %>% filter(Gene==gene1) %>% 
                                      select(Group=CancerType,SampleId,Factor=ClusterType,Count=DlCount),clusterTypes,
                                     'Cancer Type','Resolved Type',paste(gene1,': DEL Resolved Type by Cancer Type'),0,F)
    
    if(geneIndex+1 > topGenesCount)
    {
      grid.arrange(plot1, ncol=1, nrow=1, newpage = TRUE)
      break
    }
    
    gene2 = topGenesList[geneIndex+1]

    plot2 = driverFactorByGroupPlot(delResolvedTypes %>% filter(Gene==gene2) %>% 
                                      select(Group=CancerType,SampleId,Factor=ClusterType,Count=DlCount),clusterTypes,
                                    'Cancer Type','Resolved Type',paste(gene2,': DEL Resolved Type by Cancer Type'),0,F)
    
    grid.arrange(plot1, plot2, ncol=1, nrow=2, newpage = TRUE)
    geneIndex = geneIndex + 2
  }

  print(sprintf("PDF complete"))
  
  dev.off()
}

driverDELReport(delResolvedTypes,'driver_del')


## AMPLIFICATIONS
View(ampDrivers)

ampDrivers = drivers %>% filter(DriverCategory=='AMP') %>%
  mutate(AmpType=ifelse(DriverCategory!='AMP','NONE',
                 ifelse(EventType=='GAIN_ARM','GAIN_ARM',ifelse(EventType=='GAIN_CHR','GAIN_CHR',
                 ifelse(grepl('BFB',Annotations),'BFB',ifelse(grepl('DM',Annotations),'DM',ifelse(ResolvedType=='COMPLEX','COMPLEX',
                 ifelse(ResolvedType=='DUP','DUP','COMPLEX'))))))),
         CnGainBucket=ifelse(GeneMinCN>128,128,2**round(log(GeneMinCN/SamplePloidy,2))),
         GeneLoc=paste(Gene,':',Chromosome,Arm,sep=''))

# nrow(ampDrivers)

ampDrivers = ampDrivers %>% filter(!(!is.na(ClusterId)&ClusterId>=0&is.na(Foldbacks)))

View(ampDrivers %>% group_by(CnGainBucket=2**round(log(CNGain,2))) %>% count())
View(ampDrivers %>% group_by(CnGainBucket=2**round(log(GeneMinCN/SamplePloidy,2))) %>% count())

ampSummary = ampDrivers %>% group_by(CancerType,SampleId,Gene,GeneLoc,DriverCategory,CnGainBucket) %>% 
  summarise(ClusterCount=n(),
            Count=1,
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

View(ampSummary)
#View(ampSummary %>% filter(is.na(AmpMainType)))
#View(ampDrivers %>% filter(is.na(AmpType)))


driverAMPReport<-function(summaryData,fileName)
{
  outputFile = paste("~/logs/r_output/pdfs/", fileName, ".pdf", sep = "")
  print(paste("writing output to file: ", outputFile, sep=''))
  
  # DATA OUTPUT TO PDF
  pdf(file=outputFile, height = 14, width = 20)
  
  par(mar=c(1,1,1,1))
  
  print(sprintf("producing AMP driver report"))
  
  summaryData = summaryData %>% arrange(AmpMainType)
  ampTypes = unique(summaryData$AmpMainType)
  
  panCancer = summaryData %>% mutate(CancerType='PanCancer')
  panCancer = rbind(panCancer,summaryData)

  ctPlot = driverFactorByGroupPlot(panCancer %>% select(Group=CancerType,SampleId,Factor=AmpMainType,Count),ampTypes,
                                   'Cancer Type','Resolved Type','AMP Cluster Type by Cancer Type')
  
  genePlot = driverFactorByGroupPlot(summaryData %>% select(Group=GeneLoc,SampleId,Factor=AmpMainType,Count),ampTypes,
                                     'Gene','AMP Type','AMP Cluster Type Type by Gene',50)
  
  cnGainPlot = driverFactorByGroupPlot(summaryData %>% select(Group=CnGainBucket,SampleId,Factor=AmpMainType,Count),ampTypes,
                                     'Gene','AMP Type','AMP Cluster Type Type by CN Gain',0)
  
  grid.arrange(ctPlot, genePlot, cnGainPlot, ncol=1, nrow=3, newpage = TRUE)

  # top X genes by cancer type  
  topGeneCount = 20
  geneEventCounts = summaryData %>% group_by(Gene,SampleId) %>% count() %>% group_by(Gene) %>% summarise(SampleCount=n()) %>% arrange(-SampleCount)
  topGenes = head(geneEventCounts,topGeneCount)
  
  print(sprintf("plotting top %d genes", topGeneCount))
  
  topGenesList = unique(topGenes$Gene)
  topGenesCount = length(topGenesList)
  
  plotsPerPage = 2
  geneIndex = 1
  
  for(i in 1:ceiling(topGenesCount/2))
  {
    gene1 = topGenesList[geneIndex]
    
    plot1 = driverFactorByGroupPlot(summaryData %>% filter(Gene==gene1) %>% 
                                      select(Group=CancerType,SampleId,Factor=AmpMainType,Count),ampTypes,
                                    'Cancer Type','Cluster Type',paste(gene1,': AMP Cluster Type by Cancer Type'),0,F)
    
    if(geneIndex+1 > topGenesCount)
    {
      grid.arrange(plot1, ncol=1, nrow=1, newpage = TRUE)
      break
    }
    
    gene2 = topGenesList[geneIndex+1]
    
    plot2 = driverFactorByGroupPlot(summaryData %>% filter(Gene==gene2) %>% 
                                      select(Group=CancerType,SampleId,Factor=AmpMainType,Count),ampTypes,
                                    'Cancer Type','Cluster Type',paste(gene2,': AMP Cluster Type by Cancer Type'),0,F)
    
    grid.arrange(plot1, plot2, ncol=1, nrow=2, newpage = TRUE)
    geneIndex = geneIndex + 2
  }
  
  print(sprintf("PDF complete"))
  
  dev.off()
}

driverAMPReport(ampSummary,'driver_amp_new')

write.csv(ampSummary,'~/data/sv/drivers/amp_summary.csv',row.names = F, quote = F)

## DEBUG
summaryData = ampSummary %>% arrange(AmpMainType)
ampTypes = unique(summaryData$AmpMainType)
View(ampTypes)

panCancer = summaryData %>% mutate(CancerType='PanCancer')
panCancer = rbind(panCancer,summaryData)

ctPlot = driverFactorByGroupPlot(panCancer %>% select(Group=CancerType,SampleId,Factor=AmpMainType,Count),ampTypes,
                                 'Cancer Type','Resolved Type','AMP Cluster Type by Cancer Type')

print(driverFactorByGroupPlot(summaryData %>% select(Group=CnGainBucket,SampleId,Factor=AmpMainType,Count),ampTypes,
                                     'Gene','AMP Type','AMP Cluster Type Type by CN Gain',0))



## FUSIONS

fusions = read.csv('~/data/sv/fusions/LNX_FUSIONS.csv')
fusions = fusions %>% filter(Reportable=='true')

fusions = fusions %>% filter(SampleId %in% sampleCancerTypes$SampleId)
fusions = merge(fusions,sampleCancerTypes,by='SampleId',all.x=T)
nrow(fusions)

fusions = fusions %>% mutate(GenePair=paste(GeneNameUp,GeneNameDown,sep='_'),Count=1)

View(fusions)
View(fusions %>% group_by(ResolvedType) %>% count())
View(fusions %>% group_by(CancerType) %>% count())
View(fusions %>% group_by(GenePair) %>% count())

# get data from clusters
# drivers = merge(drivers,clusters %>% select(SampleId,ClusterId,Foldbacks,Annotations,SuperType),by=c('SampleId','ClusterId'),all.x=T)

fusionReport<-function(fusions,fileName)
{
  outputFile = paste("~/logs/r_output/pdfs/", fileName, ".pdf", sep = "")
  print(paste("writing output to file: ", outputFile, sep=''))
  
  # DATA OUTPUT TO PDF
  pdf(file=outputFile, height = 14, width = 20)
  
  par(mar=c(1,1,1,1))
  
  print(sprintf("producing fusion report"))
  
  fusions = fusions %>% arrange(ResolvedType)
  resolvedTypes = unique(fusions$ResolvedType)
  
  panCancer = fusions %>% mutate(CancerType='PanCancer')
  panCancer = rbind(panCancer,fusions)
  
  ctPlot = driverFactorByGroupPlot(panCancer %>% select(Group=CancerType,SampleId,Factor=ResolvedType,Count),resolvedTypes,
                                   'Cancer Type','Resolved Type','Fusion Resolved Type by Cancer Type')
  
  genePlot = driverFactorByGroupPlot(fusions %>% select(Group=GenePair,SampleId,Factor=ResolvedType,Count),resolvedTypes,
                                     'Gene-Pair','Resolved Type','Fusion Type by Gene-Pair',5)
  
  grid.arrange(ctPlot, genePlot, ncol=1, nrow=2, newpage = TRUE)
  
  # top gene-pairs by cancer type  
  topGeneCount = 10
  geneEventCounts = fusions %>% group_by(GenePair,SampleId) %>% count() %>% group_by(GenePair) %>% summarise(SampleCount=n()) %>% arrange(-SampleCount)
  topGenes = head(geneEventCounts,topGeneCount)
  
  print(sprintf("plotting top %d gene pairs", topGeneCount))
  
  topGenesList = unique(topGenes$GenePair)
  topGenesCount = length(topGenesList)
  
  plotsPerPage = 2
  geneIndex = 1
  
  for(i in 1:ceiling(topGenesCount/2))
  {
    genePair1 = topGenesList[geneIndex]
    
    plot1 = driverFactorByGroupPlot(fusions %>% filter(GenePair==genePair1) %>% 
                                      select(Group=CancerType,SampleId,Factor=ResolvedType,Count),resolvedTypes,
                                    'Cancer Type','Resolved Type',paste(genePair1,': Fusion Resolved Type by Cancer Type'),0,F)
    
    if(geneIndex+1 > topGenesCount)
    {
      grid.arrange(plot1, ncol=1, nrow=1, newpage = TRUE)
      break
    }
    
    genePair2 = topGenesList[geneIndex+1]
    
    plot2 = driverFactorByGroupPlot(fusions %>% filter(GenePair==genePair2) %>% 
                                      select(Group=CancerType,SampleId,Factor=ResolvedType,Count),resolvedTypes,
                                    'Cancer Type','Resolved Type',paste(genePair2,': Fusion Resolved Type by Cancer Type'),0,F)
    
    grid.arrange(plot1, plot2, ncol=1, nrow=2, newpage = TRUE)
    geneIndex = geneIndex + 2
  }
  
  print(sprintf("PDF complete"))
  
  dev.off()
}

fusionReport(fusions,'driver_fusion')

View(fusions)



#####
## Multiple Biopsies
mbAmpDrivers = ampDrivers %>% filter(!is.na(PatientId))
colnames(mbAmpDrivers)

nrow(mbAmpDrivers)

mbAmpDrivers = mbAmpDrivers %>% group_by(PatientId,SampleId,Gene) %>% summarise(GeneMinCN=first(GeneMinCN),
                                                                                SamplePloidy=first(SamplePloidy),
                                                                                ClusterId=first(ClusterId),
                                                                                ClusterCount=first(ClusterCount),
                                                                                ResolvedType=first(ResolvedType))

View(mbAmpDrivers %>% group_by(PatientId,Gene) %>% summarise(SampleCount=n()))
View(mbAmpDrivers %>% group_by(PatientId,Gene) %>% summarise(SampleCount=n()) %>% group_by(Gene,SampleCount) %>% count %>% spread(SampleCount,n,fill=0))

View(mbAmpDrivers %>% filter(Gene=='ERBB2'))

View(doubleBiopsy)


                                                                                                                                                                                    

## DEBUG ONLY

## VIOLIN PLOTS

ggplot(merge(sampleList,svData%>% filter(ResolvedType=='LINE') %>% group_by(SampleId) %>% tally(),by='SampleId',all.x=T) %>% replace_na(list(n=0)) %>% group_by(CancerType)  , aes(CancerType, n)) + 
  geom_violin(scale="area",fill="#6baed6")



View(drivers %>% filter(DriverCategory=='LOH'&Gene=='TP53'))
View(driverResolvedTypes %>% filter(DriverCategory=='LOH'&Gene=='TP53'))

nrow(svData)

infLohEvents = driverResolvedTypes %>% filter(DriverCategory=='LOH'&Gene=='TP53'&ResolvedType=='INF')
nrow(driverResolvedTypes %>% filter(DriverCategory=='LOH'&Gene=='TP53'))
nrow(infLohEvents)
View(infLohEvents)

infLohEvents = drivers %>% filter(DriverCategory=='LOH'&Gene=='TP53'&ResolvedType=='INF')
infLohEvents = infLohEvents %>% mutate(InfPosition=ifelse(SvIdStart!=-1,SvPosStart,SvPosEnd),
                                       InfPosBucket=1e4*round(InfPosition/1e4))

temp = infLohEvents %>% filter(InfPosBucket==18930000) %>% .$SampleId
View(svData %>% filter((ChrStart==17&PosStart<23e6)|(ChrEnd==17&PosEnd<23e6),SampleId %in% temp))
View(infLohEvents %>% group_by(InfPosBucket) %>% count)


View(infLohEvents %>% filter(InfPosition>213e5&InfPosition<218e5) %>% mutate(SvId=ifelse(SvIdStart!=-1,SvIdStart,SvIdEnd)))
write.csv(infLohEvents %>% filter(InfPosition>213e5&InfPosition<218e5) %>% mutate(SvId=ifelse(SvIdStart!=-1,SvIdStart,SvIdEnd)) %>% select(SampleId),
          '~/logs/inf_chr17.csv', row.names = F, quote = F)


tp53LohEvents = drivers %>% filter(DriverCategory=='LOH'&Gene=='TP53'&EventType=='LOH_SV_TELO')
tp53LohEvents = tp53LohEvents %>% mutate(IsINF=(ResolvedType=='INF'),
                                         SvId=ifelse(SvIdStart!=-1,SvIdStart,SvIdEnd),
                                         SvPosition=ifelse(SvIdStart!=-1,SvPosStart,SvPosEnd),
                                         SvPosBucket=1e3*round(SvPosition/1e3))

tp53LohEvents = merge(tp53LohEvents,svData %>% select(SampleId,Id,Type),by.x=c('SampleId','SvId'),by.y=c('SampleId','Id'),all.x=T)

View(tp53LohEvents %>% group_by(SvPosBucket,ResolvedType) %>% count %>% spread(ResolvedType,n))
View(tp53LohEvents %>% group_by(SvPosBucket,Type) %>% count %>% spread(Type,n))
View(tp53LohEvents %>% filter(SvPosBucket=='18900000'))
print(ggplot(lohGeneCtSummary %>% filter(Gene=='APC'), aes(x=reorder(CancerType,-EventCount), y=Count, fill = LohType))
      + geom_bar(stat = "identity", colour = "black")
      + labs(x = "", y = "Event Percent")
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))


