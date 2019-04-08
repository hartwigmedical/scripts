library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringi)


annotate_fusions<-function(fusionData)
{
  fusionData$SameSV = (fusionData$SvIdUp==fusionData$SvIdDown)
  
  # chaining info
  fusionData = fusionData %>% separate(ChainInfo,c('ChainId','ChainLinks','ChainLength','ValidTraversal'),sep = ';')
  fusionData$ChainLength = as.numeric(fusionData$ChainLength)
  fusionData$ChainLinks = as.numeric(fusionData$ChainLinks)
  fusionData$InChain = (fusionData$ChainId>=0)
  
  # chain & cluster validity
  fusionData = fusionData %>% separate(OverlapUp,c('FacingBEsUp','FacingClusterBEsUp','TotalBEsUp','FacingDistanceUp','DisruptedExonsUp','TerminatedUp'),sep = ';')
  fusionData$FacingBEsUp = as.numeric(fusionData$FacingBEsUp)
  fusionData$FacingClusterBEsUp = as.numeric(fusionData$FacingClusterBEsUp)
  fusionData$TotalBEsUp = as.numeric(fusionData$TotalBEsUp)
  fusionData$FacingDistanceUp = as.numeric(fusionData$FacingDistanceUp)
  
  fusionData = fusionData %>% separate(OverlapDown,c('FacingBEsDown','FacingClusterBEsDown','TotalBEsDown','FacingDistanceDown','DisruptedExonsDown','TerminatedDown'),sep = ';')
  fusionData$FacingBEsDown = as.numeric(fusionData$FacingBEsDown)
  fusionData$FacingClusterBEsDown = as.numeric(fusionData$FacingClusterBEsDown)
  fusionData$TotalBEsDown = as.numeric(fusionData$TotalBEsDown)
  fusionData$FacingDistanceDown = as.numeric(fusionData$FacingDistanceDown)
  
  fusionData = (fusionData 
                        %>% mutate(ValidChain=ValidTraversal=='true'&DisruptedExonsUp==0&DisruptedExonsDown==0&TerminatedUp==0&TerminatedDown==0,
                                   NonDisruptedSingle=FacingBEsUp==0&FacingBEsDown==0&DisruptedExonsUp==0&DisruptedExonsDown==0))
  return (fusionData)
}

####
## Single-SV fusions from VariantAnntotator

simpleFusions = read.csv('~/data/sv/fusions/FUSIONS.csv')
simpleFusions = read.csv('~/data/sv/fusions/VAR_ANN_FUSIONS_PROD.csv')
write.csv(reportedSimpleFusions, '~/data/sv/fusions/VAR_ANN_FUSIONS_PROD.csv', row.names = F, quote = F)
nrow(simpleFusions)
reportedSimpleFusions = simpleFusions %>% filter(Reportable=='true')
View(reportedSimpleFusions)
nrow(reportedSimpleFusions)
rm(simpleFusions)


####
## SVA Fusions from single SVs and clustered & chained SVs
svaFusions = read.csv('~/data/sv/fusions/SVA_FUSIONS.csv')
nrow(svaFusions)
reportedSvaFusions = svaFusions %>% filter(Reportable=='true')
View(reportedSvaFusions)
nrow(reportedSvaFusions)


# Annotations
reportedSvaFusions = annotate_fusions(reportedSvaFusions)


View(reportedSvaFusions)
write.csv(reportedSvaFusions, '~/data/sv/fusions/SVA_FUSIONS_ANNOT.csv', row.names = F, quote = F)

# Comparison with VariantAnnotator
reportedSvaSingleSvFusions = reportedSvaFusions %>% filter(SameSV)
View(reportedSvaSingleSvFusions)
nrow(reportedSvaSingleSvFusions)

simpleOverlapVA = merge(reportedSimpleFusions,reportedSvaSingleSvFusions, by=c('SampleId','SvIdUp'),all.x=T)
View(simpleOverlapVA)
nrow(simpleOverlapVA %>% filter(is.na(Reportable.y))) # all are accounted for by the SVA
View(simpleOverlapVA %>% filter(is.na(Reportable.y)))

simpleOverlapSVA = merge(reportedSvaSingleSvFusions,reportedSimpleFusions,by=c('SampleId','SvIdUp'),all.x=T)
View(simpleOverlapSVA)
nrow(simpleOverlapSVA %>% filter(is.na(Reportable.y)))
View(simpleOverlapSVA %>% filter(is.na(Reportable.y)))


# basic numbers
nrow(reportedSvaFusions %>% filter(SameSV)) # 585 single-SV fusions
nrow(reportedSvaFusions %>% filter(ChainLinks==0))
nrow(reportedSvaFusions %>% filter(!SameSV))
nrow(reportedSvaFusions %>% filter(!SameSV&ValidChain))

# View(reportedSvaFusions %>% filter(ChainLinks>0) %>% group_by(SameSV) %>% count())
# View(reportedSvaFusions %>% group_by(SameSV,InChain) %>% count())
# View(reportedSvaFusions %>% group_by(SameSV,InChain,ValidTraversal) %>% count())


# possibly invalid single-SV fusions
nrow(reportedSvaSingleSvFusions %>% filter(!ValidChain|!NonDisruptedSingle))
View(reportedSvaSingleSvFusions %>% filter(!ValidChain|!NonDisruptedSingle))
View(reportedSvaSingleSvFusions %>% filter(!ValidChain|!NonDisruptedSingle) %>% group_by(Clustered=ClusterCount>1) %>% count())
View(reportedSvaSingleSvFusions %>% group_by(Clustered=ClusterCount>1,InChain,IsValid=ValidChain|NonDisruptedSingle) %>% count())

reportedSvaSingleSvFusions$TotalFacingBreakends = reportedSvaSingleSvFusions$FacingBEsUp + reportedSvaSingleSvFusions$FacingBEsDown
reportedSvaSingleSvFusions$FacingBECountBucket = 2**round(log(reportedSvaSingleSvFusions$TotalFacingBreakends,2))
View(reportedSvaSingleSvFusions %>% group_by(Clustered=ClusterCount>1,InChain,FacingBECountBucket) %>% count())
nrow(reportedSvaSingleSvFusions %>% filter(!InChain&!NonDisruptedSingle))


# number of single-SV fusions with possible disruptions from other breakends


# unique valid fusions and their type
sampleFusions = (reportedSvaFusions %>% filter((InChain&ValidChain)|(!InChain&NonDisruptedSingle)) %>% group_by(SampleId,GeneUp,GeneDown) 
     %>% summarise(Count=n(),
                   SimpleSVCount=sum(SameSV&ClusterCount==1),
                   SingleSVUnchainedCount=sum(SameSV&ClusterCount>1&!InChain),
                   SingleSVChainedCount=sum(SameSV&InChain),
                   MultiSVChainedCount=sum(!SameSV&ValidChain))
     %>% mutate(FusionType=ifelse(MultiSVChainedCount==Count,'MultiSV',
                           ifelse(SimpleSVCount==Count,'SimpleSV',
                           ifelse(SingleSVUnchainedCount==Count,'SingleSVUnchained',
                           ifelse(SingleSVChainedCount==Count,'SingleSVChained','Unclear'))))))

View(sampleFusions)
View(sampleFusions %>% filter(GeneUp))

View(sampleFusions %>% group_by(FusionType) %>% count())


# number of valid chained (non-single SV) fusions without a matching single SV fusion for same genes
View(sampleFusions %>% filter(MultiSVChainedCount==Count))
nrow(sampleFusions %>% filter(MultiSVChained==Count)) # 52 valid chained fusions not also found by a single SV






# evaluation of facing breakend data
View(reportedSvaFusions %>% filter(DisruptedExonsUp>0|DisruptedExonsDown>0|TerminatedUp==1|TerminatedDown==1) %>% group_by(SameSV,InChain=ChainId!='') %>% count())
View(reportedSvaFusions %>% filter(DisruptedExonsUp>0|DisruptedExonsDown>0|TerminatedUp==1|TerminatedDown==1))
View(reportedSvaFusions %>% filter(DisruptedExonsUp>0|DisruptedExonsDown>0|TerminatedUp==1|TerminatedDown==1) %>% filter(SameSV&!InChain))
View(reportedSvaFusions %>% filter(DisruptedExonsUp==0&DisruptedExonsDown==0&TerminatedUp==0&TerminatedDown==0) %>% group_by(SameSV,InChain=ChainId!='') %>% count())
View(reportedSvaFusions %>% filter(DisruptedExonsUp==0&DisruptedExonsDown==0&TerminatedUp==0&TerminatedDown==0))



facingDistancesUp = reportedSvaFusions %>% filter(FacingDistanceUp>0)
facingDistancesUp$DistanceBucket = 2**round(log(facingDistancesUp$FacingDistanceUp,2))
facingDistancesUp$DistanceBucket = 10*round(facingDistancesUp$FacingDistanceUp/10)
facingDistancesDown = reportedSvaFusions %>% filter(FacingDistanceDown>0)
facingDistancesDown$DistanceBucket = 2**round(log(facingDistancesDown$FacingDistanceDown,2))
facingDistancesDown$DistanceBucket = 10*round(facingDistancesDown$FacingDistanceDown/10)
facingDistances = rbind(facingDistancesUp,facingDistancesDown)
facingDistances$InChain=facingDistances$ChainId!=''



plot_length_facetted(facingDistances, 'ClusterCount==1', 'DistanceBucket,InChain', 'DistanceBucket', 'InChain', 
                     'Fusion facing breakend distance', logScale=T)

plot_length_facetted(facingDistances, 'DistanceBucket<500', 'DistanceBucket', 'DistanceBucket', '', 
                     'DSB ratio to SV Count, by CN Gain or Loss', logScale=T)


View(reportedNewSimpleFusions %>% filter(FacingDistanceUp>0))

write.csv(reportedSvaFusions, '~/data/sv/SVA_FUSIONS.csv', row.names = F, quote = F)






# DEBUG ONLY
View(reportedSvaFusions %>% filter(SampleId=='CPCT02030426T'))
View(reportedSimpleFusions %>% filter(SampleId=='CPCT02030426T'))

