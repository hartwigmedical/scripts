library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringi)

annotate_fusions<-function(fusionData)
{
  fusionData = fusionData %>% mutate(SameSV = (SvIdUp==SvIdDown))
  
  # chaining info
  fusionData = fusionData %>% separate(ChainInfo,c('ChainId','ChainLinks','ChainLength','ValidTraversal','TraversalAssembled'),sep = ';')
  fusionData$ChainLength = as.numeric(fusionData$ChainLength)
  fusionData$ChainLinks = as.numeric(fusionData$ChainLinks)
  fusionData$InChain = (fusionData$ChainId>=0)
  
  # chain & cluster validity
  fusionData = fusionData %>% separate(OverlapUp,c('FacingBEsUp','AssembledLinksUp','TotalBEsUp','FacingDistanceUp','DisruptedExonsUp','TerminatedUp'),sep = ';')
  fusionData$FacingBEsUp = as.numeric(fusionData$FacingBEsUp)
  fusionData$TotalBEsUp = as.numeric(fusionData$TotalBEsUp)
  fusionData$FacingDistanceUp = as.numeric(fusionData$FacingDistanceUp)
  fusionData$TerminatedUp = !is.na(fusionData$TerminatedUp) & fusionData$TerminatedUp=='true'
  
  fusionData = fusionData %>% separate(OverlapDown,c('FacingBEsDown','AssembledLinksDown','TotalBEsDown','FacingDistanceDown','DisruptedExonsDown','TerminatedDown'),sep = ';')
  fusionData$FacingBEsDown = as.numeric(fusionData$FacingBEsDown)
  fusionData$TotalBEsDown = as.numeric(fusionData$TotalBEsDown)
  fusionData$FacingDistanceDown = as.numeric(fusionData$FacingDistanceDown)
  fusionData$TerminatedDown = !is.na(fusionData$TerminatedDown) & fusionData$TerminatedDown=='true'
  
  fusionData = (fusionData %>% 
                  mutate(ValidChain=ValidTraversal=='true'&DisruptedExonsUp==0&DisruptedExonsDown==0&TerminatedUp=='false'&TerminatedDown=='false',
                         NonDisruptedSingle=FacingBEsUp==0&FacingBEsDown==0&DisruptedExonsUp==0&DisruptedExonsDown==0,
                         BreakendDistUp=ifelse(StrandUp==1,TransStartUp-PosUp,PosUp-TransEndUp),
                         BreakendDistDown=ifelse(StrandUp==1,TransStartDown-PosDown,PosDown-TransEndDown)))

  fusionData[is.na(fusionData)] = 0
  
  return (fusionData)
}



####
## SVA Fusions from single SVs and clustered & chained SVs
svaFusions = read.csv('~/data/sv/fusions/SVA_FUSIONS.csv')
nrow(svaFusions)
reportedSvaFusions = svaFusions %>% filter(Reportable=='true')
# read.csv('~/data/sv/fusions/SVA_FUSIONS.csv')
View(reportedSvaFusionsRaw)
nrow(reportedSvaFusions)
rm(svaFusions)

# Annotations
reportedSvaFusions = annotate_fusions(reportedSvaFusions)

# check for duplicates
View(reportedSvaFusions %>% group_by(SampleId,GeneNameUp,GeneNameDown) %>% count() %>% filter(n>1))

# multiple fusions per sample
View(reportedSvaFusions %>% filter(Clustered) %>% group_by(SampleId,GeneIdUp,GeneIdDown) %>% count() %>% filter(n>1))




#######
# Comparison with previous run
reportedFusionsPrev = read.csv('~/data/sv/fusions/SVA_FUSIONS_CM_OLD.csv')
nrow(reportedFusionsPrev)

fusionComparison = merge(reportedSvaFusions,reportedFusionsPrev, by=c('SampleId','GeneIdUp','GeneIdDown'),all=T)

# new fusions
nrow(fusionComparison %>% filter(is.na(Reportable.y))) 
View(fusionComparison %>% filter(is.na(Reportable.y)))
View(fusionComparison %>% filter(is.na(Reportable.y)) %>% select(SampleId,GeneNameUp.x,GeneNameDown.x,KnownType.x,ClusterId.x,SvIdUp.x,SvIdDown.x))
View(fusionComparison %>% filter(is.na(Reportable.y)) %>% group_by(GeneNameUp.x,GeneNameDown.x,KnownType.x) %>% count())
View(fusionComparison %>% filter(is.na(Reportable.y)) %>% group_by(KnownType.x) %>% count())

# missing old fusions
nrow(fusionComparison %>% filter(is.na(Reportable.x))) 
View(fusionComparison %>% filter(is.na(Reportable.x)))
View(fusionComparison %>% filter(is.na(Reportable.x)) %>% select(SampleId,GeneNameUp.y,GeneNameDown.y,KnownType.y,ClusterId.y,SvIdUp.y,SvIdDown.y))
View(fusionComparison %>% filter(is.na(Reportable.x)) %>% group_by(GeneNameUp.y,GeneNameDown.y,KnownType.y) %>% count())
View(fusionComparison %>% filter(is.na(Reportable.x)) %>% group_by(KnownType.y) %>% count())


## Fusion comparison after disruption for known fusions change
prevFusions = read.csv('~/data/sv/fusions/LNX_FUSIONS.csv.prev')
prevFusions = prevFusions %>% filter(Reportable=='true')
nrow(prevFusions)

newFusionsAll = read.csv('~/data/sv/fusions/LNX_FUSIONS.csv')
newFusions = newFusionsAll %>% filter(Reportable=='true')
nrow(newFusions)


View(fusions %>% filter(RegionTypeUp=='Exonic'&RegionTypeDown=='Exonic') %>% group_by(KnownType) %>% count)
View(fusions %>% filter(RegionTypeUp=='Exonic'&RegionTypeDown=='Exonic') %>% select(ExactBaseUp,ExactBaseDown))

View(fusions %>% filter(RegionTypeUp=='Exonic'&RegionTypeDown=='Exonic'&CodingTypeUp=='5P_UTR'&CodingTypeDown=='5P_UTR'))


View(fusions %>% filter(RegionTypeUp=='Exonic'&RegionTypeDown=='Exonic'&PhaseMatched=='true'&CodingTypeUp=='Coding'&CodingTypeDown=='Coding') %>% 
       filter(!((ExactBaseUp==0&ExactBaseDown==1)|(ExactBaseUp==1&ExactBaseDown==2)|(ExactBaseUp==2&ExactBaseDown==0))) %>%
       select(ExactBaseUp,ExactBaseDown,everything()))




prevFusions = prevFusions %>% mutate(SampleGenePair=paste(SampleId,GeneNameUp,GeneNameDown,sep='_'))
newFusions = newFusions %>% mutate(SampleGenePair=paste(SampleId,GeneNameUp,GeneNameDown,sep='_'))

flcFusions = unfilteredFusions %>% filter(Reportable=='true') %>% mutate(SampleGenePair=paste(SampleId,GeneNameUp,GeneNameDown,sep='_'))

View(newFusions %>% filter(!(SampleGenePair %in% prevFusions$SampleGenePair)) %>% filter(RegionTypeUp=='Exonic'&RegionTypeDown=='Exonic') %>%
       select(SampleId,KnownType,GeneNameUp,GeneNameDown,RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown,TypeUp,SvIdUp,TypeDown,SvIdDown,ProteinsKept,ProteinsLost,everything()))

View(prevFusions %>% filter(!(SampleGenePair %in% newFusions$SampleGenePair)) %>% filter(RegionTypeUp=='Exonic'&RegionTypeDown=='Exonic') %>%
       select(SampleId,KnownType,GeneNameUp,GeneNameDown,RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown,TypeUp,SvIdUp,TypeDown,SvIdDown,ProteinsKept,ProteinsLost,everything()))

rm(newFusions)
rm(prevFusions)




#######
# DISRUPTIONS 

View(highestPurityCohort)

# comparison with prod
linxDisruptions = read.csv('~/data/sv/fusions/LNX_DISRUPTIONS.csv')
nrow(linxDisruptions) # 18265
linxDisruptionSummary = linxDisruptions %>% group_by(SampleId,GeneName) %>% summarise(Count=n(),
                                                                                      Standard=sum(Reportable=='true'),
                                                                                      Excluded=sum(Reportable=='false'))

View(linxDisruptions %>% group_by(ExcludedReason) %>% count())
View(linxDisruptionSummary)
View(linxDisruptions %>% filter(Reportable=='false') %>% group_by(GeneName) %>% count())
View(linxDisruptions %>% filter(Reportable=='false'&GeneName=='PTEN'))

View(linxDisruptions %>% group_by(SampleId,GeneName) %>% summarise(Count=n(),
                                                              Standard=sum(Reportable=='true'),
                                                              SameIntronNoSPA=sum(ExcludedReason=='SameIntronNoSPA'),
                                                              RemoteIntron=sum(ExcludedReason=='IntronicSection')) %>%
       filter(Count==RemoteIntron))


write.csv(linxDisruptions %>% filter(ExcludedReason=='SameIntronNoSPA') %>% group_by(SampleId) %>% count(),
          '~/logs/no_spa_disruption_samples.csv', row.names = F, quote = F)


sameIntronNoSpa = read.csv('~/logs/LNX_DISRUPTIONS.csv')
View(sameIntronNoSpa)
sameIntronNoSpa = sameIntronNoSpa %>% filter(ExcludedReason=='SameIntronNoSPA')
sameIntronNoSpa = sameIntronNoSpa %>% separate(ExtraInfo,c('ChainLinks','ChainLength'),sep='-')
sameIntronNoSpa = sameIntronNoSpa %>% mutate(ChainLinks=as.numeric(as.character(ChainLinks)),
                                             ChainLength=as.numeric(as.character(ChainLength)))
View(sameIntronNoSpa)
View(sameIntronNoSpa %>% group_by(ChainLengthBucket=2**round(log(ChainLength,2))) %>% summarise(Count=n(),
                                                                                                MedLinks=median(ChainLinks),
                                                                                                MaxLinks=max(ChainLinks)))

write.csv(sameIntronNoSpa,'~/logs/chained_non_disruptions.csv', quote = F, row.names = F)

linxDisruptions = linxDisruptions %>% mutate(UndisruptedCNBucket=ifelse(UndisruptedCN>0.1,2**round(log(UndisruptedCN,2)),0))
View(linxDisruptions %>% group_by(UndisruptedCN) %>% count())

oldDisruptions = read.csv('~/data/sv/fusions/LNX_DISRUPTIONS_OLD.csv')
nrow(oldDisruptions) # 18625
View(oldDisruptions) # 17819

oldDisruptionsummary = oldDisruptions %>% group_by(SampleId,GeneName) %>% count()
View(oldDisruptionsummary)

# samples and genes with no disruptions
mergedDisruptions = merge(oldDisruptionsummary,linxDisruptionSummary %>% filter(Standard>0),by=c('SampleId','GeneName'),all=T)

mergedDisruptions = merge(linxDisruptions %>% group_by(SampleId,GeneName) %>% count(),
                          linxDisruptionSummary %>% filter(Standard>0),by=c('SampleId','GeneName'),all=T)

View(mergedDisruptions %>% filter(is.na(Count)))


# unaffected samples & genes
View(mergedDisruptions)


# oldDisruptions = read.csv('~/logs/prod_disruptions.csv')
oldDisruptions = oldDisruptions %>% filter(SampleId %in% highestPurityCohort$sampleId)
nrow(oldDisruptions)
#oldDisruptionsummary = oldDisruptions %>% group_by(SampleId=sampleId,GeneName=gene) %>% count()

# genes entirely dropped from samples


droppedDisruptions = read.csv('~/logs/SVA_DISRUPTIONS.csv')
View(droppedDisruptions %>% filter(Reportable=='false'))
View(droppedDisruptions %>% filter(Reportable=='false') %>% group_by(ExcludedReason) %>% count())
View(droppedDisruptions %>% filter(Reportable=='false') %>% group_by(SampleId) %>% count())
View(droppedDisruptions %>% filter(Reportable=='true') %>% group_by(SampleId) %>% count())

View(droppedDisruptions %>% filter(Reportable=='false'&ExcludedReason!='SimpleSV') %>% group_by(SampleId,GeneName) %>% count())

View(droppedDisruptions %>% filter(Reportable=='true') %>% group_by(SampleId,GeneName) %>% count())

View(droppedDisruptions %>% filter(ExcludedReason!='SimpleSV'&GeneName=='PTEN'))

View(oldDisruptions %>% filter(GeneName=='PTEN'&SampleId %in% droppedDisruptions$SampleId))


View(droppedDisruptions %>% filter(ExcludedReason!='SimpleSV') %>% 
       group_by(GeneName,Reportable) %>% count() %>% spread(Reportable,n))


mergedDisruptions = merge(linxDisruptions %>% select(SampleId,GeneName,SvId),oldDisruptions %>% select(SampleId,GeneName,SvId),by=c('SampleId','GeneName'),all=T)

View(mergedDisruptions)
nrow(mergedDisruptions)
nrow(mergedDisruptions %>% filter(is.na(SvId.x)))
nrow(mergedDisruptions %>% filter(is.na(SvId.y)))
View(mergedDisruptions %>% filter(is.na(SvId.x)))

View(mergedDisruptions %>% filter(is.na(SvId.x)) %>% group_by(SampleId) %>% count())

write.csv(mergedDisruptions %>% filter(is.na(SvId.x)) %>% group_by(SampleId) %>% count() %>% select(SampleId),
          '~/logs/dropped_disruption_sample_ids.csv', row.names = F, quote = F)









load('~/data/hmf_cohort_may_2019.RData')
View(highestPurityCohort) # 3524 has multiple biopsy samples removed



nrow(reportedSvaFusions %>% filter(GeneUp=='TMPRSS2'&GeneDown=='ERG'))
nrow(reportedSimpleFusions %>% filter(GeneUp=='TMPRSS2'&GeneDown=='ERG'))
View(reportedSvaFusions %>% filter(GeneUp=='TMPRSS2'&GeneDown=='ERG') %>% group_by(SameSV) %>% count())
View(reportedSvaFusionsPrev %>% filter(SameSV&GeneUp=='TMPRSS2'&GeneDown=='ERG'))

View(reportedSvaFusions %>% group_by(RegionTypeUp,RegionTypeDown) %>% count())
View(reportedSvaFusions %>% filter(RegionTypeUp=='Intronic'&RegionTypeDown=='Exonic'))



nrow(reportedSvaFusionsPrev %>% filter(Clustered&ValidChain))
nrow(reportedSvaFusions %>% filter(Clustered&ValidChain))
nrow(reportedSvaFusionsPrev %>% filter(Clustered&ValidChain&SameSV))
nrow(reportedSvaFusions %>% filter(Clustered&ValidChain&SameSV))
nrow(reportedSvaFusionsPrev %>% filter(Clustered&SameSV)) # 585
nrow(reportedSvaFusions %>% filter(Clustered&SameSV)) # 538
nrow(reportedSvaFusionsPrev %>% filter(Clustered&ValidChain&!SameSV)) # 74
nrow(reportedSvaFusions %>% filter(Clustered&ValidChain&!SameSV)) # 90

# basic numbers - taken from HPC de-duped on 5/8/2019
nrow(reportedSvaFusions %>% filter(SameSV)) # 508
nrow(reportedSvaFusions %>% filter(SameSV&!InChain)) # 282
nrow(reportedSvaFusions %>% filter(ChainLinks==0))
nrow(reportedSvaFusions %>% filter(!SameSV)) # 82


# Summary: Unique valid fusions and their type
sampleFusions = (reportedSvaFusions %>% filter(Clustered&ValidChain) %>% group_by(SampleId,GeneUp,GeneDown) 
                 %>% summarise(Count=n(),
                               KnownType=first(KnownType),
                               SimpleSVCount=sum(SameSV&ClusterCount==1),
                               SingleSVUnchainedCount=sum(SameSV&ClusterCount>1&!InChain),
                               SingleSVChainedCount=sum(SameSV&InChain),
                               MultiSVChainedCount=sum(Clustered&!SameSV&InChain),
                               UnclusteredCount=sum(!Clustered))
                 %>% mutate(FusionType=ifelse(MultiSVChainedCount==Count,'MultiSV',ifelse(SimpleSVCount==Count,'SimpleSV',
                                       ifelse(SingleSVUnchainedCount==Count,'SingleSVUnchained',
                                       ifelse(SingleSVChainedCount==Count,'SingleSVChained','Unclear'))))))

View(sampleFusions)

View(sampleFusions %>% group_by(FusionType,KnownType) %>% count() %>% spread(KnownType,n))


# Comparison with previous fusions

# comparison with previous Linx run
reportedSvaFusionsPrev = read.csv('~/data/sv/fusions/LINX_FUSIONS_REPORTED.csv')
View(reportedSvaFusionsPrev)

knownNew = reportedSvaFusions %>% filter(KnownType=='Known')
nrow(knownNew)
knownPrev = reportedSvaFusionsPrev %>% filter(KnownType=='Known')
nrow(knownPrev)

knownByPrev = merge(knownPrev,knownNew,by=c('SvIdUp','SvIdDown'),all.x=T)
View(knownByPrev %>% filter(is.na(GeneUp.y)))
knownByNew = merge(knownNew,knownPrev,by=c('SvIdUp','SvIdDown'),all.x=T)
View(knownByNew %>% filter(is.na(GeneUp.y)))

promNew = reportedSvaFusions %>% filter(grepl('Prom',KnownType))
nrow(promNew)
promPrev = reportedSvaFusionsPrev %>% filter(grepl('Prom',KnownType))
nrow(promPrev)

promByPrev = merge(promPrev,promNew,by=c('SvIdUp','SvIdDown'),all.x=T)
View(promByPrev %>% filter(is.na(GeneUp.y)))
promByNew = merge(promNew,promPrev,by=c('SvIdUp','SvIdDown'),all.x=T)
View(promByNew %>% filter(is.na(GeneUp.y)))

promNew = reportedSvaFusions %>% filter(grepl('Prom',KnownType))
nrow(promNew)
promPrev = reportedSvaFusionsPrev %>% filter(grepl('Prom',KnownType))
nrow(promPrev)

promByPrev2 = merge(promPrev,promNew,by=c('SampleId','GeneUp','GeneDown'),all.x=T)
View(promByPrev2 %>% filter(is.na(SvIdUp.y)))
promByNew2 = merge(promNew,promPrev,by=c('SampleId','GeneUp','GeneDown'),all.x=T)
View(promByNew2 %>% filter(is.na(SvIdUp.y)))


# specific sample
specificFusions = annotate_fusions(specificFusions)
View(specificFusions)
View(specificFusions %>% filter(Reportable=='true'))

View(specificFusions %>% filter(GeneUp=='TMPRSS2'&GeneDown=='ERG'
                                &(ValidTraversal=='true'&PhaseMatched=='true'&BiotypeDown!='nonsense_mediated_decay')))




# DRUP vs DNDS TSGs
drupTsgs = read.csv('~/data/drup_genes.csv')
dndsTsgs = read.csv('~/hmf/repos/hmftools/hmf-common/src/main/resources/dnds/DndsDriverLikelihoodTsg.tsv',sep='\t')
View(dndsTsgs)
View(drupTsgs %>% filter(!(Gene %in% dndsTsgs$gene)))




# number of valid chained (non-single SV) fusions without a matching single SV fusion for same genes
View(sampleFusions %>% filter(MultiSVChainedCount==Count))
nrow(sampleFusions %>% filter(MultiSVChained==Count)) # 52 valid chained fusions not also found by a single SV


####
## Known Unclustered Fusions
knownUnclusteredFusions = svaFusions %>% filter(grepl('Unclustered',ResolvedType))
nrow(knownUnclusteredFusions)
knownUnclusteredFusions = annotate_fusions(knownUnclusteredFusions)
View(knownUnclusteredFusions)

knownUnclusteredFusionsSummary = knownUnclusteredFusions %>% group_by(SampleId,GeneUp,GeneDown) %>% 
  summarise(Count=n(),PhaseMatched=sum(PhaseMatched=='true'),Reported=sum(Reportable=='true'))

View(knownUnclusteredFusionsSummary)

knownUnclusteredFusions = read.csv('~/data/sv/fusions/SVA_KNOWN_UNCLUSTERED_FUSIONS.csv')

write.csv(knownUnclusteredFusions, '~/data/sv/fusions/SVA_KNOWN_UNCLUSTERED_FUSIONS.csv', quote=F, row.names=F)
write.csv(knownUnclusteredFusionsSummary, '~/data/sv/fusions/SVA_KNOWN_UNCLUSTERED_FUSION_SUMMARY.csv', quote=F, row.names=F)




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

rm(clusters)


#######
## Known Fusion Data

knownFusionData = read.csv('~/data/sv/sva_known_fusion_data.csv')
View(knownFusionData)
View(knownFusionData %>% group_by(SampleId,GeneUp,GeneDown) %>% summarise(Count=n(),
                                                                          Unclustered=sum(grepl('Unclustered',InvalidReasons)),
                                                                          Unchained=sum(grepl('Unchained',InvalidReasons)),
                                                                          Orientation=sum(grepl('Orientation',InvalidReasons)),
                                                                          Coding=sum(grepl('Coding',InvalidReasons)),
                                                                          Unphased=sum(grepl('Unphased',InvalidReasons))))



# DEBUG ONLY

specSampleFusions = read.csv('~/data/sv/fusions/tmp.csv')
specSampleFusions = read.csv('~/logs/SVA_FUSIONS.csv')
specSampleFusions = annotate_fusions(specSampleFusions)
View(specSampleFusions)
View(specSampleFusions %>% filter(GeneUp=='CCDC6'&GeneDown=='RET'))

View(specSampleFusions %>% filter(GeneUp=='TMPRSS2'&GeneDown=='ERG')) # &!grepl('Unclustered',ResolvedType)

View(specSampleFusions %>% filter(PhaseMatched=='false'))
View(specSampleFusions %>% filter(GeneUp=='TMPRSS2'&GeneDown=='ERG'&PhaseMatched=='false') %>% group_by(SampleId) %>% count())

newFusionsTmp = reportedSvaFusions %>% filter(Clustered&GeneUp=='TMPRSS2'&GeneDown=='ERG') %>% group_by(SameSV,InChain,ValidChain) %>% count()
prevFusionsTmp = reportedSvaFusionsPrev %>% filter(Clustered&GeneUp=='TMPRSS2'&GeneDown=='ERG') %>% group_by(SameSV,InChain,ValidChain) %>% count()


View(newFusionsTmp)
View(prevFusionsTmp)
newFusionsTmp2 = reportedSvaFusions %>% filter(Clustered&GeneUp=='TMPRSS2'&GeneDown=='ERG'&SameSV&InChain&ValidChain)
prevFusionsTmp2 = reportedSvaFusionsPrev %>% filter(Clustered&GeneUp=='TMPRSS2'&GeneDown=='ERG'&SameSV&InChain&ValidChain)
View(newFusionsTmp2)
View(prevFusionsTmp2)


