library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringi)

annotate_fusions<-function(fusionData)
{
  fusionData = fusionData %>% mutate(SameSV = (SvIdUp==SvIdDown),
                                     Clustered = ifelse(grepl('Unclustered',ResolvedType),F,T))
  
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
reportedSvaFusions = read.csv('~/data/sv/fusions/SVA_FUSIONS.csv')
reportedSvaFusionsRaw = read.csv('~/logs/SVA_FUSIONS.csv')
reportedSvaFusionsPrev2 = reportedSvaFusions
# read.csv('~/data/sv/fusions/SVA_FUSIONS.csv')
View(reportedSvaFusionsRaw)
nrow(reportedSvaFusionsRaw)
rm(svaFusions)

# Annotations
reportedSvaFusions = annotate_fusions(reportedSvaFusionsRaw)

View(reportedSvaFusions)
write.csv(reportedSvaFusions, '~/data/sv/fusions/LINX_FUSIONS_REPORTED_2_rules.csv', row.names = F, quote = F)
reportedSvaFusions = read.csv('~/data/sv/fusions/LINX_FUSIONS_REPORTED.csv')
nrow(reportedSvaFusions)
View(reportedSvaFusions)

reportedSvaFusions = read.csv('~/data/sv/fusions/LINX_FUSIONS_REPORTED_2_rules.csv')
View(reportedSvaFusions)

# exons skipped
nrow(reportedSvaFusions %>% filter(ExonsSkippedUp>0|ExonsSkippedDown>0))
View(reportedSvaFusions %>% filter(ExonsSkippedUp>0|ExonsSkippedDown>0))

# known with terminated ends
View(reportedSvaFusions %>% filter(KnownType=='Known'&(TerminatedUp|TerminatedDown)))
View(reportedSvaFusions %>% filter(TerminatedUp|TerminatedDown) %>% group_by(KnownType) %>% count())

# distance upstream from transcript for known and promicuous fusions
View(reportedSvaFusions %>% filter(BreakendDistUp>=5e4|BreakendDistDown>=5e4) %>% group_by(KnownType) %>% count())
View(reportedSvaFusions %>% filter(BreakendDistUp>=5e5|BreakendDistDown>=5e5) %>% group_by(KnownType) %>% count())
View(reportedSvaFusions %>% filter(BreakendDistUp>=1e4|BreakendDistDown>=1e4))

knownNotReported = svaFusions %>% filter(KnownType=='Known'&Reportable=='false')
knownNotReported = annotate_fusions(knownNotReported)
View(knownNotReported)
View(knownNotReported %>% filter(!grepl('Unclustered',ResolvedType)))

# check for duplicates
View(reportedSvaFusions %>% group_by(SampleId,GeneUp,GeneDown) %>% count())


nrow(reportedSvaFusions %>% filter(GeneUp=='TMPRSS2'&GeneDown=='ERG'))
nrow(reportedSimpleFusions %>% filter(GeneUp=='TMPRSS2'&GeneDown=='ERG'))
View(reportedSvaFusions %>% filter(GeneUp=='TMPRSS2'&GeneDown=='ERG') %>% group_by(SameSV) %>% count())
View(reportedSvaFusionsPrev %>% filter(SameSV&GeneUp=='TMPRSS2'&GeneDown=='ERG'))

# multiple fusions per sample
View(reportedSvaFusions %>% filter(Clustered) %>% group_by(SampleId,GeneUp,GeneDown) %>% count())
View(reportedSvaFusionsPrev %>% group_by(SampleId,GeneUp,GeneDown) %>% count())

View(reportedSvaFusions %>% filter(Clustered) %>% group_by(SampleId,GeneUp,GeneDown) %>% count())

# comparison with previous run
nrow(reportedSvaFusionsPrev %>% filter(Clustered&ValidChain))
nrow(reportedSvaFusions %>% filter(Clustered&ValidChain))
nrow(reportedSvaFusionsPrev %>% filter(Clustered&ValidChain&SameSV))
nrow(reportedSvaFusions %>% filter(Clustered&ValidChain&SameSV))
nrow(reportedSvaFusionsPrev %>% filter(Clustered&SameSV)) # 585
nrow(reportedSvaFusions %>% filter(Clustered&SameSV)) # 538
nrow(reportedSvaFusionsPrev %>% filter(Clustered&ValidChain&!SameSV)) # 74
nrow(reportedSvaFusions %>% filter(Clustered&ValidChain&!SameSV)) # 90

# basic numbers
nrow(reportedSvaFusions %>% filter(SameSV)) # 456 vs 585 single-SV fusions
nrow(reportedSvaFusions %>% filter(SameSV&!InChain)) # 274
nrow(reportedSvaFusions %>% filter(ChainLinks==0))
nrow(reportedSvaFusions %>% filter(!SameSV)) # 79


# fusion differences after using transcript
newReported = allFusions %>% filter(Reportable=='true')
nrow(newReported)

reportedSvaFusions = read.csv('~/data/sv/fusions/LINX_FUSIONS_REPORTED_2_rules.csv')
nrow(reportedSvaFusions)

View(reportedSvaFusions %>% group_by(SampleId,GeneUp,GeneDown) %>% count())
View(reportedSvaFusions %>% group_by(KnownType) %>% count())
View(newReported %>% group_by(KnownType) %>% count())

knownPrev = reportedSvaFusions %>% filter(KnownType=='Known')
knownNew = newReported %>% filter(KnownType=='Known')
knownPrevDiff = merge(knownPrev,knownNew,by.x=c('SampleId','GeneUp','GeneDown'),by.y=c('SampleId','GeneNameUp','GeneNameDown'),all.x=T)
View(knownPrevDiff %>% filter(is.na(SvIdUp.y)))
View(knownPrevDiff)


View(newReported %>% filter(Reportable=='true'&KnownType=='Known'&GeneNameUp=='TMPRSS2'&GeneNameDown=='ERG'))
View(newReported %>% filter(Reportable=='true'&KnownType=='Known'&GeneNameUp=='TMPRSS2'&GeneNameDown=='ERG') %>% group_by(TranscriptDown) %>% summarise(ProteinsKept=first))


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


# Comparison with VariantAnnotator
reportedSvaSingleSvFusions = reportedSvaFusions %>% filter(SameSV)
View(reportedSvaSingleSvFusions)
nrow(reportedSvaSingleSvFusions)

nrow(reportedSvaFusions %>% filter(KnownType=='Known'))
nrow(reportedSimpleFusions %>% filter(KnownType=='Known'))

vaSvaOverlap = merge(reportedSimpleFusions,reportedSvaFusions, by=c('SampleId','GeneUp','GeneDown'),all.x=T)
View(vaSvaOverlap)
nrow(vaSvaOverlap %>% filter(is.na(Reportable.y))) 
View(vaSvaOverlap %>% filter(is.na(Reportable.y)))
View(vaSvaOverlap %>% filter(is.na(Reportable.y)) %>% select(SampleId,GeneUp,GeneDown,KnownType.x))
View(vaSvaOverlap %>% filter(is.na(Reportable.y)) %>% group_by(GeneUp,GeneDown,KnownType.x) %>% count())
View(vaSvaOverlap %>% filter(is.na(Reportable.y)) %>% group_by(KnownType.x) %>% count())

svaVaOverlap = merge(reportedSvaFusions,reportedSimpleFusions,by=c('SampleId','GeneUp','GeneDown'),all.x=T)
View(svaVaOverlap)
nrow(svaVaOverlap %>% filter(is.na(Reportable.y)))
View(svaVaOverlap %>% filter(is.na(Reportable.y)))
View(svaVaOverlap %>% filter(is.na(Reportable.y)) %>% select(SampleId,GeneUp,GeneDown,KnownType.x))
View(svaVaOverlap %>% filter(is.na(Reportable.y)) %>% group_by(GeneUp,GeneDown,KnownType.x) %>% count())
View(svaVaOverlap %>% filter(is.na(Reportable.y)) %>% group_by(KnownType.x) %>% count())

View(svaVaOverlap %>% filter(is.na(Reportable.y)) %>% group_by(KnownType.x,SameSV) %>% count())

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


View(svData %>% filter(Id==15629227))


# specific sample
specificFusions = read.csv('~/logs/CPCT02070292T.linx.fusions_detailed.csv')
specificFusions = annotate_fusions(specificFusions)
View(specificFusions)
View(specificFusions %>% filter(Reportable=='true'))

View(specificFusions %>% filter(GeneUp=='TMPRSS2'&GeneDown=='ERG'
                                &(ValidTraversal=='true'&PhaseMatched=='true'&BiotypeDown!='nonsense_mediated_decay')))



write.csv(specificFusions%>% filter(GeneUp=='TMPRSS2'&GeneDown=='ERG'),'~/logs/CPCT02070292T_fusions.csv', row.names = F, quote = F)


# DRUP vs DNDS TSGs
drupTsgs = read.csv('~/data/drup_genes.csv')
dndsTsgs = read.csv('~/hmf/repos/hmftools/hmf-common/src/main/resources/dnds/DndsDriverLikelihoodTsg.tsv',sep='\t')
View(dndsTsgs)
View(drupTsgs %>% filter(!(Gene %in% dndsTsgs$gene)))






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

View(knownFusionData %>% filter(SampleId=='CPCT02020449T'))



########
# RNA vs SVA comparison and precision estimation

# rm(svaFusions)
rnaMatchData = load_rna_match_data('~/data/sv/rna/SVA_RNA_DATA.csv')
rnaMatchDataBothSVs = annotate_rna_both_svs(rnaMatchData)
View(rnaMatchData)

rm(rnaSvaFusions)  
rnaSvaFusions = svaFusions %>% filter(SampleId %in% rnaMatchData$SampleId)
nrow(rnaSvaFusions)
View(rnaSvaFusions)

rnaSvaFusions = annotate_fusions(rnaSvaFusions)

View(rnaSvaFusions %>% filter(ValidTraversal=='true'&TerminatedUp==0) %>% group_by(GeneUp,GeneDown,SampleId) %>% summarise(min(ChainLength)))

View(rnaSvaFusions %>% group_by(ValidTraversal,TerminatedUp,TerminatedDown) %>% count())

write.csv(rnaSvaFusions, '~/data/sv/fusions/SVA_RNA_FUSIONS_ANNOT.csv', row.names = F, quote = F)
rnaSvaFusions = read.csv('~/data/sv/fusions/SVA_RNA_FUSIONS_ANNOT.csv')
View(head(rnaSvaFusions,100))

View(rnaMatchData %>% group_by(SampleId,GeneUp,GeneDown) %>% count()) # 3048 RNA fusions by sample and gene-pair
rnaSampleData = rnaMatchData %>% group_by(SampleId,GeneUp,GeneDown) %>% count()

rnaSvaSampleData = rnaSvaFusions %>% group_by(SampleId,GeneUp,GeneDown) %>% 
  summarise(Total=n(),
            ValidCount=sum(ValidTraversal=='true'&TerminatedUp==0&TerminatedDown==0),
            SingleSVs=sum(SameSV),
            DiffSVs=sum(!SameSV&InChain),
            Reported=sum(Reportable=='true'),
            Known=sum(KnownType=='Known'))

View(rnaSvaSampleData) # 158K
nrow(rnaSvaSampleData %>% filter(ValidCount==0)) # 126K
nrow(rnaSvaSampleData %>% filter(Clustered==0)) # only 150
nrow(rnaSvaSampleData %>% filter(ValidCount>0)) # 32K

rnaSvaValidSampleData = rnaSvaSampleData %>% filter(ValidCount>0)
nrow(rnaSvaValidSampleData)

rnaSvaSampleDataMatched = merge(rnaSvaValidSampleData,rnaMatchData,by=c('SampleId','GeneUp','GeneDown'),all.x=T)
View(rnaSvaSampleDataMatched)
View(rnaSvaSampleDataMatched %>% group_by(is.na(FusionName)) %>% count())
View(rnaSvaSampleDataMatched %>% group_by(GeneUp,GeneDown,Matched=!is.na(FusionName)) %>% count())


View(rnaSvaSampleDataMatched %>% group_by(Matched=!is.na(FusionName)) %>% summarise(ValidCount=sum(ValidCount>0),
                                                                                    SingleSVs=sum(SingleSVs>0),
                                                                                    DiffSVs=sum(DiffSVs>0),
                                                                                    Reported=sum(Reported>0),
                                                                                    Known=sum(Known>0)))



# Annotations

# we already looked from point of view of the RNA fusions, now need to do the equivalent from POV of the SVA fusions
# ie how many find an RNA fusions of those which ought to (being viable ones)


validSvaRnaFusions = rnaSvaFusions %>% filter(ValidTraversal=='true')
View(validSvaRnaFusions)







# DEBUG ONLY
nrow(svaFusions %>% filter(SampleId=='CPCT02020924T'))

specSampleFusions = read.csv('~/data/sv/fusions/tmp.csv')
specSampleFusions = read.csv('~/logs/SVA_FUSIONS.csv')
specSampleFusions = annotate_fusions(specSampleFusions)
View(specSampleFusions)
View(specSampleFusions %>% filter(GeneUp=='CCDC6'&GeneDown=='RET'))

View(specSampleFusions %>% filter(GeneUp=='TMPRSS2'&GeneDown=='ERG')) # &!grepl('Unclustered',ResolvedType)
specSvData = read.csv('~/logs/CPCT02070338T_SVA.csv')
View(specSvData)

View(specSampleFusions %>% filter(PhaseMatched=='false'))
View(specSampleFusions %>% filter(GeneUp=='TMPRSS2'&GeneDown=='ERG'&PhaseMatched=='false') %>% group_by(SampleId) %>% count())

View(svaFusions %>% filter(SampleId=='CPCT02020348T'&Reportable=='true'))
View(reportedSvaFusions %>% filter(SampleId=='CPCT02030426T'))
View(reportedSimpleFusions %>% filter(SampleId=='CPCT02030426T'))


tmp = svaFusions %>% filter(SampleId=='CPCT02140052T'&SvIdUp==14975994)
tmp = annotate_fusions(tmp)
View(tmp %>% filter(TranscriptUp=='ENST00000332149'&TranscriptDown=='ENST00000417133'))
View(reportedSvaFusions %>% filter(SampleId=='CPCT02140025T'&SvIdUp==15785387&TranscriptUp=='ENST00000398585'&TranscriptDown=='ENST00000417133'))
View(reportedSvaFusionsPrev %>% filter(SampleId=='CPCT02140052T'&SvIdUp==14975994&TranscriptUp=='ENST00000398585'&TranscriptDown=='ENST00000417133'))


newFusionsTmp = reportedSvaFusions %>% filter(Clustered&GeneUp=='TMPRSS2'&GeneDown=='ERG') %>% group_by(SameSV,InChain,ValidChain) %>% count()
prevFusionsTmp = reportedSvaFusionsPrev %>% filter(Clustered&GeneUp=='TMPRSS2'&GeneDown=='ERG') %>% group_by(SameSV,InChain,ValidChain) %>% count()
View(reportedSvaFusionsPrev %>% filter(SampleId=='CPCT02011025'&Clustered&GeneUp=='TMPRSS2'&GeneDown=='ERG'))
View(reportedSvaFusionsPrev %>% filter(SampleId=='CPCT02030417T'&Clustered&GeneUp=='CCDC6'&GeneDown=='RET'))
View(reportedSvaFusionsPrev %>% filter(SampleId=='CPCT02150031T'&Clustered&GeneUp=='TMPRSS2'&GeneDown=='ETV4'))


View(reportedSvaFusions %>% filter(SampleId=='CPCT02070036T'&Clustered&GeneUp=='FGFR1'&GeneDown=='ANXA4'))

View(newFusionsTmp)
View(prevFusionsTmp)
newFusionsTmp2 = reportedSvaFusions %>% filter(Clustered&GeneUp=='TMPRSS2'&GeneDown=='ERG'&SameSV&InChain&ValidChain)
prevFusionsTmp2 = reportedSvaFusionsPrev %>% filter(Clustered&GeneUp=='TMPRSS2'&GeneDown=='ERG'&SameSV&InChain&ValidChain)
View(newFusionsTmp2)
View(prevFusionsTmp2)


