library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringi)

###############################
# LINX RNA-DNA FUSION ANALYSIS

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
                  mutate(ValidChain=ValidTraversal=='true'&DisruptedExonsUp==0&DisruptedExonsDown==0&TerminatedUp==F&TerminatedDown==F,
                         NonDisruptedSingle=FacingBEsUp==0&FacingBEsDown==0&DisruptedExonsUp==0&DisruptedExonsDown==0,
                         BreakendDistUp=ifelse(StrandUp==1,TransStartUp-PosUp,PosUp-TransEndUp),
                         BreakendDistDown=ifelse(StrandUp==1,TransStartDown-PosDown,PosDown-TransEndDown)))
  
  fusionData[is.na(fusionData)] = 0
  
  return (fusionData)
}

load_rna_match_data<-function(filename)
{
  rnaMatchData = read.csv(filename)
  
  #Annotations
  rnaMatchData = rnaMatchData %>% mutate(SvMatchUp=!is.na(SvIdUp)&TransValidLocUp=='true',
                                         SvMatchDown=!is.na(SvIdDown)&TransValidLocDown=='true',
                                         SvViableUp=TransViableUp=='true',
                                         SvViableDown=TransViableDown=='true',
                                         SvValidLocUp=TransValidLocUp=='true',
                                         SvValidLocDown=TransValidLocDown=='true',
                                         SameSV=SvIdUp==SvIdDown,
                                         SameChr=ChrUp==ChrDown,
                                         FusionDistance=ifelse(SameChr,StrandUp*(RnaPosDown-RnaPosUp),0),
                                         Proximate=FusionDistance>0&FusionDistance<5e5)
  
  rnaMatchData = rnaMatchData %>% mutate(SvMatchType=ifelse(SvMatchUp&SvMatchDown,'BothSVs',
                                                            ifelse(SvMatchUp|SvMatchDown,'SingleSV','NoSV')))
  
  # DEDUP 
  rnaMatchData = rnaMatchData %>% group_by(SampleId,FusionName) %>% arrange(-JunctionReadCount,-SpanningFragCount) %>% filter(row_number()==1) %>% ungroup() 
  rnaMatchData = rnaMatchData %>% group_by(SampleId,RnaPosUp,RnaPosDown) %>% arrange(-JunctionReadCount,-SpanningFragCount) %>% filter(row_number()==1) %>% ungroup() 
  
  # PON: Remove recurrent fusions where NONE have SV support at both ends
  # rnaMatchData = rnaMatchData %>% group_by(FusionName) %>% filter(!((sum(ifelse(SvMatchType=='NoSV',1,0))>1&(sum(ifelse(SvMatchType=='BothSVs',1,0))<1)))) %>% ungroup()
  rnaMatchData = rnaMatchData %>% group_by(FusionName) %>% 
    filter(!((sum(ifelse(SvMatchType=='NoSV',1,0))>1&sum(ifelse(SvMatchType=='BothSVs',1,0))<pmax(1,sum(ifelse(SvMatchType=='NoSV',1,0))-1)))) %>% ungroup()
  
  rnaBothSVsData = rnaMatchData %>% filter(SvMatchType=='BothSVs')
  rnaNotBothSVsData = rnaMatchData %>% filter(SvMatchType!='BothSVs')
  
  rnaBothSVsData = rnaBothSVsData %>% separate(ClusterInfoUp,c('ClusterIdUp','ClusterCountUp','ChainIdUp','ChainCountUp'),sep = ';')
  rnaBothSVsData = rnaBothSVsData %>% separate(ClusterInfoDown,c('ClusterIdDown','ClusterCountDown','ChainIdDown','ChainCountDown'),sep = ';')
  
  rnaBothSVsData = rnaBothSVsData %>% mutate(SameCluster=(ClusterIdUp==ClusterIdDown),
                                             SameChain=ifelse(SameSV,T,SameCluster&ChainIdUp==ChainIdDown))
  
  rnaNotBothSVsData = rnaNotBothSVsData %>% mutate(SameCluster=F,
                                                   SameChain=F,
                                                   ClusterIdUp=-1,ClusterCountUp=0,ChainIdUp=-1,ChainCountUp=0,
                                                   ClusterIdDown=-1,ClusterCountDown=0,ChainIdDown=-1,ChainCountDown=0)
  
  rnaNotBothSVsData = within(rnaNotBothSVsData,rm(ClusterInfoUp))
  rnaNotBothSVsData = within(rnaNotBothSVsData,rm(ClusterInfoDown))
  
  
  rnaMatchData = rbind(rnaBothSVsData,rnaNotBothSVsData)
  
  return (rnaMatchData)
}

annotate_rna_both_svs<-function(rnaMatchData)
{
  rnaBothSVsData = rnaMatchData %>% filter(SvMatchType=='BothSVs')
  
  # rnaBothSVsData = rnaBothSVsData %>% separate(ClusterInfoUp,c('ClusterIdUp','ClusterCountUp','ChainIdUp','ChainCountUp'),sep = ';')
  # rnaBothSVsData = rnaBothSVsData %>% separate(ClusterInfoDown,c('ClusterIdDown','ClusterCountDown','ChainIdDown','ChainCountDown'),sep = ';')
  
  rnaBothSVsData$ClusterCountUp = as.numeric(rnaBothSVsData$ClusterCountUp)
  rnaBothSVsData$ClusterCountDown = as.numeric(rnaBothSVsData$ClusterCountDown)
  rnaBothSVsData$ChainCountUp = as.numeric(rnaBothSVsData$ChainCountUp)
  rnaBothSVsData$ChainCountDown = as.numeric(rnaBothSVsData$ChainCountDown)
  
  rnaBothSVsData$IsChainedUp = (rnaBothSVsData$ChainCountUp>1)
  rnaBothSVsData$IsChainedDown = (rnaBothSVsData$ChainCountDown>1)
  
  rnaBothSVsData = rnaBothSVsData %>% separate(ChainInfo,c('ChainLinks','ChainLength'),sep = ';')
  rnaBothSVsData$ChainLength = as.numeric(rnaBothSVsData$ChainLength)
  rnaBothSVsData$ChainLinks = as.numeric(rnaBothSVsData$ChainLinks)
  
  rnaBothSVsData$ChainLength = as.numeric(rnaBothSVsData$ChainLength)
  rnaBothSVsData$ChainLinks = as.numeric(rnaBothSVsData$ChainLinks)
  
  rnaBothSVsData$FacingInChain = (rnaBothSVsData$ChainLength>0)
  
  return (rnaBothSVsData)
}

known_type_category<-function(knownType)
{
  newKnownType = ifelse(knownType=='Known','Known',ifelse(knownType=='5P-Prom'|knownType=='3P-Prom'|knownType=='Both-Prom','Promiscuous','Unknown'))
  return (newKnownType)
}


#####################
# RNA Summary results

# 1. RNA / DNA Sensitivity
# using RNA Match Data only

hpcDedupedSamples = read.csv('~/data/sv/hpc_non_dup_sample_ids.csv')
# View(hpcDedupedSamples)
#View(highestPurityCohort)

rnaMatchData = load_rna_match_data('~/data/sv/rna/LNX_RNA_DATA.csv')

rnaSampleIds = read.csv('~/data/sv/rna/rna_starfusion_sample_ids.csv')

nrow(rnaSampleIds %>% filter(SampleId %in% hpcDedupedSamples$sampleId))
write.csv(rnaSampleIds %>% filter(SampleId %in% hpcDedupedSamples$sampleId) %>% select(SampleId),
          '~/data/sv/rna/rna_sample_ids_hpc_dedup.csv', row.names = F, quote = F)

# filter out sample which was swapped with another during RNA sampling
rnaMatchData = rnaMatchData %>% filter(SampleId!='CPCT02330014T')

# restrict to HPC deduped cohort
rnaMatchData = rnaMatchData %>% filter(SampleId %in% hpcDedupedSamples$sampleId)

# add RNA biotype to help where no DNA is found
transExonData = read.csv('~/data/sv/ensembl_trans_exon_data.csv') 
colnames(transExonData)
transData = transExonData %>% group_by(GeneId,TransId,StableId=Trans) %>% 
  summarise(Strand=first(Strand),BioType=first(BioType),IsCanonical=first(CanonicalTranscriptId)==first(TransId),
            TransStart=first(TransStart),TransEnd=first(TransEnd),CodingStart=first(CodingStart),CodingEnd=first(CodingEnd)) %>% ungroup()

View(transData)
write.csv(transData,'~/data/sv/ensembl_trans_data.csv',row.names=F, quote=F)
transData = read.csv('~/data/sv/ensembl_trans_data.csv')

# View(transBiotypes)
rnaMatchData = merge(rnaMatchData,transData %>% mutate(RnaTransIdUp=StableId,RnaBiotypeUp=BioType) %>% select(RnaTransIdUp,RnaBiotypeUp),
                     by='RnaTransIdUp',all.x=T)
rnaMatchData = merge(rnaMatchData,transData %>% mutate(RnaTransIdDown=StableId,RnaBiotypeDown=BioType) %>% select(RnaTransIdDown,RnaBiotypeDown),
                     by='RnaTransIdDown',all.x=T)

#View(rnaMatchData)
#View(rnaMatchData %>% group_by(KnownType,DnaMatchType) %>% count %>% spread(DnaMatchType,n))

###############################
# Matching with SVA Fusion data

# load all fusions found for the 630 samples with RNA
svaRnaFusions = read.csv('~/data/sv/rna/LNX_FUSIONS.csv')
svaRnaFusions = annotate_fusions(svaRnaFusions)
svaRnaFusions = svaRnaFusions %>% filter(SampleId %in% hpcDedupedSamples$sampleId)

svaRnaFusions = svaRnaFusions %>% filter(PhaseMatched=='true')

nrow(svaRnaFusions) # 12149 -> 10760 (27-11-19), 11652 (08-12-19)
# 10275 with old priority scheme

#View(svaRnaFusions)

#View(svaRnaFusions %>% group_by(RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown) %>% 
#       summarise(Count=n(),ExonsSkippedUp=sum(ExonsSkippedUp>0),ExonsSkippedDown=sum(ExonsSkippedDown>0)))

rnaReadData = load_rna_match_data('~/data/sv/rna/read_data/SVA_RNA_READ_DATA.csv')
rnaReadData = rnaReadData %>% mutate(SampleGenePair=paste(SampleId,GeneNameUp,GeneNameDown,sep='_'))
#View(rnaReadData)

# create a combined file from the RNA and LINX fusions files
rnaCombinedData = merge(svaRnaFusions, 
                        rnaMatchData, 
                        by=c('SampleId','GeneNameUp','GeneNameDown'),all=T)

rnaCombinedData = rnaCombinedData %>% mutate(HasDnaData=!is.na(KnownType.x),
                                             HasRnaData=!is.na(KnownType.y),
                                             DnaMatchType=ifelse(HasDnaData&grepl('INVALID',DnaMatchType),'GENES',as.character(DnaMatchType)),
                                             OldHasDnaData=!is.na(KnownType.x)|DnaMatchType=='GENES'|DnaMatchType=='SV',
                                             HasSvSupport=SvMatchUp&SvMatchDown,
                                             DnaRnaMatch=HasDnaData&HasRnaData&HasSvSupport&(DnaMatchType=='GENES'|DnaMatchType=='SV'),
                                             SampleGenePair=paste(SampleId,GeneNameUp,GeneNameDown,sep='_'))

rnaCombinedData = rnaCombinedData %>% mutate(HasReadSupport=(HasDnaData&!HasRnaData&SampleGenePair %in% rnaReadData$SampleGenePair))

dnaRnaCombinedData = rnaCombinedData %>% 
  mutate(KnownType=ifelse(!is.na(KnownType.x),as.character(KnownType.x),as.character(KnownType.y)),
         SameSV=ifelse(!is.na(SameSV.x),SameSV.x,SameSV.y),
         SameCluster=ifelse(is.na(SameCluster),T,SameCluster),
         SameChain=ifelse(is.na(SameChain),T,SameChain),
         MatchCategory=ifelse(DnaRnaMatch,'DNA & RNA',ifelse(HasRnaData,'RNA Only',ifelse(HasDnaData,'DNA Only','Other'))),
         MatchType=MatchCategory,
         OldMatchCategory=ifelse(HasRnaData&!OldHasDnaData&SvMatchType!='BothSVs','RNA Only',
                  ifelse(HasReadSupport,'DNA with RNA Read Support',ifelse(OldHasDnaData&!HasRnaData,'DNA Only',
                  ifelse(OldHasDnaData,'DNA & RNA','RNA with DNA Support')))),
         OldMatchType=ifelse(MatchCategory=='DNA & RNA','DNA & RNA',ifelse(MatchCategory=='DNA Only'|MatchCategory=='DNA with RNA Read Support','DNA Only',
                                                                        ifelse(MatchCategory=='RNA with DNA Support','RNA with DNA Support','RNA Only'))),
         KnownCategory=ifelse(KnownType=='Both-Prom','Both promiscuous',ifelse(KnownType=='5P-Prom',"5' promiscuous",ifelse(KnownType=='3P-Prom',"3' promiscuous",ifelse(KnownType=='Known','Known','Unknown')))),
         SameGeneFusion=(as.character(GeneNameUp)==as.character(GeneNameDown)))

# now filter for RNA unphased with no DNA support
dnaRnaCombinedData = dnaRnaCombinedData %>% filter(!(HasRnaData&RnaPhaseMatched=='false'
                                                     &(MatchCategory=='RNA with DNA Support'|MatchCategory=='RNA Only')))

# and filter out if both RNA and DNA are unphased
dnaRnaCombinedData = dnaRnaCombinedData %>% filter(!(PhaseMatched.x=='false'&RnaPhaseMatched=='false'))

# convert RNA with DNA Support to be classified as RNA only
dnaRnaCombinedData = dnaRnaCombinedData %>% mutate(MatchType=ifelse(MatchType=='RNA with DNA Support','RNA Only',MatchType))

# filter out same-gene fusions
dnaRnaCombinedData = dnaRnaCombinedData %>% filter(!SameGeneFusion)

dnaRnaSummary = dnaRnaCombinedData %>% filter(KnownCategory!='Unknown') %>% group_by(MatchType,KnownCategory) %>% count()
#View(dnaRnaSummary)
View(dnaRnaSummary %>% spread(MatchType,n))
View(dnaRnaCombinedData %>% group_by(MatchType,KnownCategory) %>% count() %>% spread(MatchType,n))

#View(dnaRnaCombinedData)

# swap SampleIds for HMF IDs
sampleIdMapping = read.csv('~/data/sv/sample_id_mapping.csv')
View(sampleIdMapping)
dnaRnaCombinedData = merge(dnaRnaCombinedData,sampleIdMapping,by='SampleId',all.x=T)

dnaRnaCombinedData = dnaRnaCombinedData %>% mutate(Reportable=(Reportable=='true'))

dnaRnaCombinedOutputData = dnaRnaCombinedData %>% 
  mutate(KnownType=ifelse(!is.na(KnownType.x),as.character(KnownType.x),as.character(KnownType.y)),
         ChrUp=ifelse(!is.na(ChrUp.x),as.character(ChrUp.x),as.character(ChrUp.y)),
         ChrDown=ifelse(!is.na(ChrDown.x),as.character(ChrDown.x),as.character(ChrDown.y)),
         PosUp=ifelse(!is.na(PosUp.x),PosUp.x,PosUp.y),PosDown=ifelse(!is.na(PosDown.x),PosDown.x,PosDown.y),
         OrientUp=ifelse(!is.na(OrientUp.x),OrientUp.x,OrientUp.y),OrientDown=ifelse(!is.na(OrientDown.x),OrientDown.x,OrientDown.y),
         RnaPosUp,RnaPosDown,TransValidLocUp,TransViableUp,TransValidLocDown,TransViableDown,
         TranscriptIdUp=ifelse(!is.na(TranscriptUp),as.character(TranscriptUp),as.character(TransIdUp)),
         TranscriptIdDown=ifelse(!is.na(TranscriptDown),as.character(TranscriptDown),as.character(TransIdDown)),
         CodingTypeUp=ifelse(!is.na(CodingTypeUp.x),as.character(CodingTypeUp.x),as.character(CodingTypeUp.y)),
         RegionTypeUp=ifelse(!is.na(RegionTypeUp.x),as.character(RegionTypeUp.x),as.character(RegionTypeUp.y)),
         CodingTypeDown=ifelse(!is.na(CodingTypeDown.x),as.character(CodingTypeDown.x),as.character(CodingTypeDown.y)),
         RegionTypeDown=ifelse(!is.na(RegionTypeDown.x),as.character(RegionTypeDown.x),as.character(RegionTypeDown.y)),
         ExonsSkippedUp=ifelse(!is.na(ExonsSkippedUp.x),as.character(ExonsSkippedUp.x),as.character(ExonsSkippedUp.y)),
         ExonsSkippedDown=ifelse(!is.na(ExonsSkippedDown.x),as.character(ExonsSkippedDown.x),as.character(ExonsSkippedDown.y)),
         SameCluster=ifelse(is.na(SameCluster),T,SameCluster),SameChain=ifelse(is.na(SameChain),T,SameChain)) %>%
  select(HmfId,SampleId,GeneNameUp,GeneNameDown,MatchType,KnownCategory,Reportable,ClusterId,BiotypeUp,BiotypeDown,
         ChrUp,RnaPosUp,DnaPosUp=PosUp,DnaOrientUp=OrientUp,
         ChrDown,RnaPosDown,DnaPosDown=PosDown,DnaOrientDown=OrientDown,
         RnaJunctionReadCount=JunctionReadCount,RnaSpanningFragCount=SpanningFragCount,
         CodingTypeUp,CodingTypeDown,RegionTypeUp,RegionTypeDown,TranscriptIdUp,TranscriptIdDown,
         ChainLinks,ChainLength,BreakendExonUp,BreakendExonDown,ExonsSkippedUp,ExonsSkippedDown,
         RnaTransIdUp,RnaTransIdDown,RnaBiotypeUp,RnaBiotypeDown)

dnaRnaCombinedOutputData = dnaRnaCombinedOutputData %>% mutate(DnaPosUp=ifelse(DnaPosUp>0,DnaPosUp,''),
                                                               DnaPosDown=ifelse(DnaPosDown>0,DnaPosDown,''),
                                                               DnaOrientUp=ifelse(DnaOrientUp!=0,DnaOrientUp,''),
                                                               DnaOrientDown=ifelse(DnaOrientDown!=0,DnaOrientDown,''),
                                                               ChainLinks=ifelse(is.na(ChainLinks),0,ChainLinks),
                                                               ChainLength=ifelse(is.na(ChainLength),0,ChainLength))

View(dnaRnaCombinedOutputData)

# write.csv(dnaRnaCombinedOutputData,'~/data/sv/rna/LINX_dna_rna_combined_data_tmp.csv', quote = F, row.names = F)
write.csv(dnaRnaCombinedOutputData %>% select(-SampleId,-ClusterId),'~/data/sv/rna/LINX_dna_rna_combined_data_20191211.csv',quote=F,row.names=F)
write.csv(dnaRnaCombinedOutputData %>% select(-SampleId,-ClusterId),'~/data/sv/rna/LINX_dna_rna_combined_data.csv',quote=F,row.names=F)


# PLOT 1: Precision report

dnaRnaCombinedOutputData = read.csv('~/data/sv/rna/LINX_dna_rna_combined_data.csv')
View(dnaRnaCombinedOutputData)

# create a summary view to plot the precision results 
dnaRnaSummary = dnaRnaCombinedOutputData %>% filter(KnownCategory!='Unknown') %>% group_by(MatchType,KnownCategory) %>% count()

View(dnaRnaSummary)
View(dnaRnaSummary %>% spread(MatchType,n))

plotColours4 = c('blue','skyblue3','lightblue','khaki4','khaki3','sienna1')

dnaRnaSummaryPlot = (ggplot(dnaRnaSummary , aes(x=KnownCategory, y=n, fill=MatchType))
                     + geom_bar(stat = "identity", colour = "black", position = position_stack(reverse = TRUE))
                     + labs(x = "", y="Fusions", fill='Fusion Prediction', title = "")
                     + scale_fill_manual(values = plotColours4)
                     + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                     + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                     + theme(axis.text.x = element_text(angle=90, hjust=1,size=10))
                     + coord_flip())

print(dnaRnaSummaryPlot)


######
## PLOT 2: LINX Fusion Sensitivity
rnaMatchDataBothSVs = annotate_rna_both_svs(rnaMatchData)

summaryBothData = rnaMatchDataBothSVs %>% 
  mutate(SvaCategory=ifelse(SameCluster&SameChain,'Matched',ifelse(SameCluster&!SameChain,'DiffChain','DiffCluster')),
         ValidBreakends=SvViableUp&SvViableDown,
         ViableFusion=ViableFusion=='true',
         PhaseMatched=PhaseMatched=='true',
         RnaPhaseMatched=RnaPhaseMatched=='true',
         SameSV=SameSV)

summaryNotBothData = rnaMatchData %>% filter(SvMatchType!='BothSVs') %>% 
  mutate(SvaCategory = ifelse(SvMatchType=='SingleSV','SingleBEMatch','NoMatch'),
         ValidBreakends=F,
         ViableFusion=F,
         PhaseMatched=F,
         RnaPhaseMatched=RnaPhaseMatched=='true',
         SameSV=T)

summaryRnaData = rbind(summaryBothData %>% select(KnownType,SvaCategory,ValidBreakends,ViableFusion,PhaseMatched,RnaPhaseMatched,SameSV),
                       summaryNotBothData %>% select(KnownType,SvaCategory,ValidBreakends,ViableFusion,PhaseMatched,RnaPhaseMatched,SameSV))

summaryRnaData = summaryRnaData %>% 
  mutate(SvaCategory2=ifelse(SvaCategory=='Matched',ifelse(ValidBreakends,'Matched','MatchedExonsSkipped'),
                             ifelse(RnaPhaseMatched,'NotCalled','NoRnaPhasedFusion')),
         FusionType=known_type_category(KnownType))

summaryRnaData$FusionType = "AllFusions" # no split by KnownType

rnaCategorySummary1 = summaryRnaData %>% group_by(SvaCategory2,FusionType) %>% count() %>% spread(SvaCategory2,n) %>% 
  arrange(FusionType) %>% ungroup()
rnaCategorySummary1[is.na(rnaCategorySummary1)] = 0

rnaCategorySummaryData1 = rnaCategorySummary1 %>% select(FusionType,Matched,MatchedExonsSkipped,NotCalled)
rnaCategorySummaryData1 = rnaCategorySummaryData1 %>% gather('Category','Count', 2:ncol(rnaCategorySummaryData1)) 
View(rnaCategorySummaryData1)

plotColours3 = c('blue','light blue','orangered','sienna1','khaki4','khaki3','palegreen', 'seagreen')

rnaSummaryDataPlot1 = (ggplot(rnaCategorySummaryData1, aes(x=FusionType, y=Count, fill=Category))
                       + geom_bar(stat = "identity", colour = "black", position = position_stack(reverse = TRUE))
                       + labs(x='',y="Fusion Count", fill='Category', title='Fusion Sensitivity')
                       + scale_fill_manual(values = plotColours3)
                       + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                       + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                       + theme(axis.text.x = element_text(angle=90, hjust=1,size=10))
                       + coord_flip())

plot(rnaSummaryDataPlot1)




#####
## Clean Up
rm(dnaRnaCombinedData)
rm(dnaRnaCombinedOutputData)
rm(svaRnaFusions)
rm(rnaReadData)
rm(rnaMatchData)
rm(rnaMatchDataBothSVs)
rm(rnaCombinedData)
rm(dnaRnaSummaryPlot)
rm(dnaRnaSummary)


## Validation
View(dnaRnaCombinedData %>% group_by(KnownType,KnownCategory) %>% count())
View(dnaRnaCombinedData %>% group_by(KnownCategory,MatchCategory,MatchType) %>% count())


View(dnaRnaCombinedData %>% filter(!HasDnaData&SameChain&SameCluster&SvMatchUp&SvMatchDown) %>%
       group_by(PhaseMatched.y,ViableFusion,KnownCategory,SameSV.y) %>% count %>% spread(SameSV.y,n))

View(dnaRnaCombinedData %>% filter(!HasDnaData&SameChain&SameCluster&SvMatchUp&SvMatchDown) %>%
       group_by(PhaseMatched.y,ViableFusion,KnownCategory,SameSV.y,CodingTypeUp.y,CodingTypeDown.y) %>% count %>% spread(SameSV.y,n))

View(dnaRnaCombinedData %>% filter(HasDnaData&is.na(KnownType.x)) %>% group_by(KnownCategory) %>% count())

View(dnaRnaCombinedData %>% filter(!HasDnaData&SameChain&SameCluster&SvMatchUp&SvMatchDown&CodingTypeUp.y=='5P_UTR'&CodingTypeDown.y=='Coding'))

View(dnaRnaCombinedData %>% filter(MatchCategory=='RNA Only'&KnownCategory!='Unknown'))
View(dnaRnaCombinedData %>% filter(MatchCategory=='RNA with DNA Support'&KnownCategory!='Unknown'))
View(dnaRnaCombinedData %>% filter(MatchCategory=='DNA & RNA'&KnownCategory=='Known'))

View(dnaRnaCombinedData %>% filter(MatchCategory=='RNA with DNA Support') %>% 
       group_by(PhaseMatched.y,ViableFusion,SameSV.y,RegionTypeUp.y,RegionTypeDown.y,CodingTypeUp.y,CodingTypeDown.y) %>% count())

View(dnaRnaCombinedData %>% filter(MatchCategory=='RNA with DNA Support') %>% group_by(SampleId) %>% count())
View(dnaRnaCombinedData %>% filter(MatchCategory=='RNA with DNA Support'&SameCluster&SameChain&DnaMatchType=='NONE'))
View(dnaRnaCombinedData %>% filter(MatchCategory=='RNA with DNA Support'&DistancePrevDown.y<0))
View(dnaRnaCombinedData %>% filter(DistancePrevDown.y<0) %>% group_by(MatchCategory) %>% count())

View(dnaRnaCombinedData %>% filter(MatchCategory=='RNA with DNA Support') %>% 
       group_by(DnaMatchType,Unchained=(DnaMatchType=='NONE'&!SameChain),
                Unclustered=(DnaMatchType=='NONE'&!SameCluster),
                InvalidDistanceDown=DistancePrevDown.y<0) %>% count())

View(dnaRnaCombinedData %>% filter(MatchCategory=='RNA with DNA Support') %>% 
       group_by(PhaseMatched.y,ViableFusion,SameSV.y,SameCluster,SameChain) %>% count())


colnames(dnaRnaCombinedData)




