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
                  mutate(ValidChain=ValidTraversal=='true'&DisruptedExonsUp==0&DisruptedExonsDown==0&TerminatedUp=='false'&TerminatedDown=='false',
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
  
  rnaMatchData$SvMatchType = ifelse(rnaMatchData$SvMatchUp&rnaMatchData$SvMatchDown,'BothSVs',
                                    ifelse(rnaMatchData$SvMatchUp|rnaMatchData$SvMatchDown,'SingleSV','NoSV'))
  
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
View(hpcDedupedSamples)

View(highestPurityCohort)

rnaMatchData = load_rna_match_data('~/data/sv/rna/LNX_RNA_DATA.csv')

rnaSampleIds = read.csv('~/data/sv/rna/rna_starfusion_sample_ids.csv')
View(rnaSampleIds)

nrow(rnaSampleIds %>% filter(SampleId %in% hpcDedupedSamples$sampleId))
write.csv(rnaSampleIds %>% filter(SampleId %in% hpcDedupedSamples$sampleId) %>% select(SampleId),
          '~/data/sv/rna/rna_sample_ids_hpc_dedup.csv', row.names = F, quote = F)


#nrow(rnaSampleIds %>% filter(SampleId %in% cohort$sampleId))
#View(rnaMatchData %>% group_by(SampleId) %>% count())

# restrict to HPC deduped cohort
rnaMatchData = rnaMatchData %>% filter(RnaPhaseMatched=='true')

# filter out sample which was swapped with another during RNA sampling
rnaMatchData = rnaMatchData %>% filter(SampleId!='CPCT02330014T')

# filter out unphased RNA fusions for all subsequent analysis
rnaMatchData = rnaMatchData %>% filter(SampleId %in% hpcDedupedSamples$sampleId)

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

plotColours3 = c('royal blue','light blue','orangered','sienna1','khaki4','khaki3','palegreen', 'seagreen')

rnaSummaryDataPlot1 = (ggplot(rnaCategorySummaryData1, aes(x=FusionType, y=Count, fill=Category))
                       + geom_bar(stat = "identity", colour = "black", position = position_stack(reverse = TRUE))
                       + labs(x='',y="Fusion Count", fill='Category', title='Fusion Sensitivity')
                       + scale_fill_manual(values = plotColours3)
                       + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                       + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                       + theme(axis.text.x = element_text(angle=90, hjust=1,size=10))
                       + coord_flip())

## PLOT 1: LINX Fusion Sensitivity

plot(rnaSummaryDataPlot1)




###############################
# Matching with SVA Fusion data
# 'Precision' report

# load all fusions found for the 630 samples with RNA
svaRnaFusions = read.csv('~/data/sv/rna/LNX_FUSIONS.csv')
svaRnaFusions = annotate_fusions(svaRnaFusions)
svaRnaFusions = svaRnaFusions %>% filter(SampleId %in% hpcDedupedSamples$sampleId)
svaRnaFusions = svaRnaFusions %>% filter(PhaseMatched=='true')
nrow(svaRnaFusions) # 12149

View(svaRnaFusions %>% group_by(RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown) %>% 
       summarise(Count=n(),ExonsSkippedUp=sum(ExonsSkippedUp>0),ExonsSkippedDown=sum(ExonsSkippedDown>0)))

View(svaRnaFusions %>% filter(CodingTypeUp=='3P_UTR') %>% group_by(KnownType) %>% count())
View(svaRnaFusions %>% filter(CodingTypeUp=='3P_UTR'&KnownType!=''))

rnaReadData = load_rna_match_data('~/data/sv/rna/read_data/SVA_RNA_READ_DATA.csv')
rnaReadData = rnaReadData %>% mutate(SampleGenePair=paste(SampleId,GeneNameUp,GeneNameDown,sep='_'))
View(rnaReadData)

# create a combined file from the RNA and LINX fusions files
rnaCombinedData = merge(svaRnaFusions, 
                        rnaMatchData %>% filter(RnaPhaseMatched=='true'),
                        by=c('SampleId','GeneNameUp','GeneNameDown'),all=T)

rnaCombinedData = rnaCombinedData %>% mutate(HasDnaData=!is.na(KnownType.x)|DnaMatchType=='GENES'|DnaMatchType=='SV',
                                             HasRnaData=!is.na(KnownType.y),
                                             HasSvSupport=SvMatchUp&SvMatchDown,
                                             SampleGenePair=paste(SampleId,GeneNameUp,GeneNameDown,sep='_'))

rnaCombinedData = rnaCombinedData %>% mutate(HasReadSupport=(HasDnaData&!HasRnaData&SampleGenePair %in% rnaReadData$SampleGenePair))

dnaRnaCombinedData = rnaCombinedData %>% 
  mutate(KnownType=ifelse(!is.na(KnownType.x),as.character(KnownType.x),as.character(KnownType.y)),
         SameSV=ifelse(!is.na(SameSV.x),SameSV.x,SameSV.y),
         SameCluster=ifelse(is.na(SameCluster),T,SameCluster),
         SameChain=ifelse(is.na(SameChain),T,SameChain),
         MatchCategory=ifelse(HasRnaData&!HasDnaData&SvMatchType!='BothSVs','RNA Only',
                  ifelse(HasReadSupport,'DNA with RNA Read Support',ifelse(HasDnaData&!HasRnaData,'DNA Only',
                  ifelse(HasDnaData,'DNA & RNA','RNA with DNA Support')))),
         KnownCategory=ifelse(KnownType=='Both-Prom','Both promiscuous',ifelse(KnownType=='5P-Prom',"5' promiscuous",ifelse(KnownType=='3P-Prom',"3' promiscuous",ifelse(KnownType=='Known','Known','Unknown')))),
         MatchType=ifelse(MatchCategory=='DNA & RNA','DNA & RNA',ifelse(MatchCategory=='DNA Only'|MatchCategory=='DNA with RNA Read Support','DNA Only',
                   ifelse(MatchCategory=='RNA with DNA Support','RNA with DNA Support','RNA Only'))))

# validation
View(dnaRnaCombinedData %>% group_by(KnownType,KnownCategory) %>% count())
View(dnaRnaCombinedData %>% group_by(KnownCategory,MatchCategory,MatchType) %>% count())

View(dnaRnaCombinedData %>% filter(!HasDnaData&SameChain&SameCluster))
View(dnaRnaCombinedData %>% filter(!HasDnaData&SameChain&SameCluster&TransViableUp=='true'&TransViableDown=='true'))
View(dnaRnaCombinedData %>% filter(!HasDnaData&SameChain&SameCluster&TransViableUp=='true'&TransViableDown=='true'))

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

View(dnaRnaCombinedData %>% filter(MatchCategory=='RNA with DNA Support') %>% 
       group_by(DnaMatchType,Unchained=(DnaMatchType=='NONE'&!SameChain),
                Unclustered=(DnaMatchType=='NONE'&!SameCluster)) %>% count())

View(dnaRnaCombinedData %>% filter(MatchCategory=='RNA with DNA Support') %>% 
       group_by(PhaseMatched.y,ViableFusion,SameSV.y,SameCluster,SameChain) %>% count())


colnames(dnaRnaCombinedData)

View(dnaRnaCombinedData %>% group_by(KnownCategory,MatchCategory) %>% count() %>% spread(MatchCategory,n))

dnaRnaSummary = dnaRnaCombinedOutputData %>% filter(KnownCategory!='Unknown') %>% group_by(MatchType,KnownCategory) %>% count()
View(dnaRnaSummary)
View(dnaRnaSummary %>% spread(MatchType,n))


# swap SampleIds for HMF IDs
sampleIdMapping = read.csv('~/data/sv/sample_id_mapping.csv')
dnaRnaCombinedData = merge(dnaRnaCombinedData,sampleIdMapping,by='SampleId',all.x=T)

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
  select(HmfId,SampleId,GeneNameUp,GeneNameDown,MatchType,KnownCategory,ClusterId,
         ChrUp,RnaPosUp,DnaPosUp=PosUp,DnaOrientUp=OrientUp,
         ChrDown,RnaPosDown,DnaPosDown=PosDown,DnaOrientDown=OrientDown,
         RnaJunctionReadCount=JunctionReadCount,RnaSpanningFragCount=SpanningFragCount,
         CodingTypeUp,CodingTypeDown,RegionTypeUp,RegionTypeDown,TranscriptIdUp,TranscriptIdDown,
         ChainLinks,ChainLength,BreakendExonUp,BreakendExonDown,ExonsSkippedUp,ExonsSkippedDown)


View(dnaRnaCombinedOutputData)
View(dnaRnaCombinedOutputData %>% filter(is.na(HmfId)))

dnaRnaCombinedOutputData = dnaRnaCombinedOutputData %>% mutate(DnaPosUp=ifelse(DnaPosUp>0,DnaPosUp,''),
                                                               DnaPosDown=ifelse(DnaPosDown>0,DnaPosDown,''),
                                                               DnaOrientUp=ifelse(DnaOrientUp!=0,DnaOrientUp,''),
                                                               DnaOrientDown=ifelse(DnaOrientDown!=0,DnaOrientDown,''))

write.csv(dnaRnaCombinedOutputData,'~/data/sv/rna/LINX_dna_rna_combined_data_tmp.csv', quote = F, row.names = F)
write.csv(dnaRnaCombinedOutputData %>% select(-SampleId,-ClusterId),'~/data/sv/rna/LINX_dna_rna_combined_data.csv', quote = F, row.names = F)


#### 
# Plot Results

dnaRnaCombinedOutputData = read.csv('~/data/sv/rna/LINX_dna_rna_combined_data.csv')

View(dnaRnaCombinedOutputData)


# previously used dnaRnaCombinedData
# create a summary view to plot the precision results 
dnaRnaSummary = dnaRnaCombinedOutputData %>% filter(KnownCategory!='Unknown'&MatchType!='RNA with DNA Support') %>% group_by(MatchType,KnownCategory) %>% count()



View(dnaRnaSummary)
View(dnaRnaSummary %>% spread(MatchType,n))

plotColours4 = c('royal blue','skyblue3','lightblue','khaki4','khaki3','sienna1')

dnaRnaSummaryPlot = (ggplot(dnaRnaSummary , aes(x=KnownCategory, y=n, fill=MatchType))
                     + geom_bar(stat = "identity", colour = "black", position = position_stack(reverse = TRUE))
                     + labs(x = "", y="Fusions", fill='Fusion Prediction', title = "")
                     + scale_fill_manual(values = plotColours4)
                     + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                     + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                     + theme(axis.text.x = element_text(angle=90, hjust=1,size=10))
                     + coord_flip())

print(dnaRnaSummaryPlot)

