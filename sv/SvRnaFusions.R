library(data.table)
library(dplyr)
library(tidyr)
library(stringi)

###############################
# RNA-FUSION MATCHING  ANALYSIS

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
rawRnaData = read.csv('~/data/sv/rna/rna_data_all_samples.csv')
rnaSampleIds = rawRnaData %>% group_by(SampleId) %>% count() 
nrow(rnaSampleIds)
View(rnaSampleIds %>% filter(!(SampleId %in% hpcDedupedSamples$sampleId)))
View(rnaSampleIds %>% filter(SampleId %in% hpcDedupedSamples$sampleId))

rnaMatchData = load_rna_match_data('~/data/sv/rna/SVA_RNA_DATA.csv')

View(rnaMatchData)

write.csv(rnaMatchData,'~/data/sv/rna/SVA_RNA_DATA_hpc_dedup.csv',row.names = F, quote = F)

# restrict to HPC deduped cohort
rnaMatchData = rnaMatchData %>% filter(RnaPhaseMatched=='true')

# filter out unphased RNA fusions for all subsequent analysis
rnaMatchData = rnaMatchData %>% filter(SampleId %in% hpcDedupedSamples$sampleId)

rnaMatchDataBothSVs = annotate_rna_both_svs(rnaMatchData)
View(rnaMatchDataBothSVs)

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

# View(summaryRnaData)

summaryRnaData$FusionType = "AllFusions" # no split by KnownType

rnaCategorySummary1 = summaryRnaData %>% group_by(SvaCategory2,FusionType) %>% count() %>% spread(SvaCategory2,n) %>% 
  arrange(FusionType) %>% ungroup()
rnaCategorySummary1[is.na(rnaCategorySummary1)] = 0

#View(rnaCategorySummary1)
#View(rnaCategorySummary1 %>% select(FusionType,Matched,MatchedExonsSkipped,NotCalled)) # ,NoRnaPhasedFusion

rnaCategorySummaryData1 = rnaCategorySummary1 %>% select(FusionType,Matched,MatchedExonsSkipped,NotCalled)
rnaCategorySummaryData1 = rnaCategorySummaryData1 %>% gather('Category','Count', 2:ncol(rnaCategorySummaryData1)) 
#View(rnaCategorySummaryData1)

rnaCategorySummaryData1 = merge(rnaCategorySummaryData1,catData,by='Category',all.x=T)

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

rnaSampleIds = read.csv('~/data/sv/rna/rna_starfusion_sample_ids.csv')
View(rnaSampleIds)
nrow(rnaSampleIds) # 630 samples

# load all fusions found for the 630 samples with RNA
svaRnaFusions = read.csv('~/data/sv/rna/SVA_FUSIONS.csv')
svaRnaFusions = annotate_fusions(svaRnaFusions)
svaRnaFusions = svaRnaFusions %>% filter(SampleId %in% rnaSampleIds$SampleId)
svaRnaFusions = svaRnaFusions %>% filter(SampleId %in% hpcDedupedSamples$sampleId)
View(svaRnaFusions)

# limit to reported fusions
nrow(svaRnaFusions) # 12711 fusions called by LINX, 11396 in the HPC deduped set
# nrow(svaRnaFusions %>% filter(Reportable=='true'))
# View(svaRnaFusions)

rnaReadData = load_rna_match_data('~/data/sv/rna/read_data/SVA_RNA_DATA.csv')
rnaReadData = rnaReadData %>% mutate(SampleGenePair=paste(SampleId,GeneNameUp,GeneNameDown,sep='_'))


# create a combined file from the RNA and LINX fusions files
View(rnaMatchData)
rnaCombinedData = merge(svaRnaFusions, 
                         rnaMatchData %>% filter(RnaPhaseMatched=='true'),
                         by=c('SampleId','GeneNameUp','GeneNameDown'),all=T)

rnaCombinedData = rnaCombinedData %>% mutate(HasDnaData=!is.na(KnownType.x),
                                             HasRnaData=!is.na(KnownType.y),
                                             SampleGenePair=paste(SampleId,GeneNameUp,GeneNameDown,sep='_'))

rnaCombinedData = rnaCombinedData %>% mutate(HasReadSupport=(HasDnaData&!HasRnaData&SampleGenePair %in% rnaReadData$SampleGenePair))

View(rnaCombinedData)
View(rnaCombinedData %>% filter(!is.na(HasReadSupport)))
colnames(rnaCombinedData)

dnaRnaCombinedData = rnaCombinedData %>% 
  mutate(KnownType=ifelse(!is.na(KnownType.x),as.character(KnownType.x),as.character(KnownType.y)),
         ChrUp=ifelse(!is.na(ChrUp.x),ChrUp.x,ChrDown.y),ChrDown=ifelse(!is.na(ChrDown.x),ChrDown.x,ChrDown.y),
         PosUp=ifelse(!is.na(PosUp.x),PosUp.x,PosUp.y),PosDown=ifelse(!is.na(PosDown.x),PosDown.x,PosDown.y),
         OrientUp=ifelse(!is.na(OrientUp.x),OrientUp.x,OrientUp.y),OrientDown=ifelse(!is.na(OrientDown.x),OrientDown.x,OrientDown.y),
         StrandUp=ifelse(!is.na(StrandUp.x),StrandUp.x,StrandUp.y),StrandDown=ifelse(!is.na(StrandDown.x),StrandDown.x,StrandDown.y),
         RnaPosUp,RnaPosDown,TransValidLocUp,TransViableUp,TransValidLocDown,TransViableDown,
         TransIdUp=ifelse(!is.na(TranscriptUp),as.character(TranscriptUp),as.character(TransIdUp)),
         CodingTypeUp=ifelse(!is.na(CodingTypeUp.x),as.character(CodingTypeUp.x),as.character(CodingTypeUp.y)),
         RegionTypeUp=ifelse(!is.na(RegionTypeUp.x),as.character(RegionTypeUp.x),as.character(RegionTypeUp.y)),
         SameSV=ifelse(!is.na(SameSV.x),SameSV.x,SameSV.y),
         SameCluster=ifelse(is.na(SameCluster),T,SameCluster),SameChain=ifelse(is.na(SameChain),T,SameChain)) %>%
  mutate(Category=ifelse(HasRnaData&!HasDnaData&SvMatchType!='BothSVs','RNA Only',
                  ifelse(HasReadSupport,'DNA with RNA Read Support',
                  ifelse(HasDnaData&!HasRnaData,'DNA Only',
                  ifelse(HasDnaData|(SameCluster&SameChain),'DNA & RNA','RNA with DNA Support')))),
         KnownCategory=known_type_category(KnownType))

write.csv(dnaRnaCombinedData,'~/data/sv/rna/dnaRnaCombinedData_hpc_dedup.csv', quote = F, row.names = F)

View(dnaRnaCombinedData)
View(dnaRnaCombinedData %>% group_by(Category,KnownCategory) %>% count() %>% spread(Category,n))
View(dnaRnaCombinedData %>% group_by(Category,KnownType) %>% count() %>% spread(Category,n))

# create a summary view to plot the precision results
dnaRnaSummary = dnaRnaCombinedData %>% filter(KnownCategory!='Unknown') %>%
  mutate(MatchType=ifelse(Category=='DNA & RNA','DNA & RNA',
                   ifelse(Category=='DNA Only'|Category=='DNA with RNA Read Support','DNA Only','RNA Only'))) %>%
  group_by(MatchType,KnownType) %>% count()

View(dnaRnaSummary)
View(dnaRnaSummary %>% spread(MatchType,n))

# dnaRnaSummary = merge(dnaRnaSummary,catData,by='MatchType',all.x=T)

plotColours4 = c('royal blue','skyblue3','lightblue','khaki4','khaki3','sienna1')

dnaRnaSummaryPlot = (ggplot(dnaRnaSummary, aes(x=KnownType, y=n, fill=MatchType))
                        + geom_bar(stat = "identity", colour = "black", position = position_stack(reverse = TRUE))
                        + labs(x = "", y="Fusion Count", fill='Match Category', title = "DNA vs RNA Fusion Prediction")
                        + scale_fill_manual(values = plotColours4)
                        + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                        + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                        + theme(axis.text.x = element_text(angle=90, hjust=1,size=10))
                        + coord_flip())

## PLOT 2: DNA vs RNA Fusion Prediction

plot(dnaRnaSummaryPlot)



View(dnaRnaCombinedData %>% filter(HasDnaData) %>% group_by(KnownType,Category,RegionTypeUp.x,CodingTypeUp.x) %>% count())
View(dnaRnaCombinedData %>% filter(HasDnaData&(KnownCategory=='Known'|KnownCategory=='Promiscuous')&CodingTypeUp.x=='NonCoding'&HasRnaData))



#####
# previous debug and summary views


# merge SVA and RNA fusion data to check for overlap
rnaFusionMatches = merge(svaRnaFusions, 
                         rnaMatchDataBothSVs %>% filter(RnaPhaseMatched=='true'),
                         by=c('SampleId','GeneNameUp','GeneNameDown'),all.x=T)

rnaFusionMatches = rnaFusionMatches %>% mutate(SampleGenePair=paste(SampleId,GeneNameUp,GeneNameDown,sep='_'))

View(rnaFusionMatches)

unmatchedRnaFusions = rnaFusionMatches %>% filter(is.na(FusionName))
matchedRnaFusions = rnaFusionMatches %>% filter(!is.na(FusionName))
nrow(unmatchedRnaFusions) # 11.8K, 35 for reportable
nrow(matchedRnaFusions) # 822, 53 for reportable


# A. RNA only - where LINX does not find a valid Both SV match
rnaOnlyData = rnaMatchData %>% filter(RnaPhaseMatched=='true'&SvMatchType!='BothSVs') %>% # &KnownType!='None'&KnownType!=''
  select(SampleId,GeneNameUp,GeneNameDown,KnownType) %>% 
  mutate(Category='RNA Only',
         KnownType=known_type_category(KnownType),
         SameSV=F,
         RNADNAMatch='N/A')

# nrow(rnaOnlyData)
View(rnaOnlyData)
View(rnaOnlyData %>% group_by(KnownType) %>% count())

# B. DNA with RNA read support only
rnaReadData = load_rna_match_data('~/data/sv/rna/read_data/SVA_RNA_DATA.csv')
rnaReadData = rnaReadData %>% mutate(SampleGenePair=paste(SampleId,GeneNameUp,GeneNameDown,sep='_'))
View(rnaReadData)

dnaRnaReadData = rnaReadData %>% 
  select(SampleId,GeneNameUp,GeneNameDown,KnownType,SameSV) %>% 
  mutate(Category='DNA with RNA Read Support',
         KnownType=known_type_category(KnownType),
         RNADNAMatch='N/A')

View(dnaRnaReadData)
View(dnaRnaReadData %>% group_by(KnownType) %>% count()) # 6 known, 5 promiscuous

# C. DNA only, no RNA support
dnaOnlyData = unmatchedRnaFusions %>% filter(!(SampleGenePair %in% rnaReadData$SampleGenePair)) %>%
  select(SampleId,GeneNameUp,GeneNameDown,KnownType=KnownType.x,SameSV=SameSV.x) %>% 
  mutate(Category='DNA Only',
         KnownType=known_type_category(KnownType),
         RNADNAMatch='N/A')

View(dnaOnlyData) # 27 in DNA only
View(dnaOnlyData %>% group_by(KnownType) %>% count()) 

# D. DNA and RNA matched
View(matchedRnaFusions)
# rnaMatchDataBothSVs = merge(rnaMatchDataBothSVs,matchedRnaFusions %>% select(SampleGenePair,Reportable),by='SampleGenePair',all.x=T)
View(rnaMatchDataBothSVs)

View(svaRnaFusions %>% group_by(KnownType,RegionTypeUp,CodingTypeUp,ExonsSkipped=ExonsSkippedUp>0|ExonsSkippedDown>0) %>% count())
View(svaRnaFusions %>% group_by(KnownType,RegionTypeUp,CodingTypeUp) %>% count() )

View(rnaMatchDataBothSVs)

dnaRnaMatchedData = rnaMatchDataBothSVs %>% filter(SvMatchType=='BothSVs') %>%
  select(SampleId,GeneNameUp,GeneNameDown,KnownType,SameSV,SameCluster,SameChain,SvViableUp,SvViableDown,ViableFusion) %>% 
  mutate(KnownType=known_type_category(KnownType),
         ValidBreakends=SvViableUp&SvViableDown,
         RNADNAMatch=ifelse(SameCluster&SameChain&ValidBreakends,'Matched',
                     ifelse(SameCluster&SameChain&!ValidBreakends,'MatchedExonsSkipped',
                     ifelse(SameCluster&!SameChain,ifelse(ValidBreakends,'DiffChain','DiffChainExonsSkipped'),'DiffCluster'))),
         Category=ifelse(RNADNAMatch=='Matched'|RNADNAMatch=='MatchedExonsSkipped','DNA & RNA','RNA with DNA Support'))

View(dnaRnaMatchedData)
View(dnaRnaMatchedData %>% group_by(KnownType,Category,RNADNAMatch) %>% count())
View(svaRnaFusions)

dnaRnaMatchedData = dnaRnaMatchedData %>% select(SampleId,GeneNameUp,GeneNameDown,KnownType,SameSV,Category,KnownType,RNADNAMatch)
View(dnaRnaMatchData)

View(dnaRnaMatchedData %>% group_by(Category,KnownType) %>% count())

# bind together and summarise
dnaRnaMatchData = rbind(rnaOnlyData,dnaOnlyData)
dnaRnaMatchData = rbind(dnaRnaMatchData,dnaRnaReadData)
dnaRnaMatchData = rbind(dnaRnaMatchData,dnaRnaMatchedData)
View(dnaRnaMatchData)
View(dnaRnaMatchData %>% group_by(Category,KnownType) %>% count() %>% spread(Category,n))
dnaRnaMatchData$KnownType ="AllTypes"
View(dnaRnaMatchData %>% group_by(Category,KnownType) %>% count() %>% spread(Category,n))


# Plot Key Results

categories = unique(dnaRnaMatchData$Category)
print(categories)
catIndex = c(5,3,2,1,4)
catData = data.frame(cbind(catIndex,categories))
colnames(catData) = c('CatIndex','Category')
catData$CatIndex = as.numeric(catData$CatIndex)
View(catData %>% arrange(CatIndex))
View(catData)

# dnaRnaMatchData = within(dnaRnaMatchData, rm(CatIndex))
dnaRnaMatchData = merge(dnaRnaMatchData,catData,by='Category',all.x=T)

View(dnaRnaMatchData %>% group_by(KnownType,Category,CatIndex) %>% count() %>% arrange(CatIndex))

plotColours4 = c('royal blue','skyblue3','lightblue','khaki4','khaki3','sienna1')

dnaRnaMatchDataPlot2 = (ggplot(dnaRnaMatchData %>% group_by(Category,KnownType,CatIndex) %>% summarise(Count=n()), 
                               aes(x=KnownType, y=Count, fill=reorder(Category,CatIndex)))
                       + geom_bar(stat = "identity", colour = "black", position = position_stack(reverse = TRUE))
                       + labs(x = "Match Category", y="Fusion Count", fill='Match Category', title = "Reportable Fusions - DNA vs RNA Sensitiivity")
                       + scale_fill_manual(values = plotColours4)
                       + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                       + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                       + theme(axis.text.x = element_text(angle=90, hjust=1,size=10))
                       + coord_flip())

plot(dnaRnaMatchDataPlot2)


# write out all input and summary data
write.csv(dnaRnaMatchData, '~/data/sv/rna/dna_rna_match_data.csv', row.names = F, quote = F)



## RNA phase-matched fusion

rnaMatchData = load_rna_match_data('~/data/sv/rna/SVA_RNA_DATA.csv')

View(rnaMatchData %>% group_by(SpliceType,RnaPhaseMatched) %>% count())
View(rnaMatchData %>% group_by(SpliceType,ExonsFoundUp=!is.na(RnaExonRankUp),ExonsFoundDown=!is.na(RnaExonRankDown)) %>% count())
View(rnaMatchData %>% filter(is.na(RnaExonRankUp)|is.na(RnaExonRankDown)))
View(rnaMatchData %>% filter(RnaPhaseMatched=='false'))
View(svaRnaFusions)


# DEBUG only
rnaRawMatchData = read.csv('~/data/sv/rna/SVA_RNA_DATA.csv')


View(svaRnaFusions %>% filter(SampleId=='CPCT02020502T'))
View(rnaMatchData %>% filter(SampleId=='CPCT02020502T'))
View(rnaRawMatchData %>% filter(SampleId=='CPCT02020502T'))

View(rnaCombinedData %>% filter(!HasDnaData&(SameCluster&SameChain)) %>% group_by(KnownType.y) %>% count())
View(rnaCombinedData %>% filter(!HasDnaData&(SameCluster&SameChain)&KnownType.y!='Unknown'&KnownType.y!=''))

View(dnaRnaCombinedData %>% filter(is.na(ChainLinks) &KnownType.y=='Known'&MatchType=='DNA & RNA'))

View(dnaRnaCombinedOutputData)
View(dnaRnaCombinedOutputData %>% filter(is.na(ChainLinks)&MatchType=='RNA Only'))


View(dnaRnaCombinedData %>% filter(GeneNameUp=='TMPRSS2'&GeneNameDown=='ERG'))
View(dnaRnaCombinedOutputData %>% filter(GeneNameUp=='TMPRSS2'&GeneNameDown=='ERG'))
View(svaRnaFusions %>% filter(GeneNameUp=='TMPRSS2'&GeneNameDown=='ERG'))
View(rnaMatchData %>% filter(GeneNameUp=='TMPRSS2'&GeneNameDown=='ERG'))
View(rnaCombinedData %>% filter(GeneNameUp=='TMPRSS2'&GeneNameDown=='ERG'))


View(dnaRnaCombinedData %>% filter(KnownCategory=='Known'&MatchType=='DNA & RNA'))

View(dnaRnaCombinedOutputData %>% filter(KnownCategory=='Known'&MatchType=='DNA & RNA'))

View(dnaRnaCombinedData %>% filter(!(PosUp.x %in% dnaRnaCombinedOutputData$DnaPosUp)))



## previous sensitivty plots

# previous plots

summaryRnaData = summaryRnaData %>% 
  mutate(KnownCategory = ifelse(KnownType=='',ifelse(Proximate,'NoneProximate','NoneNotProximate'),as.character(KnownType)))

summaryRnaData = summaryRnaData %>% mutate(SvaCategory2=ifelse(SvaCategory=='Matched',ifelse(ValidBreakends,'Matched','MatchedExonsSkipped'),
                                                               ifelse(SvaCategory=='DiffChain','DiffChain',
                                                                      ifelse(SvaCategory=='DiffCluster',ifelse(ValidBreakends,'DiffCluster','DiffChainExonsSkipped'),SvaCategory))))

rnaCategorySummary1 = summaryRnaData %>% group_by(SvaCategory2,Proximate) %>% count() %>% spread(SvaCategory2,n) %>% arrange(Proximate) %>% ungroup()
rnaCategorySummary1 = rnaCategorySummary1 %>% mutate(Proximate=ifelse(Proximate,'Proximate','NotProximate'))
rnaCategorySummary1[is.na(rnaCategorySummary1)] = 0
View(rnaCategorySummary1 %>% select(Proximate,Matched,MatchedExonsSkipped,DiffChain,DiffCluster,SingleBEMatch,NoMatch))

rnaCategorySummaryData1 = rnaCategorySummary1 %>% select(Proximate,Matched,MatchedExonsSkipped,DiffChain,DiffCluster,SingleBEMatch,NoMatch)
rnaCategorySummaryData1 = rnaCategorySummaryData1 %>% gather('Category','Count', 2:ncol(rnaCategorySummaryData1)) 
rnaCategorySummaryData1 = rnaCategorySummaryData1 %>% mutate(CountBucket=2**round(log(Count,2)))
View(rnaCategorySummaryData1)

categories = unique(rnaCategorySummaryData1$Category)
# print(categories)
catIndex = c(1,2,3,4,5,6)
catData = data.frame(cbind(catIndex,categories))
colnames(catData) = c('CatIndex','Category')
catData$CatIndex = as.numeric(catData$CatIndex)
View(catData %>% arrange(CatIndex))

rnaCategorySummaryData1 = merge(rnaCategorySummaryData1,catData,by='Category',all.x=T)

plotColours3 = c('royal blue','light blue','orangered','sienna1','khaki4','khaki3','palegreen', 'seagreen')

rnaSummaryDataPlot1 = (ggplot(rnaCategorySummaryData1, aes(x=Proximate, y=Count, fill=reorder(Category,CatIndex)))
                       + geom_bar(stat = "identity", colour = "black", position = position_stack(reverse = TRUE))
                       + labs(x="SV Proximity", y="Fusion Count", fill='Category', title='All RNA Fusions by DNA match type')
                       + scale_fill_manual(values = plotColours3)
                       + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                       + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                       + theme(axis.text.x = element_text(angle=90, hjust=1,size=10))
                       + coord_flip())

plot(rnaSummaryDataPlot1)

# now split chained an unchained up for matched types
rnaCategorySummary2 = summaryRnaData %>% filter(SvaCategory=='Matched') %>% group_by(SameSV) %>% summarise(Count=n())
totalCount = sum(rnaCategorySummary2$Count)
rnaCategorySummary2 = rnaCategorySummary2 %>% mutate(FusionType=ifelse(SameSV,'Single SV Fusion', 'Chained Fusion'), 
                                                     Percent=round(Count/totalCount*100))

View(rnaCategorySummary2)

rnaSummaryDataPlot2 = (ggplot(rnaCategorySummary2, aes(x='', y=Count, fill=FusionType))
                       + geom_bar(stat = "identity", colour = "black", width=1)
                       + labs(x="Fusion Type", y="Fusion Count", title='DNA chained vs single SV for matched RNA fusions')
                       + scale_fill_manual(values = plotColours3)
                       + theme(axis.text.x=element_blank()) 
                       + geom_text(aes(x=1, y=cumsum(Count) - Count/2, label=paste(Percent,' %',sep='')))
                       # + geom_text(aes(y = Count/2 + c(0, cumsum(Count)[-length(Count)]), label = Percent, size=5)
                       # + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                       # + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                       #+ theme(axis.text.x = element_text(angle=90, hjust=1,size=10))
                       # + theme(axis.ticks.y = element_blank())
                       + coord_polar("y", start=0))

plot(rnaSummaryDataPlot2)

rnaSummaryDataPlot2 = (ggplot(rnaCategorySummaryData, aes(x=Type, y=Count, fill=Category))
                       + geom_bar(stat = "identity", colour = "black")
                       + labs(x = "Match Category", y="Fusion Count")
                       + scale_fill_manual(values = plotColours3)
                       # + scale_y_log10()
                       + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                       + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                       + theme(axis.text.x = element_text(angle=90, hjust=1,size=10))
                       + coord_flip())

plot(rnaSummaryDataPlot2)


rnaSummaryDataPlot2 = (ggplot(rnaCategorySummaryData, aes(x=Type, y=Count, fill=Category))
                       + geom_bar(stat = "identity", colour = "black")
                       + labs(x = "Match Category", y="Fusion Count")
                       + scale_fill_manual(values = plotColours3)
                       # + scale_y_log10()
                       + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                       + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                       + theme(axis.text.x = element_text(angle=90, hjust=1,size=10))
                       + coord_flip())

plot(rnaSummaryDataPlot2)


View(svaRnaFusions %>% filter(GeneNameUp=='TMPRSS2'&GeneNameDown=='ERG'&ClusterCount==1))





#############
## Preparation of RNA Read Count support for unmatched DNA fusions

View(rnaFusionMatches)
View(rnaFusionMatches %>% group_by(KnownType.x,Matched=!is.na(FusionName)) %>% count())
View(rnaFusionMatches %>% filter(KnownType.x=='Known'&!is.na(FusionName)))

# mismatches on phasing
View(rnaFusionMatches %>% filter(KnownType.x=='Known'&!is.na(FusionName)) 
     %>% select(SampleId.x,GeneUp.x,GeneDown.x,SvIdUp.x,SvIdDown.x,SvIdUp.y,SvIdDown.y,ViableFusion,PhaseMatched.y,
                TranscriptUp,TranscriptDown,TransIdUp,TransIdDown,SpliceType,
                RnaPosUp,RnaPosDown,CodingTypeUp.x,CodingTypeDown.x,CodingTypeUp.y,CodingTypeDown.y))

View(unmatchedRnaFusions)
View(unmatchedRnaFusions %>% group_by(KnownType.x,SameSV.x) %>% count())
nrow(unmatchedRnaFusions) # 35 SVA fusions unmatched in the RNA by StarFusion
View(unmatchedRnaFusions %>% select(SampleId.x,GeneUp.x,GeneDown.x,ChrUp.x,PosUp.x,ChrDown.x,PosDown.x))
View(unmatchedRnaFusions %>% filter(ChrUp.x!=ChrDown.x) %>% select(SampleId.x,GeneUp.x,GeneDown.x,ChrUp.x,PosUp.x,ChrDown.x,PosDown.x))
View(unmatchedRnaFusions %>% filter(as.character(GeneUp.x)==as.character(GeneDown.x)))

View(matchedRnaFusions %>% group_by(KnownType.x,SameSV.x) %>% count())


path = "/data/experiments/190502_sv_analysis_and_rna_comparison/rna_star_fusion/"
cjFile = "Chimeric.out.junction"
for(i in 1:nrow(unmatchedRnaFusions))
{
  data = unmatchedRnaFusions[i,]
  sampleId = data$SampleId.x
  sampleFile = paste(path,'runs/',sampleId,'/',sampleId,'_fusion/',cjFile,sep='')
  
  genePair = paste(data$GeneUp.x,data$GeneDown.x,sep='_')
  outputFile = paste(path,sampleId,'_',genePair,'_unmatched_read_data.tsv',sep='')
  
  chrUp = as.character(data$ChrUp.x)
  posUp = floor(data$PosUp.x / 1e6)
  chrDown = as.character(data$ChrDown.x)
  posDown = floor(data$PosDown.x / 1e6)
  
  if(chrUp < chrDown | (chrUp == chrDown && posUp <= posDown))
  {
    chrPosStr = sprintf("chr%s.%.0f.*chr%s.%.0f", chrUp, posUp, chrDown, posDown)
  }
  else
  {
    chrPosStr = sprintf("chr%s.%.0f.*chr%s.%.0f", chrDown, posDown, chrUp, posUp)
  }
  
  grepStr = sprintf("grep %s %s > %s", chrPosStr, sampleFile, outputFile)
  print(grepStr)
}

# load the results
# View(unmatchedRnaFusions)
sourcePath = '~/data/sv/rna/read_data/'
readCountData = data.frame(matrix(ncol = 9, nrow = 0))
colnames(readCountData) = c('SampleId','GeneNameUp','GeneNameDown','ChrUp','ChrDown','PosUp','PosDown','OrientUp','OrientDown')
for(i in 1:nrow(unmatchedRnaFusions))
{
  data = unmatchedRnaFusions[i,]
  sampleId = data$SampleId.x

  genePair = paste(data$GeneNameUp.x,data$GeneNameDown.x,sep='_')
  sourceFile = paste(sourcePath,sampleId,'_',genePair,'_unmatched_read_data.tsv',sep='')
  print(sourceFile)
  
  if(!file.exists(sourceFile)) 
  {
    print(paste(sampleId,' has no read count records',sep=''))
  }
  else
  {
    # load TSV
    tsvData = read.csv(sourceFile, sep='\t')
    
    chrUp = as.character(data$ChrUp.x)
    posUp = floor(data$PosUp.x / 1e6)
    chrDown = as.character(data$ChrDown.x)
    posDown = floor(data$PosDown.x / 1e6)
    orientUp = data$OrientUp.x
    orientDown = data$OrientDown.x
    switchPositions = !(chrUp < chrDown | (chrUp == chrDown && posUp <= posDown))
    
    if(nrow(tsvData) > 0)
    {
      print(paste(sourceFile,' has read count: ',nrow(tsvData),sep=''))
      
      for(j in 1:nrow(tsvData))
      {
        index = nrow(readCountData)+1
        readCountData[index,1] = as.character(sampleId)
        readCountData[index,2] = as.character(data$GeneNameUp.x)
        readCountData[index,3] = as.character(data$GeneNameDown.x)
        readCountData[index,4] = chrUp
        readCountData[index,5] = chrDown
        
        if(switchPositions)
        {
          readCountData[index,6] = tsvData[j,5]
          readCountData[index,7] = tsvData[j,2]
        }
        else
        {
          readCountData[index,6] = tsvData[j,2]
          readCountData[index,7] = tsvData[j,5]
        }
        
        readCountData[index,8] = orientUp
        readCountData[index,9] = orientDown
      }
    }
  }
}

View(readCountData)

View(readCountData %>% group_by(SampleId,GeneUp,GeneDown) %>% count())
View(readCountData %>% group_by(SampleId,GeneUp,GeneDown,ChrUp,ChrDown,PosUp,PosDown) %>% count() %>% group_by(SampleId,GeneUp,GeneDown) %>% count())
readCountsByPosition = readCountData %>% group_by(SampleId,GeneNameUp,GeneNameDown,ChrUp,ChrDown,PosUp,PosDown,OrientUp,OrientDown) %>% summarise(JunctionReadCount=n())
nrow(readCountsByPosition)

# test out read counts against SVA fusion data

rnaStarFusionData = read.csv('~/data/sv/rna/rna_data_all_samples.csv')
View(rnaStarFusionData)
colnames(rnaStarFusionData)
rm(rnaStarFusionData)

ensemblGeneData = read.csv('~/data/sv/ensembl_gene_data.csv')
View(ensemblGeneData)

readCountFusionData = merge(readCountsByPosition,ensemblGeneData %>% select(GeneName,GeneIdUp=GeneId),by.x='GeneNameUp',by.y='GeneName',all.x=T)
readCountFusionData = merge(readCountFusionData,ensemblGeneData %>% select(GeneName,GeneIdDown=GeneId),by.x='GeneNameDown',by.y='GeneName',all.x=T)
View(readCountFusionData)
View(readCountFusionData %>% filter(is.na(GeneIdUp)|is.na(GeneIdDown)))

# convert to a form similar to the rnaStartFusion input for SVA
readCountFusionData = readCountFusionData %>% 
  mutate(FusionName=paste(GeneNameUp,GeneNameDown,sep='--'),
         SpanningFragCount=0,
         SpliceType='UNKNOWN',
         JunctionReads='',SpanningFrags='',LargeAnchorSupport='',FFPM=0,LeftBreakDinuc='',LeftBreakEntropy='',RightBreakDinuc='',RightBreakEntropy='',annots='')

View(readCountFusionData)
# colnames(readCountFusionData)
View(readCountFusionData %>% group_by(SampleId,GeneNameUp,GeneNameDown,ChrUp,ChrDown,PosUp,PosDown,OrientUp,OrientDown) %>% count())

write.csv(readCountFusionData %>% select(SampleId,FusionName,JunctionReadCount,SpanningFragCount,SpliceType,
                                         GeneNameUp,GeneIdUp,ChrUp,PosUp,OrientUp,
                                         GeneNameDown,GeneIdDown,ChrDown,PosDown,OrientDown,
                                         JunctionReads,SpanningFrags,LargeAnchorSupport,FFPM,LeftBreakDinuc,LeftBreakEntropy,RightBreakDinuc,RightBreakEntropy,annots),
          '~/data/sv/rna/rna_read_data_samples.csv', row.names = F, quote=F)


rnaReadDataRaw = read.csv('~/data/sv/rna/read_data/SVA_RNA_DATA.csv', stringsAsFactors = F)
View(rnaReadDataRaw)
View(rnaReadDataRaw %>% group_by(SampleId,GeneUp,GeneDown,ChrUp,ChrDown,PosUp,PosDown,OrientUp,OrientDown) %>% count())
View(rnaReadDataRaw %>% group_by(SampleId,GeneUp,GeneDown) %>% count())
rnaReadData = load_rna_match_data('~/data/sv/rna/read_data/SVA_RNA_DATA.csv')
View(rnaReadData)


View(rnaReadDataRaw %>% filter(GeneUp==GeneDown))
View(rnaReadDataRaw %>% filter(GeneUp==GeneDown) %>% group_by(SampleId,GeneUp,GeneDown) %>% count())




#############
## Misc Analysis


# skipped exons
View(svaRnaFusions %>% filter(PhasingExonsSkipped>0) %>% select(SampleId,GeneUp,GeneDown,TranscriptUp,TranscriptDown,ExonUp,ExonDown,PhaseUp,PhaseDown,PhaseMatched,
                                                                PhasingExonsSkipped,SvIdUp,SvIdDown,ChrUp,ChrDown,PosUp,PosDown,OrientUp,OrientDown))





#############
## DEBUG ONLY

rnaClusters = read.csv('~/data/sv/rna/SVA_CLUSTERS.csv')
rnaClusters = rnaClusters %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))
nrow(rnaClusters)

rnaMatchDataBothSVs = rnaMatchDataBothSVs %>% mutate(SampleClusterIdUp=paste(SampleId,ClusterIdUp,sep='_'))
rnaMatchDataBothSVs = rnaMatchDataBothSVs %>% mutate(SampleClusterIdDown=paste(SampleId,ClusterIdDown,sep='_'))

rnaMatchDataBothSVs = merge(rnaMatchDataBothSVs,rnaClusters %>% select(SampleClusterId,ResolvedType),
                            by.x='SampleClusterIdUp',by.y='SampleClusterId',all.x=T)

rnaMatchDataBothSVs = merge(rnaMatchDataBothSVs,rnaClusters %>% select(SampleClusterId,ResolvedTypeDown=ResolvedType),
                            by.x='SampleClusterIdDown',by.y='SampleClusterId',all.x=T)

View(rnaMatchDataBothSVs)



tmp = read.csv('~/data/sv/rna/read_data/CPCT02010440T_BRAF_BRAF_unmatched_read_data.tsv',sep='\t')
View(tmp)
info = file.info('~/data/sv/rna/read_data/CPCT02010389T_AGBL1-AS1_NTRK3_unmatched_read_data.tsv')
View(info)
empty = rownames(info[info$size == 0, ])



View(svaRnaFusions %>% group_by(SampleId,GeneUp,GeneDown,KnownType) %>% count())
View(svaRnaFusions %>% group_by(SampleId,GeneUp,GeneDown,KnownType) %>% count() %>% spread(KnownType,n))

validSvaRnaFusions = svaRnaFusions %>% filter((InChain&ValidChain)|(!InChain&NonDisruptedSingle))
rm(svaRnaFusions)
nrow(validSvaRnaFusions)
nrow(rnaMatchDataBothSVs)
nrow(rnaMatchDataBothSVs %>% filter(SvValidLocUp,SvValidLocDown))
nrow(rnaMatchDataBothSVs %>% filter(SvViableUp&SvViableDown))
nrow(rnaMatchDataBothSVs %>% filter(ViableFusion=='true'))


sampleRnaSvaFusions = (validSvaRnaFusions %>% group_by(SampleId,GeneUp,GeneDown) 
                       %>% summarise(Count=n(),
                                     SimpleSVCount=sum(SameSV&ClusterCount==1),
                                     SingleSVUnchainedCount=sum(SameSV&ClusterCount>1&!InChain),
                                     SingleSVChainedCount=sum(SameSV&InChain),
                                     MultiSVChainedCount=sum(!SameSV&ValidChain))
                       %>% mutate(FusionType=ifelse(MultiSVChainedCount==Count,'MultiSV',
                                             ifelse(SimpleSVCount==Count,'SimpleSV',
                                             ifelse(SingleSVUnchainedCount==Count,'SingleSVUnchained',
                                             ifelse(SingleSVChainedCount==Count,'SingleSVChained','Unclear'))))))

nrow(sampleRnaSvaFusions)
View(sampleRnaSvaFusions)


# merge 2 data sets on gene fusion
validSvaRnaFusionsByFusion = validSvaRnaFusions %>% group_by(SampleId,GeneUp,GeneDown) %>% summarise(SvaCount=n())
rnaMatchDataByFusion = rnaMatchDataBothSVs %>% group_by(SampleId,GeneUp,GeneDown) %>% summarise(RnaCount=n())
rnaMatchVsSvaByFusion = merge(rnaMatchDataByFusion,validSvaRnaFusionsByFusion,by=c('SampleId','GeneUp','GeneDown'),all.x=T)
View(rnaMatchVsSvaByFusion %>% group_by(SvaMatched=!is.na(SvaCount)) %>% summarise(Matched=n(),Perc=round(n()/nrow(rnaMatchVsSvaByFusion),3)))
svaMatchVsRnaByFusion = merge(rnaMatchDataByFusion,validSvaRnaFusionsByFusion,by=c('SampleId','GeneUp','GeneDown'),all.y=T)
View(svaMatchVsRnaByFusion)
View(svaMatchVsRnaByFusion %>% group_by(RnaMatched=!is.na(RnaCount)) %>% summarise(Matched=n(),Perc=round(n()/nrow(svaMatchVsRnaByFusion),3)))


View(rnaMatchDataBothSVs)
sampleRnaFusions = rnaMatchDataBothSVs %>% group_by(SampleId,GeneUp,GeneDown) 



# 1. Genuine Misses out of 1294
nrow(rnaMatchData %>% filter(SvMatchType!='BothSVs',Proximate==T))  #PROXIMATE   148
nrow(rnaMatchData %>% filter(SvMatchType=='NoSV',Proximate==F) %>% group_by(SampleId) %>% count())   # 90
nrow(rnaMatchData %>% filter(SvMatchType=='SingleSV',Proximate==F))  #One SV FOUND   118

# 2.  BothSVs found - clustering and chaining summary
View(rnaMatchDataBothSVs %>% group_by(ViableBothSides=SvViableUp&SvViableDown,SameSV,SameCluster,SameChain) %>% count() %>% spread(ViableBothSides,n) %>% arrange(-`TRUE`))

# negative distances to upstream gene
View(rnaMatchDataBothSVs %>% filter(SvViableUp&SvViableDown))
View(rnaMatchDataBothSVs %>% filter(DistancePrevDown<0))

View(rnaMatchDataBothSVs %>% group_by(DownPrevDistance=ifelse(DistancePrevDown>0,'Pos','Neg'),BothViable=(SvViableUp&SvViableDown)) %>% count())


# 3. Chaining analysis
View(rnaMatchDataBothSVs %>% group_by(ClusterCountUp) %>% count())

View(rnaMatchDataBothSVs %>% filter(SvViableUp&SvViableDown&!SameSV&SameCluster&SameChain) 
     %>% group_by(FacingInChain,ChainLinks) %>% count() %>% arrange(FacingInChain,ChainLinks))

View(rnaMatchDataBothSVs %>% filter(!SvViableUp|!SvViableDown) %>% group_by(SvViableUp,SvViableDown,SvValidLocUp,SvValidLocDown,SameSV,SameCluster,SameChain) %>% count())



# merge 2 data sets on SV Id
validSvaRnaFusionsById = validSvaRnaFusions %>% group_by(SvIdUp,SvIdDown) %>% summarise(SvaCount=n())
rnaMatchDataById = rnaMatchDataBothSVs %>% group_by(SvIdUp,SvIdDown) %>% count()
rnaMatchVsSvaById = merge(rnaMatchDataById,validSvaRnaFusionsById,by=c('SvIdUp','SvIdDown'),all.x=T)
View(rnaMatchVsSvaById)
View(rnaMatchVsSvaById %>% group_by(SvaMatched=!is.na(SvaCount)) %>% summarise(Matched=n(),Perc=round(n()/nrow(rnaMatchDataById),3)))
svaMatchVsRnaById = merge(rnaMatchDataById,validSvaRnaFusionsById,by=c('SvIdUp','SvIdDown'),all.y=T)
View(svaMatchVsRnaById)
View(svaMatchVsRnaById %>% group_by(RnaMatched=!is.na(RnaCount)) %>% summarise(Matched=n(),Perc=round(n()/nrow(svaMatchVsRnaById),3)))



# 5. Analysis of RNA fusions only matching at one end
View(rnaMatchData %>% filter(SvMatchType=='SingleSV'))
svRnaData = read.csv('~/data/sv/rna/SVA_SVS.csv')

rnaSingleMatches = rnaMatchData %>% filter(SvMatchType=='SingleSV')

rnaSingleMatchUpSvData = merge(rnaSingleMatches %>% filter(!is.na(SvIdUp)), svRnaData, by.x=c('SampleId','SvIdUp'), by.y=c('SampleId','Id'), all.x=T)

View(rnaSingleMatchUpSvData %>% filter((ChrUp==ChrDown&TypeUp!='BND')|(ChrUp!=ChrDown&TypeUp=='BND')) %>% filter(TypeUp!='SGL'&TypeUp!='NONE')
     %>% select(SampleId,GeneUp,GeneDown,SvIdUp,ChrUp,PosUp,RnaPosUp,OrientUp,StrandUp,TypeUp,TransViableUp,TransValidLocUp,
                ChrDown,RnaPosDown,StrandDown,
                ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd,ClusterCount,ClusterDesc,GeneStart,GeneEnd,
                RegionTypeUp,CodingTypeUp))

rnaSingleMatchDownSvData = merge(rnaSingleMatches %>% filter(!is.na(SvIdDown)), svRnaData, by.x=c('SampleId','SvIdDown'), by.y=c('SampleId','Id'), all.x=T)

View(rnaSingleMatchDownSvData %>% filter((ChrUp==ChrDown&TypeDown!='BND')|(ChrUp!=ChrDown&TypeDown=='BND')) %>% filter(TypeDown!='SGL'&TypeDown!='NONE')
     %>% mutate(OtherChr=ifelse(PosDown==PosStart,ChrEnd,ChrStart),
                OtherPos=ifelse(PosDown==PosStart,PosEnd,PosStart),
                OtherOrient=ifelse(PosDown==PosStart,OrientEnd,OrientStart),
                Distance=abs(RnaPosUp-OtherPos),
                CorrectLoc=ifelse((StrandUp==1&OtherPos>=RnaPosUp)|(StrandUp==-1&OtherPos<=RnaPosUp),T,F))
     %>% select(SampleId,GeneUp,GeneDown,SvIdDown,ChrDown,PosDown,RnaPosDown,OrientDown,StrandDown,TypeDown,TransViableDown,TransValidLocDown,
                ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd,ClusterCount,ClusterDesc,GeneStart,GeneEnd,
                ChrUp,RnaPosUp,StrandUp,OtherChr,OtherPos,OtherOrient,Distance,CorrectLoc,DistancePrevDown))

colnames(rnaSingleMatchUpSvData)


View(svRnaData %>% filter(SampleId=='CPCT02020767T'&(ChrStart==2&ChrEnd==2)))


# 6. Summary results

# Total Found in RNA (after deduplication)	
# Called in SVA
# Clustered breakends found but not chained
# Breakends found for both partners, not clustered 	
# Breakend missing for 1 partner	
# Breakend Missing for both partners

# split by Known, 5’ Promiscuous, 3’ Promiscuous, Other Proximate, Other Non Proximate

View(rnaMatchDataBothSVs)

summaryBothData = rnaMatchDataBothSVs %>% 
  mutate(SvaCategory=ifelse(SameCluster&SameChain,'Matched',ifelse(SameCluster&!SameChain,'DiffChain','DiffCluster')),
         ValidBreakends=SvViableUp&SvViableDown,
         ViableFusion=ViableFusion=='true',
         PhaseMatched=PhaseMatched=='true')

View(summaryBothData)
View(summaryBothData %>% group_by(SvaCategory) %>% count())
View(summaryBothData %>% group_by(ViableFusion,ValidBreakends,PhaseMatched) %>% count())
View(summaryBothData %>% filter(KnownType=='Known'))

summaryNotBothData = rnaMatchData %>% filter(SvMatchType!='BothSVs') %>% 
  mutate(SvaCategory = ifelse(SvMatchType=='SingleSV','SingleBEMatch','NoMatch'),
         ValidBreakends=F,
         ViableFusion=F,
         PhaseMatched=F)

summaryRnaData = rbind(summaryBothData %>% select(KnownType,SvaCategory,Proximate,ValidBreakends,ViableFusion,PhaseMatched),
                       summaryNotBothData %>% select(KnownType,SvaCategory,Proximate,ValidBreakends,ViableFusion,PhaseMatched))


summaryRnaData = summaryRnaData %>% 
  mutate(KnownCategory = ifelse(KnownType=='',ifelse(Proximate,'NoneProximate','NoneNotProximate'),as.character(KnownType)))

# View(summaryRnaData)

summaryRnaData = summaryRnaData %>% mutate(SvaCategory2=ifelse(SvaCategory=='Matched',ifelse(ValidBreakends,'Matched','MatchedExonsSkipped'),
                                                               ifelse(SvaCategory=='DiffChain',ifelse(ValidBreakends,'DiffChain','DiffChainExonsSkipped'),
                                                                      ifelse(SvaCategory=='DiffCluster',ifelse(ValidBreakends,'DiffCluster','DiffChainExonsSkipped'),SvaCategory))))

rnaCategorySummary3 = summaryRnaData %>% group_by(SvaCategory2,KnownCategory) %>% count() %>% spread(SvaCategory2,n) %>% arrange(KnownCategory)
rnaCategorySummary3[is.na(rnaCategorySummary3)] = 0
View(rnaCategorySummary3)

rnaCategorySummary3Prev = rnaCategorySummary3
View(rnaCategorySummary3Prev)

rnaCategorySummary = summaryRnaData %>% group_by(SvaCategory,KnownCategory) %>% count() %>% spread(SvaCategory,n)
rnaCategorySummary[is.na(rnaCategorySummary)] = 0
View(rnaCategorySummary)

rnaCategorySummary2 = summaryRnaData %>% group_by(ValidBreakends,SvaCategory,KnownCategory) %>% count() %>% spread(SvaCategory,n) %>% 
  arrange(!ValidBreakends,KnownCategory)
rnaCategorySummary2[is.na(rnaCategorySummary2)] = 0
View(rnaCategorySummary2)






######################
# RNA with STAR FUSION

load("~/data/r_data/highestPurityCohortSummary.RData")
load('~/Dropbox/HMF Australia team folder/Paper_Analyses/Frozen RData/Processed/hpcFusions.RData')
View(highestPurityCohortSummary)

View(rnaSamples %>% filter(SampleId %in% highestPurityCohortSummary$sampleId))

hpcSamplesWithRna = rnaSamples %>% filter(SampleId %in% hpcFusions$sampleId)
View(hpcSamplesWithRna)

# View(rnaSamples %>% filter(SampleId %in% hpcFusions$sampleId) %>% filter(SampleId %in% prodFusions$SampleId))
write.csv(hpcSamplesWithRna, '~/data/sv/driver_paper_rna_samples.csv', row.names = F, quote = F)


prepare_star_fusion_file<-function(filename, sampleId)
{
  starFusionData = read.csv(filename, sep='\t')
  
  if(nrow(starFusionData) == 0)
    return (starFusionData)
  
  names(starFusionData)[1] = 'FusionName'
  starFusionData = starFusionData %>% separate(LeftGene,c('GeneNameUp','GeneIdUp'),sep='\\^')
  starFusionData = starFusionData %>% separate(RightGene,c('GeneNameDown','GeneIdDown'),sep='\\^')
  
  starFusionData = starFusionData %>% separate(LeftBreakpoint,c('ChrUp','PosUp','OrientUp'),remove=F)
  starFusionData$OrientUp = ifelse(grepl('-',starFusionData$LeftBreakpoint),'-1','1')
  starFusionData$ChrUp = stri_replace_all_fixed(starFusionData$ChrUp,'chr','')
  starFusionData = within(starFusionData, rm(LeftBreakpoint))
  
  starFusionData = starFusionData %>% separate(RightBreakpoint,c('ChrDown','PosDown','OrientDown'),remove=F)
  starFusionData$OrientDown = ifelse(grepl('-',starFusionData$RightBreakpoint),'-1','1')
  starFusionData$ChrDown = stri_replace_all_fixed(starFusionData$ChrDown,'chr','')
  starFusionData = within(starFusionData, rm(RightBreakpoint))
  
  starFusionData$SampleId = sampleId
  colCount = ncol(starFusionData)
  starFusionData = cbind(starFusionData[,colCount],starFusionData[,1:colCount-1])
  names(starFusionData)[1] = 'SampleId'
  
  return (starFusionData)
}

load_rna_fusions<-function(sampleIds, sourceDir)
{
  print(paste('loading ', nrow(sampleIds), ' sample files', sep=''))
  
  starFusionData = data.frame()
  
  for(i in 1:nrow(sampleIds))
  {
    # filename = "~/data/sv.log"
    sampleData = sampleIds[i,]
    sampleId = sampleData[1]
    filename = paste(sourceDir, sampleId, '_star_fusion_predictions.tsv',sep='')
    
    # print(paste('processing sample: ', sampleId, " and file: ", filename, sep=''))
    
    if(!file.exists(filename)) 
    {
      print(paste('missing file: ', filename, sep=''))
    }
    else
    {
      print(paste('loading file: ', filename, sep=''))
      
      newData = prepare_star_fusion_file(filename,sampleId)
      
      if(nrow(newData) > 0)
      {
        newData$JunctionReads = stri_replace_all_fixed(newData$JunctionReads, ',', ';')
        newData$SpanningFrags = stri_replace_all_fixed(newData$SpanningFrags, ',', ';')
        newData$annots = stri_replace_all_fixed(newData$annots, ',', ';')
        starFusionData = rbind(starFusionData,newData)
      }
    }
  }
  
  return (starFusionData)
}

samplesWithRna = read.csv('~/data/sv/rna/rna_samples_01.csv', header=F)
samplesWithRna = rbind(samplesWithRna, read.csv('~/data/sv/rna/rna_samples_02.csv', header=F))
samplesWithRna = rbind(samplesWithRna, read.csv('~/data/sv/rna/rna_samples_03.csv', header=F))
samplesWithRna = rbind(samplesWithRna, read.csv('~/data/sv/rna/rna_samples_04.csv', header=F))
samplesWithRna = rbind(samplesWithRna, read.csv('~/data/sv/rna/rna_samples_05.csv', header=F))
samplesWithRna = rbind(samplesWithRna, read.csv('~/data/sv/rna/rna_samples_06.csv', header=F))
samplesWithRna = rbind(samplesWithRna, read.csv('~/data/sv/rna/rna_samples_07.csv', header=F))
samplesWithRna = rbind(samplesWithRna, read.csv('~/data/sv/rna/rna_samples_08.csv', header=F))
samplesWithRna = rbind(samplesWithRna, read.csv('~/data/sv/rna/rna_samples_09.csv', header=F))
samplesWithRna = rbind(samplesWithRna, read.csv('~/data/sv/rna/rna_samples_10.csv', header=F))
samplesWithRna = rbind(samplesWithRna, read.csv('~/data/sv/rna/rna_samples_11.csv', header=F))
samplesWithRna = rbind(samplesWithRna, read.csv('~/data/sv/rna/rna_samples_12.csv', header=F))
samplesWithRna = rbind(samplesWithRna, read.csv('~/data/sv/rna/rna_samples_13.csv', header=F))
samplesWithRna = rbind(samplesWithRna, read.csv('~/data/sv/rna/rna_samples_14.csv', header=F))
samplesWithRna = rbind(samplesWithRna, read.csv('~/data/sv/rna/rna_samples_16.csv', header=F))
colnames(samplesWithRna) = c('SampleId')
View(samplesWithRna)

write.csv(samplesWithRna, '~/data/sv/rna/rna_starfusion_sample_ids.csv', row.names = F, quote = F)


sourceDir = '~/data/sv/rna/runs/'
View(samplesWithRna)
nrow(samplesWithRna)
allRnaData = load_rna_fusions(samplesWithRna, sourceDir)
View(allRnaData)
nrow(allRnaData)

write.csv(allRnaData, '~/data/sv/rna/rna_data_all_samples.csv', row.names = F, quote = F)

rm(allRnaData)

# adding new samples:
allRnaData = read.csv('~/data/sv/rna/rna_data_all_samples.csv')

newSamplesWithRna = read.csv('~/data/sv/rna/rna_samples_14.csv', header=F)
colnames(newSamplesWithRna) = c('SampleId')
View(newSamplesWithRna)

sourceDir = '~/data/sv/rna/runs/'
newRnaData = load_rna_fusions(newSamplesWithRna, sourceDir)
View(newRnaData)
allRnaData = rbind(allRnaData,newRnaData)
write.csv(allRnaData, '~/data/sv/rna/rna_data_all_samples.csv', row.names = F, quote = F)

rm(allRnaData)



####
# TEMP

rnaMatchData = read.csv('~/data/sv/rna/SVA_RNA_DATA.csv')

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

nrow(rnaMatchData)

nrow(rnaMatchData %>% group_by(FusionName) %>% 
       filter(!((sum(ifelse(SvMatchType=='NoSV',1,0))>1&(sum(ifelse(SvMatchType=='BothSVs',1,0))<1)))))

oldFilters = rnaMatchData %>% group_by(FusionName) %>% 
  filter(((sum(ifelse(SvMatchType=='NoSV',1,0))>1&(sum(ifelse(SvMatchType=='BothSVs',1,0))<1))))

nrow(rnaMatchData %>% group_by(FusionName) %>% 
       filter(!((sum(ifelse(SvMatchType=='NoSV',1,0))>1&sum(ifelse(SvMatchType=='BothSVs',1,0))<pmax(1,sum(ifelse(SvMatchType=='NoSV',1,0))-1)))))

newFilters = rnaMatchData %>% group_by(FusionName) %>% 
  filter(((sum(ifelse(SvMatchType=='NoSV',1,0))>1&sum(ifelse(SvMatchType=='BothSVs',1,0))<pmax(1,sum(ifelse(SvMatchType=='NoSV',1,0))-1))))

View(oldFilters %>% group_by(FusionName,SvMatchType) %>% count() %>% spread(SvMatchType,n))
View(oldFilters)
View(newFilters %>% group_by(FusionName,SvMatchType) %>% count() %>% spread(SvMatchType,n))

# PON Filtering
# PON: Remove recurrent fusions where NONE have SV support at both ends
rnaMatchData = rnaMatchData %>% group_by(FusionName) %>% 
  filter(!((sum(ifelse(SvMatchType=='NoSV',1,0))>1&(sum(ifelse(SvMatchType=='BothSVs',1,0))<1)))) %>% ungroup()

rnaMatchData = rnaMatchData %>% group_by(FusionName) %>% 
  filter(!((sum(ifelse(SvMatchType=='NoSV',1,0))>1&sum(ifelse(SvMatchType=='BothSVs',1,0))<pmax(1,sum(ifelse(SvMatchType=='NoSV',1,0))-1))))%>% ungroup()



