

annotate_dna_rna_match<-function(dnaRnaMatchData)
{
  dnaRnaMatchData = dnaRnaMatchData %>% 
    mutate(SameChr=as.character(ChrUp)==as.character(ChrDown),
           FusionDistance=ifelse(SameChr,abs(RnaPosDown-RnaPosUp),0),
           SVsMatched=ifelse(SvIdUp>0&SvIdDown>0,'BOTH',ifelse(SvIdUp>0|SvIdDown>0,'ONE','NONE')),
           DiffSVs=SvIdUp>0&SvIdDown>0&SvIdUp!=SvIdDown,
           DnaFusionMatched=(DnaFusionMatchType=='GENES'|DnaFusionMatchType=='SVS'),
           GeneType=ifelse(KnownType=='KNOWN_PAIR','KNOWN',ifelse(grepl('PROM',KnownType),'PROM','OTHER')))

  # cluster info
  dnaRnaMatchData = dnaRnaMatchData %>% separate(ClusterInfoUp,c('ClusterIdUp','ClusterCountUp','ChainIdUp','ChainCountUp'),sep = ';')
  dnaRnaMatchData = dnaRnaMatchData %>% separate(ClusterInfoDown,c('ClusterIdDown','ClusterCountDown','ChainIdDown','ChainCountDown'),sep = ';')
  dnaRnaMatchData = dnaRnaMatchData %>% mutate(DiffClusters=SvIdUp>0&SvIdDown>0&ClusterIdUp!=ClusterIdDown,
                                               DiffChains=SvIdUp>0&SvIdDown>0&ChainIdUp!=ChainIdDown)
  
  # chain info
  dnaRnaMatchData = dnaRnaMatchData %>% separate(ChainInfo,c('ChainLinks','ChainLength'),sep = ';') %>%
    mutate(ChainLength=as.numeric(ChainLength),ChainLinks=as.numeric(ChainLinks))
  
  dnaRnaMatchData = dnaRnaMatchData %>% mutate(ValidSVs=SVsMatched=='BOTH'&!DiffSVs&!DiffClusters&!DiffChains&ChainLength<1e5)
  
  return (dnaRnaMatchData)
}

## Load results from Linx's fusion DNA-RNA matching routine
# dnaRnaMatch122 = read.csv('~/data/rna/fusions/LNX_RNA_FUSION_MATCH_ISOFOX.csv')
#dnaRnaMatch = read.csv('~/data/rna/fusions/LNX_RNA_FUSION_MATCH.csv')

## Arriba and Isofox matched
dnaRnaMatch122 = read.csv('~/data/rna/fusions/LNX_RNA_FUSION_MATCH_MATCHED.csv')
dnaRnaMatch122 = annotate_dna_rna_match(dnaRnaMatch122)
View(dnaRnaMatch122)

write.csv(dnaRnaMatch122,'~/data/rna/fusions/dna_rna_fusion_validation_122.csv',row.names = F,quote = F)

dnaRnaMatch2131 = read.csv('~/data/rna/cohort/LNX_RNA_FUSION_MATCH_2131.csv')
dnaRnaMatch2131 = annotate_dna_rna_match(dnaRnaMatch2131)
View(dnaRnaMatch2131)


write.csv(dnaRnaMatch2131,'~/data/rna/fusions/dna_rna_fusion_validation_2131.csv',row.names = F,quote = F)

View(dnaRnaMatch2131 %>% filter(KnownType!='NONE'))
write.csv(dnaRnaMatch2131 %>% filter(KnownType!='NONE'),'~/data/rna/fusions/dna_rna_fusion_validation_known_2131.csv',row.names = F,quote = F)

dnaRnaMatch = dnaRnaMatch2131
# dnaRnaMatch = dnaRnaMatch122

sampleCancerTypes = read.csv('~/data/cancer_types_5200.csv')
dnaRnaMatch = merge(dnaRnaMatch,sampleCancerTypes,by='SampleId',all.x=T)

# Validation - no duplicates
View(dnaRnaMatch)

View(dnaRnaMatch %>% filter(Source=='MATCH') %>% group_by(SampleId,RnaPosUp,RnaPosDown,GeneNameUp,GeneNameDown) %>% count %>% filter(n>1))
View(dnaRnaMatch %>% group_by(FusionId,SampleId) %>% mutate(Dups=n()) %>% ungroup() %>% filter(Dups>1,FusionId!='NONE')) 

View(dnaRnaMatch %>% group_by(KnownType,GeneNameUp,GeneNameDown) %>% count) 

View(dnaRnaMatch %>% filter(KnownType!='Known',RnaJuncUp=="KNOWN",RnaJuncDown=="KNOWN",Source=='MATCH_ISF_UNFILTERED') %>%
       mutate(AF=as.numeric(substr(RnaOtherData,22,26))) %>%
       filter(!grepl('anchorDistance=1',RnaOtherData)|DiscordantFrags>0, # minAnchor pass
              RnaCohortCount<2, # passes cohort frequency
              AF>=0.005,  #passesAF
              JunctionFrags+DiscordantFrags>1 #passes minimumn support
       ) %>%
       select(SampleId,Source,GeneNameUp,GeneNameDown,RnaJuncUp,RnaJuncDown,JunctionFrags,DiscordantFrags,RnaCohortCount,RnaOtherData,AF,everything()))

# Summary of RNA support in DNA
View(dnaRnaMatch %>% filter(KnownType=='KNOWN_PAIR') %>% group_by(SampleId,FusionName) %>% mutate(TransCount=n()) %>% ungroup() %>%
       group_by(DnaFusionMatchType,Source,TransCount) %>% count %>% spread(Source,n))

View(dnaRnaMatch %>% mutate(Type=ifelse(RnaSvType=='BND','BND',ifelse(FusionDistance<1e5,paste(RnaSvType,'SHORT',sep='_'),paste(RnaSvType,'LONG',sep='_')))) %>% 
       group_by(Source,ValidSVs,GeneType,Type) %>% count() %>% spread(Type,n,fill=0))

View(dnaRnaMatch %>% mutate(Type=ifelse(RnaSvType=='BND','BND',ifelse(FusionDistance<1e5,paste(RnaSvType,'SHORT',sep='_'),paste(RnaSvType,'LONG',sep='_')))) %>% 
       group_by(Source,SVsMatched=(SVsMatched=='BOTH'),GeneType,Type) %>% count() %>% spread(Type,n,fill=0))

# level of juncton support
View(dnaRnaMatch %>% filter(Source=='MATCH'|Source=='ISOFOX_ONLY') %>% filter(JunctionFrags+DiscordantFrags>=2) %>%
       group_by(Source,KnownType,FragSupport=2**round(log(JunctionFrags+DiscordantFrags,2))) %>% count %>% spread(KnownType,n,fill=0))

View(dnaRnaMatch %>% group_by(SVsMatched,RnaPhaseMatched,DnaFusionMatched) %>% count())
View(dnaRnaMatch %>% group_by(SVsMatched,RnaPhaseMatched,DnaFusionMatched) %>% count())

View(dnaRnaMatch %>% group_by(SVsMatched,KnownType) %>% count())
View(dnaRnaMatch %>% group_by(SVsMatched,GeneType) %>% count())
View(dnaRnaMatch %>% group_by(ValidSVs,GeneType) %>% count())
View(dnaRnaMatch %>% group_by(ValidSVs,GeneType,RnaSvType) %>% count() %>% spread(RnaSvType,n,fill=0))

View(dnaRnaMatch %>% group_by(Source,ValidSVs,GeneType,RnaSvType) %>% count() %>% spread(RnaSvType,n,fill=0))

## joining with gene expression
dnaFusionGeneExp = read.csv('~/data/rna/fusions/isofox_linx_up_genes.sample_gene_perc_data.csv')
View(dnaFusionGeneExp)

dnaRnaMatch = merge(dnaRnaMatch,dnaFusionGeneExp %>% select(SampleId,GeneIdUp=GeneId,TpmUp=TPM),by=c('SampleId','GeneIdUp'),all.x=T)
View(dnaRnaMatch %>% select(SampleId,GeneIdUp,GeneNameUp,TpmUp,KnownType,RnaJuncUp,everything()))
View(dnaRnaMatch %>% filter(GeneIdUp==''))

dnaRnaMatch = dnaRnaMatch %>% mutate(TpmBucket=ifelse(is.na(TpmUp)|TpmUp<0.1,0,10**round(log(TpmUp,10))))
View(dnaRnaMatch %>% group_by(TpmBucket) %>% count)

View(dnaRnaMatch %>% group_by(KnownType,TpmBucket,DnaFusionMatched) %>% count %>% spread(DnaFusionMatched,n,fill=0))


## Investigation into RNA fusions with no DNA support
noDnaFusions = dnaRnaMatch %>% filter(KnownType=='KNOWN_PAIR'&DnaFusionMatchType=='NONE')

noDnaFusions = merge(noDnaFusions,
                     ensemblGeneData %>% mutate(GeneIdUp=GeneId,GeneStartUp=GeneStart,GeneEndUp=GeneEnd,GeneLengthUp=GeneEnd-GeneStart) %>% 
                       select(GeneIdUp,GeneStartUp,GeneEndUp,GeneLengthUp),by='GeneIdUp',all.x=T)

noDnaFusions = merge(noDnaFusions,
                     ensemblGeneData %>% mutate(GeneIdDown=GeneId,GeneStartDown=GeneStart,GeneEndDown=GeneEnd,GeneLengthDown=GeneEnd-GeneStart) %>% 
                       select(GeneIdDown,GeneStartDown,GeneEndDown,GeneLengthDown),by='GeneIdDown',all.x=T)

View(noDnaFusions %>% mutate(FusionLength=ifelse(SameChr,abs(RnaPosUp-RnaPosDown),0)) %>%
       select(SampleId,CancerType,RnaSvType,RnaPhaseMatched,FusionLength,JunctionFrags,DiscordantFrags,GeneNameUp,GeneNameDown,GeneLengthUp,GeneLengthDown,
              SVsMatched,DiffSVs,RnaJuncUp,RnaJuncDown,SvTypeUp,SvTypeDown,everything()))


# short INVs without DNA support
noDnaINVs = dnaRnaMatch %>% filter(RnaJuncUp=='KNOWN'&RnaJuncDown=='KNOWN'&DnaFusionMatchType=='NONE'&RnaSvType=='INV')
nrow(noDnaINVs)

noDnaINVs = merge(noDnaINVs,
                     ensemblGeneData %>% mutate(GeneIdUp=GeneId,GeneStartUp=GeneStart,GeneEndUp=GeneEnd,GeneLengthUp=GeneEnd-GeneStart) %>% 
                       select(GeneIdUp,GeneStartUp,GeneEndUp,GeneLengthUp),by='GeneIdUp',all.x=T)

noDnaINVs = merge(noDnaINVs,
                     ensemblGeneData %>% mutate(GeneIdDown=GeneId,GeneStartDown=GeneStart,GeneEndDown=GeneEnd,GeneLengthDown=GeneEnd-GeneStart) %>% 
                       select(GeneIdDown,GeneStartDown,GeneEndDown,GeneLengthDown),by='GeneIdDown',all.x=T)

View(noDnaINVs %>% mutate(FusionLength=ifelse(SameChr,abs(RnaPosUp-RnaPosDown),0)) %>%
       select(SampleId,CancerType,RnaSvType,RnaPhaseMatched,FusionLength,JunctionFrags,DiscordantFrags,GeneNameUp,GeneNameDown,GeneLengthUp,GeneLengthDown,
              SVsMatched,DiffSVs,RnaJuncUp,RnaJuncDown,SvTypeUp,SvTypeDown,everything()))


View(noDnaINVs %>% group_by(FusionLength=2**round(log(abs(RnaPosUp-RnaPosDown),2)),
                            JuncFrags=ifelse(JunctionFrags==0,1,2**round(log(JunctionFrags,2)))) %>% count %>% spread(JuncFrags,n,fill=0))

View(noDnaINVs %>% group_by(GeneNameUp,
                            JuncFrags=ifelse(JunctionFrags==0,1,2**round(log(JunctionFrags,2)))) %>% count %>% spread(JuncFrags,n,fill=0))

View(dnaRnaMatch %>% filter(KnownType=='KNOWN_PAIR'&(DnaFusionMatchType=='SVS'|DnaFusionMatchType=='GENES')) %>% group_by(JunctionFrags) %>% count)
View(dnaRnaMatch %>% group_by(Source) %>% count)

# chained fusions with RNA support
View(dnaRnaMatch %>% filter(DnaFusionMatchType=='SVS'|DnaFusionMatchType=='GENES') %>% 
  filter(DiffSVs==F|DiffChains==F) %>%
  group_by(KnownType,Chained=DiffSVs) %>% count)



# Summary of DNA support in RNA
View(dnaRnaMatch %>% 
       filter(DnaFusionMatched&Reportable=='true') %>% 
       group_by(Source,SVsMatched,DnaFusionMatchType,RnaPhaseMatched) %>% count())

View(dnaRnaMatch %>% 
       filter(DnaFusionMatched&Reportable=='true') %>% 
       group_by(Source,DnaFusionMatchType,KnownType) %>% count() %>% spread(KnownType,n,fill=0))

# group by gene pair
rnaGenePairSummary = dnaRnaMatch %>% filter(GeneIdUp!=''&GeneIdDown!='') %>%
  group_by(SampleId,GeneNameUp,GeneNameDown) %>% 
  summarise(RnaCount=n(),
            ExactMatches=sum(DnaFusionMatchType=='SVS'),
            GeneMatches=sum(DnaFusionMatchType=='GENES'),
            PassingRNA=sum(Source=='MATCH'|Source=='ISOFOX_ONLY'|Source=='ARRIBA_ONLY'),
            PassMatchSVs=sum((Source=='MATCH'|Source=='ISOFOX_ONLY'|Source=='ARRIBA_ONLY')&SVsMatched=='BOTH'))

nrow(rnaGenePairSummary)
View(rnaGenePairSummary)



## Load Linx's (DNA) fusions

linxFusions122 = read.csv("~/data/rna/fusions/LNX_FUSIONS.csv") 
linxFusions2131 = read.csv('~/data/rna/cohort/LNX_FUSIONS_2131.csv')
linxFusions = linxFusions2131

nrow(linxFusions)
View(linxFusions)
nrow(linxFusions %>% group_by(SampleId) %>% count)
View(linxFusions %>% group_by(SampleId) %>% count)
View(linxFusions %>% group_by(KnownType,Reportable,PhaseMatched) %>% count %>% spread(PhaseMatched,n,fill=0))

View(linxFusions %>%  filter(is.na(GeneNameUp)))

dnaGenePairSummary = linxFusions %>% 
  filter(GeneNameUp!=''&GeneNameDown!='') %>%
  group_by(SampleId,GeneNameUp=as.character(GeneNameUp),GeneNameDown=as.character(GeneNameDown),KnownType,Reportable) %>% 
  summarise(DnaCount=n())

# link to expression
dnaGenePairSummary = merge(dnaGenePairSummary,dnaFusionGeneExp %>% select(SampleId,GeneNameUp=GeneName,TpmUp=TPM),by=c('SampleId','GeneNameUp'),all.x=T)
dnaGenePairSummary = dnaGenePairSummary %>% mutate(TpmBucket=ifelse(is.na(TpmUp)|TpmUp<0.1,0,10**round(log(TpmUp,10))))
nrow(dnaGenePairSummary)
View(dnaGenePairSummary)
View(dnaGenePairSummary %>% group_by(KnownType) %>% count)
View(dnaGenePairSummary %>% group_by(SampleId) %>% count)

dnaRnaResults = merge(dnaGenePairSummary %>% select(-DnaCount),rnaGenePairSummary,by=c('SampleId','GeneNameUp','GeneNameDown'),all.x=T)
dnaRnaResults = merge(dnaRnaResults,dnaFusionGeneExp %>% select(SampleId,GeneNameUp=GeneName,TpmUp=TPM),by=c('SampleId','GeneNameUp'),all.x=T)
dnaRnaResults = dnaRnaResults %>% mutate(TpmBucket=ifelse(is.na(TpmUp),0,10**round(log(TpmUp,10)))) 
View(dnaRnaResults)
View(dnaRnaResults %>% group_by(SampleId) %>% count)
View(dnaRnaResults %>% filter(is.na(RnaCount)))

# write.csv(dnaRnaResults,'~/data/rna/fusions/dna_rna_fusion_match_results_122.csv',row.names = F,quote = F)
write.csv(dnaRnaResults,'~/data/rna/fusions/dna_rna_fusion_match_results_2131.csv',row.names = F,quote = F)

View(dnaRnaResults %>% filter(KnownType %in% c('KNOWN_PAIR','PROMISCUOUS_5','PROMISCUOUS_3','PROMISCUOUS_BOTH')) %>% 
       group_by(KnownType,Reportable) %>% 
       summarise(DnaFusions=n(),
                 RnaSupport=sum(!is.na(RnaCount)&RnaCount>0),RnaPassSupport=sum(!is.na(RnaCount)&PassingRNA>0),
                 RnaSvMatchSupport=sum(!is.na(RnaCount)&PassMatchSVs>0)) %>%
       mutate(HasRnaSupport=RnaSvMatchSupport>0))

View(dnaRnaResults %>% filter(KnownType %in% c('KNOWN_PAIR','PROMISCUOUS_5','PROMISCUOUS_3','PROMISCUOUS_BOTH')) %>% 
       group_by(KnownType,
                TpmBucket=ifelse(TpmBucket<=1,'TPM 0-1',ifelse(TpmBucket<=100,'TPM 10-100','TPM 100+')),
                HasRnaSupport=!is.na(PassMatchSVs)&PassMatchSVs>0) %>% count %>% spread(TpmBucket,n,fill=0))



View(dnaRnaResults %>% group_by(KnownType,Reportable,TpmBucket) %>% 
       summarise(DnaFusions=n(),
                 RnaSupport=sum(!is.na(RnaCount)&RnaCount>0),RnaPassSupport=sum(!is.na(RnaCount)&PassingRNA>0),
                 RnaSvMatchSupport=sum(!is.na(RnaCount)&PassMatchSVs>0)) %>%
       spread(TpmBucket,DnaFusions,fill=0))


View(dnaRnaResults %>% group_by(TpmBucket,
                                RnaPassSupport=(!is.na(RnaCount)&PassingRNA>0)) %>% count %>%
       spread(RnaPassSupport,n,fill=0))

View(dnaRnaResults %>% group_by(KnownType,TpmBucket,
                                RnaSupport=(!is.na(RnaCount)&RnaCount>0),RnaPassSupport=sum(!is.na(RnaCount)&PassingRNA>0)) %>% count %>%
       spread(RnaSupport,n,fill=0))

View(dnaRnaMatch %>% 
       filter(DnaFusionMatched&Reportable=='true') %>% 
       group_by(Source,DnaFusionMatchType,KnownType) %>% count() %>% spread(KnownType,n,fill=0))


View(linxFusions %>% filter(Reportable=='true') %>% group_by(KnownType) %>% count)






# Arriba fusions

dnaRnaArriba122 = read.csv('~/data/rna/fusions/LNX_RNA_FUSION_MATCH_ARRIBA.csv')
View(dnaRnaArriba122)
nrow(dnaRnaArriba122)

dnaRnaArriba122 = annotate_dna_rna_match(dnaRnaArriba122)

View(dnaRnaArriba122 %>% mutate(Type=ifelse(RnaSvType=='BND','BND',ifelse(FusionDistance<1e5,paste(RnaSvType,'SHORT',sep='_'),paste(RnaSvType,'LONG',sep='_')))) %>% 
       group_by(ValidSVs,GeneType,Type) %>% count() %>% spread(Type,n,fill=0))


## so few BNDs in Arriba vs Isofox
View(dnaRnaMatch122 %>% filter(RnaSvType=='BND'&KnownType=='Known'))
View(dnaRnaArriba122 %>% filter(RnaSvType=='BND'&KnownType=='Known'))
View(fusionCompareDna123)



##########
## Debug

View(ensemblTransData %>% filter(StableId=='ENST00000437593'))

# Single Sample


dnaRnaMatch = read.csv('~/logs/LNX_RNA_FUSION_MATCH_MATCHED.csv')
dnaRnaMatch = read.csv('~/logs/LNX_RNA_FUSION_MATCH_ISOFOX.csv')
#dnaRnaMatch = read.csv('~/data/rna/fusions/LNX_RNA_FUSION_MATCH.csv')

dnaRnaMatch = annotate_dna_rna_match(dnaRnaMatch)
View(dnaRnaMatch)

View(ensemblTransExonData %>% filter(ExonStart==65168338|ExonEnd==65168338))

# no support at either end
View(dnaRnaMatch %>% 
       # filter(TransValidLocUp=='false'|TransValidLocDown=='false') %>% 
       select(FusionId,SVsMatched,DnaFusionMatched,GeneType,RnaSvType,FusionDistance,ValidSVs,DiffSVs,DiffClusters,DiffChains,
              SvIdUp,ChrUp,RnaPosUp,SvPosUp,RnaOrientUp,RnaJuncUp,TransViableUp,TransValidLocUp,
              SvIdDown,ChrDown,RnaPosDown,SvPosDown,RnaOrientDown,RnaJuncDown,TransViableDown,TransValidLocDown,JunctionFrags,everything()))


# with support
View(dnaRnaMatch %>% filter(TransValidLocUp=='true'&TransValidLocDown=='true') %>% 
       select(FusionId,SvIdUp,ChrUp,RnaPosUp,SvPosUp,RnaOrientUp,RnaJuncUp,SvIdDown,ChrDown,RnaPosDown,SvPosDown,RnaOrientDown,RnaJuncDown,JunctionFrags,everything()))

sampleId = 'CPCT02020378T'
isofoxFusions = read.csv('~/data/rna/cohort/isofox_combined_fusions.csv')
View(isofoxFusions)

linxSvData = read.csv('~/logs/LNX_SVS.csv')
View(linxSvData %>% filter(SampleId==sampleId))

View(linxSvData %>% filter(SampleId==sampleId) %>% filter(ChrStart==5&ChrEnd==19&Type=='BND'))

## Working and debug


# Arriba fusions
arribaFusions = read.csv('~/data/rna/runs/CPCT02020378T.fusions.tsv',sep='\t')



View(ensemblTransExonData %>% filter(GeneId=='ENSG00000171889'))

