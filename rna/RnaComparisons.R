

cohortFusions = read.csv('~/data/rna/cohort/isofox_fusion_cohort.csv')
nrow(cohortFusions)
View(cohortFusions)
View(cohortFusions %>% group_by(Samples) %>% mutate(USC=n()) %>% ungroup() %>% filter(JuncTypeUp=='KNOWN',JuncTypeDown=='KNOWN') %>% filter(USC>1))



# passing Isofox fusions
passFusions = allRnaFusions %>% mutate(AF=(SplitFrags+RealignedFrags)/pmax(CoverageUp,CoverageDown)) %>% 
  filter(((MaxAnchorLengthUp>20&MaxAnchorLengthDown>20)|DiscordantFrags>0)&
      ((JuncTypeUp=='KNOWN'&JuncTypeDown=='KNOWN'&AF>=0.005&TotalFragments>=2)|  #splice site - splice site
         (((JuncTypeUp=='CANONICAL'&JuncTypeDown=='KNOWN')|(JuncTypeUp=='KNOWN'&JuncTypeDown=='CANONICAL'))&AF>=0.005&TotalFragments>=3)|  #canonical - splice
         (JuncTypeUp=='CANONICAL'&JuncTypeDown=='CANONICAL'&AF>=0.005&TotalFragments>=4)|  #canonical - canonical
         (AF>=0.05&TotalFragments>=10)))

nrow(passFusions)
View(passFusions)

View(passFusions %>% select(NonSupp,GeneNameUp,ChrUp,PosUp,OrientUp,GeneNameDown,ChrDown,PosDown,OrientDown,SVType,SplitFrags,RealignedFrags,DiscordantFrags,everything()))

# filter by cohort
knownPairs = read.csv('~/data/knownFusionPairs.csv')
knownPairs$KnownPairType = 'KNOWN'
passFusions = merge(passFusions,knownPairs %>% select(GeneNameUp=FiveGene,GeneNameDown=ThreeGene,KnownPairType),by=c('GeneNameUp','GeneNameDown'),all.x=T)
View(passFusions)
View(knownPairs %>% select(GeneNameUp=FiveGene,GeneNameDown=ThreeGene))

passFusions=merge(passFusions,cohortFusions %>% select(PosUp,OrientUp,PosDown,OrientDown,SampleCount),by=c('PosUp','OrientUp','PosDown','OrientDown'),all.x=T)

View(passFusions %>% filter(is.na(SampleCount)|!is.na(KnownPairType)))

isfPassFusions = read.csv('~/data/rna/runs/CPCT02020378T.isf.pass_fusions.csv')
View(isfPassFusions)


# old filter criteria
#!(grepl('HLA',GeneNameUp)),!(grepl('HLA',GeneNameDown)), #HLA issues (presumably sorted out by PON)
#!(grepl('IGH',GeneNameUp)),!(grepl('IGK',GeneNameUp)),!(grepl('IGH',GeneNameDown)),
#!(grepl('IGK',GeneNameDown)),!(grepl('IGL',GeneNameUp)),!(grepl('IGL',GeneNameDown)), #HLA issues (presumably sorted out by PON)
#(len>500000|(JuncTypeUp=='KNOWN'&JuncTypeDown=='KNOWN')|SVType=='INV'|SVType=='BND'),




View(passFusions)

passFusions=(merge(passFusions,cohortFusions %>% select(PosUp,OrientUp,PosDown,OrientDown,SampleCount),by=c('ChrUp','ChrDown','PosUp','OrientUp','PosDown','OrientDown'),all.x=T)) %>% 
  filter((GeneNameUp=='TMPRSS2'|is.na(SampleCount))) ### TMPRSS2 filter temporary


rnaAllCohortMerged=merge(allRnaFusions,cohortFusions,
                         by=c('ChrUp','ChrDown','PosUp','OrientUp','PosDown','OrientDown','GeneIdUp','GeneNameUp','GeneIdDown','GeneNameDown','SVType',
                              'JuncTypeUp','JuncTypeDown'),
                         all.x=T)

nrow(rnaAllCohortMerged %>% filter(!is.na(SampleCount)))
View(rnaAllCohortMerged %>% group_by(InCohort=!is.na(SampleCount)) %>% count)
colnames(cohortFusions)

# write for testing purposes
write.csv(rnaAllCohortMerged %>% filter(!is.na(SampleCount)) %>% 
            select(GeneIdUp,GeneNameUp,ChrUp,PosUp,OrientUp,JuncTypeUp,GeneIdDown,GeneNameDown,ChrDown,PosDown,OrientDown,JuncTypeDown,SVType,SampleCount,TotalFragments=TotalFragments.y,MaxFragments),
          '~/data/rna/cohort/isofox_cohort_fusions_subset.csv',row.names = F, quote = F)

isfPassFusions = read.csv('~/data/rna/cohort/CPCT02020378T.isf.pass_fusions.csv')
View(isfPassFusions)

isfPassFusions2 = read.csv('~/data/rna/fusions/sample_fusions/CPCT02020378T.isf.pass_fusions.csv')
View(isfPassFusions2)

# manual validation
View(isfPassFusions %>% mutate(AF=(SplitFrags+RealignedFrags)/pmax(CoverageUp,CoverageDown),
                               KnownPair=(KnownFusionTypr=='KNOWN_PAIR'),
                               TotalFragments=SplitFrags+RealignedFrags+DiscordantFrags) %>%
       filter(KnownPair|((DiscordantFrags>0|(pmin(MaxAnchorLengthUp,MaxAnchorLengthDown)>=20))&
                           ((JuncTypeUp=='KNOWN'&JuncTypeDown=='KNOWN'&AF>=0.005&TotalFragments>=2)|
                              (((JuncTypeUp=='KNOWN'&JuncTypeDown=='CANONICAL')|(JuncTypeUp=='CANONICAL'&JuncTypeDown=='KNOWN'))&AF>=0.005&TotalFragments>=3)|
                              (JuncTypeUp=='CANONICAL'&JuncTypeDown=='CANONICAL'&AF>=0.005&TotalFragments>=4)|
                              (AF>=0.05&TotalFragments>=10)))) %>% 
       select(GeneNameUp,GeneNameDown,SplitFrags,RealignedFrags,DiscordantFrags,MaxAnchorLengthUp,MaxAnchorLengthDown,AF,CoverageUp,CoverageDown,everything()))


isfAnnFusions = read.csv('~/data/rna/cohort/CPCT02020378T.isf.fusions.csv')
View(isfAnnFusions)
View(isfAnnFusions %>% group_by(InCohort=CohortCount>0,Filter) %>% count)

View(isfAnnFusions %>% filter(Filter=='SUPPORT'&JuncTypeUp=='KNOWN'&JuncTypeDown=='KNOWN') %>% 
       mutate(AF=(SplitFrags+RealignedFrags)/pmax(CoverageUp,CoverageDown)) %>%
       filter(AF>=0.005&(DiscordantFrags>0|(pmin(MaxAnchorLengthUp,MaxAnchorLengthDown)>=20))) %>% 
       select(GeneNameUp,GeneNameDown,SplitFrags,RealignedFrags,DiscordantFrags,MaxAnchorLengthUp,MaxAnchorLengthDown,AF,CoverageUp,CoverageDown,everything()))

# Arriba fusions
arribaFusions = read.csv('~/data/rna/runs/CPCT02020378T.fusions.tsv',sep='\t')

arribaFusions = arribaFusions %>% separate(breakpoint1,c('Chr1','Pos1'),sep=':') %>% separate(breakpoint2,c('Chr2','Pos2'),sep=':')
nrow(arribaFusions)
View(arribaFusions)
colnames(arribaFusions)


arribaFusions = arribaFusions %>% mutate(SampleId=sampleId)
write.csv(arribaFusions %>% mutate(JunctionFrags=split_reads1+split_reads2,
                                   OrientUp=ifelse(strand1.gene.fusion.=='+/+',1,ifelse(strand1.gene.fusion.=='-/-',-1,0)),
                                   OrientDown=ifelse(strand2.gene.fusion.=='-/-',1,ifelse(strand1.gene.fusion.=='+/+',-1,0))) %>%
            select(SampleId,GeneNameUp=X.gene1,GeneNameDown=gene2,ChrUp=Chr1,ChrDown=Chr2,PosUp=Pos1,PosDown=Pos2,OrientUp,OrientDown,JuncTypeUp=site1,JuncTypeDown=site2,
                                   JunctionFrags,DiscordantFrags=discordant_mates,CoverageUp=coverage1,CoverageDown=coverage2),
          '~/data/rna/cohort/arriba_fusions.csv', quote = F, row.names = F)

fusionCompareDna122 = read.csv('~/data/rna/cohort/isofox_ext_fusions_compare_122.csv')
View(fusionCompareDna122)
nrow(fusionCompareDna122)
View(fusionCompareDna122 %>% group_by(MatchType) %>% count)
View(fusionCompareDna122 %>% group_by(SampleId) %>% count)
View(fusionCompareDna122 %>% group_by(SampleId,MatchType) %>% count %>% spread(MatchType,n,fill=0))


fusionCompareDnaAll = read.csv('~/data/rna/cohort/isofox_ext_fusions_compare_2131.csv')
View(fusionCompareDnaAll)

write.csv(fusionCompareDnaAll %>% filter(SampleId %in% dnaSamples$SampleId),'~/data/rna/cohort/isofox_ext_fusions_compare_122b.csv',row.names = F,quote = F)

passFusionsAll = read.csv('~/data/rna/cohort/isofox_combined_fusions.csv')
nrow(passFusionsAll)
View(passFusionsAll)

fusionCompareDna122Comb = merge(fusionCompareDna122 %>% mutate(FusionId=ifelse(FusionId!='-1',paste('Id',FusionId,sep='_'),'-1')),
                                passFusionsAll %>% select(SampleId,FusionId,JuncTypeUp,JuncTypeDown,TotalFragments,SplitFrags,RealignedFrags,
                                                          CoverageUp,CoverageDown,MaxAnchorLengthUp,MaxAnchorLengthDown,CohortCount,KnownFusionType),
                                by=c('SampleId','FusionId'),all.x=T)

View(fusionCompareDna122Comb)
write.csv(fusionCompareDna122Comb,'~/data/rna/cohort/isofox_ext_fusions_compare_dna_122_ann.csv',row.names = F,quote = F)
colnames(passFusionsAll)
colnames(fusionCompareDna122)


# Arriba fusions missing from Isofox
arribaOnly = fusionCompare %>% filter(MatchType=='EXT_ONLY')
nrow(arribaOnly)

arribaIsfAll = merge(arribaFusions %>% select(ChrUp=Chr1,ChrDown=Chr2,PosUp=Pos1,PosDown=Pos2),allRnaFusions,
                     by=c('ChrUp','ChrDown','PosUp','PosDown'),all.x=T)

View(arribaIsfAll)

####################
# Isofox cohort file

dnaSamples = read.csv('~/data/rna/samples/dl_dna_fusion_samples.csv')
nrow(dnaSamples) # 122
sourceDir = '~/data/rna/fusions/sample_fusions'

isofoxFusions122 = data.frame()

for(sampleId in dnaSamples$SampleId)
{
  filename = sprintf('%s/%s.isf.pass_fusions.csv',sourceDir,sampleId)
  print(paste('loading fusions for :',filename,sep = ''))
  fusions = read.csv(filename)
  fusions = fusions %>% mutate(ChrUp=as.character(ChrUp),ChrDown=as.character(ChrDown))
  fusions$SampleId = sampleId
  isofoxFusions122 = rbind(isofoxFusions122,fusions)
}

View(isofoxFusions122)
nrow(isofoxFusions122)
nrow(isofoxFusions122 %>% group_by(SampleId) %>% count) # 122 samples
View(isofoxFusions122 %>% group_by(SampleId) %>% count)

write.csv(isofoxFusions122,'~/data/rna/fusions/isofox_pass_fusions_dna122.csv',row.names = F, quote = F)

#####################
## Arriba cohort file

arribaFusions122 = data.frame()

for(sampleId in dnaSamples$SampleId)
{
  filename = sprintf('%s/%s.fusions.tsv',sourceDir,sampleId)
  print(paste('loading fusions for :',filename,sep=''))
  fusions = read.csv(filename,sep='\t')
  fusions$SampleId = sampleId

  fusions = fusions %>% separate(breakpoint1,c('Chr1','Pos1'),sep=':') %>% separate(breakpoint2,c('Chr2','Pos2'),sep=':')
  
  fusions = fusions %>% 
    mutate(GeneNameUp=stri_replace_all_fixed(X.gene1,',',';'),
           GeneNameDown=stri_replace_all_fixed(gene2,',',';'),
           JunctionFrags=split_reads1+split_reads2,
           OrientUp=ifelse(strand1.gene.fusion.=='+/+',1,ifelse(strand1.gene.fusion.=='-/-',-1,0)),
           OrientDown=ifelse(strand2.gene.fusion.=='-/-',1,ifelse(strand1.gene.fusion.=='+/+',-1,0)),
           Filters=stri_replace_all_fixed(filters,',',';')) %>%
    select(SampleId,GeneNameUp,GeneNameDown,ChrUp=Chr1,ChrDown=Chr2,PosUp=Pos1,PosDown=Pos2,OrientUp,OrientDown,
           Type=type,Direction1=direction1,Direction2=direction2,Strand1=strand1.gene.fusion.,Strand2=strand2.gene.fusion.,
           JuncTypeUp=site1,JuncTypeDown=site2,JunctionFrags,DiscordantFrags=discordant_mates,CoverageUp=coverage1,CoverageDown=coverage2,
           Confidence=confidence,Filters)
            
  arribaFusions122 = rbind(arribaFusions122,fusions)
}

View(arribaFusions122)
nrow(arribaFusions122) # 10302
nrow(arribaFusions122 %>% group_by(SampleId) %>% count) # 122 samples
View(arribaFusions122 %>% group_by(SampleId) %>% count)
View(arribaFusions122 %>% filter(is.na(ChrUp)|is.na(ChrDown)))

arribaTmp = arribaFusions122 %>% filter(!grepl(';',GeneNameUp)&!grepl(';',GeneNameDown))
View(arribaTmp)
arribaTmp = merge(arribaTmp,ensemblGeneData %>% select(GeneNameUp=GeneName,ChrUp=Chromosome,StrandUp=Strand),by=c('GeneNameUp','ChrUp'),all.x=T)
arribaTmp = merge(arribaTmp,ensemblGeneData %>% select(GeneNameDown=GeneName,ChrDown=Chromosome,StrandDown=Strand),by=c('GeneNameDown','ChrDown'),all.x=T)
View(arribaTmp %>% select(GeneNameUp,GeneNameDown,Direction1,Direction2,Strand1,Strand2,StrandUp,StrandDown,Type,everything()))


View(arribaTmp %>% group_by(Strand1,Strand2,StrandUp,StrandDown,Type,PosUpHigher=(PosUp>PosDown)) %>% count)

# multiple genes listsed
View(arribaFusions122 %>% filter(grepl(';',GeneNameUp)|grepl(';',GeneNameDown)))


write.csv(arribaFusions122,'~/data/rna/fusions/arriba_fusions_dna122.csv',row.names = F, quote = F)


arribaFusions = read.csv('~/data/rna/runs/CPCT02020378T.fusions.tsv',sep='\t')
View(arribaFusions)
arribaFusions2 = read.csv('~/data/rna/fusions/sample_fusions/CPCT02010386T.fusions.tsv',sep='\t')
View(arribaFusions2)

View(ensemblTransExonData %>% filter(GeneId=='ENSG00000206892'))
View(ensemblGeneData %>% filter(GeneName=='MALAT1'|GeneName=='SCYL1'))
View(ensemblGeneData %>% filter(GeneName=='HRH1'|GeneName=='ATG7'))
View(ensemblGeneData %>% filter(GeneName=='MDM2'|GeneName=='RPSAP52'))

# genePair rna(RPSAP52-RP11-611O2.5) differs from dna(MDM2-HMGA2) for same SVs(1087 & 1087)
# 16:19:17.454 [main] [INFO ] genePair rna(VGLL4-HRH1) differs from dna(VGLL4-ATG7) for same SVs(41 & 41)
# genePair rna(CDC42EP2-MALAT1) differs from dna(CDC42EP2-SCYL1) for same SVs(299 & 299)

nrow(arribaFusions)
View(arribaFusions)
colnames(arribaFusions)

View(arribaFusions %>% group_by(type) %>% count)


arribaFusions = arribaFusions %>% mutate(SampleId=sampleId)
write.csv(arribaFusions %>% mutate(JunctionFrags=split_reads1+split_reads2,
                                   OrientUp=ifelse(strand1.gene.fusion.=='+/+',1,ifelse(strand1.gene.fusion.=='-/-',-1,0)),
                                   OrientDown=ifelse(strand2.gene.fusion.=='-/-',1,ifelse(strand1.gene.fusion.=='+/+',-1,0))) %>%
            select(SampleId,GeneNameUp=X.gene1,GeneNameDown=gene2,ChrUp=Chr1,ChrDown=Chr2,PosUp=Pos1,PosDown=Pos2,OrientUp,OrientDown,JuncTypeUp=site1,JuncTypeDown=site2,
                   JunctionFrags,DiscordantFrags=discordant_mates,CoverageUp=coverage1,CoverageDown=coverage2),
          '~/data/rna/cohort/arriba_fusions.csv', quote = F, row.names = F)




#########################
## StarFusion cohort file

sfFusions = read.csv('~/data/rna/fusions/star_fusion_predictions_630.csv')
nrow(sfFusions)
View(sfFusions)
nrow(sfFusions %>% group_by(SampleId) %>% count) # 581 samples

View(dnaSamples)

View(sfFusions %>% filter(SampleId %in% dnaSamples$SampleId) %>% group_by(SampleId) %>% count)


#########
## Debug

fusionCompare = read.csv('~/logs/isofox_ext_fusions_compare.csv')
View(fusionCompare)
View(fusionCompare %>% group_by(MatchType) %>% count)




tmpCohortFusions = read.csv('~/logs/isofox_fusion_cohort.csv')
nrow(tmpCohortFusions)
View(tmpCohortFusions)

sampleId = 'CPCT02020378T'
allRnaFusions = read.csv(formFilename('~/data/rna/runs/',sampleId,'all.fusions.csv'))
nrow(allRnaFusions)
View(allRnaFusions)









