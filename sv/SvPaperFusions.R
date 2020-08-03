library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringi)

## DEPRECATED
annotate_fusions<-function(fusionData)
{
  fusionData = fusionData %>% mutate(SameSV=(SvIdUp==SvIdDown),
                                     KnownType=ifelse(!is.na(KnownType)&KnownType!='',as.character(KnownType),'Unknown'))
  
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


########################
## Fusion Reference Data

newKnownPairs = read.csv('~/data/sv/known_pairs_new.csv',sep='\t')
View(newKnownPairs)
nrow(newKnownPairs)

newProms = read.csv('~/data/sv/promiscuous_genes.csv',sep='\t')
View(newProms)
igGenes = read.csv('~/data/sv/ig_fusion_genes.csv')
View(igGenes)
exonDelDups = read.csv('~/data/sv/exon_del_dup_genes.csv',sep='\t')
View(exonDelDups)
colnames(exonDelDups)

write.csv(rbind(newKnownPairs %>% select(GeneName=FiveGene),
                newKnownPairs %>% select(GeneName=ThreeGene)),'~/logs/known_pair_genes.csv',row.names = F,quote = F)

combinedKnownFusions = newKnownPairs %>% mutate(Type='KNOWN_PAIR',PubMedId=as.character(PubMedId),OtherData='') %>% 
  select(Type,FiveGene,ThreeGene,CancerTypes,PubMedId,OtherData)

View(combinedKnownFusions)

combinedKnownFusions = rbind(combinedKnownFusions,
                             newProms %>% filter(Stream=='Up') %>% mutate(Type='PROMISCUOUS_5',FiveGene=GeneName,ThreeGene='',PubMedId='',OtherData='') %>%
                               select(Type,FiveGene,ThreeGene,CancerTypes,PubMedId,OtherData))

combinedKnownFusions = rbind(combinedKnownFusions,
                             newProms %>% filter(Stream=='Down') %>% mutate(Type='PROMISCUOUS_3',FiveGene='',ThreeGene=GeneName,PubMedId='',OtherData='') %>%
                               select(Type,FiveGene,ThreeGene,CancerTypes,PubMedId,OtherData))

combinedKnownFusions = rbind(combinedKnownFusions,igGenes)

combinedKnownFusions = rbind(combinedKnownFusions,
                             exonDelDups %>% mutate(Type='EXON_DEL_DUP',PubMedId='',
                                                    FiveGene=Gene,ThreeGene=Gene,
                                                    OtherData=sprintf('%s;%d;%d;%d;%d',Transcript,MinFusedExonStart,MaxFusedExonStart,MinFusedExonEnd,MaxFusedExonEnd)) %>%
                               select(Type,FiveGene,ThreeGene,CancerTypes,PubMedId,OtherData))

View(combinedKnownFusions)
write.csv(combinedKnownFusions,'~/data/sv/known_fusion_data.csv',row.names = F, quote=F)

knownFusionData = read.csv('~/data/sv/known_fusion_data.csv')
View(knownFusionData)

# write in old-style format
write.csv(newKnownPairs %>% select(FiveGene,ThreeGene),'~/data/sv/knownFusionPairs_20200603.csv',row.names = F, quote=F)

# known and promiscuous genes for impact of expression on fusions
expGenes = combinedKnownFusions %>% select(GeneName=FiveGene)
expGenes = rbind(expGenes,combinedKnownFusions %>% select(GeneName=ThreeGene))
expGenes = merge(ensemblGeneData %>% select(GeneId,GeneName),expGenes,by='GeneName',all.y=T)

write.csv(expGenes %>% group_by(GeneId,GeneName) %>% count %>% select(GeneId,GeneName),
          '~/data/sv/fusions/known_gene_ids.csv',row.names = F, quote=F)


## check gene names for HG38
ensemblGeneData = read.csv('~/data/sv/ensembl_gene_data.csv')
View(ensemblGeneData)

ensemblGeneData38 = read.csv('~/data/ensembl_hg38/ensembl_gene_data.csv')
View(ensemblGeneData38)

View(knownFusionData %>% group_by(Type) %>% count)
knownFusionData38 = knownFusionData 
knownFusionData38 = merge(knownFusionData38,ensemblGeneData %>% select(FiveGene=GeneName,FiveGeneId=GeneId),by='FiveGene',all.x=T)
knownFusionData38 = merge(knownFusionData38,ensemblGeneData %>% select(ThreeGene=GeneName,ThreeGeneId=GeneId),by='ThreeGene',all.x=T)

knownFusionData38 = merge(knownFusionData38,ensemblGeneData38 %>% select(FiveGene=GeneName,FiveGeneId38=GeneId),by='FiveGene',all.x=T)
knownFusionData38 = merge(knownFusionData38,ensemblGeneData38 %>% select(ThreeGene=GeneName,ThreeGeneId38=GeneId),by='ThreeGene',all.x=T)

knownFusionData38 = merge(knownFusionData38,ensemblGeneData38 %>% select(FiveGene38=GeneName,FiveGeneId=GeneId),by='FiveGeneId',all.x=T)
knownFusionData38 = merge(knownFusionData38,ensemblGeneData38 %>% select(ThreeGene38=GeneName,ThreeGeneId=GeneId),by='ThreeGeneId',all.x=T)

View(knownFusionData38)
View(knownFusionData38 %>% filter(Type=='KNOWN_PAIR') %>% filter(is.na(FiveGeneId38)|is.na(ThreeGeneId38)) %>% 
       select(Type,FiveGene,FiveGene38,FiveGeneId,FiveGeneId38,ThreeGene,ThreeGene38,ThreeGeneId,ThreeGeneId38))


ampGenes = read.csv('~/hmf/repos/hmftools/hmf-common/src/main/resources/cna/AmplificationTargets.tsv',sep='\t',header = F)
colnames(ampGenes) = c('gene')

delGenes = read.csv('~/hmf/repos/hmftools/hmf-common/src/main/resources/cna/DeletionTargets.tsv',sep='\t',header = F)
colnames(delGenes) = c('gene','loca')

dndsGenesOnco = read.csv('~/hmf/repos/hmftools/hmf-common/src/main/resources/dnds/DndsDriverLikelihoodOnco.tsv',sep='\t')
dndsGenesTsg = read.csv('~/hmf/repos/hmftools/hmf-common/src/main/resources/dnds/DndsDriverLikelihoodTsg.tsv',sep='\t')

View(dndsGenesOnco)
View(dndsGenesTsg)
View(ampGenes)
View(delGenes)

allGenes = rbind(ampGenes,
                 delGenes %>% select(gene),
                 dndsGenesOnco %>% select(gene),
                 dndsGenesTsg %>% select(gene))

allGenes2 = allGenes %>% group_by(gene) %>% count %>% ungroup()
allGenes2 = allGenes2 %>% select(GeneName=gene)
View(allGenes2)

allGenes2 = merge(allGenes2,ensemblGeneData %>% select(GeneName,GeneId),by='GeneName',all.x=T)

allGenes2 = merge(allGenes2,ensemblGeneData38 %>% select(GeneName,GeneId38=GeneId),by='GeneName',all.x=T)

allGenes2 = merge(allGenes2,ensemblGeneData38 %>% select(GeneName38=GeneName,GeneId),by='GeneId',all.x=T)
View(allGenes2)
View(allGenes2 %>% filter(is.na(GeneId)))
View(allGenes2 %>% filter(is.na(GeneName38)|is.na(GeneId38)))
allGeneChanges = allGenes2 %>% filter(is.na(GeneName38)|is.na(GeneId38))

allGeneChanges = merge(allGeneChanges,ensemblGeneData %>% select(GeneId,Chromosome,GeneStart,GeneEnd,EntrezIds,Synonyms),by='GeneId',all.x=T)

allGeneChangesNameChg = merge(allGeneChanges %>% filter(!is.na(GeneName38)) %>% select(-GeneId38),
                              ensemblGeneData38 %>% select(GeneName38=GeneName,GeneId38=GeneId,Chromosome,GeneStart,GeneEnd,EntrezIds,Synonyms),by='GeneName38',all.x=T)

View(allGeneChangesNameChg %>% select(GeneName,GeneName38,GeneId,GeneId38,Chromosome.x,Chromosome.y,GeneStart.x,GeneStart.y,GeneEnd.x,GeneEnd.y,everything()))

allGeneChangesIdChg = merge(allGeneChanges %>% filter(!is.na(GeneId38)) %>% select(-GeneName38),
                              ensemblGeneData38 %>% select(GeneId38=GeneId,GeneName38=GeneName,Chromosome,GeneStart,GeneEnd,EntrezIds,Synonyms),by='GeneId38',all.x=T)

View(allGeneChangesIdChg %>% select(GeneName,GeneName38,GeneId,GeneId38,Chromosome.x,Chromosome.y,GeneStart.x,GeneStart.y,GeneEnd.x,GeneEnd.y,everything()))


######
## Fusion Comparisons

sampleIds = read.csv('~/data/sv/fusions/sample_ids_3909.csv')
nrow(sampleIds)
oldFusions = read.csv('~/data/sv/fusions/LNX_FUSIONS_old_5200.csv')
#oldFusions = oldFusions %>% filter(SampleId %in% sampleIds$SampleId)
nrow(oldFusions)
nrow(oldFusions %>% group_by(SampleId) %>% count)

# make known types comparable
View(oldFusions %>% group_by(KnownType) %>% count)
oldFusions = oldFusions %>% mutate(KnownType=ifelse(KnownType=='3P-Prom','PROMISCUOUS_3',
                                                    ifelse(KnownType=='5P-Prom','PROMISCUOUS_5',
                                                           ifelse(KnownType=='Both-Prom','PROMISCUOUS_BOTH',
                                                                  ifelse(KnownType=='Known','KNOWN_PAIR','NONE')))))
View(oldFusions %>% group_by(KnownType) %>% count)


newFusions = read.csv('~/data/sv/fusions/LNX_FUSIONS_new_5200.csv')
# newFusions = newFusionsAll %>% filter(SampleId %in% sampleIds$SampleId)
nrow(newFusions %>% group_by(SampleId) %>% count)
View(newFusions %>% group_by(KnownType) %>% count)
nrow(newFusions)
View(newFusions)
View(newFusions %>% filter(GeneNameUp=='PAX8'))

View(newFusions %>% filter(SampleId %in% c('CPCT02020618T','CPCT02010478T','CPCT02060225T','CPCT02120095T','CPCT02150021T')))


newFusions = merge(newFusions,samplesDD %>% select(SampleId,CancerType),by='SampleId',all.x=T)

igProm = newFusions %>% filter(KnownType=='IG_PROMISCUOUS')
igProm = merge(igProm,samplesDD %>% select(SampleId,CancerType),by='SampleId',all.x=T)
View(igProm %>% select(CancerType,GeneNameDown, everything()))

View(newFusions %>% filter(KnownType=='KNOWN_PAIR') %>% 
       filter(grepl('Raf-like Ras-binding',ProteinsKept)
              |grepl('Ets domain',ProteinsLost)
              |grepl('Protein kinase domain',ProteinsLost)
              |grepl('Epidermal growth factor-like domain',ProteinsLost)
              |grepl('Ankyrin repeat-containing domain',ProteinsLost)
              |grepl('Basic-leucine zipper domain',ProteinsLost)
              |grepl('High mobility group box domain',ProteinsLost)
              |grepl('Bromodomain',ProteinsLost)) %>% 
       select(SampleId,CancerType,KnownType,GeneNameUp,GeneNameDown,ProteinsKept,ProteinsLost,everything()))


# clean-up
rm(oldFusions)
rm(newFusions)
rm(newFusionsAll)
rm(cmpFusionsNew)
rm(cmpFusionsOld)

dnaFusionGeneExp = read.csv('~/data/rna/fusions/isofox_linx_up_genes.sample_gene_perc_data.csv')
View(dnaFusionGeneExp)

# difference in reportable fusions from POV of new fusions
cmpFusionsNew = merge(newFusions,oldFusions %>% select(SampleId,ReportableOld=Reportable,KnownTypeOld=KnownType,GeneIdUp,GeneIdDown),
                      by=c('SampleId','GeneIdUp','GeneIdDown'),all.x=T)

cmpFusionsNew = merge(cmpFusionsNew,dnaFusionGeneExp %>% select(SampleId,GeneIdUp=GeneId,TPM),by=c('SampleId','GeneIdUp'),all.x=T)
cmpFusionsNew = cmpFusionsNew %>% mutate(TPM=ifelse(is.na(TPM),0,round(TPM,1)))
cmpFusionsNew = cmpFusionsNew %>% mutate(TPM=round(TPM,1))

cmpFusionsNewRep = cmpFusionsNew %>% filter(Reportable=='true') %>% 
  filter(!(!is.na(ReportableOld)&ReportableOld==Reportable&!is.na(KnownTypeOld)&as.character(KnownType)==as.character(KnownTypeOld))) %>%
  select(SampleId,GeneNameUp,GeneNameDown,Reportable,ReportableOld,KnownType,KnownTypeOld,TPM,
         SvIdUp,SvIdDown,ClusterId,ResolvedType,PhaseMatched,everything())

nrow(cmpFusionsNewRep) # 137
View(cmpFusionsNewRep)
View(cmpFusionsNewRep %>% group_by(ReportableOld,KnownType,KnownTypeOld) %>% count)
View(cmpFusionsNewRep %>% group_by(GeneNameUp,GeneNameDown) %>% 
       summarise(Fusions=n(),ReportableOld=first(ReportableOld),KnownType=first(KnownType),KnownTypeOld=first(KnownTypeOld)))

View(cmpFusionsNewRep %>% filter(KnownType!='IG_KNOWN_PAIR'&KnownType!='IG_PROMISCUOUS') %>% group_by(TpmBucket=ifelse(TPM==0,0,10**round(log(TPM,10)))) %>% count)


write.csv(cmpFusionsNewRep,'~/data/sv/fusions/linx_reportable_changes_new.csv',row.names = F,quote = F)

# difference in reportable fusions from POV of old fusions
cmpFusionsOld = merge(oldFusions,newFusions %>% select(SampleId,ReportableNew=Reportable,KnownTypeNew=KnownType,GeneIdUp,GeneIdDown),
                      by=c('SampleId','GeneIdUp','GeneIdDown'),all.x=T)

cmpFusionsOld = merge(cmpFusionsOld,dnaFusionGeneExp %>% select(SampleId,GeneIdUp=GeneId,TPM),by=c('SampleId','GeneIdUp'),all.x=T)

cmpFusionsOld = cmpFusionsOld %>% mutate(TPM=ifelse(is.na(TPM),0,round(TPM,1)))

cmpFusionsOldRep = cmpFusionsOld %>% filter(Reportable=='true') %>% 
  filter(!(!is.na(ReportableNew)&ReportableNew==Reportable)) %>%
  select(SampleId,GeneNameUp,GeneNameDown,Reportable,ReportableNew,KnownType,KnownTypeNew,TPM,
         SvIdUp,SvIdDown,ClusterId,ResolvedType,PhaseMatched,everything())

View(cmpFusionsOldRep) # 126
View(cmpFusionsOldRep %>% group_by(ReportableNew,KnownType,KnownTypeNew) %>% count)
View(cmpFusionsOldRep %>% group_by(GeneNameUp,GeneNameDown) %>% 
       summarise(Fusions=n(),ReportableNew=first(ReportableNew),KnownType=first(KnownType),KnownTypeNew=first(KnownTypeNew)))

View(cmpFusionsOldRep %>% group_by(TpmBucket=ifelse(TPM==0,0,10**round(log(TPM,10)))) %>% count)

write.csv(cmpFusionsOldRep,'~/data/sv/fusions/linx_reportable_changes_old.csv',row.names = F,quote = F)


knownFusionData = read.csv('~/data/sv/known_fusion_data.csv')
View(knownFusionData)

oldKnownPairs = read.csv('~/data/knownFusionPairs.csv')
View(oldKnownPairs)
oldProm5 = read.csv('~/data/knownPromiscuousFive.csv')
oldProm3 = read.csv('~/data/knownPromiscuousThree.csv')
oldKnownFusionData = oldKnownPairs %>% mutate(Type='KNOWN_PAIR') %>% select(Type,FiveGene=fiveGene,ThreeGene=threeGene)
oldKnownFusionData = rbind(oldKnownFusionData,oldProm5 %>% mutate(Type='PROMISCUOUS_5',ThreeGene='') %>% select(Type,FiveGene=GeneName,ThreeGene))
oldKnownFusionData = rbind(oldKnownFusionData,oldProm3 %>% mutate(Type='PROMISCUOUS_3',FiveGene='') %>% select(Type,FiveGene,ThreeGene=GeneName))
View(oldKnownFusionData)

knownFusionChanges = merge(knownFusionData,oldKnownFusionData %>% mutate(IsOld=TRUE),by=c('Type','FiveGene','ThreeGene'),all=T)
knownFusionChanges = knownFusionChanges %>% mutate(RefDataChanges=ifelse(is.na(IsOld),'NEW_REF',ifelse(is.na(PubMedId),'OLD_REF','MATCH')))
View(knownFusionChanges)

# de-deup
knownFusionChanges = knownFusionChanges %>% group_by(Type,FiveGene,ThreeGene,RefDataChanges) %>% count
View(knownFusionChanges)



# link to old and new reference data
#cmpFusionsNew = merge(cmpFusionsNew,knownFusionChanges %>% select(KnownType=Type,GeneNameUp=FiveGene,GeneNameDown=ThreeGene,KnownPairRef=RefDataChanges),
#                      by=c('KnownType','GeneNameUp','GeneNameDown'),all.x=T)

View(cmpFusionsNew %>% 
       filter(!(!is.na(KnownTypeOld)&KnownType==KnownTypeOld)) %>%
       filter(!(is.na(KnownTypeOld)&KnownType=='NONE')) %>%
       group_by(GeneNameUp,GeneNameDown,Reportable,ReportableOld,KnownType,KnownTypeOld) %>% count)


View(cmpFusionsNew %>% select(SampleId,GeneNameUp,GeneNameDown,Reportable,ReportableOld,KnownPairRef,KnownType,KnownTypeOld))
nrow(cmpFusionsNew %>% filter(!is.na(ReportableOld)))
nrow(cmpFusionsNew %>% filter(!is.na(ReportableOld)&ReportableOld!=Reportable))
nrow(cmpFusionsNew %>% filter(!is.na(KnownTypeOld)&as.character(KnownTypeOld)!=as.character(KnownType)))

View(cmpFusionsNew %>% filter(Reportable=='true') %>% group_by(GeneNameUp,GeneNameDown,KnownType,TypeUp,TypeDown) %>% count)

# changes to reportable fusions
View(cmpFusionsNew %>% filter((!is.na(ReportableOld)&ReportableOld!=Reportable)|(is.na(ReportableOld)&Reportable=='true')) %>% 
       select(SampleId,GeneNameUp,GeneNameDown,Reportable,ReportableOld,RefDataChanges,KnownType,KnownTypeOld,SvIdUp,SvIdDown,ClusterId,ResolvedType,PhaseMatched,everything()))

View(cmpFusionsNew %>% filter((!is.na(ReportableOld)&ReportableOld!=Reportable)|Reportable=='true') %>% 
       select(SampleId,GeneNameUp,GeneNameDown,Reportable,ReportableOld,KnownType,KnownTypeOld,SvIdUp,SvIdDown,ClusterId,ResolvedType,PhaseMatched,everything()))

View(cmpFusionsNew %>% filter(!is.na(KnownTypeOld)&KnownTypeOld!=KnownType) %>% 
       select(SampleId,GeneNameUp,GeneNameDown,Reportable,ReportableOld,KnownType,KnownTypeOld,SvIdUp,SvIdDown,ClusterId,ResolvedType,PhaseMatched,everything()))


View(cmpFusionsNew %>% filter(!is.na(KnownTypeOld)&KnownTypeOld!=KnownType) %>% 
       select(SampleId,GeneNameUp,GeneNameDown,Reportable,ReportableOld,KnownType,KnownTypeOld,SvIdUp,SvIdDown,ClusterId,ResolvedType,PhaseMatched,everything()))



View(newFusions %>% filter(GeneNameUp=='A4GNT'&GeneNameDown=='ARMC8'))
View(oldFusions %>% filter(GeneNameUp=='A4GNT'&GeneNameDown=='ARMC8'))
View(oldFusions %>% filter(SampleId=='CPCT02011051T'&SvIdUp==111))



####
## SVA Fusions from single SVs and clustered & chained SVs
svaFusions = read.csv('~/data/sv/fusions/LNX_FUSIONS.csv')
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

newFusions = read.csv('~/data/sv/fusions/LNX_FUSIONS_20191206_new.csv') %>% filter(KnownType!='')
newFusions = newFusions %>% mutate(SampleGenePair=paste(SampleId,GeneNameUp,GeneNameDown,sep='_'))
newFusions = annotate_fusions(newFusions)
nrow(newFusions)

preChgFusions = read.csv('~/data/sv/fusions/LNX_FUSIONS_20191206_pre_change.csv') %>% filter(KnownType!='')
preChgFusions = preChgFusions %>% mutate(SampleGenePair=paste(SampleId,GeneNameUp,GeneNameDown,sep='_'))
preChgFusions = annotate_fusions(preChgFusions)
nrow(preChgFusions)

write.csv(preChgFusions,'~/data/sv/fusions/pre_chg_fusions.csv',quote = F,row.names = F)
write.csv(newFusions,'~/data/sv/fusions/new_fusions.csv',quote = F,row.names = F)

# comparisons
View(newFusions %>% filter(!(SampleGenePair %in% preChgFusions$SampleGenePair)) %>%
       select(SampleId,KnownType,GeneNameUp,GeneNameDown,Reportable,PhaseMatched,RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown,
              TypeUp,SvIdUp,TypeDown,SvIdDown,ProteinsKept,ProteinsLost,BiotypeUp,BiotypeDown,
              ChainLinks,ChainLength,TerminatedUp,TerminatedDown,ExonsSkippedUp,ExonsSkippedDown,everything()))

requiredBiotypes = c('protein_coding','retained_intron','processed_transcript','nonsense_mediated_decay','lincRNA')

View(preChgFusions %>% filter(!(SampleGenePair %in% newFusions$SampleGenePair)) %>%
       filter(BiotypeUp %in% requiredBiotypes) %>%
       filter(Reportable=='true') %>%
       select(SampleId,KnownType,GeneNameUp,GeneNameDown,Reportable,PhaseMatched,RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown,
              TypeUp,SvIdUp,TypeDown,SvIdDown,ProteinsKept,ProteinsLost,BiotypeUp,BiotypeDown,
              ChainLinks,ChainLength,TerminatedUp,TerminatedDown,ExonsSkippedUp,ExonsSkippedDown,everything()))




## Fusion Comparison
prevFusionsAll = read.csv('~/data/sv/fusions/LNX_FUSIONS_20191206_pre_change.csv')

prevFusionsAll = prevFusionsAll %>% filter(SampleId %in% changedSamples$SampleId)
prevFusionsAll = annotate_fusions(prevFusionsAll)

prevFusionsAll = prevFusionsAll %>% mutate(SampleGenePair=paste(SampleId,GeneNameUp,GeneNameDown,sep='_'),
                                           SampleIdPair=paste(SvIdUp,SvIdDown,sep='_'))
# prevFusions = prevFusions %>% filter(Reportable=='true')
prevFusions = prevFusionsAll %>% filter(KnownType!='')
nrow(prevFusions)

newFusionsAll = read.csv('~/data/sv/fusions/LNX_FUSIONS.csv')
newFusionsAll = read.csv('~/data/sv/fusions/LNX_FUSIONS_20191206_new.csv')
newFusionsAll = read.csv('~/data/sv/fusions/LNX_FUSIONS_20191206_biotypes.csv')
newFusionsAll = annotate_fusions(newFusionsAll)
newFusionsAll = newFusionsAll %>% mutate(SampleGenePair=paste(SampleId,GeneNameUp,GeneNameDown,sep='_'),
                                         SampleIdPair=paste(SvIdUp,SvIdDown,sep='_'))
# newFusions = newFusions %>% filter(Reportable=='true')
newFusions = newFusionsAll %>% filter(KnownType!='')
# newFusions = annotate_fusions(newFusions)
nrow(newFusions)

nrow(newFusionsAll %>% group_by(SampleId) %>% count)
nrow(prevFusionsAll %>% group_by(SampleId) %>% count)

View(newFusionsAll %>% group_by(SampleId,SvIdUp,SvIdDown) %>% count %>% filter(n>1))
View(prevFusionsAll %>% group_by(SampleId,SvIdUp,SvIdDown) %>% count %>% filter(n>1))

# all sample differences
samplesWithNew = newFusionsAll %>% filter(!(SampleGenePair %in% prevFusionsAll$SampleGenePair)) %>% select(SampleId,GeneNameUp,GeneNameDown)
samplesWithPrev = prevFusionsAll %>% filter(!(SampleGenePair %in% newFusionsAll$SampleGenePair)) %>% select(SampleId,GeneNameUp,GeneNameDown)

View(samplesWithNew %>% group_by(SampleId) %>% count)
View(samplesWithPrev %>% group_by(SampleId) %>% count)
View(samplesWithPrev)

changedSamples = rbind(samplesWithNew %>% group_by(SampleId) %>% count %>% ungroup(),
                       samplesWithPrev %>% group_by(SampleId) %>% count %>% ungroup()) %>% group_by(SampleId) %>% count

View(changedSamples)
write.csv(changedSamples,'~/data/sv/fusions/changed_sample_ids.csv',row.names = F, quote = F)

View(newFusions %>% filter(!(SampleGenePair %in% prevFusions$SampleGenePair)) %>%
       select(SampleId,KnownType,GeneNameUp,GeneNameDown,Reportable,PhaseMatched,RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown,
              TypeUp,SvIdUp,TypeDown,SvIdDown,ProteinsKept,ProteinsLost,BiotypeUp,BiotypeDown,
              ChainLinks,ChainLength,TerminatedUp,TerminatedDown,ExonsSkippedUp,ExonsSkippedDown,everything()))

requiredBiotypes = c('protein_coding','retained_intron','processed_transcript','nonsense_mediated_decay','lincRNA')

View(prevFusions %>% filter(!(SampleGenePair %in% newFusions$SampleGenePair)&!(SampleIdPair %in% newFusions$SampleIdPair)) %>%
       filter(BiotypeUp %in% requiredBiotypes) %>%
       # filter(Reportable=='true') %>%
       select(SampleId,KnownType,GeneNameUp,GeneNameDown,Reportable,PhaseMatched,RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown,
              TypeUp,SvIdUp,TypeDown,SvIdDown,ProteinsKept,ProteinsLost,BiotypeUp,BiotypeDown,
              ChainLinks,ChainLength,TerminatedUp,TerminatedDown,ExonsSkippedUp,ExonsSkippedDown,everything()))

View(prevFusionsAll %>% filter(!(SampleIdPair %in% newFusionsAll$SampleIdPair)&!(SampleGenePair %in% newFusionsAll$SampleGenePair)) %>%
       filter(BiotypeUp %in% requiredBiotypes) %>%
       # filter(Reportable=='true') %>%
       select(SampleId,KnownType,GeneNameUp,GeneNameDown,Reportable,PhaseMatched,RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown,
              TypeUp,SvIdUp,TypeDown,SvIdDown,ProteinsKept,ProteinsLost,BiotypeUp,BiotypeDown,
              ChainLinks,ChainLength,TerminatedUp,TerminatedDown,ExonsSkippedUp,ExonsSkippedDown,everything()))

View(newFusionsAll %>% filter(!(SampleIdPair %in% prevFusionsAll$SampleIdPair)&!(SampleGenePair %in% prevFusionsAll$SampleGenePair)) %>%
       filter(BiotypeUp %in% requiredBiotypes) %>%
       # filter(Reportable=='true') %>%
       select(SampleId,KnownType,GeneNameUp,GeneNameDown,Reportable,PhaseMatched,RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown,
              TypeUp,SvIdUp,TypeDown,SvIdDown,ProteinsKept,ProteinsLost,BiotypeUp,BiotypeDown,
              ChainLinks,ChainLength,TerminatedUp,TerminatedDown,ExonsSkippedUp,ExonsSkippedDown,everything()))

View(newFusionsAll %>% filter(SampleId=='CPCT02010240T') %>% filter(GeneNameUp=='TMEM178B'))
View(newFusionsAll %>% filter(SampleId=='CPCT02010240T') %>% filter(SvIdUp==863))
nrow(newFusionsAll %>% filter(SampleId=='CPCT02011102T'))
nrow(prevFusionsAll %>% filter(SampleId=='CPCT02010240T'))
View(newFusionsAll %>% filter(SampleId=='CPCT02011102T'))

View(newFusions %>% filter(SampleId=='CPCT02011102T') %>%
       select(SampleId,KnownType,GeneNameUp,GeneNameDown,Reportable,PhaseMatched,RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown,
              TypeUp,SvIdUp,TypeDown,SvIdDown,ProteinsKept,ProteinsLost,BiotypeUp,BiotypeDown,
              ChainLinks,ChainLength,TerminatedUp,TerminatedDown,ExonsSkippedUp,ExonsSkippedDown,everything()))

View(prevFusionsAll %>% filter(SampleId=='DRUP01010110T') %>%
       select(SampleId,KnownType,GeneNameUp,GeneNameDown,Reportable,PhaseMatched,RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown,
              TypeUp,SvIdUp,TypeDown,SvIdDown,ProteinsKept,ProteinsLost,BiotypeUp,BiotypeDown,
              ChainLinks,ChainLength,TerminatedUp,TerminatedDown,ExonsSkippedUp,ExonsSkippedDown,everything()))

View(newFusionsAll %>% filter(SampleId=='CPCT02011102T') %>%
       select(SampleId,KnownType,GeneNameUp,GeneNameDown,Reportable,PhaseMatched,RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown,
              TypeUp,SvIdUp,TypeDown,SvIdDown,ProteinsKept,ProteinsLost,BiotypeUp,BiotypeDown,
              ChainLinks,ChainLength,TerminatedUp,TerminatedDown,ExonsSkippedUp,ExonsSkippedDown,everything()))


# sample comparisons


sampleFusions = read.csv('~/logs/LNX_FUSIONS.csv')
View(sampleFusions)

sampleFusions = annotate_fusions(sampleFusions)

View(sampleFusions %>% filter(GeneNameUp=='YAP1'))
View(sampleFusions %>% group_by(SampleId) %>% count)

View(sampleFusions %>% filter(GeneNameUp=='YAP1') %>%
       select(SampleId,KnownType,GeneNameUp,GeneNameDown,Reportable,PriorityScore,PhaseMatched,TranscriptUp,TranscriptDown,RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown,
              TypeUp,SvIdUp,TypeDown,SvIdDown,ProteinsKept,ProteinsLost,BiotypeUp,BiotypeDown,CanonicalUp,CanonicalDown,
              ChainLinks,ChainLength,TerminatedUp,TerminatedDown,ExonsSkippedUp,ExonsSkippedDown,everything()))

View(sampleFusions %>% filter(GeneNameUp=='KCNB2') %>%
       select(SampleId,KnownType,GeneNameUp,GeneNameDown,Reportable,PriorityScore,PhaseMatched,TranscriptUp,TranscriptDown,RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown,
              TypeUp,SvIdUp,TypeDown,SvIdDown,ProteinsKept,ProteinsLost,BiotypeUp,BiotypeDown,CanonicalUp,CanonicalDown,
              ChainLinks,ChainLength,TerminatedUp,TerminatedDown,ExonsSkippedUp,ExonsSkippedDown,everything()))

View(sampleFusions %>% filter(GeneNameDown=='BRAF') %>%
       select(SampleId,KnownType,GeneNameUp,GeneNameDown,Reportable,PriorityScore,PhaseMatched,TranscriptUp,TranscriptDown,RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown,
                     TypeUp,SvIdUp,TypeDown,SvIdDown,ProteinsKept,ProteinsLost,BiotypeUp,BiotypeDown,
                     ChainLinks,ChainLength,TerminatedUp,TerminatedDown,ExonsSkippedUp,ExonsSkippedDown,everything()))

View(dnaRnaCombinedOutputData %>% filter(HmfId=='HMF001945A'))
View(dnaRnaCombinedOutputData %>% filter(HmfId=='HMF002774A'))

View(svaRnaFusions %>% filter(SampleId=='CPCT02010440T'))


dnaRnaCombinedDataNewPrio = dnaRnaCombinedData
dnaRnaCombinedDataOldPrio = dnaRnaCombinedData
dnaRnaCombinedDataOldPrioPMC = dnaRnaCombinedData
dnaRnaCombinedDataPreHomOffset = dnaRnaCombinedData

View(dnaRnaCombinedData %>% filter(MatchType=='DNA Only') %>% group_by(KnownCategory,
                                                                       ExonsSkipped=(ExonsSkippedUp.x>0|ExonsSkippedDown.x>0)) %>% 
       count %>% spread(KnownCategory,n))

View(dnaRnaCombinedDataNewPrio %>% filter(!is.na(Reportable)&KnownCategory!='Both promiscuous') %>% 
       group_by(MatchType,KnownCategory,ExonsSkipped=(ExonsSkippedUp.x>0|ExonsSkippedDown.x>0)) %>% 
       count %>% spread(KnownCategory,n))

View(dnaRnaCombinedDataOldPrio %>% filter(!is.na(Reportable)&KnownCategory!='Both promiscuous') %>% 
       group_by(MatchType,KnownCategory,ExonsSkipped=(ExonsSkippedUp.x>0|ExonsSkippedDown.x>0)) %>% 
       count %>% spread(knit_with_parametersparameters(file = samplesWithNew)ownCategory,n))

View(dnaRnaCombinedDataOldPrioPMC %>% filter(!is.na(Reportable)&KnownCategory!='Both promiscuous') %>% 
       group_by(MatchType,KnownCategory,ExonsSkipped=(ExonsSkippedUp.x>0|ExonsSkippedDown.x>0)) %>% 
       count %>% spread(KnownCategory,n))


View(dnaRnaCombinedData %>% filter(MatchType=='RNA Only'&KnownType=='Known'))
View(dnaRnaCombinedDataNewPrio %>% filter(MatchType=='RNA Only'&KnownType=='Known'))
View(dnaRnaCombinedDataNewPrio %>% filter(GeneNameUp=='NAB2'&GeneNameDown=='STAT6'))

View(paperDnaRna %>% filter(GeneNameUp=='NAB2'&GeneNameDown=='STAT6')) # MatchType=='RNA Only'&
colnames(paperDnaRna)

paperDnaRna = read.csv('~/data/sv/rna/rna_tmp/supptable3_LINX_dna_rna_fusion_comparison.csv')
paperDnaRna = paperDnaRna%>% mutate(HmdIdGenePair=paste(HmfId,GeneNameUp,GeneNameDown,sep='_'))

View(paperDnaRna %>% group_by(KnownCategory,MatchType) %>% count %>% spread(MatchType,n))
View(dnaRnaCombinedData %>% group_by(KnownType,MatchCategory) %>% count %>% spread(KnownType,n))

dnaRnaCombinedDataNewPrio = merge(dnaRnaCombinedDataNewPrio,sampleIdMapping,by='SampleId',all.x=T)
dnaRnaCombinedDataNewPrio = dnaRnaCombinedDataNewPrio %>% mutate(HmdIdGenePair=paste(HmfId,GeneNameUp,GeneNameDown,sep='_'))

View(dnaRnaCombinedDataNewPrio %>% filter(!(HmdIdGenePair %in% paperDnaRna$HmdIdGenePair)) %>%
       group_by(KnownType,MatchCategory,ExonsSkipped=(ExonsSkippedUp.x>0|ExonsSkippedDown.x>0)) %>% count %>% spread(ExonsSkipped,n))

View(dnaRnaCombinedDataNewPrio %>% filter(!(HmdIdGenePair %in% paperDnaRna$HmdIdGenePair)) %>%
       group_by(SampleId,HmfId) %>% count)

View(dnaRnaCombinedDataNewPrio %>% filter(!(HmdIdGenePair %in% paperDnaRna$HmdIdGenePair)&SampleId=='CPCT02010386T') %>%
       group_by(GeneNameUp,GeneNameDown) %>% count)

View(paperDnaRna %>% filter(!(HmdIdGenePair %in% dnaRnaCombinedDataNewPrio$HmdIdGenePair)) %>%
       group_by(KnownCategory,MatchType,ExonsSkipped=((!is.na(ExonsSkippedUp)&ExonsSkippedUp>0)|(!is.na(ExonsSkippedDown)&ExonsSkippedDown>0))) %>% count %>% spread(ExonsSkipped,n))

View(paperDnaRna %>% filter(is.na(ExonsSkippedUp)|is.na(ExonsSkippedDown)))

View(paperDnaRna %>% filter(!(HmdIdGenePair %in% dnaRnaCombinedDataNewPrio$HmdIdGenePair)) %>%
       group_by(HmfId) %>% count)

View(rnaMatchData %>% filter(SampleId=='CPCT02010386T') %>% group_by(GeneNameUp,GeneNameDown) %>% count)

View(dnaRnaCombinedData %>% filter(MatchType=='RNA Only'&KnownType=='Known'&GeneNameUp=='NAB2') %>%
       select(SampleId,GeneNameUp,GeneNameDown,PosUp.x,PosUp.y,RnaPosUp,OrientUp.x,OrientUp.y,StrandUp.x,TransViableUp,TransValidLocUp,
              SvIdUp.x,SvIdDown.x,SvIdUp.y,SvIdDown.y,TranscriptUp,RnaTransIdUp,SpliceType,everything()))

View(dnaRnaCombinedData %>% filter(MatchType=='RNA Only'&KnownType=='Known'&GeneNameUp=='NAB2') %>%
       select(SampleId,GeneNameUp,GeneNameDown,PosUp.x,RnaPosUp,OrientUp.x,StrandUp.x,TransValidLocUp,
              SvIdUp.x,SvIdDown.x,SpliceType,everything()))

tmpNew = read.csv('~/logs/LNX_FUSIONS.csv')
# tmpNew = tmpNew %>% filter()
View(tmpNew)

tmpOld = read.csv('~/logs/LNX_FUSIONS.csv')
View(tmpOld)

nrow(tmpOld %>% filter(PhaseMatched=='true'))
nrow(tmpNew %>% filter(PhaseMatched=='true'))

View(dnaRnaCombinedData %>% filter(is.na(ExonsSkippedUp.x)))
       
View(dnaRnaCombinedDataNewPrio %>% filter(!is.na(Reportable)&KnownCategory!='Both promiscuous') %>% 
       group_by(RegionTypeUp.x,RegionTypeDown.x) %>% count)
       

# Match type debug
View(rnaCombinedData %>% filter(HasDnaData&HasRnaData&!DnaRnaMatch))

View(rnaCombinedData %>% filter(HasDnaData&HasRnaData&!DnaRnaMatch) %>% group_by(DnaMatchType,SvMatchType,SvMatchUp,SvMatchDown) %>% count)

View(dnaRnaCombinedData %>% filter(!SameGeneFusion) %>% 
       group_by(KnownCategory,MatchType,MatchCategory,OldMatchType,OldMatchCategory) %>% count %>% spread(KnownCategory,n))


View(dnaRnaCombinedData %>% filter(OldMatchCategory=='DNA & RNA'&MatchCategory=='RNA Only'&KnownType.x=='Known'))

View(dnaRnaCombinedData %>% filter(!SameGeneFusion) %>% group_by(KnownCategory,MatchType) %>% count %>% spread(KnownCategory,n))
View(dnaRnaCombinedData %>% filter(!SameGeneFusion) %>% group_by(KnownCategory,MatchType,DnaMatchType,SvMatchType) %>% count %>% spread(KnownCategory,n))
View(dnaRnaCombinedData %>% filter(!SameGeneFusion) %>% group_by(KnownCategory,MatchType,HasDnaData,DnaMatchType,SvMatchType) %>% count %>% spread(KnownCategory,n))
View(dnaRnaCombinedData %>% filter(!SameGeneFusion) %>% group_by(KnownCategory,MatchType,DnaMatchType,SvMatchType) %>% count %>% spread(KnownCategory,n))
View(dnaRnaCombinedData %>% filter(!SameGeneFusion) %>% group_by(KnownCategory,MatchType,MatchCategory,HasDnaData,DnaMatchType,SvMatchType) %>% count %>% spread(KnownCategory,n))

colnames(svaRnaFusions)
View(svaRnaFusions %>% group_by(PhaseMatched,ValidChain) %>% count)

View(dnaRnaCombinedData %>% filter(DnaMatchType=='INVALID_Terminated',SvMatchType=='BothSVs'))
View(dnaRnaCombinedData %>% filter(is.na(MatchType)))






# -restricted_fusion_genes
# EGFR;FGD2;BRPF3

View(fusions %>% filter(RegionTypeUp=='Exonic'&RegionTypeDown=='Exonic') %>% group_by(KnownType) %>% count)
View(fusions %>% filter(RegionTypeUp=='Exonic'&RegionTypeDown=='Exonic') %>% select(ExactBaseUp,ExactBaseDown))

View(fusions %>% filter(RegionTypeUp=='Exonic'&RegionTypeDown=='Exonic'&CodingTypeUp=='5P_UTR'&CodingTypeDown=='5P_UTR'))


View(fusions %>% filter(RegionTypeUp=='Exonic'&RegionTypeDown=='Exonic'&PhaseMatched=='true'&CodingTypeUp=='Coding'&CodingTypeDown=='Coding') %>% 
       filter(!((ExactBaseUp==0&ExactBaseDown==1)|(ExactBaseUp==1&ExactBaseDown==2)|(ExactBaseUp==2&ExactBaseDown==0))) %>%
       select(ExactBaseUp,ExactBaseDown,everything()))





rm(newFusions)
rm(prevFusions)


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


