library(data.table)
library(dplyr)
library(tidyr)
library(stringi)


#####
# FUSION LIKELIHOOD

# inputs:
# - fusions from Linx - all valid fusions regardless of whether reportable or not
# - SV Data from Linx - to generate cohort counts by type and length bucket
# - fusion likelihood for proximate DELs and DUPs
# - fusion likelihood for long non-proximate DELs, DUPs and INVs, translocation and short INVs

lengthBucketMax = 4096e3
lengthBucketMin = 100
lengthLong = 5e6
longDelDupInv = 'LONG_DDI'
shortInv = 1e5

# short INVs will be given a bucket length range of 1000-2000, medium 100K-5M
# long DELs, DUPs and INVs will be given a bucket length of 5M

get_length<-function(Type,PosStart,PosEnd)
{
  length = ifelse(Type=='BND',0,abs(PosEnd-PosStart))
  return (length)
}

get_length_min<-function(Type,Length)
{
  lengthMin=ifelse(Type=='BND',0,
            ifelse(Length>lengthBucketMax,lengthLong,
            ifelse(Type=='INV',ifelse(Length<shortInv,1000,100000),
            ifelse(Length<1000,lengthBucketMin,1000*(2**floor(log(Length/1000,2))))
            )))
  
  return (lengthMin)
}

get_length_max<-function(Type,Length,LengthMin)
{
  lengthMax=ifelse(Type=='BND',0,
            ifelse(Length>lengthBucketMax,lengthLong,
            ifelse(Type=='INV',ifelse(Length<shortInv,shortInv,lengthLong),
            ifelse(Length<1000,1000,LengthMin*2)
            )))
  
  return (lengthMax)
}

# Tests
print(get_length('BND',100,200))
print(get_length('INV',100,200))
print(get_length_min('BND',0))
print(get_length_min('INV',1))
print(get_length_min('INV',2e5))
print(get_length_min('DEL',50))
print(get_length_min('DEL',1e5))
print(get_length_min('DEL',1e7))
print(get_length_max('BND',0,0))
print(get_length_max('INV',1,1000))
print(get_length_max('INV',2e5,1e5))
print(get_length_max('DEL',50,get_length_min('DEL',50)))
print(get_length_max('DEL',1e5,get_length_min('DEL',1e5)))
print(get_length_max('DEL',1e7,get_length_min('DEL',1e7)))

getBreakendsInGenes<-function(geneSvBreakends)
{
  genesListStart1 = geneSvBreakends %>% filter(GeneStart!='') %>% select(Id,Type,LengthMin,LengthMax,Orient=OrientStart,GeneStart)
  genesListStart2 = genesListStart1 %>% separate(GeneStart,sep=';',c('G1','G2','G3','G4','G5','G6','G7','G8','G9','G10'), fill = "right")
  genesListStart2[is.na(genesListStart2)] = 'NO_GENE'
  genesListStart3 = gather(genesListStart2, "GCol", "GeneName", 6:ncol(genesListStart2))
  # genesListStart3 = genesListStart3 %>% mutate(IsStart=TRUE)
  
  genesListEnd1 = geneSvBreakends %>% filter(GeneEnd!='') %>% select(Id,Type,LengthMin,LengthMax,Orient=OrientEnd,GeneEnd)
  genesListEnd2 = genesListEnd1 %>% separate(GeneEnd,sep=';',c('G1','G2','G3','G4','G5','G6','G7','G8','G9','G10'), fill = "right")
  genesListEnd2[is.na(genesListEnd2)] = 'NO_GENE'
  genesListEnd3 = gather(genesListEnd2, "GCol", "GeneName", 6:ncol(genesListEnd2))
  # genesListEnd3 = genesListEnd3 %>% mutate(IsStart=FALSE)
  
  genesList3 = rbind(genesListStart3,genesListEnd3)
  genesList3 = genesList3 %>% filter(GeneName!=''&GeneName!='NO_GENE')
  genesList3 = genesList3 %>% separate(GeneName,sep=':',c('GeneId','GeneName'))
  
  breakendsInGenes = genesList3 %>% group_by(GeneId,GeneName,Type,LengthMin,LengthMax,Orient) %>% summarise(Count=n())  
  
  return (breakendsInGenes)
}

# limit analysis to highest purity cohort and de-duped
svCohort = load('~/data/hmf_cohort_may_2019.RData')


# load actual fusions from cohort
# this includes reportable fusions and the single highest priority fusion for every other candidate gene pair and de-duped by SV Id
allFusions = read.csv('~/data/sv/fusion_like/LNX_FUSIONS.csv')
nrow(allFusions)
View(allFusions)

allFusions = allFusions %>% mutate(SameSV=(SvIdUp==SvIdDown),
                                   SameChromosome=(as.character(ChrUp)==as.character(ChrDown)),
                                   SameGene=as.character(GeneIdUp)==as.character(GeneIdDown),
                                   GenePair=paste(GeneNameUp,GeneNameDown,sep='_'),
                                   Type=ifelse(SameSV,as.character(TypeUp),
                                               ifelse(!SameChromosome,'BND',ifelse(OrientUp==OrientDown,'INV',ifelse((PosUp<PosDown)==(OrientUp<OrientDown),'DUP','DEL')))),
                                   PreTransDownDistance=ifelse(StrandDown==1&PosDown<TransStartDown,TransStartDown-PosDown,ifelse(StrandDown==-1&PosDown>TransEndDown,PosDown-TransEndDown,0)))

# apply filters to match FLC estimations
allFusions = allFusions %>% filter(SampleId %in% highestPurityCohort$sampleId) # de-dup multiple biopsy samples
allFusions = allFusions %>% filter(TypeUp!='INS'&TypeDown!='INS') # include inserts
allFusions = allFusions %>% filter(ResolvedType!='LINE'&ResolvedType!='DUP_BE') # exclude LINE and duplicate SVs
allFusions = allFusions %>% filter(SameSV&PhaseMatched=='true'&ExonsSkippedUp==0&ExonsSkippedDown==0)
nrow(allFusions)

maxPreGeneDistance=10e3
#View(allFusions %>% filter(PreTransDownDistance>maxPreGeneDistance))
#nrow(allFusions %>% filter(PreTransDownDistance==0))
#View(allFusions %>% filter(PreTransDownDistance>maxPreGeneDistance) %>% group_by(GenePair) %>% count)
allFusions = allFusions %>% filter(PreTransDownDistance<=maxPreGeneDistance)
# nrow(allFusions %>% filter((StrandDown==1&TransStartDown-PosDown>maxPreGeneDistance)|(StrandDown==-1&PosDown-TransEndDown>maxPreGeneDistance)))

allFusions = allFusions %>% mutate(Length=get_length(Type,PosUp,PosDown),
                                   LengthMin=get_length_min(Type,Length),
                                   LengthMax=get_length_max(Type,Length,LengthMin),
                                   Type=ifelse(LengthMin==lengthLong,longDelDupInv,as.character(Type)))

nrow(allFusions) # 32K


# DEBUG: filter for single chromosome
# allFusions = allFusions %>% filter(SameChromosome&ChrUp==filterChromosome)

allFusions = allFusions %>% mutate(SampleSvId=paste(SampleId,SvIdUp,sep='_'))

# de-dedup so that for each sample and type/length there is at most 1 fusion per gene-pair
allFusionsDedup = allFusions %>% group_by(SampleId,GenePair,GeneIdUp,GeneNameUp,GeneIdDown,GeneNameDown,Type,LengthMin,LengthMax) %>% count() %>%
  ungroup() %>% select(GenePair,GeneIdUp,GeneNameUp,GeneIdDown,GeneNameDown,Type,LengthMin,LengthMax)
  
nrow(allFusionsDedup) # 31814

View(allFusions %>% select(GeneNameUp,GeneNameDown,Type,Length,LengthMin,LengthMax,PosUp,PosDown,TypeUp,TypeDown))
View(allFusions %>% group_by(Type,LengthMin,LengthMax) %>% count())
View(allFusions %>% group_by(GeneNameUp,GeneNameDown) %>% count())
View(allFusions %>% filter(Type=='LONG_DDI') %>% group_by(GenePair) %>% count())

fivePrimeActuals = allFusionsDedup %>% group_by(GeneId=GeneIdUp,GeneName=GeneNameUp,Type,LengthMin,LengthMax) %>% 
  summarise(Count=n()) %>% mutate(Stream='Upstream')

threePrimeActuals = allFusionsDedup %>% group_by(GeneId=GeneIdDown,GeneName=GeneNameDown,Type,LengthMin,LengthMax) %>% 
  summarise(Count=n()) %>% mutate(Stream='Downstream')

#View(fivePrimeActuals)
#View(threePrimeActuals)

actualFusions = rbind(fivePrimeActuals,threePrimeActuals)
View(actualFusions)
nrow(actualFusions) # 43870

rm(threePrimeActuals)
rm(fivePrimeActuals)


#####
# Load Linx SV data
svData = read.csv('~/data/sv/fusions/LNX_SVS.csv')
# svData = read.csv('~/logs/SVA_SVS.csv')
# write.csv(highestPurityCohort %>% select(sampleId), '~/data/sv/hpc_non_dup_sample_ids.csv', row.names = F, quote = F)

filterChromosome=8

# apply filters
cohortSVs = svData %>% filter(Type %in% c('BND','INV','DEL','DUP'))
cohortSVs = cohortSVs %>% filter(SampleId %in% highestPurityCohort$sampleId)

# DEBUG: filter for single chromosome
cohortSVs = cohortSVs %>% filter(Type!='BND'&ChrStart==filterChromosome)

# View(cohortSVs %>% group_by(ResolvedType) %>% count)
cohortSVs = cohortSVs %>% filter(ResolvedType!='LINE'&ResolvedType!='DUP_BE')

# is this an attempt to filter out remote shards? if so why?
# cohortSVs = cohortSVs %>% filter(!(Type=='BND'&((LnkLenStart>1&LnkLenStart<1000)|(LnkLenEnd>1&LnkLenEnd<1000))))
# nrow(cohortSVs %>% filter(Type=='BND'&((LnkLenStart>1&LnkLenStart<1000)|(LnkLenEnd>1&LnkLenEnd<1000)))) # 49.7K

cohortSVs = cohortSVs %>% mutate(SampleSvId=paste(SampleId,Id,sep='_'))

cohortSVs = cohortSVs %>% mutate(InFusion=SampleSvId %in% allFusions$SampleSvId)

#View(cohortSVs %>% group_by(Type) %>% count())
cohortSVs = cohortSVs %>% mutate(Length=get_length(Type,PosStart,PosEnd),
                                 LengthMin=get_length_min(Type,Length),
                                 LengthMax=get_length_max(Type,Length,LengthMin),
                                 Type=ifelse(LengthMin==lengthLong,longDelDupInv,as.character(Type)))

View(cohortSVs %>% group_by(Type,LengthMin,LengthMax,InFusion) %>% count() %>% spread(InFusion,n))

View(cohortSVs)
nrow(cohortSVs)
#View(cohortSVs %>% filter(LengthMin==LengthMax&LengthMin>0&LengthMin<1e5))


# produce numbers of each category
cohortSvSummary = cohortSVs %>% group_by(Type,LengthMin,LengthMax) %>% 
  summarise(CountAll=n(),
            CountNoFusions=sum(!InFusion))

# decide which count to use to adjusted expected counts below
# cohortSvSummary = cohortSvSummary %>% mutate(CohortCount=CountNoPosSelection)
cohortSvSummary = cohortSvSummary %>% mutate(CohortCount=CountAll)

View(cohortSvSummary)

#write.csv(cohortSvSummary, '~/data/sv/fusion_like/cohort_bucket_counts.csv', row.names = F, quote = F)
#write.csv(cohortSVs, '~/data/sv/fusion_like/cohort_svs.csv', row.names = F, quote = F)
# cohortSvSummary = read.csv('~/data/sv/fusion_like/cohort_bucket_counts.csv')




#####
# Expected fusions for short INVs, same-arm long DELs, DUPs and INVs, and remote translocation fusions

## LOAD expected non-proxiate fusions
geneRangeData = read.csv('~/data/sv/fusion_like/GFL_GENE_DATA.csv')

# DEBUG: filter for single chromosome
geneRangeData = read.csv('~/logs/GFL_GENE_DATA.csv')
geneRangeData = geneRangeData %>% filter(Chromosome==filterChromosome)
bndCount = 0

nrow(geneRangeData)
View(geneRangeData)

# set these from the cohort summary
print(shortInv)
delDupInvCount = sum(cohortSvSummary %>% ungroup() %>% filter(Type==longDelDupInv&LengthMin==lengthLong&LengthMax==lengthLong) %>% select(CohortCount))
bndCount = sum(cohortSvSummary %>% ungroup() %>% filter(Type=='BND') %>% select(CohortCount))
shortInvCount = sum(cohortSvSummary %>% ungroup() %>% filter(Type=='INV'&LengthMin==1e3&LengthMax==shortInv) %>% select(CohortCount))
medInvCount = sum(cohortSvSummary %>% ungroup() %>% filter(Type=='INV'&LengthMin==shortInv&LengthMax==lengthLong) %>% select(CohortCount))
print(delDupInvCount)
print(bndCount)
print(medInvCount)
print(shortInvCount)
# for non-proximate rates, multiple these by the cohort counts

# convert expect rates into bucket categories
bndExpCounts = rbind(geneRangeData %>% select(GeneId,GeneName,RemoteRateUp) %>% mutate(Stream='Upstream') 
                     %>% mutate(ExpFusion=RemoteRateUp*bndCount) %>% select(GeneId,GeneName,Stream,ExpFusion),
                     geneRangeData %>% select(GeneId,GeneName,RemoteRateDown) %>% mutate(Stream='Downstream')
                     %>% mutate(ExpFusion=RemoteRateDown*bndCount) %>% select(GeneId,GeneName,Stream,ExpFusion)) %>%
  mutate(Type='BND',LengthMin=0,LengthMax=0)

# View(bndExpCounts)
sum(bndExpCounts$ExpFusion)

sameArmCounts = rbind(geneRangeData %>% select(GeneId,GeneName,LongDDIRateUp) %>% mutate(Stream='Upstream') 
                      %>% mutate(ExpFusion=LongDDIRateUp*delDupInvCount) %>% select(GeneId,GeneName,Stream,ExpFusion),
                      geneRangeData %>% select(GeneId,GeneName,LongDDIRateDown) %>% mutate(Stream='Downstream')
                      %>% mutate(ExpFusion=LongDDIRateDown*delDupInvCount) %>% select(GeneId,GeneName,Stream,ExpFusion)) %>%
  mutate(Type=longDelDupInv,LengthMin=lengthLong,LengthMax=lengthLong)

#View(sameArmCounts)
sum(sameArmCounts$ExpFusion)

shortInvCounts = rbind(geneRangeData %>% select(GeneId,GeneName,ShortInvRateUp) %>% mutate(Stream='Upstream') 
                       %>% mutate(ExpFusion=ShortInvRateUp*shortInvCount) %>% select(GeneId,GeneName,Stream,ExpFusion),
                       geneRangeData %>% select(GeneId,GeneName,ShortInvRateDown) %>% mutate(Stream='Downstream')
                       %>% mutate(ExpFusion=ShortInvRateDown*shortInvCount) %>% select(GeneId,GeneName,Stream,ExpFusion)) %>%
  mutate(Type='INV',LengthMin=1000,LengthMax=shortInv)

# View(shortInvCounts)
sum(shortInvCounts$ExpFusion)

medInvCounts = rbind(geneRangeData %>% select(GeneId,GeneName,MedInvRateUp) %>% mutate(Stream='Upstream') 
                       %>% mutate(ExpFusion=MedInvRateUp*medInvCount) %>% select(GeneId,GeneName,Stream,ExpFusion),
                       geneRangeData %>% select(GeneId,GeneName,MedInvRateDown) %>% mutate(Stream='Downstream')
                       %>% mutate(ExpFusion=MedInvRateDown*medInvCount) %>% select(GeneId,GeneName,Stream,ExpFusion)) %>%
  mutate(Type='INV',LengthMin=shortInv,LengthMax=lengthLong)

# View(medInvCounts)
sum(medInvCounts$ExpFusion)


#####
# Proximate DELs and DUP expected fusion rates

## LOAD expected non-proxiate fusions
proximateDelDups = read.csv('~/data/sv/fusion_like/GFL_DEL_DUP_PROXIMATES.csv')

# DEBUG for single chromosome
proximateDelDups = read.csv('~/logs/GFL_DEL_DUP_PROXIMATES.csv')
proximateDelDups = proximateDelDups %>% filter(Chromosome==filterChromosome)

proximateDelDups = merge(proximateDelDups,cohortSvSummary,by=c('Type','LengthMin','LengthMax'),all.x=T)
proximateDelDups = proximateDelDups %>% mutate(ExpFusion=CohortCount*ProximateRate)
View(proximateDelDups)
# View(proximateDelDups %>% group_by(Type,LengthMin,LengthMax) %>% summarise(Genes=n(),ExpFusion=sum(ExpFusion)))
# proximateDelDupSummary = proximateDelDups %>% group_by(GeneIdUp,GeneNameUp,GeneIdDown,GeneNameDown,Type) %>% summarise(Probability=sum(Count*ProximateRate))

# split the proximate fusions into their 5' and 3' genes
fivePrimeDelDupExpCounts = proximateDelDups %>% select(GeneId=GeneIdUp,GeneName=GeneNameUp,Type,LengthMin,LengthMax,ExpFusion) %>% mutate(Stream='Upstream')
threePrimeDelDupExpCounts = proximateDelDups %>% select(GeneId=GeneIdDown,GeneName=GeneNameDown,Type,LengthMin,LengthMax,ExpFusion) %>% mutate(Stream='Downstream')

# de-dep for same gene fusions
fivePrimeDelDupExpCounts = fivePrimeDelDupExpCounts %>% group_by(GeneId,GeneName,Type,LengthMin,LengthMax,Stream) %>% summarise(ExpFusion=sum(ExpFusion))
threePrimeDelDupExpCounts = threePrimeDelDupExpCounts %>% group_by(GeneId,GeneName,Type,LengthMin,LengthMax,Stream) %>% summarise(ExpFusion=sum(ExpFusion))
fivePrimeDelDupExpCounts = fivePrimeDelDupExpCounts %>% ungroup()
threePrimeDelDupExpCounts = data.frame(threePrimeDelDupExpCounts %>% ungroup())


# combine all expected fusions into a single dataset
allExpFusion = rbind(bndExpCounts,sameArmCounts)
allExpFusion = rbind(allExpFusion,shortInvCounts)
allExpFusion = rbind(allExpFusion,medInvCounts)
allExpFusion = rbind(allExpFusion,fivePrimeDelDupExpCounts)
allExpFusion = rbind(allExpFusion,threePrimeDelDupExpCounts)

rm(proximateDelDups)
rm(bndExpCounts)
rm(shortInvCounts)
rm(medInvCounts)
rm(fivePrimeDelDupExpCounts)
rm(threePrimeDelDupExpCounts)

View(allExpFusion)
# write.csv(geneExpFusions,'~/data/sv/fusions/expected_fusion_rates.csv', row.names = F, quote = F)

#####
# Calculate observed vs expected counts of breakends in genes

geneSvBreakends = cohortSVs %>% select(SampleId,Id,Type,ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd,GeneStart,GeneEnd)

geneSvBreakends = geneSvBreakends %>% mutate(Length=get_length(Type,PosStart,PosEnd),
                                             LengthMin=get_length_min(Type,Length),
                                             LengthMax=get_length_max(Type,Length,LengthMin),
                                             Type=ifelse(LengthMin==lengthLong,longDelDupInv,as.character(Type)))

breakendsInGenes = getBreakendsInGenes(geneSvBreakends)
#View(head(breakendsInGenes,1000))

View(breakendsInGenes %>% filter(GeneName=='CSMD1'&Type=='DEL'))

# use Ensembl to set strand and stream

# need to form a matrix with every gene from Ensembl with its length (and strand if required)
# filled out for every category in the cohort summary
# and populated with actual gene-breakend counts
ensemblGeneData = read.csv('~/data/sv/ensembl_gene_data.csv')
nrow(ensemblGeneData)

# DEBUG single chromosome only
ensemblGeneData = ensemblGeneData %>% filter(Chromosome==filterChromosome)
rm(ensemblGeneDataAll)

#transSpliceData = read.csv('~/data/sv/ensembl_trans_splice_data.csv')
#nrow(transSpliceData)
# View(transSpliceData %>% filter(PreSpliceAcceptorPosition<(-1)))

#preGeneDistances = transSpliceData %>% filter(PreSpliceAcceptorPosition>=0) %>% group_by(GeneId) %>% 
#  summarise(PreDistanceMax=max(Distance)) %>% mutate(PreDistanceMax=ifelse(PreDistanceMax<maxPreGeneDistance,PreDistanceMax,maxPreGeneDistance))
# View(preGeneDistances)

fullGenome=3e9
genomeLength=fullGenome

# DEBUG for single chromosome
genomeLength=146364022 # chr 8
genomeLength=81195210 # chr 17
genomeLength=48129895 # chr 21

ensemblGeneData = ensemblGeneData %>% mutate(GeneLength=GeneEnd-GeneStart,
                                             ExpBreakendRate=round(GeneLength/fullGenome,12))

geneFrequencies = merge(ensemblGeneData %>% select(GeneId,GeneName,ExpBreakendRate,Strand),cohortSvSummary,all=T)
geneFrequencies = geneFrequencies %>% select(GeneId,GeneName,ExpBreakendRate,Strand,Type,LengthMin,LengthMax,CohortCount)

geneFrequencies2 = merge(geneFrequencies,breakendsInGenes,by=c('GeneId','GeneName','Type','LengthMin','LengthMax'),all.x=T)
geneFrequencies2[is.na(geneFrequencies2)] = 0
geneFrequencies2 = geneFrequencies2 %>% mutate(Stream=ifelse((Strand==1&Orient==1)|(Strand==-1&Orient==-1),'Upstream','Downstream'))
# View(head(geneFrequencies2,1000))

View(geneFrequencies2 %>% filter(GeneName=='CSMD1'&Type=='DEL'))

geneExpAndActCounts = geneFrequencies2 %>% group_by(GeneId,GeneName,Type,LengthMin,LengthMax,Stream) %>% 
  summarise(ActBeCount=sum(Count),
            CohortCount=first(CohortCount),
            ExpBreakendRate=first(ExpBreakendRate))

rm(breakendsInGenes)
rm(geneSvBreakends)
rm(geneFrequencies)
rm(geneFrequencies2)


geneExpAndActCounts = geneExpAndActCounts %>% mutate(ExpBeCount=ExpBreakendRate*CohortCount*2) # since each SV has 2 breakends
View(geneExpAndActCounts)
View(head(geneExpAndActCounts,100))
nrow(geneExpAndActCounts)

View(geneExpAndActCounts %>% group_by(Type,LengthMin,LengthMax) %>% 
       summarise(TotalExpBeCount=round(sum(ExpBeCount)),
                 TotalActBeCount=round(sum(ActBeCount))))


# merge actual and expected fusions
allFusionData = merge(allExpFusion,actualFusions %>% ungroup() %>% select(GeneId,GeneName,Type,LengthMin,LengthMax,ActFusions=Count,Stream),
                      by=c('GeneId','GeneName','Type','LengthMin','LengthMax','Stream'),all=T)

allFusionData[is.na(allFusionData)] = 0
# View(allFusionData)

# finally combined expected & actual SV breakends in genes with expected & actual fusions
combinedGeneData = merge(geneExpAndActCounts %>% select(GeneId,GeneName,Type,LengthMin,LengthMax,Stream,ActBeCount,ExpBeCount),
                         allFusionData,by=c('GeneId','GeneName','Type','LengthMin','LengthMax','Stream'),all=T)

combinedGeneData = merge(combinedGeneData,ensemblGeneData %>% select(GeneId,GeneStart,GeneEnd,GeneLength),by='GeneId',all.x=T)
combinedGeneData[is.na(combinedGeneData)] = 0

combinedGeneData = combinedGeneData %>% 
  mutate(ExpSvCount=ifelse(ExpBeCount>ExpFusion,ExpBeCount-ExpFusion,ExpFusion),
         AdjExpFusions=ifelse(ExpBeCount>0,round((ActBeCount-ActFusions)*ExpFusion/ExpBeCount,4),ExpFusion))

combinedGeneData = combinedGeneData %>% 
  mutate(Prob=ifelse(ActFusions>AdjExpFusions&AdjExpFusions>0,1-ppois(ActFusions-1,AdjExpFusions,T),1))

View(combinedGeneData %>% filter(Prob<0.001) %>% arrange(Prob))

rm(allExpFusion)
rm(allFusionData)
rm(actualFusions)


## SAVE key data sets
write.csv(combinedGeneData, '~/data/sv/fusion_like/fusion_likelihood_summary.csv',row.names=F,quote=F)
combinedGeneData = read.csv('~/data/sv/fusion_like/fusion_likelihood_summary.csv')

write.csv(cohortSVs,'~/data/sv/fusion_like/cohort_sv_data.csv',row.names=F,quote=F)
cohortSVs = read.csv('~/data/sv/fusion_like/cohort_sv_data.csv')

write.csv(geneExpAndActCounts,'~/data/sv/fusion_like/gene_exp_act_counts.csv',row.names=F,quote=F)
geneExpAndActCounts = read.csv('~/data/sv/fusion_like/gene_exp_act_counts.csv')

# write.csv(cohortSVs,'~/data/sv/fusion_like/', row.names = F, quote = F)

# View(combinedGeneData %>% filter(is.na(Prob)))
#View(combinedGeneData %>% filter(Prob<0.001))



#####
# SUMMARY ANALYSIS

# Promiscous partners
knownPairs = read.csv('~/data/knownFusionPairs.csv')
knownPairs = knownPairs %>% mutate(GenePair=paste(fiveGene,threeGene,sep='_'))
threePrimeProm = read.csv('~/data/knownPromiscuousThree.csv')
fivePrimeProm = read.csv('~/data/knownPromiscuousFive.csv')
View(knownPairs)
View(threePrimeProm)
View(fivePrimeProm)

# annotate fragile sites
fragileSiteGenes = read.csv('~/data/sv/drivers/fragile_site_genes.csv')
View(fragileSiteGenes)
topTsgGenes = read.csv('~/data/sv/drivers/KnownTSGs.csv')
View(topTsgGenes)

combinedGeneData = combinedGeneData %>% mutate(FragileSite=GeneName %in% fragileSiteGenes$GeneName)
combinedGeneData = combinedGeneData %>% mutate(KnownTSG=GeneName %in% topTsgGenes$GeneName)
combinedGeneData = combinedGeneData %>% mutate(ThreePrimeProm=Stream=='Downstream'&GeneName %in% threePrimeProm$GeneName)
combinedGeneData = combinedGeneData %>% mutate(FivePrimeProm=Stream=='Upstream'&GeneName %in% fivePrimeProm$GeneName)


# calculate a probability across all buckets for each gene
geneSummary = combinedGeneData %>% group_by(GeneName,Stream,GeneLength,FragileSite,KnownTSG,ThreePrimeProm,FivePrimeProm) %>% 
  summarise(ActFusions=sum(ActFusions),
            AdjExpFusions=round(sum(AdjExpFusions),4)) %>%
  mutate(Prob=ifelse(ActFusions>AdjExpFusions&AdjExpFusions>0,1-ppois(ActFusions-1,AdjExpFusions,T),1))

View(geneSummary %>% filter(Prob<0.001))
View(geneSummary %>% filter(!FragileSite&Prob<0.001) %>% arrange(Prob))
View(geneSummary %>% filter(FivePrimeProm&ActFusions>0) %>% arrange(Prob))
View(geneSummary %>% filter(ThreePrimeProm&ActFusions>0) %>% arrange(Prob))


View(combinedGeneData %>% filter(Prob<0.001&!FragileSite&!ThreePrimeProm&!FivePrimeProm&!KnownTSG) %>% arrange(Prob))






## DEBUG ONLY

View(combinedGeneData %>% filter(GeneName=='CSMD1'&Type=='DEL'))
View(combinedGeneData %>% filter(GeneName=='PTEN'&Type=='DEL'))
View(combinedGeneData %>% filter(GeneName=='ANKRD11'&Type=='DEL'))
View(combinedGeneData %>% filter(GeneName=='STAT6'&Type=='INV'))
View(combinedGeneData %>% filter(GeneName=='ERBB4'&Type=='LONG_DDI'))

View(allFusionData %>% filter(GeneName=='STAT6'))
View(allFusions %>% filter(GeneNameDown=='STAT6'))

View(allFusionData %>% filter(Type=='INV'&LengthMin==1e3))
View(geneRangeData %>% filter(GeneName=='STAT6'))
View(geneExpAndActCounts %>% filter(GeneName=='STAT6'))
View(geneExpAndActCounts)





View(combinedGeneData)
View(head(combinedGeneData,100))

# expectedGeneFreq = expectedGeneFreq %>% mutate(PoisProb=1-ppois(Count-1,ExpCount))

View(combinedGeneData %>% group_by(Type) %>% summarise(TotalActFusions=sum(ActFusions),
                                                       TotalExpFusions=sum(ExpFusion),
                                                       TotalAdjExpFusions=sum(AdjExpFusion)))

View(combinedGeneData %>% filter(Stream=='Downstream') %>% 
       group_by(Type,LengthMin,LengthMax) %>% 
       summarise(TotalActFusions=sum(ActFusions),
                 TotalExpFusions=round(sum(ExpFusion)),
                 TotalAdjExpFusions=round(sum(AdjExpFusion)),
                 TotalExpSvCount=round(sum(ExpSvCount)),
                 TotalActSvCount=round(sum(ActSvCount))))


# excluding genes connected with fragile sites and line elements
View(combinedGeneData %>% filter(Stream=='Downstream') %>% 
       filter(!(GeneId %in% repeatingNonFusionGenes$GeneId)) %>%
       filter(!(GeneId %in% fragileSiteGenes$GeneId)&!(GeneId %in% lineElementGenes$GeneId)) %>%
       group_by(Type,LengthMin,LengthMax) %>% 
       summarise(TotalActFusions=sum(ActFusions),
                 TotalExpFusions=round(sum(ExpFusion)),
                 TotalAdjExpFusions=round(sum(AdjExpFusion)),
                 TotalExpSvCount=round(sum(ExpSvCount)),
                 TotalActSvCount=round(sum(ActSvCount))))


# analysis of genes with excess expected or actuals
View(combinedGeneData %>% 
       # filter(Type=='BND') %>%
       filter(Type=='DEL'&LengthMin==100) %>% 
       #filter(!(GeneId %in% fragileSiteGenes$GeneId)&!(GeneId %in% lineElementGenes$GeneId)) %>%
       #filter(!(GeneId %in% repeatingNonFusionGenes$GeneId)) %>%
       filter((AdjExpFusion>ActFusions+5)|(ActFusions>0.5&AdjExpFusion>ActFusions*2)) %>%
       select(Stream,GeneId,GeneName,GeneLength,Type,LengthMin,LengthMax,ExpSvCount,ActSvCount,ExpFusion,TotalExpFusion,ActFusions,AdjExpFusion))

View(combinedGeneData %>% filter(Stream=='Downstream') %>% 
       filter(!(GeneId %in% repeatingNonFusionGenes$GeneId)) %>%
       filter(!(GeneId %in% fragileSiteGenes$GeneId)&!(GeneId %in% lineElementGenes$GeneId)) %>%
       filter((AdjExpFusion>ActFusions+5)|(ActFusions>0.5&AdjExpFusion>ActFusions*2)) %>%
       group_by(Type,LengthMin,LengthMax) %>% 
       summarise(TotalActFusions=sum(ActFusions),
                 TotalExpFusions=round(sum(ExpFusion)),
                 TotalAdjExpFusions=round(sum(AdjExpFusion)),
                 TotalExpSvCount=round(sum(ExpSvCount)),
                 TotalActSvCount=round(sum(ActSvCount))))

