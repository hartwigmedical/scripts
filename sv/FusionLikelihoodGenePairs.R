## SPECIFIC PAIR FUSIONS

# GENE PAIR REPORT
# - GenePair
# - Type and Length bucket
# - Gene locations and length
# - Expected vs Actual fusions
# - act SVs in 1M based window of either gene, actual SVs touching both genes in length bucket
# - connection with a known fusion or promiscuous partner
# - connection with a fragile site or driver AMP or DEL
# - any recurrent fused exons or contexts
# - any recurrent cancer types
# - any repeated positions (eg from contamination)
# - any repeated cluster types

annotate_drivers<-function(genePairData,driverPositions)
{
  genePairData$NearDriver = FALSE
  genePairData$DriverGenes = ''
  colIndex = ncol(genePairData)
  
  for(i in 1:nrow(genePairData))
  {
    gpData = genePairData[i,]
    genePair = gpData$GenePair
    chrUp = gpData$ChrUp
    chrDown = gpData$ChrDown
    geneStartUp = gpData$GeneStartUp
    geneEndUp = gpData$GeneEndUp
    geneStartDown = gpData$GeneStartDown
    geneEndDown = gpData$GeneEndDown
    
    print(paste(i,': ',genePair,' - adding driver annotation',sep=''))
    
    drivers = driverPositions %>% filter((as.character(Chromosome)==as.character(chrUp)&!(geneStartUp>GeneEnd+1e6|geneEndUp<GeneStart-1e6))
                                         |(as.character(Chromosome)==as.character(chrDown)&!(geneStartDown>GeneEnd+1e6|geneEndDown<GeneStart-1e6))
                                         |(as.character(Chromosome)==as.character(chrUp)&as.character(Chromosome)==as.character(chrDown)&
                                             ((geneEndDown<GeneStart&geneStartUp>GeneEnd)|(geneEndUp<GeneStart&geneStartDown>GeneEnd))))
    
    if(nrow(drivers)>0)
    {
      genePairData[i,which(colnames(genePairData)=='NearDriver')] = TRUE

      driverGenes = ''
      for(j in 1:pmin(nrow(drivers),3))
      {
        dgData = drivers[j,]
        gene = as.character(dgData$Gene)
        driverGenes = ifelse(driverGenes=='',gene,paste(driverGenes,gene,sep=';'))
      }
      
      genePairData[i,which(colnames(genePairData)=='DriverGenes')] = driverGenes
    }
  }
  
  return (genePairData)
}

# tmp = annotate_drivers(head(lowProbPairs,10),ampDriverPositions)
#View(tmp)

annotate_proximate_svs<-function(genePairData,cohortSVs,log=F)
{
  genePairData$LocalSvCount = 0
  genePairData$GenePairSvCount = 0
  genePairData$ExpGenePairCount=0
  genePairData$ExpLocalFusions=0
  
  localSvThreshold=5e5

  typeLengths = genePairData %>% group_by(Type,LengthMin,LengthMax) %>% count
  genePairResults = data.frame(matrix(ncol=ncol(genePairData), nrow=0))
  index = 1
  
  for(j in 1:nrow(typeLengths))
  {
    tlData = typeLengths[j,]
    type = tlData$Type
    lengthMin = tlData$LengthMin
    lengthMax = tlData$LengthMax
    
    gpSubset = genePairData %>% filter(Type==type&LengthMin==lengthMin&LengthMax==lengthMax)
    cohortSubset = cohortSVs %>% filter(Type==type&LengthMin==lengthMin&LengthMax==lengthMax)

    for(i in 1:nrow(gpSubset))
    {
      gpData = gpSubset[i,]
      genePair = gpData$GenePair
      chrUp = gpData$ChrUp
      chrDown = gpData$ChrDown
      geneIdUp = gpData$GeneIdUp
      geneIdDown = gpData$GeneIdDown
      geneNameUp = gpData$GeneNameUp
      geneNameDown = gpData$GeneNameDown
      geneStartUp = gpData$GeneStartUp
      geneEndUp = gpData$GeneEndUp
      geneStartDown = gpData$GeneStartDown
      geneEndDown = gpData$GeneEndDown
      type = as.character(gpData$Type)
      lengthMin = gpData$LengthMin
      lengthMax = gpData$LengthMax
      
      genePairCount = 0
      localCount = 0
      
      print(paste(index,': ',genePair,' - adding local SV info',sep=''))
      index = index + 1
      
      # count number of SVs hitting both these genes
      # matchingSVs = cohortSubset %>% filter(Type==type&LengthMin==lengthMin&LengthMax==lengthMax)
      
      genePairCount = nrow(cohortSubset %>% 
                             filter((grepl(geneIdUp,GeneStart)&grepl(geneIdDown,GeneEnd))|(grepl(geneIdDown,GeneStart)&grepl(geneIdUp,GeneEnd))))

      geneMidpointUp = (geneStartUp+geneEndUp)/2
      geneMidpointDown = (geneStartDown+geneEndDown)/2

      rangeStartUp = pmin(geneStartUp,geneMidpointUp-localSvThreshold)
      rangeEndUp = pmax(geneEndUp,geneMidpointUp+localSvThreshold)
      rangeStartDown = pmin(geneStartDown,geneMidpointDown-localSvThreshold)
      rangeEndDown = pmax(geneEndDown,geneMidpointDown+localSvThreshold)
      
      # count number of SVs in vacinity
      cohortSubset = cohortSubset %>% mutate(StartInGeneUp=ChrStart==chrUp&PosStart>rangeStartUp&PosStart<rangeEndUp,
                                           StartInGeneDown=ChrStart==chrDown&PosStart>rangeStartDown&PosStart<rangeEndDown,
                                           EndInGeneUp=ChrEnd==chrUp&PosEnd>rangeStartUp&PosEnd<rangeEndUp,
                                           EndInGeneDown=ChrEnd==chrDown&PosEnd>rangeStartDown&PosEnd<rangeEndDown)
      
      localCount = nrow(cohortSubset %>% filter((StartInGeneUp|StartInGeneDown)&(EndInGeneUp|EndInGeneDown)))
  
      gpSubset[i,which(colnames(gpSubset)=='LocalSvCount')] = localCount
      gpSubset[i,which(colnames(gpSubset)=='GenePairSvCount')] = genePairCount
      
      sameChromosome = (as.character(chrUp)==as.character(chrDown))
      
      distanceBetweenRanges = ifelse(sameChromosome,pmax(rangeStartUp,rangeStartDown)-pmin(rangeEndUp,rangeEndDown),0)
      
      localSvDistance = ifelse(!sameChromosome|distanceBetweenRanges>0,
                               (rangeEndUp-rangeStartUp)+(rangeEndDown-rangeStartDown),
                               pmax(rangeEndUp,rangeEndDown)-pmin(rangeStartUp,rangeStartDown))
      
      if(log)
      {
        print(sprintf("geneUp(%s %d-> %d len=%d) geneDown(%s %d-> %d len=%d) distance(%d) localSvDistance(%d)", 
                      geneNameUp, geneStartUp, geneEndUp, geneEndUp-geneStartUp, 
                      geneNameDown, geneStartDown, geneEndDown, geneEndDown-geneStartDown, 
                      distanceBetweenRanges, localSvDistance))
      }
      
      gpSubset = gpSubset %>% mutate(ExpGenePairCount=GenePairRate*fullGenome/localSvDistance*LocalSvCount,
                                             ExpLocalFusions=FusionRate*fullGenome/localSvDistance*LocalSvCount)
    }
    
    genePairResults = rbind(genePairResults,gpSubset)
  }

  return (genePairResults)
}

# tmp = annotate_proximate_svs(lowProbPairs %>% filter(GenePair %in% c('TMPRSS2_ERG','AUTS2_AUTS2')),cohortSVs,T)
# View(tmp)


annotate_fusion_data<-function(genePairData,unfilteredFusions,byTyeAndLength=T)
{
  genePairData$CancerTypes = ""
  genePairData$FusedExons = ""
  genePairData$FusionContexts = ""
  genePairData$FusionClusters = ""
  genePairData$ExonicCount = 0
  genePairData$RepeatedPositions = F
  
  colIndex = ncol(genePairData)
  
  for(i in 1:nrow(genePairData))
  {
    gpData = genePairData[i,]
    genePair = gpData$GenePair
    
    print(paste(i,': ',genePair,' - adding fusion annotations',sep=''))

    if(byTyeAndLength)
    {
      type = as.character(gpData$Type)
      lengthMin = gpData$LengthMin
      lengthMax = gpData$LengthMax
      gpFusions = unfilteredFusions %>% filter(GenePair==genePair&Type==type&LengthMin==lengthMin&LengthMax==lengthMax)
    }
    else
    {
      gpFusions = unfilteredFusions %>% filter(GenePair==genePair)
    }
    
    fusionCount = nrow(gpFusions)
    
    if(fusionCount>0)
    {
      # cancer types
      fusionsByType = gpFusions %>% group_by(CancerType) %>% count() %>% arrange(-n)
      
      ctInfo = ""
      
      for(j in 1:nrow(fusionsByType))
      {
        cancerType = fusionsByType[j,1]
        count = fusionsByType[j,2]
        if(count > 0.25*fusionCount)
        {
          ctStr = paste(cancerType,'=',round(count/fusionCount*100),'%',sep='')
          ctInfo = ifelse(j==1,ctStr,paste(ctInfo,ctStr,sep=';'))
        }
        else
        {
          break
        }
      }
      
      genePairData[i,which(colnames(genePairData)=='CancerTypes')] = ctInfo
      
      # exonic fusions 
      genePairData[i,which(colnames(genePairData)=='ExonicCount')] = nrow(gpFusions %>% filter(RegionTypeUp=='Exonic'|RegionTypeDown=='Exonic'))
  
      # fused exons
      fusionsByType = gpFusions %>% group_by(FusedExonUp,FusedExonDown) %>% count() %>% arrange(-n)
      
      exonInfo = ""
      
      for(j in 1:nrow(fusionsByType))
      {
        fbt = fusionsByType[j,]
        count = fbt$n
        
        if(count > 0.25*fusionCount)
        {
          exonStr = paste(fbt$FusedExonUp,'->',fbt$FusedExonDown,'=',round(count/fusionCount*100),'%',sep='')
          exonInfo = ifelse(j==1,exonStr,paste(exonInfo,exonStr,sep=';'))
        }
        else
        {
          break
        }
      }
      
      genePairData[i,which(colnames(genePairData)=='FusedExons')] = exonInfo
      
      # coding context
      fusionsByType = gpFusions %>% group_by(CodingTypeUp,CodingTypeDown) %>% count() %>% arrange(-n)
      
      ctInfo = ""
      
      for(j in 1:nrow(fusionsByType))
      {
        fbt = fusionsByType[j,]
        count = fbt$n
        
        if(count > 0.25*fusionCount)
        {
          str = paste(fbt$CodingTypeUp,'->',fbt$CodingTypeDown,'=',round(count/fusionCount*100),'%',sep='')
          ctInfo = ifelse(j==1,str,paste(ctInfo,str,sep=';'))
        }
        else
        {
          break
        }
      }
      
      genePairData[i,which(colnames(genePairData)=='FusionContexts')] = ctInfo
      
      # resolved type
      fusionsByType = gpFusions %>% group_by(ResolvedType) %>% count() %>% arrange(-n)
      
      rtInfo = ""
      
      for(j in 1:nrow(fusionsByType))
      {
        fbt = fusionsByType[j,]
        resolvedType = fbt$ResolvedType
        count = fbt$n
        
        if(count > 0.25*fusionCount)
        {
          str = paste(resolvedType,'=',round(count/fusionCount*100),'%',sep='')
          rtInfo = ifelse(j==1,str,paste(rtInfo,str,sep=';'))
        }
        else
        {
          break
        }
      }
      
      genePairData[i,which(colnames(genePairData)=='FusionClusters')] = rtInfo
      
      colIndex = which(colnames(genePairData)=='RepeatedPositions')
      
      # check for repeated positions
      fusionsByType = gpFusions %>% group_by(PosUp=round(PosUp,-1),PosDown=round(PosDown,-1)) %>% count() %>% filter(n>1) %>% arrange(-n)
      
      if(nrow(fusionsByType) > 0)
      {
        for(j in 1:nrow(fusionsByType))
        {
          fbt = fusionsByType[j,]
          count = fbt$n
          
          if(count > 0.25*fusionCount)
          {
            genePairData[i,colIndex] = TRUE
            break
          }
        }
      }
    }
  }
  
  return (genePairData)
}


# 1. Get counts of all fusion, including those usually filtered out for FLC analysis
unfilteredFusions = read.csv('~/data/sv/fusion_like/LNX_FUSIONS.csv')

unfilteredFusions = unfilteredFusions %>% mutate(SameSV=(SvIdUp==SvIdDown),
                                                 SameChromosome=(as.character(ChrUp)==as.character(ChrDown)),
                                                 SameGene=as.character(GeneIdUp)==as.character(GeneIdDown),
                                                 GenePair=paste(GeneNameUp,GeneNameDown,sep='_'),
                                                 Type=ifelse(SameSV,as.character(TypeUp),
                                                             ifelse(!SameChromosome,'BND',ifelse(OrientUp==OrientDown,'INV',ifelse((PosUp<PosDown)==(OrientUp<OrientDown),'DUP','DEL')))),
                                                 PreTransDownDistance=ifelse(StrandDown==1&PosDown<TransStartDown,TransStartDown-PosDown,ifelse(StrandDown==-1&PosDown>TransEndDown,PosDown-TransEndDown,0)))

unfilteredFusions = unfilteredFusions %>% filter(SampleId %in% highestPurityCohort$sampleId) # de-dup multiple biopsy samples

unfilteredFusions = unfilteredFusions %>% mutate(Length=get_length(Type,PosUp,PosDown),
                                                 LengthMin=get_length_min(Type,Length),
                                                 LengthMax=get_length_max(Type,Length,LengthMin),
                                                 Type=ifelse(LengthMin==lengthLong,longDelDupInv,as.character(Type)))

unfilteredFusions = unfilteredFusions %>% filter(TypeUp!='INS'&TypeDown!='INS') # include inserts
unfilteredFusions = unfilteredFusions %>% filter(ResolvedType!='LINE'&ResolvedType!='DUP_BE') # exclude LINE and duplicate SVs
nrow(unfilteredFusions)

maxPreGeneDistance=10e3
unfilteredFusions = unfilteredFusions %>% mutate(ExceedsDownPreGeneLimit=((StrandDown==1&TransStartDown-PosDown>maxPreGeneDistance)
                                                                          |(StrandDown==-1&PosDown-TransEndDown>maxPreGeneDistance)),
                                                 Filtered=PhaseMatched=='false'|ExonsSkippedUp>0|ExonsSkippedDown>0|ExceedsDownPreGeneLimit)

sampleCancerTypes = highestPurityCohort %>% select(SampleId=sampleId,CancerType=cancerType)

unfilteredFusions = merge(unfilteredFusions,sampleCancerTypes,by='SampleId',all.x=T)
View(unfilteredFusions %>% group_by(CancerType) %>% count)
View(unfilteredFusions %>% group_by(Type,LengthMin,LengthMax) %>% count)
View(unfilteredFusions %>% filter(PhaseMatched=='true') %>% group_by(KnownType,ExonsSkipped=(ExonsSkippedUp>0|ExonsSkippedDown>0)|ExceedsDownPreGeneLimit) %>% count)
View(unfilteredFusions %>% filter(Filtered) %>% group_by(KnownType,Filtered) %>% count)
View(unfilteredFusions %>% filter(!Filtered) %>% group_by(SameSV) %>% count)
View(unfilteredFusions)
# View(unfilteredFusions %>% filter(KnownType==''&(ExonsSkippedUp>0|ExonsSkippedDown>0)))

View(unfilteredFusions %>% group_by(GenePair,SampleId) %>% count %>% filter(n>1))


# 1b
# Produce a set of gene-pairs to calculate likelihood for
# - known-fusion pairs
# - any fusion in our cohort occurring 2+ times

cohortGenePairSet = unfilteredFusions %>% group_by(SampleId,GenePair,GeneIdUp,GeneIdDown,KnownType) %>% count() %>% 
  group_by(GenePair,GeneIdUp,GeneIdDown,KnownType) %>% summarise(SampleCount=n()) %>% filter(KnownType=='Known'|SampleCount>=2)
cohortGenePairSet = cohortGenePairSet %>% ungroup()

View(cohortGenePairSet)
nrow(cohortGenePairSet)

# manually add in the 6 cosmic blacklisted  fusions
blacklistedCosmic = read.csv('~/data/sv/fusion_like/blacklisted_cosmic_fusions.csv')
blacklistedCosmic = merge(blacklistedCosmic,ensemblGeneData %>% select(GeneIdUp=GeneId,GeneName),by.x='GeneNameUp',by.y='GeneName',all.x=T)
blacklistedCosmic = merge(blacklistedCosmic,ensemblGeneData %>% select(GeneIdDown=GeneId,GeneName),by.x='GeneNameDown',by.y='GeneName',all.x=T)
blacklistedCosmic = blacklistedCosmic %>% mutate(GenePair=paste(GeneNameUp,GeneNameDown,sep='_'),
                                                 KnownType='Blacklist',
                                                 SampleCount=0)

cohortGenePairSet = cohortGenePairSet %>% filter(!(GenePair %in% blacklistedCosmic$GenePair))

# blacklistedCosmic = blacklistedCosmic %>% filter(!(GenePair %in% cohortGenePairSet$GenePair))
View(blacklistedCosmic)


genePairSet = rbind(cohortGenePairSet,blacklistedCosmic %>% select(GenePair,GeneIdUp,GeneIdDown,KnownType,SampleCount))
View(genePairSet)


# known pairs - for now only including known pairs which feature at least once in our cohort
knownPairs = read.csv('~/data/knownFusionPairs.csv')
knownFusionPairs = merge(knownPairs,ensemblGeneData %>% select(GeneIdUp=GeneId,GeneName),by.x='FiveGene',by.y='GeneName',all.x=T)
knownFusionPairs = merge(knownFusionPairs,ensemblGeneData %>% select(GeneIdDown=GeneId,GeneName),by.x='ThreeGene',by.y='GeneName',all.x=T)
knownFusionPairs = knownFusionPairs %>% filter(!(is.na(GeneIdUp)|is.na(GeneIdDown)))
knownFusionPairs = knownFusionPairs %>% mutate(GenePair=paste(FiveGene,ThreeGene,sep='_'))
View(knownFusionPairs)
write.csv(knownFusionPairs %>% select(GeneIdUp,GeneIdDown),'~/data/sv/fusion_like/knownFusionPairs.csv',row.names = F, quote = F)

knownFusionPairs = knownFusionPairs %>% select(GenePair,GeneIdUp,GeneIdDown) %>% mutate(KnownType='Known',SampleCount=1)
nrow(knownFusionPairs)
View(knownFusionPairs %>% filter(!(GenePair %in% cohortGenePairSet$GenePair)))

genePairSet = rbind(genePairSet,knownFusionPairs %>% filter(!(GenePair %in% cohortGenePairSet$GenePair)))
nrow(genePairSet)
View(genePairSet)

write.csv(genePairSet %>% group_by(GeneIdUp,GeneIdDown) %>% count() %>% select(GeneIdUp,GeneIdDown), # check for duplicates
          '~/data/sv/fusion_like/cohort_gene_pairs.csv',row.names = F, quote = F)

View(genePairSet %>% group_by(KnownType) %>% count())

View(genePairSet %>% group_by(GeneIdUp,GeneIdDown) %>% count())


# 2. Load FLC calc data - expected fusion and gene-pair SV rates
flcSpecificPairs = read.csv('~/data/sv/fusion_like/GFL_FUSION_PAIR_LIKELIHOOD.csv')
# flcSpecificPairs = read.csv('~/logs/GFL_FUSION_PAIR_LIKELIHOOD.csv')
View(flcSpecificPairs)
nrow(flcSpecificPairs)

flcSpecificPairs = flcSpecificPairs %>% mutate(GenePair=paste(GeneNameUp,GeneNameDown,sep='_'),
                                               FusionVsGenePairRateRatio=ifelse(GenePairRate>0,round(FusionRate/GenePairRate,4),0))

# work out actual fusion counts per pair
# if a sample has more than 1 fusion for a given pair, take the shortest
gpActFusions = unfilteredFusions %>% filter(GenePair %in% genePairSet$GenePair) %>% 
  filter(!Filtered) %>% 
  group_by(SampleId,GenePair,Type,LengthMin,LengthMax) %>% count() %>% 
  arrange(SampleId,GenePair,LengthMin) %>% group_by(GenePair,SampleId) %>% 
  summarise(Type=first(Type),LengthMin=first(LengthMin),LengthMax=first(LengthMax)) %>%
  group_by(GenePair,Type,LengthMin,LengthMax) %>% summarise(ActFusions=n())

View(gpActFusions)

gpActFusionsUnfiltered = unfilteredFusions %>% filter(GenePair %in% genePairSet$GenePair) %>% 
  group_by(SampleId,GenePair,Type,LengthMin,LengthMax) %>% count() %>% 
  arrange(SampleId,GenePair,LengthMin) %>% group_by(GenePair,SampleId) %>% 
  summarise(Type=first(Type),LengthMin=first(LengthMin),LengthMax=first(LengthMax)) %>%
  group_by(GenePair,Type,LengthMin,LengthMax) %>% summarise(ActFusionsNF=n())

#gpActFusionsUnfiltered = unfilteredFusions %>% filter(GenePair %in% genePairSet$GenePair) %>% 
#  group_by(SampleId,GenePair,Type,LengthMin,LengthMax) %>% count() %>% group_by(GenePair,Type,LengthMin,LengthMax) %>% summarise(ActFusionsNF=n())

gpActFusions = merge(gpActFusionsUnfiltered,gpActFusions,,by=c('GenePair','Type','LengthMin','LengthMax'),all.x=T)
gpActFusions[is.na(gpActFusions)] = 0

View(gpActFusions)
View(gpActFusions %>% filter(ActFusions>1))

lowProbPairs = within(lowProbPairs,rm(ActFusions))
lowProbPairs = within(lowProbPairs,rm(ActFusionsNF))
lowProbPairs = merge(lowProbPairs,gpActFusions,by=c('GenePair','Type','LengthMin','LengthMax'),all.x=T)
View(lowProbPairs)

flcSpecificPairs = merge(flcSpecificPairs,gpActFusions,by=c('GenePair','Type','LengthMin','LengthMax'),all.x=T)

# merge in actual observed counts by type
flcSpecificPairs = merge(flcSpecificPairs,cohortSvSummary,by=c('Type','LengthMin','LengthMax'),all.x=T)
flcSpecificPairs[is.na(flcSpecificPairs)] = 0

# View(flcSpecificPairs %>% group_by(Type,LengthMin,LengthMax) %>% count())

flcSpecificPairs = flcSpecificPairs %>% 
  mutate(ExpFusions=FusionRate*CohortCount,
         Prob=1-ppois(ActFusions-1,ExpFusions,T))

View(flcSpecificPairs %>% filter(Prob<0.01) %>% arrange(Prob))
View(flcSpecificPairs)

# append gene info - location and length
flcSpecificPairs = merge(flcSpecificPairs,ensemblGeneData %>% select(GeneIdUp=GeneId,ChrUp=Chromosome,GeneStartUp=GeneStart,GeneEndUp=GeneEnd,GeneLengthUp=GeneLength),by='GeneIdUp',all.x=T)
flcSpecificPairs = merge(flcSpecificPairs,ensemblGeneData %>% select(GeneIdDown=GeneId,ChrDown=Chromosome,GeneStartDown=GeneStart,GeneEndDown=GeneEnd,GeneLengthDown=GeneLength),by='GeneIdDown',all.x=T)

# filter out long same-gene fusions
# flcSpecificPairs = flcSpecificPairs %>% filter((as.character(GeneIdUp)!=as.character(GeneIdDown))|(GeneLengthUp<500000))

flcSpecificPairs = flcSpecificPairs %>% mutate(ReverseGenePair=paste(GeneNameDown,GeneNameUp,sep='_'))

oncoGenes = read.csv('~/data/sv/fusion_like/reportable_oncogenes.csv')
View(oncoGenes)

#  annotate fusion with known gene and fusion data
flcSpecificPairs = flcSpecificPairs %>% mutate(
  KnownPair=GenePair %in% knownFusionPairs$GenePair, 
  BlacklistCosmic=GenePair %in% blacklistedCosmic$GenePair, 
  RevKnownPair=ReverseGenePair %in% knownFusionPairs$GenePair, 
  ThreePrimeProm=GeneNameDown %in% threePrimeProm$GeneName,
  FivePrimeProm=GeneNameUp %in% fivePrimeProm$GeneName,
  FragileSite=(GeneNameUp %in% fragileSiteGenes$GeneName)|(GeneNameDown %in% fragileSiteGenes$GeneName),
  ThreePrimeOncogene=(GeneNameDown %in% oncoGenes$GeneName),
  KnownTSG=(GeneNameUp %in% topTsgGenes$GeneName)|(GeneNameDown %in% topTsgGenes$GeneName))

View(flcSpecificPairs)

# limit to those observed 1+ times
flcSpecificPairs = flcSpecificPairs %>% filter(ActFusionsNF>=1|KnownPair|BlacklistCosmic)
nrow(flcSpecificPairs) # 10509

# for now remove the probability filter
lowProbPairs = flcSpecificPairs %>% filter(ActFusionsNF>=2|KnownPair|BlacklistCosmic)
nrow(lowProbPairs) # 4633
View(lowProbPairs)

# to lower the set being evaluated, ignore any unknown gene-pairs with high probability
# View(lowProbPairs %>% filter(!KnownPair&!ThreePrimeProm&!FivePrimeProm) %>% filter(Prob>0.01) %>% arrange(Prob))
lowProbPairs = lowProbPairs %>% filter(KnownPair|ThreePrimeProm|FivePrimeProm|FragileSite|BlacklistCosmic|Prob<0.01)
nrow(lowProbPairs) # 1446
View(lowProbPairs)

# check driver catalog
drivers = read.csv('~/data/sv/drivers/LNX_DRIVERS.csv')
ampDrivers = drivers %>% filter((DriverType=='AMP'|DriverType=='DEL')&ClusterId>=0)
ampDrivers = merge(ampDrivers,ensemblGeneData %>% select(Gene=GeneName,GeneStart,GeneEnd),by='Gene',all.x=T)
ampDriverPositions = ampDrivers %>% group_by(Gene,Chromosome,GeneStart,GeneEnd) %>% summarise(DriverCount=n())
View(ampDriverPositions)
rm(drivers)

lowProbPairs = annotate_drivers(lowProbPairs,ampDriverPositions)

# add in SV counts around and inside the genes
lowProbPairs = annotate_proximate_svs(lowProbPairs,cohortSVs)
nrow(lowProbPairs)
View(lowProbPairs)

# add in any repeated exon-fusings or dominant cancer types

# - any recurrent fused exons or contexts
# - any recurrent cancer types
# - any repeated cluster types
lowProbPairs = annotate_fusion_data(lowProbPairs,unfilteredFusions)

lowProbPairs = lowProbPairs %>% mutate(ActFusions=ifelse(is.na(ActFusions),0,ActFusions),
                                       ActFusionsNF=ifelse(is.na(ActFusionsNF),0,ActFusionsNF),
                                       ExpFusions=FusionRate*CohortCount,
                                       Prob=ifelse(ActFusions>0,1-ppois(ActFusions-1,ExpFusions,T),1),
                                       LocalProb=ifelse(ActFusions>0,1-ppois(ActFusions-1,ExpLocalFusions,T),1),
                                       GenePairProb=1-ppois(GenePairSvCount-1,ExpGenePairCount,T))



View(lowProbPairs)
View(lowProbPairs %>% filter(is.na(ActFusions)))
# colnames(lowProbPairs)
write.csv(lowProbPairs,'~/data/sv/fusion_like/fusion_like_gene_pairs_by_length.csv', quote = F, row.names = F)

View(lowProbPairs %>% select(GenePair,Type,LengthMin,LengthMax,Prob,LocalProb,ActFusions,ExpFusions,ExpLocalFusions,CancerTypes,FusedExons,FusionClusters,everything())
     %>% arrange(LocalProb))

View(lowProbPairs %>% filter(!KnownPair&!RevKnownPair&!FragileSite&!NearDriver) %>% arrange(LocalProb) %>%
       mutate(Length=paste(LengthMin/1000,'-',LengthMax/1000,'K',sep='')) %>%
       select(GenePair,Type,Length,Prob,LocalProb,ActFusions,CancerTypes,FusedExons,FusionClusters,everything()))

googleSheetData = lowProbPairs %>% filter(!KnownPair&!RevKnownPair&!FragileSite) %>% filter(LocalProb<1e-4) %>% arrange(LocalProb) %>%
       mutate(Length=paste(LengthMin/1000,'-',LengthMax/1000,'K',sep='')) %>%
       select(GenePair,Type,Length,Prob,LocalProb,ActFusions,CancerTypes,FusedExons,FusionClusters,DriverGenes,everything())

View(googleSheetData)
write.csv(googleSheetData %>% select(GenePair,Type,Length,Prob,LocalProb,ActFusions,CancerTypes,FusedExons,FusionClusters,DriverGenes),
          '~/data/sv/fusion_like/gene_pair_google_sheets.csv', quote = F, row.names = F)


# Summary across all buckets

colnames(lowProbPairs)
nrow(genePairSummary)

genePairSummary = lowProbPairs %>% group_by(GenePair) %>% 
  summarise(ActFusions=sum(ActFusions),
            ActFusionsNF=sum(ActFusionsNF),
            ExpFusions=sum(FusionRate*CohortCount),
            ExpLocalFusions=sum(FusionRate*LocalSvCount),
            GenePairSvCount=sum(GenePairSvCount),
            KnownPair=first(KnownPair|RevKnownPair),
            Blacklisted=first(BlacklistCosmic),
            ThreePrimeProm=first(ThreePrimeProm),
            FivePrimeProm=first(FivePrimeProm),
            FragileSite=first(FragileSite),
            NearDriver=first(NearDriver),
            DriverGenes=first(DriverGenes),
            ThreePrimeOncogene=first(ThreePrimeOncogene),
            KnownTSG=first(KnownTSG),
            GeneStartUp=first(GeneStartUp),
            GeneEndUp=first(GeneEndUp),
            GeneLengthUp=first(GeneLengthUp),
            GeneLengthDown=first(GeneLengthDown),
            GeneStartDown=first(GeneStartDown),
            GeneEndDown=first(GeneEndDown)) %>%
  mutate(Prob=1-ppois(ActFusions-1,ExpFusions,T),
         LocalProb=1-ppois(ActFusions-1,ExpLocalFusions,T))

genePairSummary = annotate_fusion_data(genePairSummary,unfilteredFusions,F)


View(genePairSummary %>% filter(ActFusionsNF>=2|(KnownPair&ActFusionsNF>=1)|((KnownPair|ThreePrimeProm|FivePrimeProm)&ActFusionsNF>=1)) %>%
       arrange(LocalProb))

View(genePairSummary)

write.csv(genePairSummary %>% filter(ActFusionsNF>=2|(KnownPair&ActFusionsNF>=1)|((KnownPair|ThreePrimeProm|FivePrimeProm)&ActFusionsNF>=1)) %>%
       arrange(LocalProb),
       '~/data/sv/fusion_like/fusion_like_gene_pairs.csv', row.names = F, quote = )

View(knownFusionPairs)


View(lowProbPairs)
View(lowProbPairs %>% filter(!KnownPair) %>% filter(ActFusions!=ActFusionsNF))

# Investigate indications that fusions are passengers
# - expect a high rate of SV activity around the genes
# - show that SVs that fall into the gene are likely to cause an inframe fusion by virtue of the intron/exon ratios
# - how many SVs in the genes are disruptive? 
# - how localised are the fusions - eg always same exons being duped or deleted
# - does fusion occur in context of a complex cluster which further amplifies or disrupts the gene? disruptions should be picked up already

sameGeneLowProbs = lowProbPairs %>% filter(as.character(GeneIdUp)==as.character(GeneIdDown)&LocalProb<0.001&!RepeatedPositions)

View(sameGeneLowProbs %>% 
       select(GenePair,Type,LengthMin,LengthMax,Prob,LocalProb,ActFusions,ExpFusions,ExpLocalFusions,LocalSvCount,GenePairSvCount,ExpGenePairCount,GenePairProb,CohortCount,CancerTypes,FusedExons,FusionClusters,FusionClusters,FragileSite,NearDriver,everything())
     %>% arrange(LocalProb))

View(drivers)

View(lowProbPairs %>% filter(as.character(GeneIdUp)==as.character(GeneIdDown)&LocalProb>0.01&!RepeatedPositions) %>% 
       select(GenePair,Type,LengthMin,LengthMax,Prob,LocalProb,ActFusions,ExpFusions,ExpLocalFusions,LocalSvCount,GenePairSvCount,ExpGenePairCount,CohortCount,CancerTypes,FusedExons,FusionClusters,FusionClusters,FragileSite,NearDriver,everything())
     %>% arrange(LocalProb))



# manual validation of SV and gene-pair rates
tmpSVs = cohortSVs %>% filter(Type=='DEL'&LengthMin==100&LengthMax==1000&ChrStart==2)
nrow(tmpSVs)
nrow(tmpSVs %>% filter(grepl('ENSG00000155657',GeneStart)))
nrow(tmpSVs %>% filter(grepl('ENSG00000155657',GeneEnd)))
nrow(tmpSVs %>% filter(grepl('ENSG00000155657',GeneEnd)&grepl('ENSG00000155657',GeneStart)))

View(unfilteredFusions %>% filter(GenePair=='OLFM2_OLFM2'&Type=='DEL'&LengthMin==100&LengthMax==1000) %>%
       select(BreakendExonUp,BreakendExonDown,RegionTypeUp,RegionTypeDown,Length,PosUp,PosDown,ProteinsKept,ProteinsLost,everything()))


View(unfilteredFusions %>% filter(GenePair=='OLFM2_OLFM2'))

View(unfilteredFusions %>% filter(Type=='DUP'&FusedExonUp==1&FusedExonDown==2&as.character(GeneIdUp)==as.character(GeneIdDown)
                                  &ExonsSkippedUp==0&ExonsSkippedDown==0&PhaseMatched=='true'&RegionTypeUp=='Intronic'))


View(unfilteredFusions %>% filter(GenePair=='AHR_AHR'&Type=='DEL'&LengthMin==2000&LengthMax==4000))
View(unfilteredFusions %>% filter(GenePair=='AHR_AHR'&Type=='DEL'&LengthMin==2000&LengthMax==4000) %>% group_by(ResolvedType) %>% count() %>% arrange(-n))
View(unfilteredFusions %>% filter(GenePair=='AHR_AHR'&Type=='DEL'&LengthMin==2000&LengthMax==4000) %>% group_by(PosUp=round(PosUp,-1),PosDown=round(PosDown,-1)) %>% count() %>% arrange(-n))
fusionsByType = gpFusions %>% group_by(PosUp=round(PosUp,-1),PosDown=round(PosDown,-1)) %>% count() %>% arrange(-n)

ensemblTransExonData = read.csv('~/data/sv/ensembl_trans_exon_data.csv')
ensemblTransExonData = ensemblTransExonData %>% filter(Trans %in% unfilteredFusions$TranscriptUp|Trans %in% unfilteredFusions$TranscriptDown)
nrow(ensemblTransExonData)
rm(ensemblTransExonData)


View(lowProbPairs %>% filter(as.character(GeneIdUp)==as.character(GeneIdDown)) %>% group_by(GeneId=GeneIdUp,Type) %>% count())
View(unfilteredFusions %>% filter(SameGene))



generate_intron_data<-function(exonData)
{
  intronExonData = data.frame(matrix(ncol=5, nrow=0))
  colnames(intronExonData) = c('Type','PosStart','PosEnd','Phase','ExonRank')
  exonCount = nrow(exonData)-1 
  
  transLength = max(exonData$ExonEnd) - min(exonData$ExonStart)
  exonCount = nrow(exonData)
  minExonWidth = round(transLength*0.01)
  
  for(i in 1:exonCount)
  {
    exonRow = i*2 - 1
    exonPrev=exonData[i,]
    exonNext=exonData[i+1,]

    intronRow = exonRow+1
    exonMidpoint = (exonPrev$ExonStart+exonPrev$ExonEnd)/2
    intronExonData[exonRow,1] = 'Exon'
    intronExonData[exonRow,2] = exonPrev$ExonStart # pmin(exonMidpoint-minExonWidth,exonPrev$ExonStart)
    intronExonData[exonRow,3] = exonPrev$ExonEnd # pmax(exonMidpoint+minExonWidth,exonPrev$ExonEnd)
    intronExonData[exonRow,4] = ''
    intronExonData[exonRow,5] = exonPrev$ExonRank
    
    if(i < nrow(exonData))
    {
      intronExonData[intronRow,1] = 'Intron'
      intronExonData[intronRow,2] = exonPrev$ExonEnd
      intronExonData[intronRow,3] = exonNext$ExonStart
      intronExonData[intronRow,4] = ifelse(exonPrev$Strand==1,exonPrev$ExonEndPhase,exonNext$ExonEndPhase)
      intronExonData[intronRow,5] = ''
    }
  }
  
  return (intronExonData)
}

exonData = ensemblTransExonData %>% filter(Trans=='ENST00000300305')
View(exonData)
intronExonData = generate_intron_data(exonData)
View(intronExonData %>% mutate(Length=PosEnd-PosStart))


exonData = ensemblTransExonData %>% filter(Trans=='ENST00000275493')
print(generate_fusion_intron_plot('ENST00000275493','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # EGFR
print(generate_fusion_intron_plot('ENST00000300305','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # RUNX1


generate_fusion_intron_plot<-function(transId,svType,unfilteredFusions,ensemblTransExonData,cohortSVs)
{
  fusions = unfilteredFusions %>% filter(TranscriptUp==transId&TranscriptDown==transId&Type==svType&!Filtered) %>% 
    arrange(SampleId,Length) %>% group_by(SampleId,CancerType,GeneIdUp,GeneNameUp,ChrUp) %>%
    summarise(PosUp=first(PosUp),PosDown=first(PosDown),Length=first(Length)) %>%
    arrange(-Length) %>% ungroup()
    
  if(nrow(fusions) == 0)
  {
    print("no fusions found")
    return(fusions)
  }
  
  print(sprintf("%s type=%s: %d fusions found: ",transId, svType, nrow(fusions)))
  
  rowIndex = data.frame(as.numeric(as.character(rownames(fusions))))
  colnames(rowIndex) <- c("FusionIndex")
  fusions = cbind(fusions,rowIndex)
  
  geneId = fusions[1,which(colnames(fusions)=='GeneIdUp')]
  geneName = fusions[1,which(colnames(fusions)=='GeneNameUp')]
  geneChr = fusions[1,which(colnames(fusions)=='ChrUp')]

  print(sprintf("gene(%s:%s) chr=%s",geneId, geneName, geneChr))
  
  exonData = ensemblTransExonData %>% filter(Trans==transId)
  intronExonData = generate_intron_data(exonData)
  
  intronExonData = intronExonData %>% mutate(BlockColour=ifelse(Type=='Exon','black',ifelse(Phase==-1,'grey',
                                                        ifelse(Phase==0,'yellow3',
                                                               ifelse(Phase==1,'darkolivegreen3','skyblue2'))))) # tomato2
  
  transStart = min(exonData$ExonStart)
  transEnd = max(exonData$ExonEnd)
  transLength = transEnd-transStart
  
  geneSVs = cohortSVs %>% filter(Type==svType&grepl(geneId,GeneStart)&grepl(geneId,GeneEnd)&!(SampleId %in% fusions$SampleId))
  geneSVs = geneSVs %>% filter(PosStart>=transStart&PosEnd<=transEnd) %>% arrange(Length)
  
  print(paste('SVs in gene-pair=', nrow(geneSVs),sep=''))
  
  rowIndex = data.frame(as.numeric(as.character(rownames(geneSVs))))
  colnames(rowIndex) <- c("SvIndex")
  geneSVs = cbind(geneSVs,rowIndex)
  
  phaseColours = c('Red','Yellow','Green','Grey')
  exonHeight = 5
  sampleBuffer=transLength*0.02

  plot = (ggplot() 
          + geom_point(data=fusions, aes(x=PosDown,y=exonHeight*(1.1+FusionIndex),colour='Fusion Downstream'))
          + geom_point(data=fusions, aes(x=PosUp,y=exonHeight*(1.1+FusionIndex),colour='Fusion Upstream'))
          + geom_text(data=fusions, aes(x=pmax(PosDown,PosUp)+sampleBuffer,y=exonHeight*(1.1+FusionIndex),label=paste(SampleId,CancerType,sep=' '),hjust=0, vjust=0))
          + geom_rect(data=intronExonData, aes(xmin=PosStart, xmax=PosEnd, ymin=0, ymax=exonHeight), fill=intronExonData$BlockColour)
          + geom_text(data=intronExonData, aes(x=(PosStart+PosEnd)/2,y=-exonHeight/2,label=ExonRank))
          + labs(y='SampleId')
          + theme(axis.title.y = element_blank()) + theme(axis.ticks.y = element_blank()) + theme(axis.text.y = element_blank())
          + scale_color_hue("SV Positions") 
          # + scale_fill_manual("CI horizontal line", values=rep(1,4),guide=guide_legend(override.aes = list(colour=c("orange", "darkred"))),
          #                    labels=c("CI of 95%", "CI of 99%"))
          # + scale_fill_identity(name='Breakends',guide='legend') 
          # + scale_colour_manual(name='Phasing', values=phaseLegend)
          + ggtitle(sprintf("Same Gene Fusion: %s-%s type=%s length=%d", geneId, geneName, svType, transLength)))
  
  if(nrow(geneSVs) > 0)
  {
    plot = (plot 
            + geom_point(data=geneSVs, aes(x=PosStart,y=-exonHeight-SvIndex/nrow(geneSVs)*exonHeight,colour='Non-fusion SV Start'))
            + geom_point(data=geneSVs, aes(x=PosEnd,y=-exonHeight-SvIndex/nrow(geneSVs)*exonHeight,colour='Non-fusion SV End')))
  }
  
  return (plot)
}

print(generate_fusion_intron_plot('ENST00000242057','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # AHR
print(generate_fusion_intron_plot('ENST00000275493','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # EGFR
print(generate_fusion_intron_plot('ENST00000349496','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # CTNNB1
print(generate_fusion_intron_plot('ENST00000275493','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # BRAF

#print(generate_fusion_intron_plot('ENST00000300305','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # RUNX1
print(generate_fusion_intron_plot('ENST00000288135','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # KIT

View(unfilteredFusions %>% filter(GenePair=='EGFR_EGFR') %>% group_by(SampleId,Type,LengthMin,LengthMax) %>% count())
View(unfilteredFusions %>% filter(GenePair=='KIT_KIT'))

View(unfilteredFusions %>% filter(TranscriptUp=='ENST00000242057'&TranscriptDown=='ENST00000242057'&Type=='DEL'&!Filtered) %>% 
  arrange(SampleId,Length) %>% group_by(SampleId,GeneIdUp,GeneNameUp,ChrUp) %>%
  summarise(PosUp=first(PosUp),PosDown=first(PosDown),Length=first(Length)) %>%
  arrange(Length))

View(unfilteredFusions %>% filter(GenePair=='EGFR_EGFR') %>% arrange(SampleId,Length))


View(unfilteredFusions %>% filter(GenePair=='EGFR_EGFR') %>% arrange(SampleId,Length) %>%
       group_by(SampleId,GeneIdUp,GeneNameUp,ChrUp) %>%
       summarise(Count=n(),PosUp=first(PosUp),PosDown=first(PosDown),Length=first(Length)))

View(ensemblTransExonData %>% filter(ExonRank==1&Strand==1&ExonPhase==0))
View(ensemblTransExonData %>% filter(Trans=='ENST00000428131'))
View(ensemblTransExonData %>% filter(ExonRank==1&Strand==1&CodingEnd==ExonEnd))
View(ensemblTransExonData %>% filter(Strand==1&CodingEnd==ExonEnd))


print(generate_fusion_intron('ENST00000269305','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # TP53
print(generate_fusion_intron('ENST00000349496','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # CTNNB1
print(generate_fusion_intron('ENST00000426105','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # PUM1
print(generate_fusion_intron('ENST00000269701','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # AKAP8

print(generate_fusion_intron_plot('ENST00000288602','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # BRAF
print(generate_fusion_intron_plot('ENST00000288602','DUP',unfilteredFusions,ensemblTransExonData,cohortSVs)) # BRAF

print(generate_fusion_intron('ENST00000342788','DUP',unfilteredFusions,ensemblTransExonData,cohortSVs)) # ERBB4_ERBB4

print(generate_fusion_intron('ENST00000378933','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # EP300

print(generate_fusion_intron_plot('ENST00000433211','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # RBFOX1

View(unfilteredFusions %>% filter(GenePair=='KIT_KIT'))

sameGeneFusions = lowProbPairs %>% filter(as.character(GeneIdUp)==as.character(GeneIdDown)&!FragileSite&!RepeatedPositions) %>% 
  group_by(GeneId=GeneIdUp,Type,GeneNameUp) %>% summarise(ActFusions=sum(ActFusions),
                                                          LowestLocalProb=min(LocalProb),
                                                          GeneLength=first(GeneLengthUp)) %>% arrange(LowestLocalProb)

nrow(sameGeneFusions)
View(sameGeneFusions)

topN = 300
outputFile = '~/data/sv/fusion_like/same_gene_plots_top100.pdf'
pdf(file=outputFile, height = 14, width = 20)
par(mar=c(1,1,1,1))

for(i in 1:pmin(nrow(sameGeneFusions),topN))
{
  gpData = sameGeneFusions[i,]
  geneId = gpData$GeneId
  type= as.character(gpData$Type)
  
  gpFusions = head(unfilteredFusions %>% filter(as.character(GeneIdUp)==geneId&as.character(GeneIdDown)==geneId
                                                &CanonicalUp=='true'&!Filtered),1)
  transId = as.character(gpFusions[1,which(colnames(gpFusions)=='TranscriptUp')])

  print(sprintf("plotting %d: gene(%s) trans(%s)", i, geneId, transId))
  
  print(generate_fusion_intron_plot(transId,type,unfilteredFusions,ensemblTransExonData,cohortSVs))
}

dev.off()


rnaReadCounts = read.csv('~/data/sv/fusion_like/same_gene_rna_reads.csv')
View(rnaReadCounts)
View(rnaReadCounts %>% group_by(Gene,ExonJunction) %>% summarise(ReadCount=sum(ReadCount)))

tp53SVs = cohortSVs %>% filter(Type=='DEL'&grepl('TP53',GeneStart)&grepl('TP53',GeneEnd))
tp53SVs = cohortSVs %>% filter(Type=='DEL'&grepl('TP53',GeneStart)&grepl('TP53',GeneEnd)&!(SampleId %in% tp53Fusions$SampleId))
tp53SVs = tp53SVs %>% filter(PosStart>=7571720&PosEnd<=7590856) %>% arrange(Length)
View(tp53SVs)

rowIndex = data.frame(as.numeric(as.character(rownames(tp53SVs))))
colnames(rowIndex) <- c("SvIndex")
tp53SVs = cbind(tp53SVs,rowIndex)

tp53ExonData = ensemblTransExonData %>% filter(Trans=='ENST00000269305')
View(tp53ExonData)

View(ensemblTransExonData %>% filter(Trans=='ENST00000316626'))

tmp = unfilteredFusions %>% filter(GenePair=='TP53_TP53'&Type=='DEL') %>% 
  select(SampleId,ChrUp,PosUp,PosDown,Length,RegionTypeUp,RegionTypeDown,BreakendExonUp,BreakendExonDown,TranscriptUp,CancerType) %>%
  mutate(SvLength=PosUp-PosDown) %>% arrange(SvLength)

View(lowProbPairs %>% filter(GenePair=='NTRK2_NTRK2') %>% select(GenePair,Type,LengthMin,ActFusions,FusedExons,CancerTypes,everything()))

View(unfilteredFusions %>% filter(GenePair %in% c('KIT_KIT','GSK3B_GSK3B','ESR1_ESR1','BRAF_BRAF','RUNX1_RUNX1','PUM1_PUM1','EP300_EP300','PREX2_PREX2','NTRK2_NTRK2'))
     %>% select(SampleId,GenePair,BreakendExonUp,BreakendExonDown,RegionTypeUp,RegionTypeDown,LengthMin,ChrUp,TranscriptUp,TranscriptDown,everything()))

View(ensemblTransExonData %>% filter(Trans %in% c('ENST00000316626','ENST00000288602','ENST00000263253','ENST00000440973','ENST00000376214','ENST00000300305',
                                                  'ENST00000426105','ENST00000288368','ENST00000288135')) %>%
       select(TransId,Trans,Strand,ExonRank,ExonStart,ExonEnd,everything()))

View(ensemblGeneData)


View(unfilteredFusions %>% filter(RegionTypeUp=='Exonic'&RegionTypeDown=='Intronic'))

View(unfilteredFusions %>% filter(SameSV&PhaseMatched=='true') %>% filter(RegionTypeUp=='Exonic'&RegionTypeDown=='Intronic') %>%
       select(SampleId,GenePair,RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown,BreakendExonUp,BreakendExonDown,FusedExonUp,FusedExonDown,
              ExonsSkippedUp,ExonsSkippedDown,everything()))

View(unfilteredFusions %>% filter(SameSV&PhaseMatched=='true') %>% filter(RegionTypeUp=='Intronic'&RegionTypeDown=='Exonic') %>%
       select(SampleId,GenePair,RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown,BreakendExonUp,BreakendExonDown,FusedExonUp,FusedExonDown,
              ExonsSkippedUp,ExonsSkippedDown,everything()))

View(unfilteredFusions %>% filter(SameSV&PhaseMatched=='true') %>% filter(RegionTypeUp=='Exonic'&RegionTypeDown=='Intronic') %>%
       # filter(FusedExonUp==1&ExonsSkippedUp==0) %>%
       mutate(AdjExonSkippedUp=ifelse(RegionTypeUp=='Exonic'&RegionTypeDown=='Intronic'&ExonsSkippedUp==0,1,ExonsSkippedUp),
              AdjFusedExonUp=pmax(BreakendExonUp-AdjExonSkippedUp,1)) %>%
       select(SampleId,GenePair,RegionTypeUp,RegionTypeDown,CodingTypeUp,BreakendExonUp,FusedExonUp,
              ExonsSkippedUp,AdjExonSkippedUp,AdjFusedExonUp,CodingTypeDown,BreakendExonDown,FusedExonDown,ExonsSkippedDown,everything()))

View(unfilteredFusions %>% filter(SameSV&PhaseMatched=='true') %>% filter(RegionTypeUp=='Exonic'&RegionTypeDown=='Intronic') %>%
       filter(FusedExonUp==0&ExonsSkippedUp==0) %>%
       select(SampleId,GenePair,RegionTypeUp,RegionTypeDown,CodingTypeUp,CodingTypeDown,BreakendExonUp,BreakendExonDown,FusedExonUp,FusedExonDown,
              ExonsSkippedUp,ExonsSkippedDown,everything()))

## DEBUG


ensemblGeneData = read.csv('~/data/sv/ensembl_gene_data.csv')
View(ensemblData)

# 3' must have coding bases after the first exon
# 5' if coding must have coding bases after the first exon
ensemblTransData = read.csv('~/data/sv/ensembl_trans_exon_data.csv')
View(ensemblTransData)
View(head(ensemblTransData,100))
ensemblTransData = ensemblTransData %>% filter(TransId==CanonicalTranscriptId)
nrow(ensemblTransData)

transData = ensemblTransData %>% group_by(GeneId,TransId,Trans,Strand,CodingStart,CodingEnd,TransStart,TransEnd) %>%
  summarise(ExonCount=n(),
            LowExonStart=min(ExonStart),
            HighExonStart=max(ExonStart),
            LowExonEnd=min(ExonEnd),
            HighExonEnd=max(ExonEnd)) %>% ungroup()

transData = transData %>% mutate(CodingStart=ifelse(CodingStart=='NULL',-1,as.numeric(as.character(CodingStart))),
                                 CodingEnd=ifelse(CodingEnd=='NULL',-1,as.numeric(as.character(CodingEnd))),
                                 IsCoding=CodingStart>0&CodingEnd>0,
                                 CodingPastFirstExon=ifelse(Strand==1,IsCoding&ExonCount>1&CodingEnd>LowExonEnd,
                                                            IsCoding&ExonCount>1&CodingStart<HighExonStart))

View(transData)
View(transData %>% filter(is.na(CodingEnd)))
View(ensemblTransData %>% filter(Trans=='ENST00000442233'))
View(transData %>% group_by(GeneId) %>% count)

nrow(transData %>% filter(CodingPastFirstExon)) # 3' possibles
nrow(transData %>% filter(CodingPastFirstExon|(!IsCoding&ExonCount>1))) # 5' possibles
nrow(transData %>% filter(!CodingPastFirstExon&IsCoding))
nrow(transData %>% filter(ExonCount==1))

possibleTrans = transData %>% mutate(Possible3P=CodingPastFirstExon,
                                     Possible5P=CodingPastFirstExon|(!IsCoding&ExonCount>1))

possibleTrans = merge(possibleTrans,ensemblGeneData %>% select(GeneId,Chromosome),by='GeneId',all.x=T)
View(possibleTrans %>% filter(is.na(Chromosome)))
possibleTrans = possibleTrans %>% filter(!is.na(Chromosome))

View(possibleTrans %>% group_by(Chromosome) %>% summarise(Possible5P=sum(Possible5P),Possible3P=sum(Possible3P)))

rm(possibleTrans)
rm(transData)
rm(ensemblTransData)



tp53Introns = generate_intron_data(tp53ExonData)
View(tp53Introns)

phaseColours = c('Red','Yellow','Green','Grey')

tp53Introns = tp53Introns %>% mutate(PhaseColour=ifelse(IntronPhase==-1,'grey',
                                                 ifelse(IntronPhase==0,'yellow3',
                                                 ifelse(IntronPhase==1,'skyblue2','tomato2'))))

plot = (ggplot() 
        + geom_point(data=tp53Fusions, aes(x=PosDown,y=FusionIndex*15,colour='PosDown'))
        + geom_point(data=tp53Fusions, aes(x=PosUp,y=FusionIndex*15,colour='PosUp'))
        + geom_text(data=tp53Fusions, aes(x=PosUp+500,y=FusionIndex*15,label=paste(SampleId,CancerType,sep=' '),hjust=0, vjust=0))
        + geom_point(data=tp53SVs, aes(x=PosStart,y=-SvIndex/nrow(tp53SVs)*20,colour='PosStart'))
        + geom_point(data=tp53SVs, aes(x=PosEnd,y=-SvIndex/nrow(tp53SVs)*20,colour='PosEnd'))
        + geom_rect(data=tp53Introns, aes(xmin=IntronStart, xmax=IntronEnd, ymin=0, ymax=10), fill=tp53Introns$PhaseColour)
        + geom_text(data=tp53Introns, aes(x=(IntronStart+IntronEnd)/2,y=5,label=IntronRank))
        + labs(y='SampleId')
        # + geom_text(data=tp53ExonData, aes(x=(ExonStart+ExonEnd)/2,y=20,label=ExonRank))
        )

print(plot)



lowProbPairs = lowProbPairs %>% mutate(LocalProb=ifelse(ActFusions>ExpLocalFusions&ExpLocalFusions>0,1-ppois(ActFusions-1,ExpLocalFusions,T),1))






## INVESTIGATIONS

# same-genes
View(lowProbsByTypeLength %>% filter(as.character(GeneIdUp)==as.character(GeneIdDown)) %>%
       select(GenePair,Type,LengthMin,LengthMax,ActFusions,ActFusionsNF,ExpFusions,Prob,
              FivePrimeProm,ThreePrimeProm,NearDriver,FragileSite,everything()) %>%
       arrange(Prob))

# 2-gene DELs and DUPs
View(lowProbsByTypeLength %>% filter(as.character(GeneIdUp)!=as.character(GeneIdDown)) %>%
       select(GenePair,Type,LengthMin,LengthMax,ActFusions,ActFusionsNF,ExpFusions,Prob,
              FivePrimeProm,ThreePrimeProm,NearDriver,FragileSite,everything()) %>%
       arrange(Prob))

# other types
View(lowProbsByTypeLength %>% filter(Type!='DEL'&Type!='DUP') %>%
       select(GenePair,Type,LengthMin,LengthMax,ActFusions,ActFusionsNF,ExpFusions,Prob,
              FivePrimeProm,ThreePrimeProm,NearDriver,FragileSite,everything()) %>%
       arrange(Prob))

View(unfilteredFusions %>% filter(GenePair %in% lowProbsByTypeLength$GenePair) %>% 
       select(SampleId,CancerType,GenePair,Type,LengthMin,LengthMax,SameSV,CodingTypeUp,CodingTypeDown,FusedExonUp,FusedExonDown,KnownType,everything()))

View(unfilteredFusions)
View(unfilteredFusions %>% filter(GenePair %in% lowProbsByType$GenePair) %>% 
       group_by(SampleId,GenePair,FusedExonUp,FusedExonDown) %>% count() %>% group_by(GenePair,FusedExonUp,FusedExonDown) %>% 
       summarise(ActFusions=n()) %>% filter(ActFusions>1))

View(allFusions %>% filter(GenePair=='STAG1_MSL2'))


# factor in actual SVs which touch these 2 genes or are near by

proximatePairs = lowProbsByTypeLength %>% filter(Type=='DEL'|Type=='DUP')
nrow(proximatePairs)





# RNA Read Counts for K-mers

View(unfilteredFusions %>% filter(SameGene&PhaseUp==-1))



ahrExonReadCounts = read.csv('~/data/sv/fusion_like/ahr_exon_junction_counts.csv')
View(ahrExonReadCounts)

ctnnb1ExonReadCounts = read.csv('~/data/sv/fusion_like/ctnnb1_exon_junction_counts.csv')
View(ctnnb1ExonReadCounts)

ahrExonReadCounts = ahrExonReadCounts %>% mutate(ExonDesc=paste(ifelse(Aberrant=='true','Aberrant','Standard'),ExonJunction,sep=' '),
                                                 SampleType=ifelse(HasFusion=='true','Fusion','No Fusion')) %>%  arrange(Gene,ExonDesc)

print(ggplot(ahrExonReadCounts %>% filter(Gene=='CTNNB1'), aes(ExonDesc, ReadCount)) #  %>% mutate(ReadCount=pmin(ReadCount,50))
      + facet_wrap(~SampleType)
      + scale_y_log10()
      + theme(axis.text.x = element_text(angle = 90, hjust=1,size=10))
      + geom_boxplot(varwidth=T, fill="plum")
      + ggtitle(sprintf("SameGene Fusion Exon-Junction Read Counts: CTNNB1")))

print(ggplot(ctnnb1ExonReadCounts, aes(ExonJunction, ReadCount))
      + facet_wrap(~HasFusion)
      + scale_y_log10()
      + theme(axis.text.x = element_text(angle = 90, hjust=1,size=10))
      + geom_boxplot(varwidth=T, fill="plum"))

print(ggplot(ahrExonReadCounts %>% group_by(ExonJunction,HasFusion) %>% count, aes(ExonJunction, n))
      + facet_wrap(~HasFusion)
      + theme(axis.text.x = element_text(angle = 90, hjust=1,size=10))
      + geom_boxplot(varwidth=T, fill="plum"))

View(cohortSVs %>% filter(Type=='DEL'&LengthMin==2048000&LengthMax==4096000&grepl('ENSG00000157554',GeneStart)&grepl('ENSG00000184012',GeneEnd)))




View(unfilteredFusions %>% filter(GeneNameDown=='MYC'))
