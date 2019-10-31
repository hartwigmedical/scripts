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
  
  localSvThreshold=1e6
  
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
      
      # count number of SVs in vacinity
      cohortSubset = cohortSubset %>% mutate(StartInGeneUp=ChrStart==chrUp&PosStart>geneStartUp-localSvThreshold&PosStart<geneEndUp+localSvThreshold,
                                           StartInGeneDown=ChrStart==chrDown&PosStart>geneStartDown-localSvThreshold&PosStart<geneEndDown+localSvThreshold,
                                           EndInGeneUp=ChrEnd==chrUp&PosEnd>geneStartUp-localSvThreshold&PosEnd<geneEndUp+localSvThreshold,
                                           EndInGeneDown=ChrEnd==chrDown&PosEnd>geneStartDown-localSvThreshold&PosEnd<geneEndDown+localSvThreshold)
      
      localCount = nrow(cohortSubset %>% filter((StartInGeneUp|StartInGeneDown)&(EndInGeneUp|EndInGeneDown)))
  
      gpSubset[i,which(colnames(gpSubset)=='LocalSvCount')] = localCount
      gpSubset[i,which(colnames(gpSubset)=='GenePairSvCount')] = genePairCount
      
      sameChromosome = (as.character(chrUp)==as.character(chrDown))
      
      distanceBetweenGenes = ifelse(sameChromosome,pmax(geneStartUp,geneStartDown)-pmin(geneEndUp,geneEndDown),0)
      
      localSvDistance = ifelse(!sameChromosome|distanceBetweenGenes>localSvThreshold*2,
                               (geneEndUp-geneStartUp)+(geneEndDown-geneStartDown)+localSvThreshold*4,
                               pmax(geneEndUp,geneEndDown)-pmin(geneStartUp,geneStartDown)+localSvThreshold*2)
      
      if(log)
      {
        print(sprintf("geneUp(%s %d-> %d) geneDown(%s %d-> %d) distance(%d) localSvDistance(%d)", 
                      geneNameUp, geneStartUp, geneEndUp, geneNameDown, geneStartDown, geneEndDown, distanceBetweenGenes, localSvDistance))
      }
      
      gpSubset = gpSubset %>% mutate(ExpGenePairCount=GenePairRate*fullGenome/localSvDistance*LocalSvCount,
                                             ExpLocalFusions=FusionRate*fullGenome/localSvDistance*LocalSvCount)
    }
    
    genePairResults = rbind(genePairResults,gpSubset)
  }

  return (genePairResults)
}

# tmp = annotate_proximate_svs(head(lowProbPairs,10),cohortSVs,F)
# View(tmp)


annotate_fusion_data<-function(genePairData,unfilteredFusions)
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
    type = as.character(gpData$Type)
    lengthMin = gpData$LengthMin
    lengthMax = gpData$LengthMax
    
    print(paste(i,': ',genePair,' - adding fusion annotations',sep=''))
    
    gpFusions = unfilteredFusions %>% filter(GenePair==genePair&Type==type&LengthMin==lengthMin&LengthMax==lengthMax)
    fusionCount = nrow(gpFusions)
    
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
  
  return (genePairData)
}

generate_intron_data<-function(tp53ExonData)
{
  intronData = data.frame(matrix(ncol=4, nrow=0))
  colnames(intronData) = c('IntronStart','IntronEnd','IntronPhase','IntronRank')
  intronCount = nrow(tp53ExonData)-1 
  
  for(i in 1:intronCount)
  {
    exonPrev=tp53ExonData[i,]
    exonNext=tp53ExonData[i+1,]
    
    intronData[i,1] = exonPrev$ExonEnd
    intronData[i,2] = exonNext$ExonStart
    intronData[i,3] = ifelse(exonPrev$Strand==1,exonPrev$ExonEndPhase,exonNext$ExonEndPhase)
    intronData[i,4] = exonPrev$ExonRank-1
  }
  
  return (intronData)
}

generate_fusion_intron_plot<-function(transId,svType,unfilteredFusions,ensemblTransExonData,cohortSVs)
{
  fusions = unfilteredFusions %>% filter(TranscriptUp==transId&TranscriptDown==transId&Type==svType) %>% 
    mutate(SvLength=PosUp-PosDown) %>% arrange(SvLength)
  
  if(nrow(fusions) == 0)
  {
    print("no fusions found")
    return(fusions)
  }
  
  print(paste("fusions found: ",nrow(fusions),sep=''))
  
  rowIndex = data.frame(as.numeric(as.character(rownames(fusions))))
  colnames(rowIndex) <- c("FusionIndex")
  fusions = cbind(fusions,rowIndex)
  
  geneId = fusions[1,which(colnames(fusions)=='GeneIdUp')]
  geneName = fusions[1,which(colnames(fusions)=='GeneNameUp')]
  # print(geneId)
  
  exonData = ensemblTransExonData %>% filter(Trans==transId)
  intronData = generate_intron_data(exonData)
  
  intronData = intronData %>% mutate(PhaseColour=ifelse(IntronPhase==-1,'grey',
                                                        ifelse(IntronPhase==0,'yellow3',
                                                               ifelse(IntronPhase==1,'skyblue2','tomato2'))))
  
  transStart = min(exonData$ExonStart)
  transEnd = max(exonData$ExonEnd)
  
  geneSVs = cohortSVs %>% filter(Type==svType&grepl(geneId,GeneStart)&grepl(geneId,GeneEnd)&!(SampleId %in% fusions$SampleId))
  geneSVs = geneSVs %>% filter(PosStart>=transStart&PosEnd<=transEnd) %>% arrange(Length)
  
  print(paste('SVs in gene-pair=', nrow(geneSVs),sep=''))
  
  rowIndex = data.frame(as.numeric(as.character(rownames(geneSVs))))
  colnames(rowIndex) <- c("SvIndex")
  geneSVs = cbind(geneSVs,rowIndex)
  
  phaseColours = c('Red','Yellow','Green','Grey')
  
  plot = (ggplot() 
          + geom_point(data=fusions, aes(x=PosDown,y=FusionIndex*15,colour='PosDown'))
          + geom_point(data=fusions, aes(x=PosUp,y=FusionIndex*15,colour='PosUp'))
          + geom_text(data=fusions, aes(x=pmax(PosDown,PosUp)+500,y=FusionIndex*15,label=paste(SampleId,CancerType,sep=' '),hjust=0, vjust=0))
          + geom_rect(data=intronData, aes(xmin=IntronStart, xmax=IntronEnd, ymin=0, ymax=10), fill=intronData$PhaseColour)
          + geom_text(data=intronData, aes(x=(IntronStart+IntronEnd)/2,y=5,label=IntronRank))
          + labs(y='SampleId')
          + ggtitle(sprintf("Same Gene Fusion: %s-%s %s", geneId, geneName, svType)))
  
  if(nrow(geneSVs) > 0)
  {
    plot = (plot 
            + geom_point(data=geneSVs, aes(x=PosStart,y=-SvIndex/nrow(geneSVs)*20,colour='PosStart'))
            + geom_point(data=geneSVs, aes(x=PosEnd,y=-SvIndex/nrow(geneSVs)*20,colour='PosEnd')))
  }
  
  return (plot)
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

unfilteredFusions = merge(unfilteredFusions,sampleCancerTypes,by='SampleId',all.x=T)
View(unfilteredFusions %>% group_by(CancerType) %>% count)
View(unfilteredFusions %>% filter(PhaseMatched=='true') %>% group_by(KnownType,ExonsSkipped=(ExonsSkippedUp>0|ExonsSkippedDown>0)|ExceedsDownPreGeneLimit) %>% count)
View(unfilteredFusions %>% filter(Filtered) %>% group_by(KnownType,Filtered) %>% count)
View(unfilteredFusions %>% filter(!Filtered) %>% group_by(SameSV) %>% count)
View(unfilteredFusions)
# View(unfilteredFusions %>% filter(KnownType==''&(ExonsSkippedUp>0|ExonsSkippedDown>0)))

# 1b
# Produce a set of gene-pairs to calculate likelihood for
# - known-fusion pairs
# - any fusion in our cohort occurring 2+ times

cohortGenePairSet = unfilteredFusions %>% group_by(SampleId,GenePair,GeneIdUp,GeneIdDown,KnownType) %>% count() %>% 
  group_by(GenePair,GeneIdUp,GeneIdDown,KnownType) %>% summarise(SampleCount=n()) %>% filter(KnownType=='Known'|SampleCount>=2)
View(cohortGenePairSet)
nrow(cohortGenePairSet)

# manually add in the 6 cosmic blacklisted  fusions
blacklistedCosmic = read.csv('~/data/sv/fusion_like/blacklisted_cosmic_fusions.csv')
blacklistedCosmic = merge(blacklistedCosmic,ensemblGeneData %>% select(GeneIdUp=GeneId,GeneName),by.x='GeneNameUp',by.y='GeneName',all.x=T)
blacklistedCosmic = merge(blacklistedCosmic,ensemblGeneData %>% select(GeneIdDown=GeneId,GeneName),by.x='GeneNameDown',by.y='GeneName',all.x=T)
blacklistedCosmic = blacklistedCosmic %>% mutate(GenePair=paste(GeneNameUp,GeneNameDown,sep='_'),
                                                 KnownType='',
                                                 SampleCount=0)

blacklistedCosmic = blacklistedCosmic %>% filter(!(GenePair %in% cohortGenePairSet$GenePair))
View(blacklistedCosmic)



genePairSet = rbind(cohortGenePairSet %>% ungroup(),blacklistedCosmic %>% select(GenePair,GeneIdUp,GeneIdDown,KnownType,SampleCount))
View(genePairSet)
write.csv(genePairSet %>% group_by(GeneIdUp,GeneIdDown) %>% count() %>% select(GeneIdUp,GeneIdDown), # check for duplicates
          '~/data/sv/fusion_like/cohort_gene_pairs.csv',row.names = F, quote = F)

View(genePairSet %>% group_by(GeneIdUp,GeneIdDown) %>% count())


# known pairs - for now only including known pairs which feature at least once in our cohort
knownPairs = read.csv('~/data/knownFusionPairs.csv')

View(knownPairs)
knownFusionPairs = merge(knownPairs,ensemblGeneData %>% select(GeneIdUp=GeneId,GeneName),by.x='FiveGene',by.y='GeneName',all.x=T)
knownFusionPairs = merge(knownFusionPairs,ensemblGeneData %>% select(GeneIdDown=GeneId,GeneName),by.x='ThreeGene',by.y='GeneName',all.x=T)
knownFusionPairs = knownFusionPairs %>% filter(!(is.na(GeneIdUp)|is.na(GeneIdDown)))
knownFusionPairs = knownFusionPairs %>% mutate(GenePair=paste(FiveGene,ThreeGene,sep='_'))
View(knownFusionPairs)
write.csv(knownFusionPairs %>% select(GeneIdUp,GeneIdDown),'~/data/sv/fusion_like/knownFusionPairs.csv',row.names = F, quote = F)

knownFusionPairs = knownFusionPairs %>% select(GenePair,GeneIdUp,GeneIdDown) %>% mutate(SampleCount=1)
nrow(knownFusionPairs)
View(knownFusionPairs %>% filter(!(GenePair %in% cohortGenePairSet$GenePair)))

genePairSet = rbind(cohortGenePairSet %>% ungroup(),knownFusionPairs %>% filter(!(GenePair %in% cohortGenePairSet$GenePair)))
nrow(genePairSet)
View(genePairSet)

write.csv(genePairSet %>% select(GeneIdUp,GeneIdDown),'~/data/sv/fusion_like/cohort_known_gene_pairs.csv',row.names = F, quote = F)


# 2. Load FLC calc data - expected fusion and gene-pair SV rates
flcSpecificPairs = read.csv('~/data/sv/fusion_like/GFL_FUSION_PAIR_LIKELIHOOD.csv')
# flcSpecificPairs = read.csv('~/logs/GFL_FUSION_PAIR_LIKELIHOOD.csv')
View(flcSpecificPairs)
nrow(flcSpecificPairs)

flcSpecificPairs = flcSpecificPairs %>% mutate(GenePair=paste(GeneNameUp,GeneNameDown,sep='_'),
                                               FusionVsGenePairRateRatio=ifelse(GenePairRate>0,round(FusionRate/GenePairRate,4),0))

# for now limit to those genes-pairs actually observed in our cohort
# flcSpecificPairs = flcSpecificPairs %>% filter(GenePair %in% gpActFusions$GenePair)

gpActFusions = unfilteredFusions %>% filter(GenePair %in% flcSpecificPairs$GenePair) %>% 
  filter(!Filtered) %>% 
  group_by(SampleId,GenePair,Type,LengthMin,LengthMax) %>% count() %>% group_by(GenePair,Type,LengthMin,LengthMax) %>% summarise(ActFusions=n())

gpActFusionsUnfiltered = unfilteredFusions %>% filter(GenePair %in% flcSpecificPairs$GenePair) %>% 
  group_by(SampleId,GenePair,Type,LengthMin,LengthMax) %>% count() %>% group_by(GenePair,Type,LengthMin,LengthMax) %>% summarise(ActFusionsNF=n())

gpActFusions = merge(gpActFusions,gpActFusionsUnfiltered,by=c('GenePair','Type','LengthMin','LengthMax'),all.x=T)

View(gpActFusions)
View(gpActFusions %>% filter(ActFusions>1))


flcSpecificPairs = merge(flcSpecificPairs,gpActFusions,by=c('GenePair','Type','LengthMin','LengthMax'),all.x=T)

# merge in actual observed counts by type
flcSpecificPairs = merge(flcSpecificPairs,cohortSvSummary %>% select(-CountAll),by=c('Type','LengthMin','LengthMax'),all.x=T)
flcSpecificPairs[is.na(flcSpecificPairs)] = 0

# View(flcSpecificPairs)

# limit to those observed 1+ times
flcSpecificPairs = flcSpecificPairs %>% filter(ActFusionsNF>=1)
nrow(flcSpecificPairs) # 10121

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

#  annotate fusion with known gene and fusion data
flcSpecificPairs = flcSpecificPairs %>% mutate(
  KnownPair=GenePair %in% knownFusionPairs$GenePair, 
  RevKnownPair=ReverseGenePair %in% knownFusionPairs$GenePair, 
  ThreePrimeProm=GeneNameDown %in% threePrimeProm$GeneName,
  FivePrimeProm=GeneNameUp %in% fivePrimeProm$GeneName,
  FragileSite=(GeneNameUp %in% fragileSiteGenes$GeneName)|(GeneNameDown %in% fragileSiteGenes$GeneName),
  KnownTSG=(GeneNameUp %in% topTsgGenes$GeneName)|(GeneNameDown %in% topTsgGenes$GeneName))

View(flcSpecificPairs)

# lowProbPairs = flcSpecificPairs %>% filter(ActFusionsNF>=2&Prob<0.001)

# for now remove the probability filter
lowProbPairs = flcSpecificPairs %>% filter(ActFusionsNF>=2|KnownPair)
nrow(lowProbPairs) # 4243
View(lowProbPairs)

# to lower the set being evaluated, ignore any unknown gene-pairs with high probability
# View(lowProbPairs %>% filter(!KnownPair&!ThreePrimeProm&!FivePrimeProm) %>% filter(Prob>0.01) %>% arrange(Prob))
lowProbPairs = lowProbPairs %>% filter(KnownPair|ThreePrimeProm|FivePrimeProm|FragileSite|Prob<0.01|(GenePair %in% blacklistedCosmic$GenePair))
nrow(lowProbPairs)
View(lowProbPairs)

# check driver catalog
drivers = read.csv('~/data/sv/drivers/LNX_DRIVERS.csv')
ampDrivers = drivers %>% filter((DriverType=='AMP'|DriverType=='DEL')&ClusterId>=0)
ampDrivers = merge(ampDrivers,ensemblGeneData %>% select(Gene=GeneName,GeneStart,GeneEnd),by='Gene',all.x=T)
ampDriverPositions = ampDrivers %>% group_by(Gene,Chromosome,GeneStart,GeneEnd) %>% summarise(DriverCount=n())
View(ampDriverPositions)

lowProbPairs = annotate_drivers(lowProbPairs,ampDriverPositions)

# add in SV counts around and inside the genes
lowProbPairs = annotate_proximate_svs(lowProbPairs,cohortSVs)
nrow(lowProbPairs)

# add in any repeated exon-fusings or dominant cancer types

# - any recurrent fused exons or contexts
# - any recurrent cancer types
# - any repeated cluster types
lowProbPairs = annotate_fusion_data(lowProbPairs,unfilteredFusions)

lowProbPairs = lowProbPairs %>% mutate(LocalProb=1-ppois(ActFusions-1,ExpLocalFusions,T),
                                       GenePairProb=1-ppois(GenePairSvCount-1,ExpGenePairCount,T))



View(lowProbPairs)
# colnames(lowProbPairs)
write.csv(lowProbPairs,'~/data/sv/fusion_like/gene_pair_summary.csv', quote = F, row.names = F)

View(lowProbPairs %>% select(GenePair,Type,LengthMin,LengthMax,Prob,LocalProb,ActFusions,ExpFusions,ExpLocalFusions,CancerTypes,FusedExons,FusionClusters,everything())
     %>% arrange(LocalProb))

View(lowProbPairs %>% filter(!KnownPair&!RevKnownPair&!FragileSite&!NearDriver) %>% arrange(LocalProb) %>%
       mutate(Length=paste(LengthMin/1000,'-',LengthMax/1000,'K',sep='')) %>%
       select(GenePair,Type,Length,Prob,LocalProb,ActFusions,CancerTypes,FusedExons,FusionClusters,everything()))

View(lowProbPairs %>% filter(!KnownPair&!RevKnownPair&!FragileSite&!NearDriver) %>% arrange(LocalProb) %>%
       mutate(Length=paste(LengthMin/1000,'-',LengthMax/1000,'K',sep='')) %>%
       select(GenePair,Type,Length,Prob,LocalProb,ActFusions,CancerTypes,FusedExons,FusionClusters,everything()))

write.csv(lowProbPairs %>% filter(!KnownPair&!RevKnownPair&!FragileSite&!NearDriver&!RepeatedPositions) %>% arrange(LocalProb) %>%
       mutate(Length=paste(LengthMin/1000,'-',LengthMax/1000,'K',sep='')) %>%
       select(GenePair,Type,Length,Prob,LocalProb,ActFusions,CancerTypes,FusedExons,FusionClusters,FusionContexts),
       '~/data/sv/fusion_like/gene_pair_google_sheets.csv', quote = F, row.names = F)


# Summary across all buckets

genePairSummary = lowProbPairs %>% group_by(GenePair) %>% 
  summarise(ActFusions=sum(ActFusions),
            ActFusionsNF=sum(ActFusionsNF),
            ExpFusions=sum(FusionRate*CohortCount),
            ExpLocalFusions=sum(FusionRate*LocalSvCount),
            GenePairSvCount=sum(GenePairSvCount),
            KnownPair=first(KnownPair|RevKnownPair),
            ThreePrimeProm=first(ThreePrimeProm),
            FivePrimeProm=first(FivePrimeProm),
            FragileSite=first(FragileSite),
            NearDriver=first(NearDriver),
            DriverGenes=first(DriverGenes),
            KnownTSG=first(KnownTSG)) %>%
  mutate(Prob=1-ppois(ActFusions-1,ExpFusions,T),
         LocalProb=1-ppois(ActFusions-1,ExpLocalFusions,T))

View(genePairSummary %>% filter(ActFusionsNF>=2|(KnownPair&ActFusionsNF>=1)|((KnownPair|ThreePrimeProm|FivePrimeProm)&ActFusionsNF>=1)) %>%
       arrange(LocalProb))

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

View(unfilteredFusions %>% filter(GenePair=='TTN_TTN'&Type=='DEL'&LengthMin==100&LengthMax==1000) %>%
       select(BreakendExonUp,BreakendExonDown,RegionTypeUp,RegionTypeDown,Length,PosUp,PosDown,ProteinsKept,ProteinsLost,everything()))


nrow(tmpSVs %>% filter(PosStart>179390716-1e6&PosEnd<179695529+1e6))
geneLength=179695529-179390716
print(geneLength/(geneLength+2e6))




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

print(generate_fusion_intron('ENST00000269305','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # TP53
print(generate_fusion_intron('ENST00000242057','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # AHR
print(generate_fusion_intron('ENST00000349496','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # CTNNB1
print(generate_fusion_intron('ENST00000426105','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # PUM1
print(generate_fusion_intron('ENST00000269701','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # AKAP8

print(generate_fusion_intron('ENST00000288602','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # BRAF
print(generate_fusion_intron('ENST00000288602','DUP',unfilteredFusions,ensemblTransExonData,cohortSVs)) # BRAF

print(generate_fusion_intron('ENST00000342788','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # ERBB4_ERBB4
print(generate_fusion_intron('ENST00000342788','DUP',unfilteredFusions,ensemblTransExonData,cohortSVs)) # ERBB4_ERBB4

print(generate_fusion_intron('ENST00000262189','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # KMT2C
print(generate_fusion_intron('ENST00000378933','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # EP300
print(generate_fusion_intron('ENST00000275493','DEL',unfilteredFusions,ensemblTransExonData,cohortSVs)) # EGFR
View(unfilteredFusions)

sameGeneFusions = lowProbPairs %>% filter(as.character(GeneIdUp)==as.character(GeneIdDown)) %>% 
  group_by(GeneId=GeneIdUp,Type,GeneNameUp) %>% summarise(LowestLocalProb=min(LocalProb)) %>% arrange(LowestLocalProb)

nrow(sameGeneFusions)
View(sameGeneFusions)

topN = 100
outputFile = '~/data/sv/fusion_like/same_gene_plots_top100.pdf'
pdf(file=outputFile, height = 14, width = 20)
par(mar=c(1,1,1,1))

for(i in 1:pmin(nrow(sameGeneFusions),topN))
{
  gpData = sameGeneFusions[i,]
  geneId = gpData$GeneId
  type= as.character(gpData$Type)
  
  gpFusions = head(unfilteredFusions %>% filter(as.character(GeneIdUp)==geneId&as.character(GeneIdDown)==geneId&CanonicalUp=='true'),1)
  transId = as.character(gpFusions[1,which(colnames(gpFusions)=='TranscriptUp')])

  print(sprintf("plotting %d: gene(%s) trans(%s)", i, geneId, transId))
  
  print(generate_fusion_intron(transId,type,unfilteredFusions,ensemblTransExonData,cohortSVs))
}

dev.off()



tp53SVs = cohortSVs %>% filter(Type=='DEL'&grepl('TP53',GeneStart)&grepl('TP53',GeneEnd))
tp53SVs = cohortSVs %>% filter(Type=='DEL'&grepl('TP53',GeneStart)&grepl('TP53',GeneEnd)&!(SampleId %in% tp53Fusions$SampleId))
tp53SVs = tp53SVs %>% filter(PosStart>=7571720&PosEnd<=7590856) %>% arrange(Length)
View(tp53SVs)

rowIndex = data.frame(as.numeric(as.character(rownames(tp53SVs))))
colnames(rowIndex) <- c("SvIndex")
tp53SVs = cbind(tp53SVs,rowIndex)

tp53ExonData = ensemblTransExonData %>% filter(Trans=='ENST00000269305')
View(tp53ExonData)

tp53Fusions = unfilteredFusions %>% filter(GenePair=='TP53_TP53'&Type=='DEL') %>% 
  select(SampleId,ChrUp,PosUp,PosDown,Length,RegionTypeUp,RegionTypeDown,BreakendExonUp,BreakendExonDown,TranscriptUp,CancerType) %>%
  mutate(SvLength=PosUp-PosDown) %>% arrange(SvLength)

rowIndex = data.frame(as.numeric(as.character(rownames(tp53Fusions))))
colnames(rowIndex) <- c("FusionIndex")
tp53Fusions = cbind(tp53Fusions,rowIndex)

View(tp53Fusions)

View(unfilteredFusions %>% filter(as.character(GeneIdUp)==as.character(GeneIdDown)&Type=='DUP'))
View(unfilteredFusions %>% filter(as.character(GeneIdUp)==as.character(GeneIdDown)&Type=='DUP') %>% group_by(GenePair) %>% count()
     %>% filter(n>1) %>% arrange(-n))




## DEBUG


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



# manual checks
matchingSVs = cohortSVs %>% filter(Type=='DEL'&LengthMin==2000&LengthMax==4000)

localCount = matchingSVs %>% filter(ChrStart==7) %>%
  filter(PosStart>pmin(17338246)-localSvThreshold&PosEnd<pmax(17385776)+localSvThreshold)
nrow(localCount)     

genePairCount = matchingSVs %>% filter((grepl('ENSG00000106546',GeneStart)&grepl('ENSG00000106546',GeneEnd))
                                       |(grepl('ENSG00000106546',GeneStart)&grepl('ENSG00000106546',GeneEnd)))
nrow(genePairCount)
View(genePairCount)



# RNA Read Counts for K-mers

View(unfilteredFusions %>% filter(SameGene&PhaseUp==-1))



ahrExonReadCounts = read.csv('~/data/sv/fusion_like/ahr_exon_junction_counts.csv')
View(ahrExonReadCounts)

ctnnb1ExonReadCounts = read.csv('~/data/sv/fusion_like/ctnnb1_exon_junction_counts.csv')
View(ctnnb1ExonReadCounts)

print(ggplot(ahrExonReadCounts, aes(ExonJunction, ReadCount)) #  %>% mutate(ReadCount=pmin(ReadCount,50))
      + facet_wrap(~HasFusion)
      + scale_y_log10()
      + theme(axis.text.x = element_text(angle = 90, hjust=1,size=10))
      + geom_boxplot(varwidth=T, fill="plum"))

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

genePairData = merge(flcSpecificPairs,svFusionBreakends,by=c('GeneIdUp','GeneIdDown','Type','LengthMin','LengthMax'),all.x=T)

genePairData = genePairData %>% 
  mutate(ExpFusions=FusionRate*CountAll,
         #ExpOffFusions=ifelse(ExpBreakendFusions>ExpFusion,ExpBreakendFusions-ExpFusion,ExpFusion),
         #AdjExpFusions=ifelse(ExpOffFusions>0,round((ActBeCount-ActFusions)*ExpFusion/ExpOffFusions,4),ExpFusion),
         Prob=ifelse(ActFusions>ExpFusions&ExpFusions>0,1-ppois(ActFusions-1,ExpFusions,T),1))

View(genePairData)


## DEBUG

genePairData = lowProbPairs
View(genePairData)

genePairData = within(genePairData,rm(LocalSvCount))
genePairData = within(genePairData,rm(LocalSvCount))
genePairData = within(genePairData,rm(GenePairSvCount))

genePairData$LocalSvCount = 0
genePairData$GenePairSvCount = 0
localSvThreshold=1e6
colIndex = ncol(genePairData)

gpData = head(genePairData %>% filter(GeneNameUp=='AHR'&GeneNameDown=='AHR'),1)
nrow(gpData)
genePair = gpData$GenePair
chrUp = gpData$ChrUp
chrDown = gpData$ChrDown
geneIdUp = gpData$GeneIdUp
geneIdDown = gpData$GeneIdDown
geneStartUp = gpData$GeneStartUp
geneEndUp = gpData$GeneEndUp
geneStartDown = gpData$GeneStartDown
geneEndDown = gpData$GeneEndDown
type = as.character(gpData$Type)
lengthMin = gpData$LengthMin
lengthMax = gpData$LengthMax


# count number of SVs hitting both these genes
matchingSVs = cohortSVs %>% filter(Type==type&LengthMin==lengthMin&LengthMax==lengthMax)

genePairCount = matchingSVs %>% filter((grepl(geneIdUp,GeneStart)&grepl(geneIdDown,GeneEnd))|(grepl(geneIdDown,GeneStart)&grepl(geneIdUp,GeneEnd)))

# count number of SVs in vacinity
localCount = matchingSVs %>% filter(as.character(ChrStart)==as.character(chrUp)) %>%
  filter((PosStart>geneStartUp-localSvThreshold&PosEnd<geneEndUp+localSvThreshold)
         |(PosStart>geneStartDown-localSvThreshold&PosEnd<geneEndDown+localSvThreshold))

View(localCount)

gpData$LocalSvCount = nrow(localCount)
gpData$GenePairSvCount = nrow(genePairCount)

distanceBetweenGenes = pmax(geneStartUp,geneStartDown)-pmin(geneEndUp,geneEndDown)

localSvDistance = ifelse(distanceBetweenGenes>localSvThreshold*2,
                         (geneEndUp-geneStartUp)+(geneEndDown-geneStartDown)+localSvThreshold*4,
                         pmax(geneEndUp,geneEndDown)-pmin(geneStartUp,geneStartDown)+localSvThreshold*2)

gpData = gpData %>% mutate(ExpGenePairCount=round(GenePairRate*fullGenome/localSvDistance*LocalSvCount,2),
                                       ExpLocalFusions=round(FusionRate*fullGenome/localSvDistance*LocalSvCount,2))
  
View(gpData)





