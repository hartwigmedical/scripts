library(data.table)
library(dplyr)
library(tidyr)
library(stringi)


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
colnames(samplesWithRna) = c('SampleId')
View(samplesWithRna)


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


###############################
# RNA-FUSION MATCHING  ANALYSIS

load_rna_match_data<-function(filename)
{
  rnaMatchData = read.csv(filename)
  
  #Annotations
  rnaMatchData = rnaMatchData %>% mutate(SvMatchUp=!is.na(SvIdUp),
                                          SvMatchDown=!is.na(SvIdDown),
                                          SvViableUp=TransViableUp=='true',
                                          SvViableDown=TransViableDown=='true',
                                          SvValidLocUp=TransValidLocUp=='true'|SvViableUp,
                                          SvValidLocDown=TransValidLocDown=='true'|SvViableDown,
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
  rnaMatchData = rnaMatchData %>% group_by(FusionName) %>% filter(!((sum(ifelse(SvMatchType=='NoSV',1,0))>1&(sum(ifelse(SvMatchType=='BothSVs',1,0))<1)))) %>% ungroup()
  
  # rnaMatchData$FragmentSupport = (rnaMatchData$SpanningFragCount>0)
  return (rnaMatchData)
}

annotate_rna_both_svs<-function(rnaMatchData)
{
  rnaBothSVsData = rnaMatchData %>% filter(SvMatchType=='BothSVs')
  
  rnaBothSVsData = rnaBothSVsData %>% separate(ClusterInfoUp,c('ClusterIdUp','ClusterCountUp','ChainIdUp','ChainCountUp'),sep = ';')
  rnaBothSVsData = rnaBothSVsData %>% separate(ClusterInfoDown,c('ClusterIdDown','ClusterCountDown','ChainIdDown','ChainCountDown'),sep = ';')
  
  rnaBothSVsData$ClusterCountUp = as.numeric(rnaBothSVsData$ClusterCountUp)
  rnaBothSVsData$ClusterCountDown = as.numeric(rnaBothSVsData$ClusterCountDown)
  rnaBothSVsData$ChainCountUp = as.numeric(rnaBothSVsData$ChainCountUp)
  rnaBothSVsData$ChainCountDown = as.numeric(rnaBothSVsData$ChainCountDown)
  
  rnaBothSVsData$SameCluster = (rnaBothSVsData$ClusterIdUp==rnaBothSVsData$ClusterIdDown)
  
  rnaBothSVsData$IsChainedUp = (rnaBothSVsData$ChainCountUp>1)
  rnaBothSVsData$IsChainedDown = (rnaBothSVsData$ChainCountDown>1)
  rnaBothSVsData$SameChain = ifelse(rnaBothSVsData$SameSV,T,rnaBothSVsData$SameCluster&rnaBothSVsData$ChainIdUp==rnaBothSVsData$ChainIdDown)
  rnaBothSVsData = rnaBothSVsData %>% separate(ChainInfo,c('ChainLinks','ChainLength'),sep = ';')
  rnaBothSVsData$ChainLength = as.numeric(rnaBothSVsData$ChainLength)
  rnaBothSVsData$ChainLinks = as.numeric(rnaBothSVsData$ChainLinks)
  
  rnaBothSVsData$ChainLength = as.numeric(rnaBothSVsData$ChainLength)
  rnaBothSVsData$ChainLinks = as.numeric(rnaBothSVsData$ChainLinks)
  
  rnaBothSVsData$FacingInChain = (rnaBothSVsData$ChainLength>0)
  
  return (rnaBothSVsData)
}

View(allRnaData)
# write.csv(allRnaData, '~/data/sv/rna/rna_data_all_samples.csv', row.names = F, quote = F)


rnaMatchData = load_rna_match_data('~/data/sv/rna/SVA_RNA_DATA.csv')
View(rnaMatchData)
rnaMatchDataBothSVs = annotate_rna_both_svs(rnaMatchData)
View(rnaMatchDataBothSVs)
View(rnaMatchDataBothSVs %>% group_by(SameSV,SameCluster,SameChain) %>% count())
View(rnaMatchData %>% group_by(ViableFusion,TransViableUp,TransViableDown) %>% count())
View(rnaMatchData)
View(rnaMatchData %>% group_by(SvMatchUp,SvMatchDown,SvViableUp,SvViableDown,SpliceType,ViableFusion) %>% count())

# 1. Genuine Misses out of 1294
nrow(rnaMatchData %>% filter(SvMatchType!='BothSVs',Proximate==T))  #PROXIMATE   110 => 103
nrow(rnaMatchData %>% filter(SvMatchType=='NoSV',Proximate==F) %>% group_by(SampleId) %>% count())   #NO SV FOUND  675 => 87. # TO DO - Relax rules
nrow(rnaMatchData %>% filter(SvMatchType=='SingleSV',Proximate==F))  #One SV FOUND   154 => 161

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


# 4. Match with standard fusions

svaRnaFusions = read.csv('~/data/sv/rna/SVA_FUSIONS.csv')
nrow(svaRnaFusions)
nrow(svaFusions)
reportedSvaRnaFusions = svaRnaFusions %>% filter(Reportable=='true')
nrow(reportedSvaRnaFusions)

svaRnaFusions = annotate_fusions(svaRnaFusions)
View(svaRnaFusions)

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

# merge 2 data sets on SV Id
validSvaRnaFusionsById = validSvaRnaFusions %>% group_by(SvIdUp,SvIdDown) %>% summarise(SvaCount=n())
rnaMatchDataById = rnaMatchDataBothSVs %>% group_by(SvIdUp,SvIdDown) %>% count()
rnaMatchVsSvaById = merge(rnaMatchDataById,validSvaRnaFusionsById,by=c('SvIdUp','SvIdDown'),all.x=T)
View(rnaMatchVsSvaById)
View(rnaMatchVsSvaById %>% group_by(SvaMatched=!is.na(SvaCount)) %>% summarise(Matched=n(),Perc=round(n()/nrow(rnaMatchDataById),3)))
svaMatchVsRnaById = merge(rnaMatchDataById,validSvaRnaFusionsById,by=c('SvIdUp','SvIdDown'),all.y=T)
View(svaMatchVsRnaById)
View(svaMatchVsRnaById %>% group_by(RnaMatched=!is.na(RnaCount)) %>% summarise(Matched=n(),Perc=round(n()/nrow(svaMatchVsRnaById),3)))

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


## TIDY UP
rm(sampleRnaSvaFusions)
rm(validSvaRnaFusions)

