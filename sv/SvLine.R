

# load Linx data
svData = read.csv('~/data/sv/cohort/LNX_SVS.csv')
svData = svData %>% filter(SampleId %in% samplesDD$SampleId)
svData = svData %>% filter(!(ResolvedType %in% c('LOW_VAF','DUP_BE'))) # remove artefect types



# just line clusters and elements

lineSvData = svData %>% filter(ResolvedType=='LINE')
lineSvData = lineSvData %>% mutate(IsSource=LEStart!='NONE'|LEEnd!='NONE')


#####
## LINE cohort data


# Write input data for SvTools cohort analysis
lineElementClusters = lineSvData %>% group_by(SampleId,ClusterId,ClusterCount) %>% summarise(SourceCount=sum(IsSource))
lineSvData = merge(lineSvData,lineElementClusters %>% select(-ClusterCount),by=c('SampleId','ClusterId'),all.x=T)
View(lineSvData %>% filter(SourceCount>=1) %>% select(SourceCount,everything()))
View(lineSvData %>% filter(SourceCount==0) %>% select(SourceCount,everything()))
View(lineElementClusters)
View(lineElementClusters %>% group_by(HasSource=SourceCount>0) %>% count)

write.csv(lineSvData %>% filter(SourceCount>=1) %>% arrange(SampleId,ClusterId,-IsSource) %>% select(SampleId,ClusterId,ChrStart,ChrEnd,PosStart,PosEnd,LEStart,LEEnd),
          '~/data/sv/cohort/sv_line_data.csv',row.names = F,quote = F)


## common locations

# load cohort data
lineLocations = read.csv('~/data/sv/cohort/LINE_SOURCE_DATA.csv')
# lineLocations = lineLocations %>% mutate(Length=PosEnd-PosStart)
View(lineLocations)
View(lineLocations %>% filter(Length>0))
View(lineLocations %>% filter(Type=='SUSPECT'))
View(lineLocations %>% filter(SampleSourceLocations==1))
View(lineLocations %>% filter(RmId>0&KnownPosStart==-1))
View(lineLocations %>% filter(RmPosStart<&KnownPosStart==-1))
View(lineLocations %>% filter((KnownPosStart>0&KnownPosEnd<0)|(RmPosStart>0&RmPosEnd<0)|(PolymorphPosStart>0&PolymorphPosEnd<0))) # always a pair

# non-known, repeat-masker locations
View(lineLocations %>% filter(RmId>0&KnownPosStart==-1) %>% group_by(LineId,RmId,Chromosome,RmPosStart,RmPosEnd,RmStrand,PosStart,PosEnd,SampleCount,PcawgSampleCount) %>% 
       summarise(Distance=first(ifelse(RmStrand==1,abs((RmPosStart+RmPosEnd)/2-(PosStart+PosEnd)/2),abs((PosStart+PosEnd)/2-(RmPosStart+RmPosEnd)/2)))))

# location summary data
locationSummary = lineLocations %>% group_by(LineId,Type,Chromosome,PosStart,PosEnd,SampleCount,PcawgSampleCount,
                                             KnownPosStart,KnownPosEnd,RmId,RmStrand,RmPosStart,RmPosEnd,PolymorphPosStart,PolymorphPosEnd) %>% 
  summarise(MaxSampleInserts=max(SampleInserts),
            TotalInserts=first(TotalInserts))

View(locationSummary %>% mutate(HasKnown=KnownPosStart>0,
                                HasRepeatMasker=RmPosStart>0,
                                HasPolymorphic=PolymorphPosStart>0))

write.csv(locationSummary %>% mutate(HasKnown=KnownPosStart>0,HasRepeatMasker=RmPosStart>0,HasPolymorphic=PolymorphPosStart>0),'~/data/sv/cohort/line_locations.csv',row.names = F,quote=F)


# multiple locations per cluster

View(lineLocations)

View(lineSvData %>% filter(SampleId=='CPCT02030318T'&ClusterId==906))
View(lineSvData %>% filter(SampleId=='CPCT02040242T'&ClusterId==20))

View(lineSvData %>% filter(SampleId=='CPCT02020353T'&ClusterId==46) %>% select(SampleId,ClusterId,ClusterCount,ClusterDesc,Id,Type,ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd,
                                                                               LEStart,LEEnd,DBLenStart,DBLenEnd,InsertSeq,InsSeqAlignments,everything()))






## Repeat Masker Line Annotations
rmLineData = read.csv('~/data/sv/repeat_masker_line.csv',sep='\t')
View(rmLineData)
write.csv(rmLineData %>% select(RmId,Chromosome,PosStart,PosEnd,Strand),'~/data/sv/line_repeat_mask_data.csv', row.names = F,quote = F)


# Chaining
tmpLineChains = read.csv('~/logs/LNX_LINE_CHAINS.csv')
View(tmpLineChains)

lineChains = read.csv('~/data/sv/cohort/LNX_LINE_CHAINS.csv')
View(lineChains)
View(lineChains %>% group_by(ChainDesc,ChainComplete) %>% count %>% spread(ChainComplete,n,fill=0))
View(lineChains %>% filter(ChainDesc=='BND=2'&ChainComplete=='false'))

View(lineSvData %>% filter(SampleId=='CPCT02040156T'&ClusterId==105))


# Line activity around Intact Known or Suspect Sites
lineLocSubset = lineLocations %>% filter(SampleCount>=5) %>% mutate(HasKnown=KnownPosStart>0,
                                                                    HasRepeatMasker=RmPosStart>0,
                                                                    HasPolymorphic=PolymorphPosStart>0) %>%
  filter(HasKnown|HasRepeatMasker|HasPolymorphic) %>% 
  mutate(RefPosStart=ifelse(HasKnown,KnownPosStart,ifelse(HasPolymorphic,PolymorphPosStart,RmPosStart)),
         RefPosEnd=ifelse(HasKnown,KnownPosEnd,ifelse(HasPolymorphic,PolymorphPosEnd,RmPosEnd)),
         LineLength=RefPosEnd-RefPosStart)

View(lineLocSubset)

# determine strand for the polymorphic line elements
lineLocSubset = lineLocSubset %>% mutate(SrcPreDistance=RefPosStart-PosStart,
                                         SrcPostDistance=PosEnd-RefPosEnd,
                                         RmStrand=ifelse(RmStrand==0,ifelse(SrcPreDistance>SrcPostDistance,-1,1),RmStrand))

# View(lineLocSubset %>% filter(RmStrand==0) %>% select(PolyStrand,RefPosStart,RefPosEnd,PosStart,PosEnd,SrcPreDistance,SrcPostDistance,everything()))

lineLocChains = merge(lineChains,lineLocSubset %>% select(LineId,RmId,SampleId,ClusterId,Chromosome,RefPosStart,RefPosEnd,RmStrand,SampleCount,PcawgSampleCount),
                         by=c('SampleId','ClusterId'),all.y=T)

lineLocChains = lineLocChains %>% filter(!is.na(ChainId))

lineChainTypes = c('BND=2','BND=1_SGL=1','BND=2_INV=1','BND=1_INV=1','SGL','BND')
View(lineLocChains %>% group_by(ChainDesc,ChainSvCount) %>% count)
lineLocChains = lineLocChains %>% filter(ChainDesc %in% c('BND=2','BND=1_SGL=1','BND=2_INV=1','BND=1_INV=1','SGL','BND'))

View(lineLocChains)
nrow(lineLocChains)
View(lineLocChains %>% filter(LineId==25))
View(lineLocChains %>% filter(LineId==25) %>% filter(SourceChr=='X'&SourcePosStart>RefPosStart-5e3&SourcePosStart<RefPosEnd+5e3))
colnames(lineLocChains)


lineLocChains2 = lineLocChains %>% filter(as.character(Chromosome)==as.character(SourceChr)&SourcePosStart>RefPosStart-1e4&SourcePosStart<RefPosEnd+1e4&SourcePosEnd<RefPosEnd+1e4)
nrow(lineLocChains2)

lineLocChains2 = lineLocChains2 %>% mutate(HasBothSource=(SourcePosStart>0&SourcePosEnd>0),
                                           IsReverseSrcStart=(HasBothSource&(RmStrand==SourceOrientStart)),
                                           RelSrcPosStartRaw=ifelse(RmStrand==1,SourcePosStart-RefPosEnd,RefPosStart-SourcePosStart),
                                           RelSrcPosStart=ifelse(HasBothSource,ifelse(!IsReverseSrcStart,RelSrcPosStartRaw,NA),NA),
                                           RelSrcPosReverseStart=ifelse(HasBothSource,ifelse(IsReverseSrcStart,RelSrcPosStartRaw,NA),NA),
                                           RelSrcSolo=ifelse(!HasBothSource,RelSrcPosStartRaw,NA),
                                           RelSrcPosEnd=ifelse(SourcePosEnd>0,ifelse(RmStrand==1,SourcePosEnd-RefPosEnd,RefPosStart-SourcePosEnd),NA),
                                           RelInvPosStart=ifelse(SourceInvPosStart>0,ifelse(RmStrand==1,SourceInvPosStart-RefPosEnd,RefPosStart-SourceInvPosStart),NA),
                                           RelInvPosEnd=ifelse(SourceInvPosEnd>0,ifelse(RmStrand==1,SourceInvPosEnd-RefPosEnd,RefPosStart-SourceInvPosEnd),NA))

# validation
View(lineLocChains2 %>% filter(LineId==25) %>% select(ChainDesc,ChainSvCount,RmStrand,SourceOrientStart,RelSrcPosStart,RelSrcPosReverseStart,RelSrcSolo,
                                                      HasBothSource,IsReverseSrcStart,SourcePosStart,RefPosStart,RelSrcPosEnd,SourcePosEnd,RefPosEnd,
                                                      SourceInvPosStart,SourceInvPosEnd,RelInvPosStart,RelInvPosEnd,everything()))


View(lineLocChains2)
write.csv(lineLocChains2,'~/data/sv/cohort/line_chain_locations.csv',row.names=F,quote=F)

lineLocPlotData = lineLocChains2 %>% filter(LineId==25) # 3

View(lineLocPlotData %>% filter(!HasBothSource) %>% group_by(ChainDesc) %>% count)

lineLocPlotData = lineLocPlotData %>% mutate(SortPos=ifelse(!is.na(RelSrcPosStart),RelSrcPosStart,RelSrcPosReverseStart)) %>% arrange(SortPos,RelSrcSolo)
                                             
lineLocPlotData$PosIndex = row.names(lineLocPlotData)
lineLocPlotData = lineLocPlotData %>% mutate(PosIndex=as.numeric(PosIndex))
View(lineLocPlotData)
View(lineLocPlotData %>% group_by(ChainDesc,ChainSvCount) %>% count)
View(lineLocPlotData)

View(lineChainPlotData)
View(lineLocPlotData %>% filter(!is.na(RelSrcPosStart)))
     
print(ggplot(lineLocPlotData %>% filter(),aes(x=PosIndex)) + 
        geom_point(aes(y=RelSrcPosStart,color='red'),shape=17) +
        geom_point(aes(y=RelSrcPosReverseStart,color='red4'),shape=6) +
        geom_point(aes(y=RelSrcPosEnd,color='red'),shape=6) +
        geom_point(aes(y=RelInvPosStart,color='blue'),shape=17) +
        geom_point(aes(y=RelInvPosEnd,color='blue'),shape=6) +
        geom_point(aes(y=RelSrcSolo,color='green'),shape=6) +
        ylim(-7000,7000))

print(ggplot(lineLocPlotData %>% filter(!is.na(RelSrcPosStart)),aes(x=PosIndex)) + 
        geom_point(aes(y=RelSrcPosStart,color='RelSrcPosStart'),shape=17) +
        geom_point(aes(y=RelSrcPosEnd,color='RelSrcPosEnd'),shape=6) +
        geom_point(aes(y=RelInvPosStart,color='RelInvPosStart'),shape=17) +
        geom_point(aes(y=RelInvPosEnd,color='RelInvPosEnd'),shape=6) +
        # scale_color_manual(labels = c('RelSrcPosStart','RelSrcPosEnd','RelInvPosStart','RelInvPosEnd'), values = c('red','red','blue','blue')) +
        ylim(-7000,7000))


# trianges: 2 no-fill up, 6 no-fill down, 17 fill up



lineChainPlotData = lineChainPlotData %>% arrange(RelPosition)
lineChainPlotData$PosIndex = row.names(lineChainPlotData)
lineChainPlotData = lineChainPlotData %>% mutate(PosIndex=as.numeric(PosIndex))
View(lineChainPlotData)
View(lineChainPlotData %>% group_by(ChainDesc,ChainSvCount) %>% count)

lineChainPlotData = lineChainPlotData %>% mutate(OrientType=ifelse(SvOrient==RmStrand,'FORWARD','REVERSE'),
                                                 PosType=sprintf('%s_%s_%s',ChainDesc,SvType,OrientType))

View(lineChainPlotData)

lineChainPlotData2 = lineChainPlotData %>% select(LineId,SampleId,ClusterId,ChainSvCount,ChainDesc,PosType,PosIndex,RelPosition) %>%
  gather('PosType',)

print(ggplot(lineChainPos,aes(x=PosIndex)) + 
        geom_point(aes(y=RelSrcPosStart,color='RelSrcPosStart')) +
        geom_point(aes(y=RelSrcPosEnd,color='RelSrcPosEnd')) +
        geom_point(aes(y=RelInvPosStart,color='RelInvPosStart')) +
        geom_point(aes(y=RelInvPosEnd,color='RelInvPosEnd')) +
        ylim(-5000,5000))



View(lineLocChains %>% filter(Chromosome=='X'&RefPosStart>11725366-5e3&RefPosStart<11731400+5e3) %>%
       select(Chromosome,RefPosStart,RefPosEnd,SourceChr,SourcePosStart,SourcePosEnd,SampleId,ClusterId,ChainDesc,everything())
     ) # %>% 
  filter(ChainDesc=='BND=2') # |ChainDesc=='BND=1_SGL=1'

write.csv(intactLineChains,'~/data/sv/cohort/line_location_data.csv',row.names = F,quote = F)


lineLocSubset = lineLocSubset %>% filter(SourceChr=='X'&SourcePosStart>11725366-5e3&SourcePosStart<11731400+5e3&RefPosStart>11725366-5e3&RefPosStart<11731400+5e3) %>% 
  filter(ChainDesc=='BND=2') # |ChainDesc=='BND=1_SGL=1'

lineLocPlotData = lineLocSubset %>% filter(RmId==1428876) %>% filter(ChainDesc=='BND=2_INV=1') # |ChainDesc=='BND=1_SGL=1'


View(intactLinePlotData)
intactLinePlotData = intactLinePlotData %>% filter(RefPosStart>0) %>% arrange(RelSrcPosStart)
intactLinePlotData$PosIndex = row.names(intactLinePlotData)
intactLinePlotData = intactLinePlotData %>% mutate(PosIndex=as.numeric(PosIndex))

View(intactLinePlotData)

print(ggplot(intactLinePlotData,aes(x=PosIndex)) + 
        geom_point(aes(y=RelSrcPosStart,color='RelSrcPosStart')) +
        geom_point(aes(y=RelSrcPosEnd,color='RelSrcPosEnd')) +
        geom_point(aes(y=RelInvPosStart,color='RelInvPosStart')) +
        geom_point(aes(y=RelInvPosEnd,color='RelInvPosEnd')) +
        ylim(-5000,5000))



intactLocSummary = intactLocations %>% group_by(RmId,Chromosome,RefPosStart,RefPosEnd,LineLength,SampleCount,PcawgSampleCount) %>% count %>% arrange(-SampleCount)
topIntactLocations = head(intactLocSummary,20)
View(intactLocSummary)

View(intactLocations %>% group_by(LineId,Type,Chromosome,SampleCount,PcawgSampleCount,RmId,RefPosStart,RefPosEnd,LineLength) %>% 
       summarise(MaxSampleInserts=max(SampleInserts),
                 TotalInserts=first(TotalInserts)))


intactLineChains = merge(lineChains,intactLocations %>% select(LineId,RmId,SampleId,ClusterId,Chromosome,RefPosStart,RefPosEnd,RmStrand,SampleCount,PcawgSampleCount),
                         by=c('SampleId','ClusterId'),all.y=T)

nrow(intactLineChains)

intactLineChains = intactLineChains %>% filter(!is.na(ChainId)) # remove any sites which didn't map to a chained LINE cluster
nrow(intactLineChains)

# remove any cluster data not relating to the intact line location
nrow(intactLineChains %>% filter(as.character(SourceChr)==as.character(Chromosome)&abs(SourcePosStart-RefPosStart<1e4)&abs(SourcePosEnd-RefPosEnd<1e4)))
nrow(intactLineChains %>% filter(!(as.character(SourceChr)==as.character(Chromosome)&abs(SourcePosStart-RefPosStart<1e4)&abs(SourcePosEnd-RefPosEnd<1e4))))
intactLineChains = intactLineChains %>% filter(as.character(SourceChr)==as.character(Chromosome)&abs(SourcePosStart-RefPosStart)<1e4&abs(SourcePosEnd-RefPosEnd)<1e4)

View(intactLineChains)
nrow(intactLineChains)

intactLineChains = intactLineChains %>% mutate(RelSrcPosStart=ifelse(RmStrand==1,SourcePosStart-RefPosEnd,RefPosStart-SourcePosStart),
                                             RelSrcPosEnd=ifelse(RmStrand==1,SourcePosEnd-RefPosEnd,RefPosStart-SourcePosEnd),
                                             RelInvPosStart=ifelse(SourceInvPosStart>0,ifelse(RmStrand==1,SourceInvPosStart-RefPosEnd,RefPosStart-SourceInvPosStart),0),
                                             RelInvPosEnd=ifelse(SourceInvPosEnd>0,ifelse(RmStrand==1,SourceInvPosEnd-RefPosEnd,RefPosStart-SourceInvPosEnd),0))

write.csv(intactLineChains,'~/data/sv/cohort/line_location_data.csv',row.names = F,quote = F)


intactLinePlotData = intactLineChains %>% filter(SourceChr=='X'&SourcePosStart>11725366-5e3&SourcePosStart<11731400+5e3&RefPosStart>11725366-5e3&RefPosStart<11731400+5e3) %>% 
  filter(ChainDesc=='BND=2') # |ChainDesc=='BND=1_SGL=1'

intactLinePlotData = intactLineChains %>% filter(RmId==1428876) %>% filter(ChainDesc=='BND=2_INV=1') # |ChainDesc=='BND=1_SGL=1'


View(intactLinePlotData)
intactLinePlotData = intactLinePlotData %>% filter(RefPosStart>0) %>% arrange(RelSrcPosStart)
intactLinePlotData$PosIndex = row.names(intactLinePlotData)
intactLinePlotData = intactLinePlotData %>% mutate(PosIndex=as.numeric(PosIndex))

View(intactLinePlotData)

print(ggplot(intactLinePlotData,aes(x=PosIndex)) + 
        geom_point(aes(y=RelSrcPosStart,color='RelSrcPosStart')) +
        geom_point(aes(y=RelSrcPosEnd,color='RelSrcPosEnd')) +
        geom_point(aes(y=RelInvPosStart,color='RelInvPosStart')) +
        geom_point(aes(y=RelInvPosEnd,color='RelInvPosEnd')) +
        ylim(-5000,5000))

plotIndex = 1
linePlots = list()

for(i in 1:nrow(head(intactLocSummary,10)))
{
  locData = intactLocSummary[i,]
  
  print(sprintf('generating plot for %d: chr%s: %d - %d', locData$RmId,locData$Chromosome,locData$RefPosStart,locData$RefPosEnd))
  
  plotData = intactLineChains %>% filter(RmId==locData$RmId) %>% filter(ChainDesc=='BND=2_INV=1') # |ChainDesc=='BND=1_SGL=1'
  # plotData = intactLineChains %>% filter(RmId==locData$RmId) %>% filter(ChainDesc=='BND=2') # |ChainDesc=='BND=1_SGL=1'
  
  if(nrow(plotData)>0)
  {
    plotData = plotData %>% filter(RefPosStart>0) %>% arrange(RelSrcPosStart)
    plotData$PosIndex = row.names(plotData)
    plotData = plotData %>% mutate(PosIndex=as.numeric(PosIndex))
    
    plot = ggplot(plotData,aes(x=PosIndex)) + 
      geom_point(aes(y=RelSrcPosStart,color='RelSrcPosStart')) +
      geom_point(aes(y=RelSrcPosEnd,color='RelSrcPosEnd')) +
      geom_point(aes(y=RelInvPosStart,color='RelInvPosStart')) +
      geom_point(aes(y=RelInvPosEnd,color='RelInvPosEnd')) +
      ylim(-5000,5000) +
      labs(title=sprintf('%d: chr%s: %d - %d', locData$RmId,locData$Chromosome,locData$RefPosStart,locData$RefPosEnd))
  
    linePlots[[plotIndex]] = plot
    plotIndex = plotIndex + 1
  }
}

plot_grid(linePlots[[1]],linePlots[[2]],linePlots[[3]],linePlots[[4]],linePlots[[5]],linePlots[[6]],linePlots[[7]],linePlots[[8]],nrow=4,ncol=2)

# DATA OUTPUT TO PDF
library(grid)
library(gridExtra)

outputFile = '~/data/sv/cohort/line_locations_with_invs.pdf'
pdf(file=outputFile, height = 14, width = 20)

par(mar=c(1,1,1,1))

for(i in 1:3)
{
  grid.arrange(linePlots[[i*2-1]],linePlots[[i*2]], ncol=1, nrow=2, newpage=TRUE)
}


dev.off()


View(lineChains %>% filter(ChainDesc=='BND=2_INV=1') %>% group_by(SourceOrientStart,SourceOrientEnd,SourceInvOrient) %>% count)
View(lineChains %>% filter(ChainDesc=='BND=2'&SourceChr==12) %>% group_by(SourceOrientStart,SourceOrientEnd) %>% count)


## Deletion Bridges
lineDBs = rbind(lineSvData %>% filter(LEStart=='NONE'&DBLenStart>-1000&DBLenStart<5000) %>% select(SampleId,ClusterId,ClusterCount,Id,Type,DbLength=DBLenStart,InsertSeq,ChainId),
                lineSvData %>% filter(LEEnd=='NONE'&DBLenEnd>-1000&DBLenEnd<5000) %>% select(SampleId,ClusterId,ClusterCount,Id,Type,DbLength=DBLenEnd,InsertSeq,ChainId))

lineDBs = lineDBs %>% mutate(DbLengthBucket=ifelse(DbLength==0,0,ifelse(DbLength<0,-2**round(log(abs(DbLength),2)),2**round(log(DbLength,2)))))
View(lineDBs)

View(lineDBs)

print(ggplot(lineDBs %>% group_by(DbLengthBucket) %>% count,aes(x=DbLengthBucket,y=n))
      + geom_line()
      # + scale_x_log10()
      + labs(x='Deletion Bridge Length',y='',title='LINE Insertion Deletion Bridge Lengths')
      + xlim(-100,200)
      #+ facet_wrap(~TeloSglCount)
      )

lineDBs = lineDBs %>% mutate(InsSeqLength=stri_length(InsertSeq),
                             PolyAs=stri_count_fixed(InsertSeq,'A'),
                             PolyTs=stri_count_fixed(InsertSeq,'T'),
                             PolyATPerc=ifelse(InsSeqLength>0,pmax(PolyAs,PolyTs)/InsSeqLength,0),
                             PolyATPercBucket=round(PolyATPerc,1),
                             PolyATOnly=PolyATPerc>=0.9)

View(lineDBs)

print(ggplot(lineDBs %>% filter(Type=='BND') %>% group_by(DbLength,PolyATPercBucket) %>% count,aes(x=DbLength,y=n))
      + geom_line()
      # + scale_x_log10()
      + labs(x='Deletion Bridge Length',y='',title='LINE Insertion Deletion Bridge Lengths by Poly-AT % (only BNDs)')
      + xlim(-25,50)
      + facet_wrap(~PolyATPercBucket))

lineDBByLoc = merge(lineDBs,lineLocations %>% select(SampleId,ClusterId,LineId,Type,SampleCount,Chromosome,SrcPosStart=PosStart,SrcPosEnd=PosEnd),by=c('SampleId','ClusterId'))
lineDBByLoc = lineDBByLoc %>% filter(!is.na(LineId)) %>% mutate(LineLoc=sprintf('%d: %s %d-%d sc=%d',LineId,Chromosome,SrcPosStart,SrcPosEnd,SampleCount))

topLocs = head(lineLocations %>% group_by(LineId,SampleCount) %>% count %>% arrange(-n),20)

print(ggplot(lineDBByLoc %>% filter(LineId %in% topLocs$LineId) %>% group_by(DbLength,LineLoc) %>% count,aes(x=DbLength,y=n))
      + geom_line()
      # + scale_x_log10()
      + labs(x='Deletion Bridge Length',y='',title='LINE Insertion Deletion Bridge Lengths')
      + xlim(-25,50)
      + facet_wrap(~LineLoc))

View(lineChains %>% filter(SourceInvPosStart>0&SourceInvPosEnd>0))

lineDBChained = merge(lineDBs,lineChains %>% filter(ChainSvCount>1) %>% select(SampleId,ClusterId,ChainId,ChainDesc,ChainSvCount,SourcePosStart,SourcePosEnd),
                    by=c('SampleId','ClusterId','ChainId'),all.x=T)
lineDBChained = lineDBChained %>% mutate(InChain=!is.na(ChainDesc))
View(lineDBChained)

print(ggplot(lineDBChained %>% filter(InChain & ChainDesc %in% c('BND=2','BND=2_INV=1','BND=1_SGL=1')) %>% group_by(DbLength,ChainDesc) %>% count,aes(x=DbLength,y=n))
      + geom_line()
      # + scale_x_log10()
      + labs(x='Deletion Bridge Length',y='',title='LINE Insertion Deletion Bridge Lengths by Chain Type')
      + xlim(-25,50)
      + facet_wrap(~ChainDesc))

lineDBChained = merge(lineDBChained,lineLocations %>% select(SampleId,ClusterId,LineId,Type,SampleCount,Chromosome,SrcPosStart=PosStart,SrcPosEnd=PosEnd,KnownPosStart,KnownPosEnd),
                      by=c('SampleId','ClusterId'))

lineDBChained = lineDBChained %>% filter(!is.na(LineId))
lineDBChained = lineDBChained %>% mutate(LineLoc=sprintf('%d: %s %d-%d sc=%d',LineId,Chromosome,SrcPosStart,SrcPosEnd,SampleCount))

print(ggplot(lineDBChained %>% filter(InChain&ChainDesc=='BND=2') %>% filter(LineId==3|LineId==25) %>% group_by(DbLength,LineId) %>% count,
             aes(x=DbLength,y=n))
      + geom_line()
      # + scale_x_log10()
      + labs(x='Deletion Bridge Length',y='',title='LINE Insertion Deletion Bridge Lengths by Chain Type')
      + xlim(-25,50)
      + facet_wrap(~LineLoc))


View(lineDBChained %>% filter(InChain&ChainDesc=='BND=2') %>% filter(LineId==3) %>% group_by(SourcePosEnd=100*round(SourcePosEnd/100)) %>% count)

lineDBChainedByLoc = lineDBChained %>% filter(InChain&ChainDesc=='BND=2') %>% filter(LineId==3|LineId==25) %>% 
  mutate(RelSourceEnd=SourcePosEnd-KnownPosEnd) %>% filter(RelSourceEnd>0)
  
View(lineDBChainedByLoc)

print(ggplot(lineDBChainedByLoc %>% filter(LineId==3) %>% group_by(DbLength,RelSourceEnd=100*round(RelSourceEnd/100)) %>% count %>% filter(n>=10),
             aes(x=DbLength,y=n))
      + geom_line()
      # + scale_x_log10()
      + labs(x='Deletion Bridge Length',y='',title='LINE Insertion Deletion Bridge Lengths by Source Pos End')
      + xlim(-25,50)
      + facet_wrap(~RelSourceEnd))


View(lineSvData %>% filter(SampleId=='CPCT02020277TII'&ClusterId==587))



#####
## PCAWG insert comparisions
pcawgSampleIds = read.csv('~/data/sv/pcawg/pcawg_sample_ids.csv')
linxData = read.csv('~/data/sv/pcawg/LNX_LINE_CHAINS_PCAWG.csv')
nrow(linxData)
View(linxData %>% group_by(SampleId) %>% count)
View(linxData %>% group_by(ChainDesc,ChainComplete) %>% count %>% spread(ChainComplete,n,fill=0))
View(linxData)

pcawgData = read.csv('~/data/sv/pcawg/pcawg_insert_sites.csv')
pcawgData = pcawgData %>% mutate(SampleId=paste(SampleId,'T',sep=''))
View(pcawgData)
pcawgData2 = pcawgData %>% filter(SampleId %in% pcawgSampleIds$SampleId)
nrow(pcawgData2)

View(pcawgData2 %>% group_by(SampleId) %>% count)
write.csv(pcawgData2,'~/data/sv/pcawg/pcawg_line_insert_data.csv',row.names = F,quote = F)


compareResults = read.csv('~/data/sv/pcawg/LINE_INSERT_COMPARE.csv')
compareResults = read.csv('~/data/sv/pcawg/LINE_INSERT_COMPARE_1KQS.csv')
View(compareResults)
View(compareResults %>% group_by(SampleId,MatchType) %>% count %>% spread(MatchType,n,fill=0))

compareResults = compareResults %>% filter(SampleId %in% pcawgLineSampleIds$SampleId)
View(compareResults %>% group_by(SampleId) %>% count)

lineCompareSamples = compareResults %>% group_by(SampleId,MatchType) %>% count
lineCompareSamples = merge(lineCompareSamples,compareResults %>% group_by(SampleId) %>% summarise(Total=n()),by='SampleId',all.x=T)
View(lineCompareSamples)

## Key Results Plot
print(ggplot(lineCompareSamples,aes(x=reorder(SampleId,-Total),y=n,fill=MatchType))
      + geom_bar(stat='identity') 
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + labs(title='LINE Insertion Location Comparison',x='',y=''))


## comparing QS sum of 400 vs 1000
compareResultsLowQs = compareResults
# compareResultsLowQs = read.csv('~/data/sv/pcawg/LINE_INSERT_COMPARE.csv')
View(compareResultsLowQs)
View(compareResultsLowQs %>% filter(ChainId==-1) %>% group_by(InsertType,HasInversion) %>% count)

View(compareResultsPrev %>% group_by(SampleId,MatchType) %>% count %>% spread(MatchType,n,fill=0))

compareSamples = merge(compareResults %>% group_by(SampleId,MatchType) %>% count %>% spread(MatchType,n,fill=0),
                       compareResultsPrev %>% group_by(SampleId,MatchType) %>% count %>% spread(MatchType,n,fill=0) %>% select(SampleId,MATCH_1KQS=MATCH,LINX_1KQS=LINX,PCAWG_1KQS=PCAWG),
                       by='SampleId',all.x=T)

View(compareSamples %>% select(SampleId,MATCH,MATCH_1KQS,LINX,LINX_1KQS,PCAWG,PCAWG_1KQS))
View(compareResults %>% group_by(SampleId,MatchType) %>% count %>% spread(MatchType,n,fill=0))

View(pcawgSvData %>% filter(SampleId=='DO224638T'&ResolvedType=='LINE') %>% group_by(ClusterId) %>% mutate(clusterQual=sum(QualScore)) %>% ungroup() %>% 
       select(SampleId,ClusterId,ClusterCount,ClusterDesc,ResolvedType,Id,Type,ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd,
              LEStart,LEEnd,InsertSeq,RepeatClass,QualScore,clusterQual,everything()))

pcawgLineSampleIds = read.csv('~/data/sv/pcawg/pcawg_line_sample_ids.csv')
View(pcawgLineSampleIds)



pcawgSvData = read.csv('~/data/sv/pcawg/LNX_SVS_PCAWG_1KQS.csv')
pcawgSvData = read.csv('~/data/sv/pcawg/LNX_SVS.csv')
nrow(pcawgSvData)
View(pcawgSvData)
colnames(pcawgSvData)

pcawgSvData = pcawgSvData %>% filter(SampleId %in% pcawgLineSampleIds$SampleId)
View(pcawgSvData %>% group_by(SampleId) %>% count)

# low qual-score but rescued due to in DB and having poly A/T tail
View(pcawgSvData %>% filter(Type=='SGL'&QualScore<=1000) %>% group_by(ResolvedType,HasDB=DBLenStart>-1000) %>% count)

pcawgSvData = merge(pcawgSvData,compareResults %>% select(SampleId,ClusterId,MatchType),by=c('SampleId','ClusterId'),all.x=T)

#View(pcawgSvData %>% filter(Type=='SGL'&QualScore<1000) %>% 
#       group_by(ResolvedType,HasDB=DBLenStart>-1000,PcawgMatch=ifelse(is.na(MatchType)|MatchType=='LINE','LINX_ONLY','MATCHED')) %>% count)


pcawgSvData = pcawgSvData %>% mutate(PosStartRnd=round(PosStart,-3),PosEndRnd=round(PosEnd,-3))

unmatchedSites = compareResults %>% filter(ChainId==-1) %>% mutate(InsertPositionBucket=round(InsertPosStart,-3))

unmatchedSites$UnmatchedId = row.names(unmatchedSites)

View(unmatchedSites)

unmatchedSvData = merge(pcawgSvData,unmatchedSites %>% select(SampleId,UnmatchedId,ChrStart=InsertChr,InsertPosStart,InsertPosEnd,PosStartRnd=InsertPositionBucket,InsertType),
             by=c('SampleId','ChrStart','PosStartRnd'),all.x=T) 

unmatchedSvData = merge(unmatchedSvData,unmatchedSites %>% 
                          select(SampleId,UnmatchedId2=UnmatchedId,ChrEnd=InsertChr,InsertPosStart2=InsertPosStart,InsertPosEnd2=InsertPosEnd,PosEndRnd=InsertPositionBucket,InsertType2=InsertType),
                        by=c('SampleId','ChrEnd','PosEndRnd'),all.x=T) 

View(unmatchedSvData)
View(unmatchedSvData %>% filter(!is.na(UnmatchedId)|!is.na(UnmatchedId2)))

noMatchSvData = unmatchedSvData %>% filter(is.na(UnmatchedId)&is.na(UnmatchedId2))
View(noMatchSvData)

View(unmatchedSites %>% filter(!(UnmatchedId %in% unmatchedSvData$UnmatchedId)&!(UnmatchedId %in% unmatchedSvData$UnmatchedId2)))

View(unmatchedSites %>% filter(!(UnmatchedId %in% unmatchedSvData$UnmatchedId)&!(UnmatchedId %in% unmatchedSvData$UnmatchedId2)) %>%
       group_by(SampleId,InsertType) %>% count %>% spread(InsertType,n,fill=0))

View(unmatchedSites %>% filter(!(UnmatchedId %in% unmatchedSvData$UnmatchedId)&!(UnmatchedId %in% unmatchedSvData$UnmatchedId2)) %>%
       filter(SampleId=='DO50388T'))

manMatchSvData = unmatchedSvData %>% filter(!is.na(UnmatchedId)|!is.na(UnmatchedId2))

manMatchSvData = manMatchSvData %>% mutate(UnmatchedId=ifelse(!is.na(UnmatchedId),UnmatchedId,UnmatchedId2),
                                           InsertType=ifelse(!is.na(InsertType),as.character(InsertType),as.character(InsertType2)),
                                           InsertPosStart=ifelse(!is.na(InsertPosStart),InsertPosStart,InsertPosStart2),
                                           InsertPosEnd=ifelse(!is.na(InsertPosEnd),InsertPosEnd,InsertPosEnd2))

View(manMatchSvData)
nrow(manMatchSvData)

View(manMatchSvData %>% select(SampleId,ClusterId,ClusterCount,ClusterDesc,ResolvedType,Id,Type,InsertType,LEStart,LEEnd,InsertSeq,RepeatClass,QualScore,DBLenStart,DBLenEnd,
                               CNChgStart,CNStart,UnmatchedId,ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd,InsertPosStart,InsertPosEnd,everything()))

View(manMatchSvData %>% group_by(ClusterDesc,ResolvedType,InsertType) %>% count)

View(manMatchSvData %>% filter(SampleId=='DO50388T') %>% group_by(ClusterId,ClusterDesc,ResolvedType,InsertType) %>% count)




# manual linking for X:11.7M
x117LineChains = lineChains %>% filter(SourceChr=='X'&SourcePosStart>11725366-5e3&SourcePosStart<11731400+5000) %>%
  mutate(KnownPosStart=11725366,KnownPosEnd=11731400,Strand='+')

x117LineChains = x117LineChains %>% filter(SourcePosEnd<0|abs(SourcePosEnd-SourcePosStart)<1e4)

x117LineChains = x117LineChains %>% mutate(RelSourcePosStart=ifelse(Strand=='+',SourcePosStart-KnownPosEnd,KnownPosStart-SourcePosStart),
                                           RelSourcePosEnd=ifelse(SourcePosEnd>0,ifelse(Strand=='+',SourcePosEnd-KnownPosEnd,KnownPosStart-SourcePosEnd),0),
                                           RelInvPosStart=ifelse(SourceInvPosStart>0,ifelse(Strand=='+',SourceInvPosStart-KnownPosEnd,KnownPosStart-SourceInvPosStart),0),
                                           RelInvPosEnd=ifelse(SourceInvPosEnd>0,ifelse(Strand=='+',SourceInvPosEnd-KnownPosEnd,KnownPosStart-SourceInvPosEnd),0))

View(x117LineChains)
View(x117LineChains %>% filter(ChainSvCount>1))
View(x117LineChains %>% group_by(ChainSvCount) %>% count)
View(x117LineChains %>% group_by(ChainDesc) %>% count)
View(x117LineChains %>% filter(SourceInvPosStart>11725366-5e3&SourceInvPosStart<11731400+5000&SourceInvPosEnd>11725366-5e3&SourceInvPosEnd<11731400+5000))
View(x117LineChains %>% filter(ChainDesc=='DEL'))

write.csv(x117LineChains,'~/logs/known_line_chains_x117.csv',row.names = F,quote = F)

relPosPlotData = x117LineChains %>% filter(ChainDesc=='BND=2'|ChainDesc=='BND=1_SGL=1') %>% arrange(SourcePosStart)
relPosPlotData = x117LineChains %>% filter(SourceInvPosStart>0&SourceInvPosEnd>0&SourcePosEnd>0) %>% arrange(SourcePosStart)
relPosPlotData$PosIndex = row.names(relPosPlotData)
relPosPlotData = relPosPlotData %>% mutate(PosIndex=as.numeric(PosIndex))
View(relPosPlotData)

print(ggplot(relPosPlotData,aes(x=PosIndex)) + 
        geom_point(aes(y=RelSourcePosStart,color='RelSourcePosStart')) +
        geom_point(aes(y=RelSourcePosEnd,color='RelSourcePosEnd')) +
        geom_point(aes(y=RelInvPosStart,color='RelInvPosStart')) + 
        geom_point(aes(y=RelInvPosEnd,color='RelInvPosEnd')) +
        ylim(-5000,5000))

relPosPlotData2 = x117LineChains %>% filter(ChainSvCount==1&SourcePosStart>0) %>% arrange(SourcePosStart)
relPosPlotData2$PosIndex = row.names(relPosPlotData2)
relPosPlotData2 = relPosPlotData2 %>% mutate(PosIndex=as.numeric(PosIndex))
View(relPosPlotData2)

print(ggplot(relPosPlotData2,aes(x=PosIndex)) + 
        geom_point(aes(y=RelSourcePosStart,color='RelSourcePosStart')) +
        facet_wrap(~LoneBeOrient))







######
## DEBUG



## New Suspect Rule Validation
svData = read.csv('~/data/sv/cohort/LNX_SVS.csv')
svData = svData %>% filter(SampleId %in% samplesDD$SampleId)
svData = svData %>% filter(!(ResolvedType %in% c('LOW_VAF','DUP_BE'))) # remove artefect types

# lineSvData = svData %>% filter(ResolvedType=='LINE'|LEStart!='NONE'|LEEnd!='NONE')
lineSvData = svData %>% filter(ResolvedType=='LINE')
View(lineSvData %>% filter(ClusterDesc=='BND=1_SGL=1') %>% select(SampleId,ClusterId,Type,LEStart,LEEnd,InsertSeq,Annotations,everything()))
nrow(lineSvData)
View(lineSvData %>% filter(LEStart=='SUSPECT'|LEEnd=='SUSPECT') %>% group_by(Annotations) %>% count)
View(lineSvData %>% filter(ClusterCount<=3) %>% group_by(SampleId,ClusterId,ClusterDesc) %>% count %>% group_by(ClusterDesc) %>% count)

rm(svOldData)
rm(lineOldSvData)
rm(lineCompare)
svOldData = read.csv('~/logs/LNX_SVS_pre_line_changes.csv')
# svOldData = read.csv('~/data/sv/cohort/LNX_SVS_PROD.csv')
svOldData = svOldData %>% filter(SampleId %in% samplesDD$SampleId)
svOldData = svOldData %>% filter(!(ResolvedType %in% c('LOW_VAF','DUP_BE'))) # remove artefect types

svOldData = svOldData %>% mutate(LEStart=ifelse(LEStart=='None','NONE',ifelse(LEStart=='Known','KNOWN',ifelse(LEStart=='Suspect','SUSPECT','KNOWN;SUSPECT'))),
                                 LEEnd=ifelse(LEEnd=='None','NONE',ifelse(LEEnd=='Known','KNOWN',ifelse(LEEnd=='Suspect','SUSPECT','KNOWN;SUSPECT'))))

#lineOldSvData = svOldData %>% filter(ResolvedType=='LINE'|LEStart!='None'|LEEnd!='None')
#nrow(lineOldSvData)
#View(lineOldSvData)

View(lineOldSvData %>% group_by(LEStart,LEEnd) %>% count)
View(lineOldSvData %>% group_by(SampleId) %>% count)

lineCompare = merge(svData,svOldData %>% select(SampleId,Id,OldClusterDesc=ClusterDesc,OldResolvedType=ResolvedType,OldLEStart=LEStart,OldLEEnd=LEEnd),
                    by=c('SampleId','Id'),all=T)
nrow(lineCompare)

# old take records which were or are now LINE
lineCompare = lineCompare %>% filter(ResolvedType=='LINE'|LEStart!='NONE'|LEEnd!='NONE'|OldResolvedType=='LINE'|OldLEStart!='NONE'|OldLEEnd!='NONE')
nrow(lineCompare)

write.csv(lineCompare %>% filter(as.character(ClusterDesc)!=as.character(OldClusterDesc)|as.character(ResolvedType)!=as.character(OldResolvedType)
                                 |as.character(LEStart)!=as.character(OldLEStart)|as.character(LEEnd)!=as.character(OldLEEnd)),
          '~/data/sv/cohort/sv_data_line_changes.csv',row.names = F,quote = F)


View(lineCompare %>% filter(as.character(ClusterDesc)!=as.character(OldClusterDesc)|as.character(ResolvedType)!=as.character(OldResolvedType)
                            |as.character(LEStart)!=as.character(OldLEStart)|as.character(LEEnd)!=as.character(OldLEEnd)))
View(lineCompare)
nrow(lineCompare %>% filter(as.character(ClusterDesc)==as.character(OldClusterDesc)&as.character(ResolvedType)==as.character(OldResolvedType)
                            &as.character(LEStart)==as.character(OldLEStart)&as.character(LEEnd)==as.character(OldLEEnd)))

# any difference
View(lineCompare %>% filter(as.character(ClusterDesc)!=as.character(OldClusterDesc)|as.character(ResolvedType)!=as.character(OldResolvedType)
                            |as.character(LEStart)!=as.character(OldLEStart)|as.character(LEEnd)!=as.character(OldLEEnd)) %>%
       select(SampleId,Id,OldClusterDesc,ClusterDesc,OldResolvedType,ResolvedType,OldLEStart,LEStart,OldLEEnd,LEEnd,Annotations))

# cluster resovled diffs
View(lineCompare %>% filter(as.character(ClusterDesc)!=as.character(OldClusterDesc)|as.character(ResolvedType)!=as.character(OldResolvedType)) %>%
       select(SampleId,Id,OldClusterDesc,ClusterDesc,OldResolvedType,ResolvedType,OldLEStart,LEStart,OldLEEnd,LEEnd,Annotations))

View(lineCompare %>% filter(as.character(ClusterDesc)!=as.character(OldClusterDesc)|as.character(ResolvedType)!=as.character(OldResolvedType)) %>%
       group_by(OldResolvedType,ResolvedType) %>% count)

View(lineCompare %>% filter(as.character(ClusterDesc)!=as.character(OldClusterDesc)|as.character(ResolvedType)!=as.character(OldResolvedType)) %>%
       group_by(SampleId,OldResolvedType,ResolvedType) %>% count)

View(lineCompare %>% filter(ResolvedType!=OldResolvedType&ResolvedType=='COMPLEX'&OldResolvedType=='LINE') %>% group_by(ClusterDesc,OldClusterDesc) %>% count)

View(lineCompare %>% filter(ResolvedType!=OldResolvedType&ResolvedType=='COMPLEX'&OldResolvedType=='LINE') %>% filter(ClusterCount<=4) %>%
       select(SampleId,Id,Type,ClusterId,ClusterCount,ClusterDesc,OldClusterDesc,ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd,InsertSeq,LEStart,LEEnd,RepeatClass,Annotations,everything()) %>% 
       arrange(SampleId,ClusterId))



View(lineCompare %>% filter(ClusterDesc=='BND=39_INV=3_SGL=5') %>% select(LEStart,OldLEStart,LEEnd,OldLEEnd,SampleId,Id,Type,ClusterId,ClusterCount,ClusterDesc,OldClusterDesc,ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd,InsertSeq,LEStart,LEEnd,RepeatClass,everything()))

View(svData %>% filter(SampleId=='CPCT02330132T'&ResolvedType=='SGL_PAIR_INS') %>% select(SampleId,Id,OrientStart,InsertSeq,Annotations,everything()))
View(svData %>% filter(SampleId=='CPCT02010022T'&ClusterId==155) %>% select(SampleId,Id,Type,ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd,InsertSeq,Annotations,everything()))

View(lineCompare %>% filter(as.character(ClusterDesc)!=as.character(OldClusterDesc)|as.character(ResolvedType)!=as.character(OldResolvedType)) %>%
       group_by(SampleId,Id,OldResolvedType,ResolvedType) %>% count %>% group_by(OldResolvedType,ResolvedType) %>% count)

View(svData %>% filter(Type!='SGL'&Type!='INF') %>% filter(LEStart!='NONE'|LEEnd!='NONE') %>% select(SampleId,Id,OrientStart,InsertSeq,Annotations,everything()))






# View(lineCompare %>% group_by(Match=ifelse()))

# load cohort data
lineLocations = read.csv('~/data/sv/cohort/LINE_SOURCE_DATA.csv')
lineLocations = lineLocations %>% mutate(Length=PosEnd-PosStart)
View(lineLocations)
View(lineLocations %>% filter(Length>0))
View(lineLocations %>% filter(Type=='SUSPECT'))
View(lineLocations %>% filter(SampleSourceLocations==1))
View(lineLocations %>% filter(RmId>0&KnownPosStart==-1))

View(lineLocations %>% filter(RmId>0&KnownPosStart==-1) %>% group_by(LineId,RmId,Chromosome,RmPosStart,RmPosEnd,RmStrand,PosStart,PosEnd,SampleCount,PcawgSampleCount) %>% 
       summarise(Distance=first(ifelse(RmStrand==1,abs((RmPosStart+RmPosEnd)/2-(PosStart+PosEnd)/2),abs((PosStart+PosEnd)/2-(RmPosStart+RmPosEnd)/2)))))

# location summary data
locationSummary = lineLocations %>% group_by(LineId,Type,Chromosome,PosStart,PosEnd,SampleCount,PcawgSampleCount,KnownPosStart,KnownPosEnd,RmId,Length) %>% 
  summarise(MaxSampleInserts=max(SampleInserts),
            TotalInserts=first(TotalInserts))
View(locationSummary)

View(locationSummary %>% mutate(KnownSiteLength=KnownPosEnd-KnownPosStart))

knownSites = locationSummary %>% filter(KnownPosEnd>0) %>% group_by(LineId,Chromosome,KnownPosStart,KnownPosEnd,RmId=LowerRmId) %>% 
  summarise(SampleCount=first(SampleCount),
            PcawgSampleCount=first(PcawgSampleCount)) %>% mutate(KnownSiteLength=KnownPosEnd-KnownPosStart)
View(knownSites)

knownSites = merge(knownSites,
                   rmLineData %>% select(RmId,RmPosStart=PosStart,RmPosEnd=PosEnd,MatchingRepeat,SWScore,Strand,RepeatBegin,RepeatEnd,DevPerc,DelPerc,InsPerc,ClassFamily),
                   by='RmId',all.x=T)

View(knownSites)
View(knownSites %>% mutate(RmLength=RmPosEnd-RmPosStart))
View(rmLineData %>% filter(RmId>874560&RmId<874580))


## Location of Line SVs relative to full line sites
knownSiteData = lineLocations %>% filter(KnownPosEnd-KnownPosStart>5000)

knownSiteData = merge(knownSiteData,rmLineData %>% select(LowerRmId=RmId,RmPosStart=PosStart,RmPosEnd=PosEnd,Strand),by='LowerRmId',all.x=T)
knownSiteData = knownSiteData %>% select(LineId,RmId=LowerRmId,Chromosome,Strand,KnownPosStart,KnownPosEnd,RmPosStart,RmPosEnd,SampleCount,SampleId,ClusterId,SamplePosStart,SamplePosEnd)
View(knownSiteData)

linePosData = rbind(lineSvData %>% filter(LEStart!='NONE') %>% mutate(StartEnd='Start') %>% 
                      select(SampleId,ClusterId,Id,Type,Chr=ChrStart,Pos=PosStart,Orient=OrientStart,StartEnd,LineType=LEStart,RemoteDbLen=DBLenEnd),
                    lineSvData %>% filter(LEEnd!='NONE'&Type!='SGL'&Type!='INF') %>% mutate(StartEnd='End') %>% 
                      select(SampleId,ClusterId,Id,Type,Chr=ChrEnd,Pos=PosEnd,Orient=OrientEnd,StartEnd,LineType=LEEnd,RemoteDbLen=DBLenStart))

View(linePosData)

svKnownData = merge(linePosData,knownSiteData,by=c('SampleId','ClusterId'),all.x=T)
svKnownData = svKnownData %>% filter(!is.na(RmId)) %>% filter(as.character(Chr)==as.character(Chromosome))
View(svKnownData)

svKnownData = svKnownData %>% mutate(SvDownstreamDist=ifelse(Strand=='+',Pos-RmPosEnd,RmPosStart-Pos),
                                     Orientation=ifelse((Strand=='+')==(Orient==1),'Same','Reverse'),
                                     Location=ifelse(Pos>=RmPosStart&Pos<=RmPosEnd,'INSIDE',ifelse(SvDownstreamDist>0,'DOWNSTREAM','UPSTREAM')))

svKnownData = svKnownData %>% filter(abs(SvDownstreamDist)<1e4)

View(svKnownData %>% group_by(Location,DistBucket=100*round(SvDownstreamDist/100),Type) %>% count %>% spread(Type,n,fill=0))

svKnownData = svKnownData %>% mutate(LineLabel=sprintf('%d: chr%s:%d-%d',LineId,Chromosome,RmPosStart,RmPosEnd),
                                     TypeInfo=ifelse(Type=='BND',sprintf('%s-%s',Type,Orientation),ifelse(Type=='INV',sprintf('%s-%s-%s',Type,StartEnd,Orientation),'Other')))

View(svKnownData)
write.csv(svKnownData %>% group_by(SampleId) %>% count %>% select(SampleId),'~/logs/known_line_samples.csv',row.names = F,quote = F)

View(svKnownData %>% group_by(Type,Strand,Orient,StartEnd,TypeInfo) %>% count)
View(svKnownData %>% group_by(Chromosome,Type,Strand,Orient,StartEnd,TypeInfo) %>% count)

print(ggplot(svKnownData %>% filter(Location=='DOWNSTREAM') %>% group_by(Location,DistBucket=100*round(SvDownstreamDist/100),Type) %>% count %>% spread(Type,n,fill=0),
             aes(x=DistBucket))
      + geom_line(aes(y=BND,color='BND'))
      + geom_line(aes(y=INV,color='INV')))

print(ggplot(svKnownData %>% group_by(Location,DistBucket=100*round(SvDownstreamDist/100),Type) %>% count %>% spread(Type,n,fill=0),
             aes(x=DistBucket))
      + geom_line(aes(y=BND,color='BND'))
      + geom_line(aes(y=INV,color='INV')))

svKnownDataLocs = svKnownData %>% group_by(LineLabel,SampleCount,Location,DistBucket=100*round(SvDownstreamDist/100),Type) %>% count

print(ggplot(svKnownDataLocs %>% filter(SampleCount>100) %>% spread(Type,n,fill=0),
             aes(x=DistBucket))
      + geom_line(aes(y=BND,color='BND'))
      + geom_line(aes(y=INV,color='INV'))
      + scale_y_log10()
      + xlim(-5000,5000)
      + facet_wrap(~LineLabel)
      + labs(title='SV breakend distance from end of Known Line Site'))

print(ggplot(svKnownData %>% filter(SampleCount>100&(Type=='BND'|Type=='INV')) %>%
               group_by(LineLabel,DistBucket=100*round(SvDownstreamDist/100),TypeInfo) %>% count %>% spread(TypeInfo,n,fill=0),
             aes(x=DistBucket))
      + geom_line(aes(y=`BND-Same`,color='BND-Same'))
      + geom_line(aes(y=`BND-Reverse`,color='BND-Reverse'))
      + scale_y_log10()
      + xlim(-5000,5000)
      + facet_wrap(~LineLabel)
      + labs(title='SV breakend distance from end of Known Line Site'))

print(ggplot(svKnownData %>% filter(SampleCount>100&(Type=='BND'|Type=='INV')) %>%
               group_by(LineLabel,DistBucket=100*round(SvDownstreamDist/100),TypeInfo) %>% count %>% spread(TypeInfo,n,fill=0),
             aes(x=DistBucket))
      + geom_line(aes(y=`INV-Start-Same`,color='INV-Start-Same'))
      + geom_line(aes(y=`INV-End-Same`,color='INV-End-Same'))
      + geom_line(aes(y=`INV-Start-Reverse`,color='INV-Start-Reverse'))
      + geom_line(aes(y=`INV-End-Reverse`,color='INV-End-Reverse'))
      + scale_y_log10()
      + xlim(-5000,5000)
      + facet_wrap(~LineLabel)
      + labs(title='SV breakend distance from end of Known Line Site'))


## linked to repeat masker
rmLineLocations = lineLocations %>% group_by(LineId,Type,Chromosome,PosStart,PosEnd,SampleCount,PcawgSampleCount,LowerRmId,UpperRmId,Length) %>% 
  summarise(MaxSampleInserts=max(SampleInserts),
            TotalInserts=first(TotalInserts))
View(rmLineLocations)
nrow(rmLineLocations)

rmLineLocations = merge(rmLineLocations,
                        rmLineData %>% select(LowerRmId=RmId,LowerRmPosStart=PosStart,LowerRmPosEnd=PosEnd,
                                              LowerMatchingRepeat=MatchingRepeat,LowerSWScore=SWScore,LowerStrand=Strand,LowerRepeatBegin=RepeatBegin,LowerRepeatEnd=RepeatEnd),
                        by='LowerRmId',all.x=T)

rmLineLocations = merge(rmLineLocations,
                        rmLineData %>% select(UpperRmId=RmId,UpperRmPosStart=PosStart,UpperRmPosEnd=PosEnd,
                                              UpperMatchingRepeat=MatchingRepeat,UpperSWScore=SWScore,UpperStrand=Strand,UpperRepeatBegin=RepeatBegin,UpperRepeatEnd=RepeatEnd),
                        by='UpperRmId',all.x=T)

View(rmLineLocations)


## Specific LINE site investigation
llxSvData = svData %>% filter(ResolvedType=='LINE') %>% filter((ChrStart=='X'&PosStart>=11725368&PosStart<=11736399)|(ChrEnd=='X'&PosEnd>=11725368&PosEnd<=11736399)
                                                               |grepl('X:117',InsSeqAlignments))
nrow(llxSvData)

llxPosData = rbind(llxSvData %>% mutate(StartEnd='Start') %>%
                     select(SampleId,Id,ClusterId,ClusterDesc,ClusterCount,Type,Chr=ChrStart,Pos=PosStart,Orient=OrientStart,StartEnd,LineType=LEStart,
                            DbLen=DBLenStart,InsSeqAlignments,LnkLen=LnkLenStart),
                   llxSvData %>% filter(Type!='SGL'&Type!='INF') %>% mutate(StartEnd='End') %>%
                     select(SampleId,Id,ClusterId,ClusterDesc,ClusterCount,Type,Chr=ChrEnd,Pos=PosEnd,Orient=OrientEnd,StartEnd,LineType=LEEnd,
                            DbLen=DBLenEnd,InsSeqAlignments,LnkLen=LnkLenEnd))

View(llxPosData)
View(llxPosData %>% filter(ClusterDesc=='SGL=2'))
View(llxPosData %>% filter(ClusterDesc=='BND=1_SGL=1'))

View(svData %>% filter(SampleId=='CPCT02010969T'&ClusterId==30))

llxPosData = llxPosData %>% mutate(SglMapDataPos=ifelse(Type=='SGL',regexpr('X:117',InsSeqAlignments),-1),
                                   SglMapPosition=ifelse(SglMapDataPos>0,stri_sub(InsSeqAlignments,SglMapDataPos+2,10),''),
                                   SglMapOrient=ifelse(SglMapDataPos>0,stri_sub(InsSeqAlignments,SglMapDataPos+11,12),''),
                                   SglMapOrient=ifelse(SglMapOrient=='+',1,ifelse(SglMapOrient=='-',-1,0)))

llxPosData = llxPosData %>% mutate(IsSource=LineType!='NONE')
View(llxPosData)
View(llxPosData %>% filter(ClusterDesc=='BND=1_SGL=1') %>% arrange(SampleId,ClusterId))
nrow(llxPosData)


llxClusterData = llxPosData %>% arrange(SampleId,ClusterId,-IsSource) %>% group_by(SampleId,ClusterId,ClusterCount,ClusterDesc) %>% 
  summarise(Breakends=n(),
            SourceCount=sum(IsSource),
            SourceChr=first(ifelse(IsSource,as.character(Chr),'')),
            SourcePosStart=min(ifelse(IsSource,Pos,ifelse(SglMapOrient!='0',as.numeric(SglMapPosition),1e9))),
            SourcePosEnd=max(ifelse(IsSource,Pos,ifelse(SglMapOrient!='0',as.numeric(SglMapPosition),-1))),
            InsertChr=last(ifelse(!IsSource,as.character(Chr),'')),
            InsertPosStart=min(ifelse(!IsSource,Pos,1e9)),
            InsertPosEnd=max(ifelse(!IsSource,Pos,-1)),
            SourceLinks=sum(LnkLen>0),
            InsertDBs=sum(DbLen>-1000),
            SglMapped=sum(SglMapOrient!='0'))

View(llxClusterData)
View(llxClusterData %>% filter(SampleId=='CPCT02020938T'&ClusterId==115))
View(llxClusterData %>% filter(ClusterDesc=='BND=2'&Breakends==4))
View(llxClusterData %>% filter(ClusterDesc=='BND=1_SGL=1'&Breakends==3))
View(llxClusterData %>% filter(ClusterDesc=='SGL=2'&Breakends==2))



## Repeat Masker Line Annotations
rmLineData = read.csv('~/data/sv/repeat_masker_line.csv',sep='\t')
View(rmLineData)
write.csv(rmLineData %>% select(RmId,Chromosome,PosStart,PosEnd,Strand),'~/data/sv/line_repeat_mask_data.csv', row.names = F,quote = F)


# Chaining
tmpLineChains = read.csv('~/logs/LNX_LINE_CHAINS.csv')
View(tmpLineChains)

lineChains = read.csv('~/data/sv/cohort/LNX_LINE_CHAINS.csv')
View(lineChains)
View(lineChains %>% group_by(ChainDesc) %>% count)

knownLineChains = merge(lineChains,svKnownData %>% group_by(SampleId,ClusterId,RmId,KnownPosStart,KnownPosEnd,Strand) %>% count %>% select(-n),by=c('SampleId','ClusterId'),all.x=T)
View(knownLineChains %>% filter(!is.na(RmId)))
knownLineChains = knownLineChains %>% filter(!is.na(RmId))
knownLineChains = knownLineChains %>% mutate(RelSourcePosStart=ifelse(Strand=='+',SourcePosStart-KnownPosEnd,KnownPosStart-SourcePosStart),
                                             RelSourcePosEnd=ifelse(Strand=='+',SourcePosEnd-KnownPosEnd,KnownPosStart-SourcePosEnd),
                                             RelInvPosStart=ifelse(SourceInvPosStart>0,ifelse(Strand=='+',SourceInvPosStart-KnownPosEnd,KnownPosStart-SourceInvPosStart),0),
                                             RelInvPosEnd=ifelse(SourceInvPosEnd>0,ifelse(Strand=='+',SourceInvPosEnd-KnownPosEnd,KnownPosStart-SourceInvPosEnd),0))


# manual linking for X:11.7M
x117LineChains = lineChains %>% filter(SourceChr=='X'&SourcePosStart>11725366-5e3&SourcePosStart<11731400+5000) %>%
  mutate(KnownPosStart=11725366,KnownPosEnd=11731400,Strand='+')

x117LineChains = x117LineChains %>% filter(SourcePosEnd<0|abs(SourcePosEnd-SourcePosStart)<1e4)

x117LineChains = x117LineChains %>% mutate(RelSourcePosStart=ifelse(Strand=='+',SourcePosStart-KnownPosEnd,KnownPosStart-SourcePosStart),
                                           RelSourcePosEnd=ifelse(SourcePosEnd>0,ifelse(Strand=='+',SourcePosEnd-KnownPosEnd,KnownPosStart-SourcePosEnd),0),
                                           RelInvPosStart=ifelse(SourceInvPosStart>0,ifelse(Strand=='+',SourceInvPosStart-KnownPosEnd,KnownPosStart-SourceInvPosStart),0),
                                           RelInvPosEnd=ifelse(SourceInvPosEnd>0,ifelse(Strand=='+',SourceInvPosEnd-KnownPosEnd,KnownPosStart-SourceInvPosEnd),0))

View(x117LineChains)
View(x117LineChains %>% filter(ChainSvCount>1))
View(x117LineChains %>% group_by(ChainSvCount) %>% count)
View(x117LineChains %>% group_by(ChainDesc) %>% count)
View(x117LineChains %>% filter(SampleId=='CPCT02040196T'&ClusterId==544))
View(x117LineChains %>% filter(SourceInvPosStart>11725366-5e3&SourceInvPosStart<11731400+5000&SourceInvPosEnd>11725366-5e3&SourceInvPosEnd<11731400+5000))
View(x117LineChains %>% filter(ChainDesc=='DEL'))

write.csv(x117LineChains,'~/logs/known_line_chains_x117.csv',row.names = F,quote = F)

relPosPlotData = x117LineChains %>% filter(ChainDesc=='BND=2'|ChainDesc=='BND=1_SGL=1') %>% arrange(SourcePosStart)
relPosPlotData = x117LineChains %>% filter(SourceInvPosStart>0&SourceInvPosEnd>0&SourcePosEnd>0) %>% arrange(SourcePosStart)
relPosPlotData$PosIndex = row.names(relPosPlotData)
relPosPlotData = relPosPlotData %>% mutate(PosIndex=as.numeric(PosIndex))
View(relPosPlotData)

print(ggplot(relPosPlotData,aes(x=PosIndex)) + 
        geom_point(aes(y=RelSourcePosStart,color='RelSourcePosStart')) +
        geom_point(aes(y=RelSourcePosEnd,color='RelSourcePosEnd')) +
        geom_point(aes(y=RelInvPosStart,color='RelInvPosStart')) + 
        geom_point(aes(y=RelInvPosEnd,color='RelInvPosEnd')) +
        ylim(-5000,5000))

relPosPlotData2 = x117LineChains %>% filter(ChainSvCount==1&SourcePosStart>0) %>% arrange(SourcePosStart)
relPosPlotData2$PosIndex = row.names(relPosPlotData2)
relPosPlotData2 = relPosPlotData2 %>% mutate(PosIndex=as.numeric(PosIndex))
View(relPosPlotData2)

print(ggplot(relPosPlotData2,aes(x=PosIndex)) + 
        geom_point(aes(y=RelSourcePosStart,color='RelSourcePosStart')) +
        facet_wrap(~LoneBeOrient))


print(ggplot(x117LineChains %>% filter(ChainDesc=='BND=2'), aes(x=RelSourcePosStart,y=cumsum(RelSourcePosStart))) + geom_point())

chr22LineChains = lineChains %>% filter(SourceChr==22&SourcePosStart>29059272-5e3&SourcePosStart<29065304+5000) %>%
  mutate(KnownPosStart=29059272,KnownPosEnd=29065304,Strand='+')

chr22LineChains = chr22LineChains %>% filter(SourcePosEnd<0|abs(SourcePosEnd-SourcePosStart)<1e4)

chr22LineChains = chr22LineChains %>% mutate(RelSourcePosStart=ifelse(Strand=='+',SourcePosStart-KnownPosEnd,KnownPosStart-SourcePosStart),
                                             RelSourcePosEnd=ifelse(SourcePosEnd>0,ifelse(Strand=='+',SourcePosEnd-KnownPosEnd,KnownPosStart-SourcePosEnd),0),
                                             RelInvPosStart=ifelse(SourceInvPosStart>0,ifelse(Strand=='+',SourceInvPosStart-KnownPosEnd,KnownPosStart-SourceInvPosStart),0),
                                             RelInvPosEnd=ifelse(SourceInvPosEnd>0,ifelse(Strand=='+',SourceInvPosEnd-KnownPosEnd,KnownPosStart-SourceInvPosEnd),0))

View(chr22LineChains)
View(chr22LineChains %>% group_by(ChainSvCount) %>% count)
View(chr22LineChains %>% group_by(ChainDesc) %>% count)



View(chr22LineChains %>% filter(SampleId=='CPCT02040196T'&ClusterId==544))

write.csv(chr22LineChains,'~/logs/known_line_chains_chr22.csv',row.names = F,quote = F)


knownLineChains = knownLineChains %>% mutate(LineLabel=sprintf('chr%s:%d-%d',SourceChr,KnownPosStart,KnownPosEnd))
View(knownLineChains)
View(knownLineChains %>% filter(SourceChr=='X'))

write.csv(knownLineChains %>% filter(SourceChr=='X'&RmId==1428876),'~/logs/known_line_chains_x117.csv',row.names = F,quote = F)

print(ggplot(knownLineChains %>% filter(ChainDesc=='BND=2'),
             aes(x=DistBucket))
      + geom_point(aes(y=RelSourcePosStart,color='RelSourcePosStart'))
      + geom_line(aes(y=INV,color='INV')))

print(ggplot(svKnownData %>% group_by(Location,DistBucket=100*round(SvDownstreamDist/100),Type) %>% count %>% spread(Type,n,fill=0),
             aes(x=DistBucket))
      + geom_line(aes(y=BND,color='BND'))
      + geom_line(aes(y=INV,color='INV')))

svKnownDataLocs = svKnownData %>% group_by(LineLabel,SampleCount,Location,DistBucket=100*round(SvDownstreamDist/100),Type) %>% count

print(ggplot(svKnownDataLocs %>% filter(SampleCount>100) %>% spread(Type,n,fill=0),
             aes(x=DistBucket))
      + geom_line(aes(y=BND,color='BND'))
      + geom_line(aes(y=INV,color='INV'))
      + scale_y_log10()
      + xlim(-5000,5000)
      + facet_wrap(~LineLabel)
      + labs(title='SV breakend distance from end of Known Line Site'))


View(svKnownData)

View(lineSvData %>% filter(SampleId=='CPCT02040156T'&ClusterId==50))

tmpLinks = read.csv('~/logs/LNX_LINKS.csv')
View(tmpLinks %>% filter(ResolvedType=='LINE'))

linkLinks = links %>% filter(ResolvedType=='LINE')
View(linkLinks)
View(linkLinks %>% group_by(ChainCount,IsAssembled) %>% count %>% spread(IsAssembled,n,fill=0))
View(linkLinks %>% filter(ChainCount==4))

print(ggplot(linkLinks %>% filter(TILength<5e3) %>% group_by(IsAssembled,Length=2**round(log(TILength,2))) %>% count,aes(x=Length,y=n))
      + geom_line() 
      + scale_x_log10()
      + facet_wrap(~IsAssembled)
      + labs(title='LINE Source Location TI Lengths by IsAssembled',x='',y=''))

View(svData %>% filter(SampleId=='CPCT02330049T'&ClusterId==531&ChainId==0) %>% select(Id,ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd,LEStart,LEEnd,AsmbStart,AsmbEnd,LnkLenStart,LnkLenEnd,everything()))

View(lineSvData %>% filter(ChainCount==4) %>% 
       select(SampleId,ClusterId,ClusterCount,Id,ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd,LEStart,LEEnd,AsmbStart,AsmbEnd,LnkLenStart,LnkLenEnd,everything()))

lineElements = read.csv('~/data/line_elements.csv')
View(lineElements %>% mutate(Length=PosEnd-PosStart))
View(lineElements %>% mutate(Length=PosEnd-PosStart) %>% group_by(IsLong=Length<5000) %>% count)

## most common locations
lineLocSummary = lineLocations %>% group_by(LineId,Type,Chromosome,PosStart,PosEnd,SampleCount,PcawgSampleCount,Length) %>% 
  summarise(MaxSampleInserts=max(SampleInserts),
            TotalInserts=first(TotalInserts))
View(lineLocSummary)

View(rmLineData %>% filter(RmId==98249|RmId==98250))
View(rmLineData %>% filter(Chromosome==1&PosStart>199430000&PosEnd<199450000))

# 29058968	29073261

View(lineSvData %>% filter(SampleId=='CPCT02310016T'&ClusterId==118) %>% select(SampleId,ClusterId,Id,Type,ChrStart,ChrEnd,PosStart,PosEnd,LEStart,LEEnd,everything()))
View(lineSvData %>% filter(SampleId=='CPCT02330112T'&ClusterId==55) %>% select(SampleId,ClusterId,Id,Type,ChrStart,ChrEnd,PosStart,PosEnd,LEStart,LEEnd,everything()))
View(lineSvData %>% filter(SampleId=='CPCT02050163T'&ClusterId==319) %>% select(SampleId,ClusterId,Id,Type,ChrStart,ChrEnd,PosStart,PosEnd,LEStart,LEEnd,everything()))



tmpSVs = read.csv('~/logs/LNX_SVS.csv')
View(tmpSVs %>% filter(ResolvedType=='LINE'))

View(tmpSVs %>% filter(SampleId=='DRUP01010108T'&ClusterId==559) %>% select(SampleId,ClusterId,Id,Type,ChrStart,ChrEnd,PosStart,PosEnd,LEStart,LEEnd,everything()))
View(tmpSVs %>% filter(SampleId=='CPCT02330112T'&ClusterId==55) %>% select(SampleId,ClusterId,Id,Type,ChrStart,ChrEnd,PosStart,PosEnd,LEStart,LEEnd,everything()))
View(tmpSVs %>% filter(SampleId=='CPCT02010022T'&ClusterId==132) %>% select(SampleId,ClusterId,Id,Type,ChrStart,ChrEnd,PosStart,PosEnd,LEStart,LEEnd,everything()))

View(tmpSVs %>% filter(SampleId=='CPCT02010278T'&Id %in% c(650,651)) %>% select(SampleId,ClusterId,Id,Type,ChrStart,ChrEnd,PosStart,PosEnd,LEStart,LEEnd,everything()))


# CPCT02330112T
# 


View(lineSvData %>% filter(ClusterDesc=='SGL=2'&LEStart!='None'&LEEnd!='None') %>% group_by(SampleId,ClusterId) %>% summarise(Count=n(),
                                                                                                                              OrientSum=first(OrientStart)+last(OrientStart)))

tmp23 = lineLocations %>% filter(LineId==23)

tmp23 = merge(tmp23,lineSvData %>% select(SampleId,ClusterId,Id,Type,ChrStart,ChrEnd,PosStart,PosEnd,LEStart,LEEnd),by=c('SampleId','ClusterId'),all.x=T)
View(tmp23 %>% filter(ChrEnd==12))

View(lineLocations %>% filter(LineId==23))


View(links %>% filter(ResolvedType=='LINE'))
View(links %>% filter(ResolvedType=='LINE') %>% group_by(LinkReason,IsAssembled) %>% count)



pcawgLineData = read.csv('~/data/sv/cohort/pcawg_line_element_data.csv')
nrow(pcawgLineData)
colnames(pcawgLineData)

print(ggplot(lineLocations %>% filter(Length>0) %>% group_by(Type,Length=2**round(log(Length,2))) %>% count,aes(x=Length,y=n))
      + geom_line() 
      + scale_x_log10()
      + facet_wrap(~Type)
      + labs(title='LINE Source Location Width by Type',x='',y=''))


View(lineSvData %>% filter(ClusterCount==4) %>% select(SampleId,ClusterId,Id,Type,ChrStart,ChrEnd,PosStart,PosEnd,LEStart,LEEnd,everything()) %>% arrange(SampleId,ClusterId))


View(lineSvData %>% filter(SampleId=='CPCT02330112T') %>% select(SampleId,ClusterId,Id,Type,ChrStart,ChrEnd,PosStart,PosEnd,LEStart,LEEnd,everything()) %>% arrange(SampleId,ClusterId))
View(lineSvData %>% filter(SampleId=='CPCT02030318T') %>% select(SampleId,ClusterId,Id,Type,ChrStart,ChrEnd,PosStart,PosEnd,LEStart,LEEnd,everything()) %>% arrange(SampleId,ClusterId))
View(lineSvData %>% filter(SampleId=='CPCT02330112T') %>% select(SampleId,ClusterId,Id,Type,ChrStart,ChrEnd,PosStart,PosEnd,LEStart,LEEnd,everything()) %>% arrange(SampleId,ClusterId))



nrow(svData %>% filter(ResolvedType=='LINE'))

View(lineSvData %>% filter(LEStart!='None'&LEEnd!='None') %>% group_by(Type,LEStart,LEEnd) %>% count %>% spread(Type,n,fill=0))
View(lineSvData %>% filter(LEStart!='None'&LEEnd!='None') %>% group_by(Type,LEStart,LEEnd,IsShort=(Type!='BND'&PosEnd-PosStart<1e4)) %>% count %>% spread(Type,n,fill=0))

linePosData = rbind(lineSvData %>% filter(LEStart!='None') %>% select(SampleId,Id,ClusterId,ClusterCount,Chr=ChrStart,Pos=PosStart,LineType=LEStart),
                    lineSvData %>% filter(LEEnd!='None'&Type!='SGL'&Type!='INF') %>% select(SampleId,Id,ClusterId,ClusterCount,Chr=ChrEnd,Pos=PosEnd,LineType=LEEnd))

View(linePosData %>% filter(SampleId=='CPCT02030318T'&ClusterId==1601))

# width of clusters
bucketSize=1e5
lineClusters = linePosData %>% group_by(SampleId,ClusterId,ClusterCount,Chr,PosBucket=bucketSize*round(Pos/bucketSize)) %>% 
  summarise(PosSvCount=n(),PosLow=min(Pos),PosHigh=max(Pos)) %>% mutate(RegionLength=PosHigh-PosLow)
View(lineClusters)
nrow(lineClusters)
View(lineClusters %>% group_by(RegionLength=2**round(log(RegionLength,2))) %>% count)

View(svData %>% filter(ResolvedType=='LINE',SampleId =='CPCT02010022T'))

bucketSize=1e5

linePosData = merge(linePosData,samplesDD %>% select(SampleId,CancerType),by='SampleId',all.x=T)
View(linePosData %>% filter(LineType!='None') %>% group_by(SampleId,CancerType,Chr,PosBucket=bucketSize*round(Pos/bucketSize)) %>% count)
View(linePosData  %>% group_by(SampleId,LineType,CancerType,Chr,PosBucket=bucketSize*round(Pos/bucketSize)) %>% count %>% group_by(Chr,PosBucket,LineType) %>% count %>% spread(LineType,n) )
View(linePosData %>% filter(LineType!='None') %>% group_by(SampleId,CancerType,Chr,PosBucket=bucketSize*round(Pos/bucketSize)) %>% count %>% group_by(SampleId,CancerType) %>% count)
View(linePosData %>% filter(LineType!='None') %>% group_by(SampleId,CancerType,Chr,PosBucket=bucketSize*round(Pos/bucketSize),LineType) %>% count %>% group_by(SampleId,CancerType,LineType) %>% count %>% spread(LineType,n))
temp=(linePosData %>% filter(LineType!='None') %>% group_by(SampleId,CancerType,Chr,PosBucket=bucketSize*round(Pos/bucketSize),LineType) %>% count %>% group_by(SampleId,CancerType,LineType) %>% summarise(n=sum(n)) %>% spread(LineType,n))
ggplot(temp,aes(`Known;Suspect`,Suspect)) + geom_point()
bucketSize=1e5
linePosSummary = linePosData %>% group_by(Chr,PosBucket=bucketSize*round(Pos/bucketSize),
                                          LineType=ifelse(grepl('Known',LineType),'KNOWN',ifelse(grepl('Suspect',LineType),'SUSPECT','NONE'))) %>% count
View(linePosSummary %>% spread(LineType,n))
View(linePosSummary %>% group_by(LineType) %>% summarise(Total=sum(n)))
View(linePosSummary %>% group_by(Chr,PosBucket) %>% count)

print(ggplot(linePosSummary 
             # %>% filter(Chr %in% c(1,2,3,4))
             ,aes(x=PosBucket,y=n,fill=LineType))
      + geom_bar(stat='identity') 
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
      + scale_y_log10()
      + facet_wrap(~Chr)
      + labs(title='LINE Source and Insertion Locations',x='',y=''))

print(ggplot(linePosSummary %>% spread(LineType,n,fill=0),aes(x=PosBucket)
             + geom_point(aes()) 
             + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))
             + scale_y_log10()
             + facet_wrap(~Chr)
             + labs(title='LINE Source and Insertion Locations',x='',y=''))
      
      
      View(drivers %>% group_by(Gene,ResolvedType) %>% count %>% spread(ResolvedType,n,fill=0))
      
      
      
