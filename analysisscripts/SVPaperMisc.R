################################################
####### FOLDBACK ANALYSIS #######################
#################################################

get_foldback_chain_links<-function(foldbacks)
{
  colCount = ncol(foldbacks)  
  for(i in 1:nrow(foldbacks))
  {
    data = foldbacks[i,]
    foldbacks[i,colCount-2] = as.numeric(stri_split_fixed(data$FoldbackLinkInfo,';')[[1]][1])
    foldbacks[i,colCount-1] =  as.numeric(stri_split_fixed(data$FoldbackLinkInfo,';')[[1]][2])
    foldbacks[i,colCount] = as.numeric(stri_split_fixed(data$FoldbackLinkInfo,';')[[1]][3])
  }
  
  return (foldbacks)
}


svData = read.csv('~/Dropbox/HMF Australia team folder/Structural Variant Analysis/CLUSTER_GRIDSS.csv', header = T, stringsAsFactors = F)
svData = sv_load_and_prepare('~/data/sv/CLUSTER.csv')

# Interpetation of chaining for foldbacks
foldbacksStart = svData %>% filter(FoldbackLenStart>=0)
foldbacksStart$FoldbackLength = foldbacksStart$FoldbackLenStart
foldbacksStart$FoldbackLinkInfo = foldbacksStart$FoldbackLinkInfoStart
foldbacksEnd = svData %>% filter(FoldbackLenEnd>=0)
foldbacksEnd$FoldbackLength = foldbacksEnd$FoldbackLenEnd
foldbacksEnd$FoldbackLinkInfo = foldbacksEnd$FoldbackLinkInfoEnd
foldbacks = rbind(foldbacksStart,foldbacksEnd)

foldbacks$HasLinkInfo = ifelse(grepl(';', foldbacks$FoldbackLinkInfo),T,F)
foldbacks$ChainLinks = 0
foldbacks$AssemblyLinks = 0
foldbacks$ChainLength = 0

foldbacksNoInfo = foldbacks %>% filter(!HasLinkInfo)
foldbacksNoInfo = within(foldbacksNoInfo, rm(FoldbackLinkInfo))

foldbacksInfo = foldbacks %>% filter(HasLinkInfo)

foldbacksInfo = foldbacksInfo %>% separate(FoldbackLinkInfo,c('ChainLinks','AssemblyLinks','ChainLength'),sep = ';')
View(foldbacksInfo)
# foldbacksInfo = get_foldback_chain_links(foldbacksInfo)
foldbacks = rbind(foldbacksNoInfo,foldbacksInfo)


foldbacks$ChainLinks = as.numeric(foldbacks$ChainLinks)
foldbacks$AssemblyLinks = as.numeric(foldbacks$AssemblyLinks)
foldbacks$ChainLength = as.numeric(foldbacks$ChainLength)
foldbacks$AvgLinkLength = ifelse(foldbacks$ChainLength>0,foldbacks$ChainLength/foldbacks$ChainLinks,0)

foldbacks$FoldbackType = ifelse(!is.na(foldbacks$FoldbackLnkStart)&!is.na(foldbacks$FoldbackLnkEnd)&foldbacks$FoldbackLnkStart==foldbacks$FoldbackLnkEnd&foldbacks$Type=='INV','INV','Combo')
foldbacks$FoldbackLenBucket = 2**round(log(foldbacks$FoldbackLength,2))
foldbacks$FoldbackAsmbPercent = ifelse(foldbacks$ChainLinks>0,round(foldbacks$AssemblyLinks/foldbacks$ChainLinks/0.5)*0.5,0)
foldbacks$ChainSize = ifelse(foldbacks$ChainLinks==0,'None',ifelse(foldbacks$ChainLinks==1,'Single',ifelse(foldbacks$ChainLinks<=3,'Small','Long')))

View(foldbacks)
View(foldbacks %>% filter(ChainLinks==0&FoldbackType=='Combo'))
View(foldbacks %>% group_by(ChainSize) %>% count())
View(foldbacks %>% group_by(FoldbackAsmbPercent) %>% count())
View(foldbacks %>% filter(ChainLinks>3&AssemblyLinks>=0.5*ChainLinks))

View(foldbacks %>% group_by(FoldbackType, ChainSize, FoldbackAsmbPercent) %>% count())


#1.Foldback by Chain Count and Asm percent
View(foldbacks %>% group_by(ChainLinksBucket,FoldbackLenBucket,FoldbackAsmbPercent) %>% count() %>% spread(FoldbackAsmbPercent,n))

View(foldbacks %>% group_by(SampleId,ClusterId,ResolvedType) %>% summarise(FBcount=n()/2,
                                                                           ClusterCount=first(ClusterCount),CNMax=max(pmax(AdjCNStart,AdjCNEnd)),FBPloidyMax=max(pmax(AdjCNChgStart,AdjCNChgEnd,Ploidy))))

#2.Simple + foldback length distribution for mostly assembled combos
foldbackLenSummary = foldbacks %>% filter(FoldbackType=="INV"|FoldbackAsmbPercent>0.5) %>% group_by(FoldbackLenBucket,FoldbackType) %>% summarise(Count=n()) %>% spread(FoldbackType,Count)
View(foldbackLenSummary)
print(ggplot(data = foldbackLenSummary, aes(x=FoldbackLenBucket, y=Count))
                     + geom_line(aes(y=INV, colour='INV'))
                     + geom_line(aes(y=Combo, colour='Combo > 50% assembled'))
                     + scale_x_log10()
                                   + labs(title = "Foldback Length Distribution"))

print(folbackLengthPlot)


# all combos
foldbacks$Category = paste(foldbacks$FoldbackType,"_CS=",foldbacks$ChainSize,"_ASMPerc=",foldbacks$FoldbackAsmbPercent,sep='')

# limited by chain length
foldbacks$ChainLengthGroup = ifelse(foldbacks$ChainLength<=5e3,'ShortChain','LongChain')
plotData = foldbacks %>% filter(FoldbackType!='INV') %>% group_by(Category,FoldbackLenBucket,ChainLengthGroup) %>% summarise(Count=n()) %>% spread(ChainLengthGroup,Count)
print(ggplot(data = plotData, aes(x=FoldbackLenBucket, y=Count))
      + geom_line(aes(y=ShortChain,color='ShortChain'))
      + geom_line(aes(y=LongChain,color='LongChain'))
      + scale_x_log10()
      + facet_wrap(~Category)
      + labs(title = "Foldback Length Distribution"))


foldbacks$AvgLinkLenBucket = 2**round(log(foldbacks$AvgLinkLength,2))

View(foldbacks %>% filter(ChainLinks>0))


plotData = foldbacks %>% filter(FoldbackType!='INV'&FoldbackLength>=5e2&FoldbackLength<=8e3&ChainLength<=8e3) %>% group_by(AvgLinkLenBucket) %>% summarise(Count=n())
View(plotData)

print(ggplot(data = plotData, aes(x=AvgLinkLenBucket, y=Count))
      + geom_line()
      + scale_x_log10()
      + labs(title = "Foldback Chained TI Length Distribution"))



# DEL_Ext_TI with foldbacks

delExtTIFoldbacks = foldbacks %>% filter(ResolvedType=='DEL_Ext_TI')
View(delExtTIFoldbacks)
View(delExtTIFoldbacks %>% filter((OrientStart==-1&ArmStart=='P')|(OrientStart==1&ArmStart=='Q')))

invFoldbacks = svData %>% filter(Type=='INV'&FoldbackLenStart>=0&FoldbackLnkStart==FoldbackLnkEnd)

View(invFoldbacks %>% filter(ResolvedType=='DEL_Ext_TI'&((OrientStart==-1&ArmStart=='P')|(OrientStart==1&ArmStart=='Q'))))

View(invFoldbacks)
  ifelse(!is.na(foldbacks$FoldbackLnkStart)&!is.na(foldbacks$FoldbackLnkEnd)&foldbacks$FoldbackLnkStart==foldbacks$FoldbackLnkEnd&foldbacks$Type=='INV','INV','Combo')


#3. Anlaysis of combo foldbacks by ASM Percent
foldbackComboLenSummary = foldbacks %>% filter(FoldbackType=="Combo") %>% group_by(FoldbackLenBucket,FoldbackAsmbPercent,FoldbackType) %>% summarise(Count=n()) %>% spread(FoldbackType,Count)
print(ggplot(data = foldbackComboLenSummary, aes(x=FoldbackLenBucket, y=Count))
      + geom_line(aes(y=Combo, colour='Combo'))
      + scale_x_log10()
      + facet_wrap(~FoldbackAsmbPercent)
      + labs(title = "Combo Foldback Length Distribution by Assembled PCT"))                     

#4. Anlaysis of combo foldbacks by Chainlinks
foldbackComboLenSummary = foldbacks %>% filter(FoldbackType=="Combo") %>% group_by(FoldbackLenBucket,ChainLinksBucket,FoldbackType) %>% summarise(Count=n()) %>% spread(FoldbackType,Count)
print(ggplot(data = foldbackComboLenSummary, aes(x=FoldbackLenBucket, y=Count))
      + geom_line(aes(y=Combo, colour='Combo'))
      + scale_x_log10()
      + facet_wrap(~ChainLinksBucket)
      + labs(title = "Combo Foldback Length Distribution by ChainLinksBucket"))      




nrow(svData %>% filter(ResolvedType=='SglPair_INS')) 
nrow(svData %>% filter(ResolvedType=='SglPair_INS'&(grepl('AAAAAAAA',InsertSeq)|grepl('TTTTTTTT',InsertSeq)))) 

View(svData %>% filter(ResolvedType=='SglPair_DEL')) 
nrow(svData %>% filter(ResolvedType=='SglPair_DEL'&(grepl('AAAAAAAA',InsertSeq)|grepl('TTTTTTTT',InsertSeq)))) 

sampleClusterSummary = (svData %>% group_by(SampleId,ClusterId)
                        %>% summarise(ClusterCount=first(ClusterCount),
                                      IsResolved=first(IsResolved),
                                      ResolvedType=first(ResolvedType),
                                      ArmCount=first(ArmCount),
                                      DBCount=sum(DBCount),
                                      TICount=sum(TICount),
                                      ShortTICount=sum(ShortTICount),
                                      AsmbLinkCount=sum(AsmbTICount),
                                      InferLinkCount=sum(InferTICount),
                                      ChainedCount=sum(IsChained),
                                      RepeatedChainLinkCount=sum(RepeatedChainLink),
                                      FoldbackCount=sum(FoldbackCount),
                                      Consistent=first(IsConsistent),
                                      DupBECount=sum(DoubleDupBE),
                                      DelCount=sum(Type=='DEL'),
                                      DupCount=sum(Type=='DUP'),
                                      InsCount=sum(Type=='INS'),
                                      InvCount=sum(Type=='INV'),
                                      BndCount=sum(Type=='BND'),
                                      SglCount=sum(Type=='SGL'),
                                      NoneCount=sum(Type=='NONE'),
                                      LineCount=sum(IsLINE),
                                      KnownLineCount=sum(LEStart=='Known'|LEEnd=='Known'|LEStart=='Ident'|LEEnd=='Ident'),
                                      SuspectLineCount=sum(LEStart=='Suspect'|LEEnd=='Suspect'),
                                      PolyAorTCount=sum(grepl('AAAAAAAA',InsertSeq)|grepl('TTTTTTTT',InsertSeq)),
                                      FragileSiteCount=sum(IsFS)
                        )
                        %>% arrange(SampleId,ClusterId))


View(sampleClusterSummary %>% filter(ClusterCount==2&BndCount==0&LineCount>0))
View(sampleClusterSummary %>% filter(LineCount>0))
View(sampleClusterSummary %>% filter(ResolvedType=='Line'))
View(sampleClusterSummary %>% filter(LineCount>0&BndCount>0))
nrow(sampleClusterSummary %>% filter(ResolvedType=='Line'&KnownLineCount<ClusterCount&SuspectLineCount==0&PolyAorTCount==0))








