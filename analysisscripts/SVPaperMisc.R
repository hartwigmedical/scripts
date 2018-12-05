svData = sv_load_and_prepare('~/data/sv/CLUSTER.csv')


# FOLDBACK LENGTHS and TYPES

foldbacksStart = svData %>% filter(FoldbackLenStart>=0)
foldbacksStart$FoldbackLength = foldbacksStart$FoldbackLenStart
foldbacksStart$FoldbackLinkInfo = foldbacksStart$FoldbackLinkInfoStart

foldbacksEnd = svData %>% filter(FoldbackLenEnd>=0)
foldbacksEnd$FoldbackLength = foldbacksEnd$FoldbackLenEnd
foldbacksEnd$FoldbackLinkInfo = foldbacksEnd$FoldbackLinkInfoEnd
foldbacks = rbind(foldbacksStart,foldbacksEnd)
View(foldbacks)

foldbacks$HasLinkInfo = ifelse(foldbacks$FoldbackLinkInfo!='-1;-1;-1'&grepl(';', foldbacks$FoldbackLinkInfo),T,F)
foldbacks$ChainLinks = 0
foldbacks$AssemblyLinks = 0
foldbacks$ChainLength = 0

foldbacksNoInfo = foldbacks %>% filter(!HasLinkInfo)
nrow(foldbacksNoInfo)

foldbacksInfo = foldbacks %>% filter(HasLinkInfo)
nrow(foldbacksInfo)
foldbacksInfo = get_foldback_chain_links(foldbacksInfo)
View(foldbacksInfo)

foldbacks = rbind(foldbacksNoInfo,foldbacksInfo)

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

foldbacks$ChainLinks = as.numeric(foldbacks$ChainLinks)
foldbacks$AssemblyLinks = as.numeric(foldbacks$AssemblyLinks)
foldbacks$ChainLength = as.numeric(foldbacks$ChainLength)

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


foldbackLenSummary = foldbacks %>% group_by(FoldbackLenBucket,FoldbackType) %>% summarise(Count=n()) %>% spread(FoldbackType,Count)
View(foldbackLenSummary)

# View(sampleClusterSummary %>% filter(AsmbLinkCount==ClusterCount-1&ClusterCount>=2))
# View(sampleClusterSummary %>% filter(AsmbLinkCount==ClusterCount-1&ClusterCount>=2&FoldbackCount>0))

print(ggplot(data = foldbackLenSummary, aes(x=FoldbackLenBucket, y=Count))
                     + geom_line(aes(y=INV, colour='INV'))
                     + geom_line(aes(y=Combo, colour='Combo'))
                     + scale_x_log10()
                     + labs(title = "Foldback Length Distribution"))

folbackLengthPlot = (ggplot(data = foldbackLenSummary, aes(x=FoldbackLenBucket, y=Count))
                     + geom_line(aes(y=INV, colour='INV'))
                     + geom_line(aes(y=AsmbCombo, colour='AsmbCombo'))
                     + geom_line(aes(y=OtherCombo, colour='OtherCombo'))
                     + scale_x_log10()
                     + labs(title = "Foldback Length Distribution"))

print(folbackLengthPlot)


# all combos
foldbacks$Category = paste(foldbacks$FoldbackType,"_CS=",foldbacks$ChainSize,"_ASMPerc=",foldbacks$FoldbackAsmbPercent,sep='')

print(ggplot(data = foldbacks %>% filter(FoldbackType!='INV') %>% group_by(Category,FoldbackLenBucket) %>% summarise(Count=n()), aes(x=FoldbackLenBucket, y=Count))
                     + geom_line()
                     + scale_x_log10()
                     + facet_wrap(~Category)
                     + labs(title = "Foldback Length Distribution"))

# limited by chain length
foldbacks$ChainLengthGroup = ifelse(foldbacks$ChainLength<=1e4,'ShortChain','LongChain')
print(ggplot(data = foldbacks %>% filter(FoldbackType!='INV'&ChainLength<=1e4) %>% group_by(Category,FoldbackLenBucket) %>% summarise(Count=n()), aes(x=FoldbackLenBucket, y=Count))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~Category)
      + labs(title = "Foldback Length Distribution"))

plotData = foldbacks %>% filter(FoldbackType!='INV') %>% group_by(Category,FoldbackLenBucket,ChainLengthGroup) %>% summarise(Count=n()) %>% spread(ChainLengthGroup,Count)
View(plotData)
print(ggplot(data = plotData, aes(x=FoldbackLenBucket, y=Count))
      + geom_line(aes(y=ShortChain,color='ShortChain'))
      + geom_line(aes(y=LongChain,color='LongChain'))
      + scale_x_log10()
      + facet_wrap(~Category)
      + labs(title = "Foldback Length Distribution"))


folbackLengthPlot = (ggplot(data = foldbackLenSummary, aes(x=FoldbackLenBucket, y=Count))
                     + geom_line(aes(y=INV, colour='INV'))
                     + geom_line(aes(y=Combo, colour='Combo'))
                     + scale_x_log10()
                     + labs(title = "Foldback Length Distribution"))

print(folbackLengthPlot)


# LINE ELEMENTS


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








