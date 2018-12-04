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

foldbacks$HasLinkInfo = ifelse(foldbacks$FoldbackLinkInfo!='-1;-1'&grepl(';', foldbacks$FoldbackLinkInfo),T,F)
foldbacks$ChainLinks = 0
foldbacks$AssemblyLinks = 0

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
    foldbacks[i,colCount-1] = as.numeric(stri_split_fixed(data$FoldbackLinkInfo,';')[[1]][1])
    foldbacks[i,colCount] =  as.numeric(stri_split_fixed(data$FoldbackLinkInfo,';')[[1]][2])
  }
  
  return (foldbacks)
}

foldbacks$ChainLinks = as.numeric(foldbacks$ChainLinks)
foldbacks$AssemblyLinks = as.numeric(foldbacks$AssemblyLinks)

foldbacks$FoldbackType = ifelse(!is.na(foldbacks$FoldbackLnkStart)&!is.na(foldbacks$FoldbackLnkEnd)&foldbacks$FoldbackLnkStart==foldbacks$FoldbackLnkEnd&foldbacks$Type=='INV','INV','Combo')
foldbacks$FoldbackLenBucket = 2**round(log(foldbacks$FoldbackLength,2))
foldbacks$FoldbackAsmbPercent = ifelse(foldbacks$ChainLinks>0,round(foldbacks$AssemblyLinks/foldbacks$ChainLinks/0.2)*0.2,0)
foldbacks$ChainSize = ifelse(foldbacks$ChainLinks==0,'None',ifelse(foldbacks$ChainLinks==1,'Single',ifelse(foldbacks$ChainLinks<=3,'Small','Long')))

View(foldbacks)
View(foldbacks %>% filter(ChainLinks==0&FoldbackType=='Combo'))
View(foldbacks %>% group_by(ChainSize) %>% count())
View(foldbacks %>% group_by(FoldbackAsmbPercent) %>% count())
View(foldbacks %>% filter(ChainLinks>3&AssemblyLinks>=0.5*ChainLinks))

foldbacks$Category = paste(foldbacks$FoldbackType,"_CS=",foldbacks$ChainSize,"_ASMPerc=",foldbacks$FoldbackAsmbPercent,sep='')

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


print(ggplot(data = foldbacks %>% filter(FoldbackType!='INV') %>% group_by(Category,FoldbackLenBucket) %>% summarise(Count=n()), aes(x=FoldbackLenBucket, y=Count))
                     + geom_line()
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












