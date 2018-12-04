svData = read.csv('~/Dropbox/HMF Australia team folder/Structural Variant Analysis/CLUSTER_GRIDSS.csv', header = T, stringsAsFactors = F)

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

##############################

# Interpetation of chaining for foldbacks
foldbacksStart = svData %>% filter(FoldbackLenStart>=0)
foldbacksStart$FoldbackLength = foldbacksStart$FoldbackLenStart
foldbacksStart$FoldbackLinkInfo = foldbacksStart$FoldbackLinkInfoStart
foldbacksEnd = svData %>% filter(FoldbackLenEnd>=0)
foldbacksEnd$FoldbackLength = foldbacksEnd$FoldbackLenEnd
foldbacksEnd$FoldbackLinkInfo = foldbacksEnd$FoldbackLinkInfoEnd
foldbacks = rbind(foldbacksStart,foldbacksEnd)

foldbacks$HasLinkInfo = ifelse(foldbacks$FoldbackLinkInfo!='-1;-1'&grepl(';', foldbacks$FoldbackLinkInfo),T,F)
foldbacks$ChainLinks = 0
foldbacks$AssemblyLinks = 0

foldbacksNoInfo = foldbacks %>% filter(!HasLinkInfo)
foldbacksInfo = foldbacks %>% filter(HasLinkInfo)
foldbacksInfo = get_foldback_chain_links(foldbacksInfo)
foldbacks = rbind(foldbacksNoInfo,foldbacksInfo)

# ANNOTATIONS
foldbacks$FoldbackType = ifelse(!is.na(foldbacks$FoldbackLnkStart)&!is.na(foldbacks$FoldbackLnkEnd)&foldbacks$FoldbackLnkStart==foldbacks$FoldbackLnkEnd&foldbacks$Type=='INV','INV','Combo')
foldbacks$FoldbackLenBucket = 2**round(log(foldbacks$FoldbackLength,2))
foldbacks$ChainLinksBucket = 2**round(log(foldbacks$ChainLinks,2))
foldbacks$FoldbackAsmbPercent = ifelse(foldbacks$ChainLinks>0,round(foldbacks$AssemblyLinks/foldbacks$ChainLinks/0.2)*0.2,0)
foldbacks$ChainSize = ifelse(foldbacks$ChainLinks==0,'None',ifelse(foldbacks$ChainLinks==1,'Single',ifelse(foldbacks$ChainLinks<=3,'Small','Long')))
foldbacks$Category = paste(foldbacks$FoldbackType,"_CS=",foldbacks$ChainSize,"_ASMPerc=",foldbacks$FoldbackAsmbPercent,sep='')

#### TO DO: Would be good to know a total chain length

################################################
####### FOLDBACK ANALYSIS #######################
#################################################

#1.Foldback by Chain Count and Asm percent
View(foldbacks %>% group_by(ChainLinksBucket,FoldbackLenBucket,FoldbackAsmbPercent) %>% count() %>% spread(FoldbackAsmbPercent,n))


#2.Simple + foldback length distribution for mostly assembled combos
foldbackLenSummary = foldbacks %>% filter(FoldbackType=="INV"|FoldbackAsmbPercent>0.5) %>% group_by(FoldbackLenBucket,FoldbackType) %>% summarise(Count=n()) %>% spread(FoldbackType,Count)
print(ggplot(data = foldbackLenSummary, aes(x=FoldbackLenBucket, y=Count))
                     + geom_line(aes(y=INV, colour='INV'))
                     + geom_line(aes(y=Combo, colour='Combo > 50% assembled'))
                     + scale_x_log10()
                     + labs(title = "Foldback Length Distribution"))

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














