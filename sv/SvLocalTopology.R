
svData = read.csv('~/logs/SVA_SVS.csv')
clusters = read.csv('~/logs/SVA_CLUSTERS.csv')

#####
# Validation of topology types


# gather up SV data by local topology

get_local_toplogy<-function(data)
{
localTopStart = data %>% 
  select(SampleId,ClusterId,ClusterCount,ResolvedType,Id,Type,Chromosome=ChrStart,Position=PosStart,Orient=OrientStart,
         LocTopId=LocTopIdStart,LocTopType=LocTopTypeStart,LocTopTICount=LocTopTIStart,
         FbLink=FoldbackLnkStart,FbLength=FoldbackLenStart,DBLength=DBLenStart,TILength=LnkLenStart,Assembled=AsmbMatchStart)

localTopEnd = data %>% filter(Type!='INF'&Type!='SGL') %>% 
  select(SampleId,ClusterId,ClusterCount,ResolvedType,Id,Type,Chromosome=ChrEnd,Position=PosEnd,Orient=OrientEnd,
         LocTopId=LocTopIdEnd,LocTopType=LocTopTypeEnd,LocTopTICount=LocTopTIEnd,
         FbLink=FoldbackLnkEnd,FbLength=FoldbackLenEnd,DBLength=DBLenEnd,TILength=LnkLenEnd,Assembled=AsmbMatchEnd)


ltData = rbind(localTopStart,localTopEnd)
return (ltData)
}

# example usage
View(localTopData %>% filter(ClusterCount>1) %>% arrange(SampleId,ClusterId,LocTopId))

localTopSummary = localTopData %>% group_by(SampleId,ClusterId,ClusterCount,LocTopId,LocTopType,LocTopTICount) %>% 
  group_by(SampleId,LocTopType) %>% count()



################
# Analysis of TI and DB only length distributions
tiDsbData = svData %>% filter(ClusterCount>1&ResolvedType!='SIMPLE'&ResolvedType!='LINE'&Type!='SGL'&Type!='INF') %>%
  filter(LocTopTypeStart=='TI_ONLY'|LocTopTypeEnd=='TI_ONLY'|(LocTopTypeStart=='DSB'&LocTopTIStart==0)|(LocTopTypeEnd=='DSB'&LocTopTIEnd==0))

# clean-up
rm(tiDsbLtPosData)
rm(tiDsbLtData)
rm(tiDsbData)

nrow(tiDsbData)

tiDsbLtData = get_local_toplogy(tiDsbData)
View(tiDsbLtData)

View(tiDsbLtData %>% filter(LocTopType=='TI_ONLY'|LocTopType=='DSB') %>% group_by(SampleId,ClusterId,LocTopId,LocTopType) %>% count())
View(tiDsbLtData %>% filter(LocTopType=='DSB') %>% group_by(SampleId,ClusterId,LocTopId,LocTopType) %>% count() %>% group_by(n) %>% count())

tiDsbLtPosData = tiDsbLtData %>% filter(LocTopType=='TI_ONLY'|(LocTopType=='DSB'&LocTopTICount==0)) %>% 
  group_by(SampleId,ClusterId,LocTopId,LocTopType) %>% 
  summarise(BreakendCount=n(),
            ResolvedType=first(ResolvedType),
            Chromosome=first(Chromosome),
            PosStart=ifelse(first(Position)<=last(Position),first(Position),last(Position)),
            PosEnd=ifelse(first(Position)>last(Position),first(Position),last(Position)),
            OrientStart=ifelse(first(Position)<=last(Position),first(Orient),last(Orient)),
            OrientEnd=ifelse(first(Position)>last(Position),first(Orient),last(Orient)),
            DBLenFirst=first(DBLength),
            DBLenLast=last(DBLength),
            TILenFirst=first(TILength),
            TILenLast=last(TILength),
            AsmbStart=first(Assembled=='MATCH'),
            AsmbEnd=last(Assembled=='MATCH')) %>%
  filter(BreakendCount==2) 

tiDsbLtPosData = tiDsbLtPosData %>%
  mutate(Length=ifelse(LocTopType=='DSB'&OrientStart==-1,PosStart-PosEnd,PosEnd-PosStart),
         LengthBucket=round(Length/5)*5)


# summary length stats
tiDsbStats = tiDsbLtPosData %>% group_by(LocTopType,Assembled=AsmbStart,LengthBucket) %>% count() %>% ungroup()
tiDsbStats = tiDsbStats %>% mutate(Category=ifelse(LocTopType=='DSB','Double-Stranded Break',ifelse(Assembled,'Assembled Short TI','Inferred Short TI')),
                                   LengthBucket=ifelse(LocTopType=='DSB',-LengthBucket,LengthBucket))

View(tiDsbStats)

print(ggplot(tiDsbStats %>% filter(LengthBucket>20&LengthBucket<2e3), aes(x=LengthBucket, y=n, fill=Category))
      + geom_bar(stat = "identity", colour = "black")
      + labs(x='Length Bucket', y='Count')
      #+ geom_bar(position = "fill",stat = "identity", colour = "black")
      #+ scale_y_continuous(labels = scales::percent_format())
      + scale_fill_manual(values = c('blue','red','green')))


# show short DUPs and DSBs together
shortDupData = shortDups %>% group_by(LengthBucket) %>% count() %>% mutate(LocTopType='DUP')
View(shortDupData)

comData = bind_rows(shortDupData %>% select(LocTopType,LengthBucket,n),
                tiDsbStats %>% filter(LengthBucket>0) %>% select(LocTopType,LengthBucket,n)) %>% 
  filter(LocTopType=='DSB'|LocTopType=='DUP')

View(comData %>% group_by(LocTopType) %>% count())

print(ggplot(comData, aes(x=LengthBucket, y=n, fill=LocTopType))
      + geom_bar(stat = "identity", colour = "black")
      + labs(x='Length Bucket', y='Count')
      + scale_fill_manual(values = c('blue','red')))


comDataSummary = comData %>% filter(LengthBucket<=1e3) %>% group_by(LengthBucket,LocTopType) %>% summarise(n=sum(n)) %>% spread(LocTopType,n)
comDataSummary[is.na(comDataSummary)] <- 0
View(comDataSummary)

print(ggplot(comDataSummary, aes(x=LengthBucket, y=n))
      + geom_line(aes(y=DUP,color='DUP'))
      + geom_line(aes(y=DSB,color='DSB'))
      + labs(x='Length Bucket', y='Count')
      + scale_fill_manual(values = c('blue','red')))


View(tiDsbLtPosData)
View(head(tiDsbLtPosData,200))
nrow(tiDsbLtPosData %>% filter(PosStart>PosEnd))
nrow(tiDsbLtPosData %>% filter(AsmbStart!=AsmbEnd))
View(tiDsbLtPosData %>% filter(LocTopType=='DSB') %>% group_by(Facing=OrientStart==-1&OrientEnd==1) %>% count())
View(tiDsbLtPosData %>% filter(LocTopType=='TI_ONLY') %>% group_by(Facing=OrientStart==-1&OrientEnd==1) %>% count())

View(tiDsbLtPosData %>% filter(LocTopType=='TI_ONLY') %>% group_by(TiLenDiff=TILenFirst-Length) %>% count())

View(tiDsbLtPosData %>% filter(LocTopType=='DSB'&ResolvedType=='COMPLEX'&Length<4500))
View(tiDsbLtPosData %>% filter(LocTopType=='DSB') %>% group_by(ResolvedType,NoDB=DBLenFirst==-1000) %>% count())
View(tiDsbLtPosData %>% filter(LocTopType=='DSB') %>% group_by(UnequDB=DBLenFirst!=DBLenLast) %>% count())
View(tiDsbLtPosData %>% filter(LocTopType=='DSB') %>% filter(DBLenFirst!=DBLenLast))


# validation
view_cluster_sv('CPCT02010003T',2)


# TI only data
tiOnly = svData %>% filter(LocTopTypeStart=='TI_ONLY'|LocTopTypeEnd=='TI_ONLY')
nrow(tiOnly)

tiLtData = get_local_toplogy(tiOnly)
View(tiLtData)

tiLtPosData = tiLtData %>% filter(LocTopType=='TI_ONLY') %>% group_by(SampleId,ClusterId,LocTopId) %>% 
       summarise(ResolvedType=first(ResolvedType),
                 Chromosome=first(Chromosome),
                 PosStart=min(Position),
                 PosEnd=max(Position)) %>%
       mutate(TILength=PosEnd-PosStart,TIMidPoint=(PosStart+PosEnd)/2)

write.csv(tiLtPosData,'~/logs/short_ti_positions.csv', row.names = F, quote = F)

view_cluster_sv('COLO829T',63)



# TI_ONLY topologies with a DB less than 5K
View(svData %>% filter(LocTopTypeEnd=='TI_ONLY'&DBLenEnd>-1000&DBLenEnd<5e3))

view_chromosome_sv('CPCT02220026T',3,3)
view_cluster_sv('DRUP01010092T',212)




View(localTopData %>% group_by(SampleId,ClusterId,LocTopType) %>% group_by(SampleId,LocTopType) %>% count())
