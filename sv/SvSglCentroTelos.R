#########
## SGLs going to Centromeres or Telomeres


tcColours = c('powderblue','steelblue2','steelblue3','steelblue','steelblue4','blue')

svData = read.csv('~/data/sv/drivers/LNX_SVS.csv')

sampleCancerTypes = load_cancer_types('~/data/hpc_sample_cancer_types.csv',F,10)

centroLengths = read.csv('~/data/sv/chr_cento_lengths.csv')


## Centromeric SGLs
sglCentros = svData %>% filter(Type=='SGL'&(RepeatClass=='Satellite/centr'|RepeatType=='HSATII'))  
View(sglCentros)

sglCentros = sglCentros %>% mutate(CentroType=ifelse(RepeatClass=='Satellite/centr','Centro','PeriC'),
                                   InsertSeqLength=stri_length(InsertSeq),
                                   SampleClusterId=paste(SampleId,ClusterId,sep='_'))

View(sglCentros %>% filter(CentroType=='Centro') %>% select(SampleId,Id,OrientStart,InsertSeqLength,InsertSeq,everything()))

write.csv(sglCentros %>% filter(CentroType=='Centro') %>% select(SampleId,Id,InsertSeq) %>% arrange(SampleId,Id),
          '~/logs/sgl_centro_ins_seq.csv',row.names = F, quote = F)

View(sglCentros %>% filter(CentroType=='Centro') %>% group_by(OrientStart,ArmStart,InsertSeq=10*round(InsertSeqLength/10)) %>% count)
nrow(sglCentros %>% filter(CentroType=='Centro'))


print(ggplot(data = sglCentros %>% filter(CentroType=='Centro') %>% 
               mutate(Context=paste(ArmStart,OrientStart,sep='_'),InsertSeqLength=10*round(InsertSeqLength/10)) %>% 
               group_by(Context,InsertSeqLength) %>% count, 
             aes(x=InsertSeqLength,y=n))
      + geom_line()
      + facet_wrap(~Context))


# cluster-2s
View(sglCentros %>% filter(ClusterCount==2) %>% group_by(ResolvedType) %>% count)
View(sglCentros %>% group_by(ClusterCount=pmin(ClusterCount,10)) %>% count)


# proximity to centromere
View(sglCentros %>% group_by(ChrStart,round(CentroDistPerc,1)) %>% count %>% spread(ChrStart,n))
View(sglCentros %>% group_by(ChrStart,round(CentroDistPerc,1)) %>% count %>% spread(ChrStart,n))
View(sglCentros %>% group_by(ChrStart,CentroDistBucket) %>% count %>% spread(ChrStart,n))
View(sglCentros %>% group_by(OrientStart,CentroDistBucket) %>% count %>% spread(OrientStart,n))
View(sglCentros %>% group_by(CentroType,ChrStart,CentroDistBucket) %>% count %>% spread(ChrStart,n))
View(sglCentros %>% group_by(SglCentroCount=2**round(log(SglCentroCount,2)),round(CentroDistPerc,1)) %>% count %>% spread(SglCentroCount,n))

# matching with centromeric gain/loss

View(cnChrData)

sampleCentroChanges = cnChrData %>% filter(!Chromosome %in% c(13,14,15,21,22,'X','Y')) %>% 
  group_by(SampleId,CNChgBucket=round(CentroCNChg)) %>% summarise(CentroChgCount=n())

View(sampleCentroChanges)

sampleSglPloidies = sglCentros %>% filter(CentroType=='Centro') %>% group_by(SampleId,CNChgBucket=round(Ploidy)) %>% summarise(SglChgCount=n())

# only take 1 SGL per cluster
sampleSglPloidies = sglCentros %>% filter(CentroType=='Centro') %>% group_by(SampleId,ClusterId) %>% summarise(CNChgBucket=round(mean(Ploidy))) %>%
  group_by(SampleId,CNChgBucket) %>% summarise(SglChgCount=n())

# only take 1 SGL per cluster and eliminate opposing SGLs on the same arm
sampleSglPloidies = sglCentros %>% filter(CentroType=='Centro') %>% group_by(SampleId,ClusterId,ChrStart,ArmStart) %>% 
  summarise(ArmSglCount=n(),
            ArmSglNetOrient=sum(OrientStart),
            CNChgBucket=round(mean(Ploidy))) %>% filter(ArmSglNetOrient!=0&ArmSglCount==1) %>%
  group_by(SampleId,ClusterId) %>% summarise(CNChgBucket=first(CNChgBucket)) %>%
  group_by(SampleId,CNChgBucket) %>% summarise(SglChgCount=n())

sampleSglPloidies = sampleSglPloidies%>% ungroup() %>% mutate(CNChgBucket=pmin(CNChgBucket,5))
View(sampleSglPloidies)

sampleCentroChanges = sampleCentroChanges%>% ungroup() %>% mutate(CNChgBucket=pmin(CNChgBucket,5))

sampleSglCentroData = merge(sampleCentroChanges,sampleSglPloidies,by.x=c('SampleId','CNChgBucket'),all=T)
sampleSglCentroData[is.na(sampleSglCentroData)]=0
View(sampleSglCentroData)

print(ggplot(data = sampleSglCentroData %>% filter(CNChgBucket>0&SglChgCount>0), aes(x=CentroChgCount,y=SglChgCount))
      + geom_point(position="jitter")
      + geom_smooth(,method=lm,se=FALSE, fullrange=F)
      + facet_wrap(~CNChgBucket)
      + labs(x='Centromere CN Changes', y='Centromeric SGLs', title = "Faceted by Size of Copy Number Change"))



write.csv(sglCentros %>% filter(CentroType=='Centro') %>% mutate(SampleSvId=paste(SampleId,Id,sep='_')) %>% select(SampleSvId,InsertSeq),
          '~/logs/sgl_centro_ins_seq.fa',row.names=F, quote=F)


#####
## BLAT results on insert sequences from HG38

load_blat_results<-function(blatFile)
{
  blatRawData = read.csv(blatFile)

  filteredData = blatRawData %>% filter(!grepl('chrUn',TName)&!grepl('_random',TName)&!grepl('_alt',TName)&!grepl('HLA',TName))
  filteredData = filteredData %>% filter(QName!='')
    
  filteredData = filteredData %>% mutate(QSizeMult=pmax(9-floor(QSize/100),1),
                                                       Score=round((1000-(QSizeMult*Mismatch+QGapCount+TGapCount))*pmin(Match/QSize,1)))
  
  # cull uninteresting fields
  filteredData = filteredData %>% select(QName,TName,QSize,Score,Strand,TStart,TEnd,Match,Mismatch,QGapCount,TGapCount,QSizeMult)
  
  # only take the highest score for each chromosome
  topResultByChr = filteredData %>% arrange(QName,TName,-Score) %>% group_by(QName,TName) %>% 
         summarise(QSize=first(QSize),TopScore=first(Score),SecondScore=max(ifelse(Score==TopScore,0,Score)),
                   Strand=first(Strand),TStart=first(TStart),TEnd=first(TEnd),Match=first(Match),
                   Mismatch=first(Mismatch),QGapCount=first(QGapCount),TGapCount=first(TGapCount)) %>% ungroup()
  
  # rename TopScore to Score
  topResultByChr = topResultByChr %>% mutate(Score=TopScore) %>% select(-TopScore)
  
  return (topResultByChr)
}

tmp = load_blat_results('~/dev/blat/output.csv')
View(tmp)

filteredBlatResults_00 = load_blat_results('~/data/sv/sgl_centro_results_00.csv')
filteredBlatResults_01 = load_blat_results('~/data/sv/sgl_centro_results_01.csv')
# nrow(filteredBlatResults_01)
# View(filteredBlatResults_01)
# some of the first results overlapped with 01
filteredBlatResults_00 = filteredBlatResults_00 %>% filter(!(QName %in% filteredBlatResults_01$QName))
filteredBlatResults = rbind(filteredBlatResults_00,filteredBlatResults_01)
filteredBlatResults = rbind(filteredBlatResults,load_blat_results('~/data/sv/sgl_centro_results_02.csv'))
filteredBlatResults = rbind(filteredBlatResults,load_blat_results('~/data/sv/sgl_centro_results_03.csv'))
filteredBlatResults = rbind(filteredBlatResults,load_blat_results('~/data/sv/sgl_centro_results_04.csv'))
filteredBlatResults = rbind(filteredBlatResults,load_blat_results('~/data/sv/sgl_centro_results_05.csv'))
filteredBlatResults = rbind(filteredBlatResults,load_blat_results('~/data/sv/sgl_centro_results_06.csv'))
filteredBlatResults = rbind(filteredBlatResults,load_blat_results('~/data/sv/sgl_centro_results_07.csv'))
filteredBlatResults = rbind(filteredBlatResults,load_blat_results('~/data/sv/sgl_centro_results_08.csv'))
filteredBlatResults = rbind(filteredBlatResults,load_blat_results('~/data/sv/sgl_centro_results_09.csv'))
filteredBlatResults = rbind(filteredBlatResults,load_blat_results('~/data/sv/sgl_centro_results_10.csv'))
nrow(filteredBlatResults) # 220K

# check for duplicates
nrow(filteredBlatResults %>% group_by(QName,TName) %>% count %>% filter(n>1))
write.csv(filteredBlatResults,'~/data/sv/sgl_centro_blat_results.csv',row.names = F, quote = F)

View(filteredBlatResults)

# unmapped SGLs

unmapBlatResults = load_blat_results('~/logs/sgl_unmap_results_01.csv')
unmapBlatResults = rbind(unmapBlatResults,load_blat_results('~/logs/sgl_unmap_results_02.csv'))
unmapBlatResults = rbind(unmapBlatResults,load_blat_results('~/logs/sgl_unmap_results_03.csv'))
unmapBlatResults = rbind(unmapBlatResults,load_blat_results('~/logs/sgl_unmap_results_04.csv'))
unmapBlatResults = rbind(unmapBlatResults,load_blat_results('~/logs/sgl_unmap_results_05.csv'))
View(unmapBlatResults)
nrow(unmapBlatResults)
View(unmapBlatResults %>% group_by(TName,Score=round(Score,-2)) %>% count %>% spread(Score,n))


   
sglCentroData = filteredBlatResults %>% arrange(QName,-Score) %>% group_by(QName) %>% 
  summarise(CentroChr=first(TName),CentroStrand=first(Strand),CentroPosStart=first(TStart),CentroPosEnd=first(TEnd),
            TopScore=first(Score),NextScore=first(SecondScore),NextChrScore=max(ifelse(Score==TopScore,0,Score)),
            InsSeqLength=first(QSize),Match=first(Match),Mismatch=first(Mismatch)) %>% ungroup()

View(sglCentroData)

sglCentroData = sglCentroData %>% mutate(HighConf=TopScore>900&(TopScore-NextChrScore>25))
sglCentroData = sglCentroData %>% separate(QName,c('SampleId','Id'),sep='_')
sglCentroData$CentroChr = stri_replace_all_fixed(sglCentroData$CentroChr,'chr','')

sglCentroData = sglCentroData %>% filter(SampleId %in% sampleCancerTypes$SampleId)

nrow(sglCentros)
View(sglCentros %>% group_by(CentroType) %>% count)
View(sglCentroData %>% group_by(HighConf) %>% count)

# T Gap bases - gaps in insert sequence
diffData = sglCentroData %>% filter(HighConf&TopScore<1000) %>% mutate(CentroPosLength=CentroPosEnd-CentroPosStart,QDiff=CentroPosLength-InsSeqLength)
View(diffData %>% filter(QDiff>10) %>% group_by(QDiff=2**round(log(QDiff,2))) %>% count)
nrow(sglCentroData %>% filter(HighConf))
nrow(diffData %>% filter(QDiff>10))


# merge with SGL SV data and centromere change data
View(sglCentroData)
sglCentroSvData = merge(sglCentroData,
                        sglCentros %>% select(SampleId,Id,ChrStart,PosStart,OrientStart,ArmStart,ClusterId,ClusterCount,ResolvedType,
                                              Ploidy,CNChgStart,CNStart,GeneStart),
                        by=c('SampleId','Id'),all.x=T)

View(sglCentroSvData)

# merge in centromere CN changes
sglCentroSvData = merge(sglCentroSvData,cnChrData %>% select(SampleId,Chromosome,CentroCnP,CentroCnQ,PCentroGain,QCentroGain,CentroCNChg),
                        by.x=c('SampleId','CentroChr'),by.y=c('SampleId','Chromosome'),all.x=T)


# destination of centromeric SGLs
sglDestChrs = sglCentroSvData %>% group_by(Chromosome=CentroChr) %>% count
sglDestChrs = merge(sglDestChrs,chrIndex,by='Chromosome',all.x=T)

print(ggplot(sglDestChrs, aes(x=reorder(Chromosome,ChrIndex), y=n))
      + geom_bar(stat = "identity", colour = "black")
      + labs(x='Chromosome', y='# Centromeric SGLs', title='Desintation for Centromeric SGLs'))


# cohort plot

View(sglCentroSvData)
View(sglCentroSvData %>% filter(is.na(QCentroGain)))
View(cnChrData %>% filter(Chromosome==1) %>% group_by(QCentroGain,PCentroGain,CentroCNChg=pmin(round(CentroCNChg),4)) %>% count %>% spread(CentroCNChg,n,fill=0))

View(cnChrData %>% filter(SampleId=='CPCT02011069T'))

# orientations
View(sglCentroSvData %>% group_by(ArmStart,OrientStart,OrientEnd) %>% count %>% spread(ArmStart,n))
View(sglCentroSvData %>% filter(HighConf&ChrStart==CentroChr) %>%group_by(OrientStart,CentroStrand) %>% count %>% spread(CentroStrand,n))

# inferred SV types
sglCentroSvData = sglCentroSvData %>% mutate(OrientEnd=ifelse((OrientStart==1&CentroStrand=='+')|(OrientStart==-1&CentroStrand=='-'),-1,1),
                                             PosEnd=ifelse(OrientEnd==-1,CentroPosStart,CentroPosEnd),
                                             SglType=ifelse(ChrStart!=CentroChr,'BND',
                                                      ifelse(OrientStart==OrientEnd,'INV',
                                                      ifelse(OrientStart==1&OrientEnd==-1,'DEL','DUP'))),
                                             SglLength=ifelse(SglType=='BND',0,abs(PosStart-(CentroPosStart+CentroPosEnd)*0.5)))

View(sglCentroSvData %>% group_by(SglType,LengthBucket=ifelse(SglType=='BND',0,2**round(log(SglLength,2)))) %>% count %>% spread(SglType,n))
View(sglCentroSvData %>% group_by(SglType,HighConf) %>% count %>% spread(HighConf,n))

# proportion of local vs translocation SGLs
View(sglCentroSvData %>% filter(HighConf) %>% group_by(Local=SglType!='BND') %>% count)
nrow(sglCentroSvData %>% filter(HighConf&SglType!='BND'))/nrow(sglCentroSvData %>% filter(HighConf))

print(ggplot(data = sglCentroSvData %>% filter(HighConf&SglType!='BND') %>% 
               group_by(SglType,LengthBucket=ifelse(SglType=='BND',0,2**round(log(SglLength,2)))) %>% count, 
             aes(x=LengthBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~SglType)
      + labs(title = "SGL Centromeric Lengths"))

View(sglCentroSvData %>% filter(HighConf) %>% group_by(SglType,ChrStart) %>% count %>% spread(ChrStart,n))

# clustering and types of 'simple' SVs
View(sglCentroSvData %>% group_by(ClusterCount=2**round(log(ClusterCount,2))) %>% count)

sglCentroSvData = merge(sglCentroSvData,sampleCancerTypes,by='SampleId',all.x=T)

# colnames(sglCentroSvData)

# SGLs going to the same centromere
sampleChrData = sglCentroSvData %>% filter(HighConf) %>% group_by(SampleId,CancerType,CentroChr) %>% 
  summarise(SglCount=n(),
            ClusterCount=n_distinct(ClusterId),
            SourceChrCount=n_distinct(ChrStart),
            NetOrient=sum(OrientEnd),
            NetPloidy=sum(Ploidy*OrientEnd),
            CentroChange=first(CentroCNChg),
            PCentroGain=first(PCentroGain),
            QCentroGain=first(QCentroGain),
            FirstSourceChr=first(ChrStart),
            FirstCentroOrient=first(OrientEnd))

View(sampleChrData)

# frequency of SVs originating on each chromosome and arm by whether a local or translocation SV
View(sampleChrData %>% filter(SglCount==1) %>% group_by(FirstSourceChr,CentroChr) %>% count %>% spread(CentroChr,n))

View(sglCentroSvData %>% filter(HighConf) %>% group_by(ChrStart,ArmStart,SvType=ifelse(SglType=='BND','Translocation','Local')) 
     %>% count %>% spread(SvType,n) %>% mutate(LocalPercent=round(Local/(Translocation+Local),3)))

data0 = sglCentroSvData %>% filter(HighConf) %>% group_by(Source=paste(ChrStart,ArmStart,sep='_')) %>% 
  summarise(Total=n(),
            Translocation=sum(SglType=='BND')/n(),
            Local=sum(SglType!='BND')/n()) %>% gather('SvType','Percent',3:4) %>% mutate(SourceText=sprintf('%s (%d)',Source,Total))

print(ggplot(data0, aes(x=reorder(SourceText,-Total), y=Percent, fill=SvType))
      + geom_bar(stat = "identity", colour = "black")
      + labs(x='SGL Chromosome', y='# Centromeric SGLs')
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10)))

# frequency of SVs ending on each chromosome and arm by whether a local or translocation SV
data1 = sglCentroSvData %>% filter(HighConf) %>% group_by(CentroChr) %>% 
  summarise(Total=n(),
            Translocation=sum(SglType=='BND')/n(),
            Local=sum(SglType!='BND')/n()) %>% gather('SvType','Percent',3:4) %>%
  mutate(CentroChrText=sprintf('%s (%d)',CentroChr,Total))

## cohort plot
print(ggplot(data1, aes(x=reorder(CentroChrText,-Total), y=Percent, fill=SvType))
      + geom_bar(stat = "identity", colour = "black")
      + labs(x='Centromere Chromosome (total SGLs)', y='% of Centromeric SGLs per Chromosome',title='Proportion of Local vs Translocation Centromeric SGLs by Chromosome')
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10)))

data2 = sglCentroSvData %>% filter(HighConf) %>% group_by(CentroChr) %>% 
  summarise(Total=n(),
            TransToP=sum(SglType=='BND'&OrientEnd==1)/n(),TransToQ=sum(SglType=='BND'&OrientEnd==-1)/n(),
            LocalToP=sum(SglType!='BND'&OrientEnd==1)/n(),LocalToQ=sum(SglType!='BND'&OrientEnd==-1)/n()) %>% 
  gather('SvType','Percent',3:6) %>%
  mutate(CentroChrText=sprintf('%s (%d)',CentroChr,Total))

print(ggplot(data2, aes(x=reorder(CentroChrText,-Total), y=Percent, fill=SvType))
      + geom_bar(stat = "identity", colour = "black")
      + labs(x='Centromere Chromosome', y='# Centromeric SGLs')
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10)))

# count of samples with links between chromosome pairs by cancer type
linksByChromomsome = sglCentroSvData %>% filter(HighConf) %>% group_by(CancerType,SampleId,SglChr=ChrStart,CentroChr) %>% count %>%
  group_by(CancerType,SglChr,CentroChr) %>% summarise(Count=n())

# View(sglCentroSvData %>% filter(CancerType=='Breast'&HighConf&ChrStart==8&CentroChr==8) %>% group_by(SampleId) %>% count)

View(linksByChromomsome)
# View(sglCentroSvData %>% filter(HighConf) %>% group_by(ChrStart,CentroChr) %>% count %>% spread(CentroChr,n))

# locations of SGLs ending on chr 1
View(sglCentroSvData %>% filter(CentroChr==1) %>% group_by(QCentroGain,Position=round(PosEnd,-5)) %>% count %>% spread(QCentroGain,n,fill=0))


# from point of view of chromosomal data
centroGainData = cnChrData %>% filter(PCentroGain|QCentroGain)

View(sampleChrData)
centroGainData = merge(centroGainData,sampleChrData %>% ungroup() %>% select(-PCentroGain,-QCentroGain,-CancerType),by.x=c('SampleId','Chromosome'),by.y=c('SampleId','CentroChr'),all.x=T)

View(centroGainData %>% group_by(Chromosome,QCentroGain.x,HasSgls=!is.na(SglCount)) %>% count %>% spread(HasSgls,n))
View(centroGainData %>% group_by(Chromosome,QCentroGain.x) %>% count)
View(centroGainData)
colnames(centroGainData)
View(sampleChrData %>% filter(CentroChr==1&SglCount==1) %>% group_by(SglPloidy=pmax(pmin(round(NetPloidy),3),-3),
                                                                     CentroCnChg=pmin(round(CentroChange),3)) %>% count)

##### 
## Paired SGLs to same centromere, possible DBs and TIs

pairedSgls = sglCentroSvData %>% filter(HighConf) %>% arrange(SampleId,CentroChr,PosEnd) %>% group_by(SampleId,CentroChr) %>% 
  summarise(SglCount=n(),
            ClusterCount=n_distinct(ClusterId),
            SourceChrCount=n_distinct(ChrStart),
            NetOrient=sum(OrientEnd),
            NetPloidy=sum(OrientEnd*Ploidy),
            NetScore=sum(TopScore),
            CentroLowerPos=first(PosEnd),
            CentroUpperPos=last(PosEnd),
            CentroLowerOrient=first(OrientEnd),
            CentroUpperOrient=last(OrientEnd)) %>% filter(SglCount==2) %>% 
  mutate(LinkLength=CentroUpperPos-CentroLowerPos,
         LinkType=ifelse(NetOrient!=0,'None',ifelse(CentroLowerOrient==-1&CentroUpperOrient==1&LinkLength>=30,'TI','DB')),
         PairType=ifelse(SourceChrCount==2,'Translocation','Local'))
            
View(pairedSgls %>% filter(NetOrient==0))  
View(pairedSgls %>% filter(NetOrient==0) %>% group_by(NetPloidy=pmax(pmin(round(NetPloidy),3),-3)) %>% count)  

print(ggplot(pairedSgls %>% filter(LinkLength<3e6&round(NetPloidy)==0&NetScore>=1990) %>% 
               group_by(LinkType,LengthBucket=2**round(log(LinkLength,2))) %>% count,
             aes(x=LengthBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~LinkType)
      + labs(title = "Centromeric Link Lengths"))

print(ggplot(pairedSgls %>% filter(LinkLength<3e6) %>% group_by(Context=paste(PairType,LinkType,sep='_'),
                                                                LengthBucket=2**round(log(LinkLength,2))) %>% count,
             aes(x=LengthBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~Context)
      + labs(title = "Centromeric Link Lengths"))


#####
## Enrichment of Chromosome 1 Q-Arm CN gain
# - bias of SGLs by orientation - ie are they heading out to towards the Q telomere
# - position distribution of SGL centromereic end-points
# - clustering and grouping of centromeric singles

# 1600 samples have chr 1 Q-gain, 1950 have neither, 220 have P-gain
View(cnChrData %>% filter(Chromosome==1) %>% group_by(PCentroGain,QCentroGain) %>% count)

View(sampleChrData %>% filter(CentroChr==1) %>% group_by(SglCount=pmin(SglCount,5),NeutralOrient=(NetOrient==0)) %>% count %>% 
       spread(NeutralOrient,n,fill=0))

View(sampleChrData %>% filter(CentroChr==1&QCentroGain) %>% group_by(SglCount=pmin(SglCount,5),
                                                                     NetPloidy=pmax(pmin(round(NetPloidy),3),-3)) %>% count) # NetOrient=pmax(pmin(round(NetOrient),3),-3)

# moderate excess of SGLs facing Q arm for samples with a Q-arm gain
View(sampleChrData %>% filter(CentroChr==1&QCentroGain) %>% group_by(NetPloidy=pmax(pmin(round(NetPloidy),3),-3)) %>% count)


#####
## Fisher exact test for co-occurene of SGLs and centromeric CN gain
centroGainSglData = merge(cnChrData %>% mutate(CentroGain=PCentroGain|QCentroGain) %>% select(SampleId,Chromosome,CentroGain,CentroCNChg,PCentroGain,QCentroGain,CentroCnP,CentroCnQ),
                          sampleChrData %>% ungroup() %>% mutate(Chromosome=CentroChr) %>% select(SampleId,Chromosome,SglCount,NetPloidy),
                          by=c('SampleId','Chromosome'),all=T)

centroGainSglData = centroGainSglData %>% filter(!(Chromosome %in% c(13,14,15,21,22)))
centroGainSglData = centroGainSglData %>% mutate(HasSGLs=!is.na(SglCount),
                                                 HasCentroGain=!is.na(CentroGain)&CentroGain)

View(centroGainSglData)
View(centroGainSglData %>% group_by(HasSGLs,HasCentroGain) %>% count)

results = data.frame(matrix(ncol = 10, nrow = 0))
colnames(results) = colnames(tmp)

for(chr in unique(centroGainSglData$Chromosome))
{
  print(sprintf('testing chromosome: %s', chr))
  
  testStr = sprintf('Chromosome_%s',chr)
  tmp = calc_fisher_et(centroGainSglData %>% filter(Chromosome==chr),
                 centroGainSglData %>% filter(Chromosome==chr&HasSGLs),
                 centroGainSglData %>% filter(Chromosome==chr&HasCentroGain),'HasSGLs','HasCentroGain',log=F,T,chr)
  
  results = rbind(results,tmp)
}

results = results %>% mutate(WithSglRate=round(Both/HasSGLs,3),
                             NoSglRate=round(NoHasSGLsWithHasCentroGain/(All-HasSGLs),3))

View(results %>% select(WithSglRate,NoSglRate,everything()))

resultsPlot = results %>% select(Chromosome=Test,WithSGLs=WithSglRate,WithoutSGLs=NoSglRate)
resultsPlot = resultsPlot %>% gather('RateOfCentromericGain','Rate',2:3)
resultsPlot = merge(resultsPlot,chrIndex,by='Chromosome',all.x=T)
View(resultsPlot)

## cohort plot
print(ggplot(resultsPlot, aes(x=reorder(Chromosome,ChrIndex),y=Rate,fill=RateOfCentromericGain))
      + geom_bar(position="dodge",stat='identity')
      + labs(x='Centromere Chromosome', y='% of Samples with Centromeric Gain',title='Rate of Centromeric CN Gain with and without Centromeric SGLs'))




# Location of local SGLs by chromosome
# cohort plot

View(centroLengths)
chrPosBuckets = data.frame(matrix(ncol=2,nrow=0))
colnames(chrPosBuckets) = c('Chromosome','PosBucket')

for(i in 1:nrow(centroLengths))
{
  data = centroLengths[i,]
  length = data$Length
  chr = data$Chromosome
  buckets = floor(length/1e6)
  
  print(sprintf('chr(%s) length(%d) buckets(%d)', chr, length, buckets))
  
  for(j in 1:buckets)
  {
    index = nrow(chrPosBuckets)+1
    chrPosBuckets[index,1] = as.character(chr)
    chrPosBuckets[index,2] = j
    
    # print(sprintf('chr(%s) index(%d) j(%d)', chr, index, j))
  }
}

View(chrPosBuckets)
View(chrPosBuckets %>% group_by(Chromosome) %>% count)

data4 = sglCentroSvData %>% filter(CentroChr==ChrStart) %>% mutate(PosBucket=round(PosStart/1e6)) %>% 
  group_by(ChrStart,PosBucket) %>% count

data5 = merge(chrPosBuckets,data4,by.x=c('Chromosome','PosBucket'),by.y=c('ChrStart','PosBucket'),all.x=T)
data5[is.na(data5)] = 0

View(data4)
data5 = data5 %>% mutate(SourcePosLabel=sprintf('%.0f MB',PosBucket))
View(data5)
View(data5 %>% filter(Chromosome==17))

print(ggplot(data5 %>% filter(Chromosome!=''), aes(x=reorder(SourcePosLabel,PosBucket), y=n))
      + geom_bar(stat = "identity", colour = "black")
      + facet_wrap(~Chromosome)
      + labs(x='Chromosomal Position', y='# SGLs', title = "Originating position of Local Centromeric SGLs by Chromosome")
      + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()))

colnames(sglCentroSvData)


#####
## Chromosome 17 investigations

#  550 have both SGLs and gain, 551 have gain without any SGLs, 128 have SGLs without gain
View(centroGainSglData %>% filter(Chromosome==17) %>% group_by(HasSGLs,HasCentroGain) %>% count)

colnames(sglCentroSvData)
colnames(sampleChrData)
both17 = sglCentroSvData %>% filter(HighConf&CentroChr==17&(PCentroGain|QCentroGain)) %>% filter(ResolvedType!='LOW_VAF'&ResolvedType!='DUP_BE')
sample17 = both17 %>% group_by(SampleId,CancerType,CentroChr) %>% 
  summarise(SglCount=n(),
            ClusterCount=n_distinct(ClusterId),
            SourceChrCount=n_distinct(ChrStart),
            NetOrient=sum(OrientEnd),
            NetPloidy=sum(Ploidy*OrientEnd),
            CentroChange=first(CentroCNChg),
            PCentroGain=first(PCentroGain),
            QCentroGain=first(QCentroGain),
            FirstSourceChr=first(ChrStart),
            FirstCentroOrient=first(OrientEnd))



sample17 = sample17 %>% mutate(CorrectOrient=(NetOrient>0&PCentroGain)|(NetOrient<0&QCentroGain),
                               ArmGained=ifelse(PCentroGain,'P','Q'),
                               NetPloidyRnd=pmin(pmax(round(NetPloidy),-3),3),
                               CentroChgRnd=pmin(pmax(round(CentroChange),-3),3))

# matching by orientation - correct 75% of the time
View(both17 %>% group_by(CorrectOrient=(OrientEnd==1&PCentroGain)|(OrientEnd==-1&QCentroGain),ArmGained=ifelse(PCentroGain,'P','Q')) %>% count)
View(sample17 %>% group_by(CorrectOrient=(NetOrient>0&PCentroGain)|(NetOrient<0&QCentroGain),ArmGained=ifelse(PCentroGain,'P','Q')) %>% count)
View(sample17 %>% filter(NetOrient!=0) %>% group_by(CorrectOrient,ArmGained) %>% count)

# matching by orientation and copy number - high concordance
View(sample17 %>% group_by(ArmGained,CorrectOrient,CnChgMatch=abs(NetPloidyRnd)==CentroChgRnd) %>% count %>% spread(ArmGained,n))
View(sample17 %>% filter(CorrectOrient) %>% group_by(ArmGained,NetPloidy=abs(NetPloidyRnd),CentroChg=CentroChgRnd) %>% count %>% spread(ArmGained,n))

# source of the SGLs - chromosome and position
View(both17 %>% group_by(SglType,SampleId) %>% count %>% spread(SglType,n)) # 540 BNDs, 2:1:1 ratio for others totally 450
View(both17 %>% group_by(SampleId) %>% count %>% group_by(n) %>% count)
View(both17 %>% filter(SampleId=='CPCT02070352T') %>% select(OrientEnd,everything()))
View(both17 %>% group_by(SglType) %>% count) # 540 BNDs, 2:1:1 ratio for others totally 450
View(both17 %>% group_by(ChrStart) %>% count) # chr 8 most common source
View(both17 %>% filter(ChrStart==8))
View(both17 %>% group_by(ChrStart,SglType) %>% count)


# centromere is 23763006, total length is 81M, Q arm is 2x as long as P (P=23,763,006 Q=57,432,204)
# TP53 is at 7,565,097-7,590,856

# not clear enrichment from remote locations
View(both17 %>% filter(ChrStart!=17) %>% group_by(ChrStart,SourcePos=round(PosStart,-6)) %>% count)

# lots of local activity around the centromere
 View(both17 %>% filter(ChrStart==17) %>% group_by(ChrStart,SourcePos=round(PosStart,-6)) %>% count)

print(ggplot(both17 %>% filter(ChrStart==17) %>% group_by(SglType,SourcePos=round(PosStart,-6)) %>% count, 
             aes(x=SourcePos, y=n))
      + geom_bar(stat = "identity", colour = "black")
      + facet_wrap(~SglType)
      + labs(title = "SGLs to chromosome 17"))

# clustering - mostly complex 760 (total SGLs=) vs 160 SGLs
View(both17 %>% group_by(SampleId,ClusterId,Id,ResolvedType) %>% summarise(SGLs=n()) %>% 
       group_by(SampleId,ClusterId,ResolvedType) %>% summarise(Clusters=1,SglCount=sum(SGLs)) %>%
       group_by(SampleId,ResolvedType) %>% summarise(Clusters=sum(Clusters),SglCount=sum(SglCount)) %>%
       group_by(ResolvedType) %>% summarise(SglCount=sum(SglCount),Clusters=sum(Clusters)))

# some very large clusters
View(both17 %>% filter(ResolvedType=='COMPLEX') %>% group_by(SampleId,ClusterId,ClusterCount=2**round(log(ClusterCount,2))) %>% count 
     %>% group_by(ClusterCount) %>% summarise(Count=n()))

bfbClusters = svData %>% filter(ClusterCount>=3&((FoldbackLenStart>0|FoldbackLenEnd>0)|Type=='SGL')) %>% group_by(SampleId,ClusterId) %>% 
  summarise(Foldbacks=sum(FoldbackLenStart>0|FoldbackLenEnd>0),SglCount=sum(Type=='SGL')) %>% ungroup()
rm(bfbClusters)

# numbers of SGLs per cluster
sglClusters = both17 %>% filter(ResolvedType=='COMPLEX') %>% group_by(SampleId,ClusterId) %>% summarise(SglCentroCount=n())
sglClusters = merge(sglClusters,bfbClusters,by=c('SampleId','ClusterId'),all.x=T)
sglClusters = sglClusters %>% mutate(HasBFB=!is.na(Foldbacks)&Foldbacks>0)
View(sglClusters %>% group_by(SglCount=pmin(SglCentroCount,5),HasBFB) %>% count %>% spread(HasBFB,n))
View(sglClusters %>% group_by(SglCount=pmin(SglCount,5),HasBFB) %>% count %>% spread(HasBFB,n))
View(sglClusters)

sglAllClusters = merge(bfbClusters,sglClusters %>% select(-Foldbacks,-SglCount),by=c('SampleId','ClusterId'),all.x=T)
View(sglAllClusters)
rm(sglAllClusters)
sglAllClusters = sglAllClusters %>% mutate(HasBFB=!is.na(Foldbacks)&Foldbacks>0,
                                           HasCentroSgls=!is.na(SglCentroCount)&SglCentroCount>0)

View(sglAllClusters %>% filter(SglCount>0) %>% group_by(SglCount=pmin(SglCount,5),HasCentroSgls,HasBFB) %>% count %>% spread(HasBFB,n,fill=0))


# off-setting pairs - low evidence for shards, looks more likely to be incomplete (or incorrect orientation)
View(sample17 %>% filter(SglCount==2) %>% group_by(NetOrient) %>% count)

# LOH on the other arm
View(cnChrData %>% filter(Chromosome==17&(PCentroGain|QCentroGain)) %>% group_by(PCentroGain,QCentroGain,LohP,LohQ) %>% count)
colnames(cnChrData)

centroGainData = centroGainData %>% mutate(HasCentroSGLs=!is.na(SglCount)&SglCount>0)
View(centroGainData %>% filter(Chromosome==17&QCentroGain) %>% group_by(HasCentroSGLs,QCentroGain,LohP) %>% count)
View(centroGainData %>% filter(Chromosome==17&QCentroGain) %>% group_by(HasCentroSGLs,QCentroGain,LohP,QArmGain) %>% count %>% spread(QArmGain,n,fill=0))
View(centroGainData %>% filter(Chromosome==17&QCentroGain&!HasCentroSGLs&LohP=='true'))
View(centroGainData)

View(sglCentroSvData %>% group_by(HighConf) %>% count)
View(sglUnmapSvData %>% group_by(HighConf) %>% count)


#####
## Centromeric gains matched by other centromeres where no SGLs exist to explain them
noSglCentroData = centroGainSglData %>% filter(!HasSGLs&HasCentroGain)
View(noSglCentroData)

# chromosomes 1,5,8,19 and 17 have most non-SGL-related gain
View(noSglCentroData %>% group_by(Chromosome,PCentroGain,QCentroGain) %>% summarise(Count=n(),
                                                                                    MedPloidy=round(median(CentroCNChg),1),
                                                                                    AvgPloidy=round(mean(CentroCNChg),1)))

View(noSglCentroData %>% group_by(SampleId) %>% summarise(ChrCount=n(),MedPloidy=round(median(CentroCNChg),1)) %>% 
       group_by(ChrCount) %>% summarise(Count=n(),MedPloidy=median(MedPloidy)))

chrArmNoSgls = noSglCentroData %>% group_by(SampleId,ChrArm=paste(Chromosome,ifelse(PCentroGain,'P','Q'),sep='_'),
                                            CentroCnGain=pmin(pmax(round(CentroCNChg),1),4)) %>% summarise(Count=n()) %>%
  group_by(ChrArm,CentroCnGain) %>% summarise(Count=n())

chrArmNoSgls = chrArmNoSgls %>% spread(CentroCnGain,Count,fill=0) %>% gather('CentroCnGain','Count',2:5) # regig to get zeros

print(ggplot(chrArmNoSgls,aes(x=ChrArm,y=CentroCnGain)) 
      + geom_tile(aes(fill=Count),colour="white",stat = "identity",position="identity") 
      + geom_text(aes(label=ifelse(!is.na(Count),Count,0)),size=3)
      + scale_fill_gradient(low="steelblue1",high="royalblue4")
      + theme(axis.text.x = element_text(angle=90,hjust=1,size=10)))


# faceted by the other arm's CN - ie the baseline CN for the chromosome
chrArmNoSgls2 = noSglCentroData %>% group_by(SampleId,ChrArm=paste(Chromosome,ifelse(PCentroGain,'P','Q'),sep='_'),
                                            CentroCnGain=pmin(pmax(round(CentroCNChg),1),4),
                                            BaselineCn=pmin(pmax(round(pmin(CentroCnP,CentroCnQ)),1),4)) %>% summarise(Count=n()) %>%
  group_by(ChrArm,BaselineCn,CentroCnGain) %>% summarise(Count=n())

# View(chrArmNoSgls2)

chrArmNoSgls2 = chrArmNoSgls2 %>% spread(BaselineCn,Count,fill=0) %>% gather('BaselineCn','Count',3:6)
chrArmNoSgls2 = chrArmNoSgls2 %>% spread(CentroCnGain,Count,fill=0) %>% gather('CentroCnGain','Count',3:6)
View(chrArmNoSgls2)

print(ggplot(chrArmNoSgls2,aes(x=ChrArm,y=CentroCnGain)) 
      + geom_tile(aes(fill=Count),colour="white",stat = "identity",position="identity") 
      + geom_text(aes(label=ifelse(!is.na(Count),Count,0)),size=3)
      + scale_fill_gradient(low="steelblue1",high="royalblue4")
      + facet_wrap(~BaselineCn)
      + theme(axis.text.x = element_text(angle=90,hjust=1,size=10)))

View(noSglCentroData %>% filter(Chromosome==1&round(CentroCNChg)==1) %>% group_by(BaselineCn=pmin(pmax(round(pmin(CentroCnP,CentroCnQ)),1),4)) %>% count)

# look at samples with a single chromosome with gain - are these explained by foldbacks?
sampleChrCounts = noSglCentroData %>% group_by(SampleId) %>% summarise(ChrGainCount=n())

singleChrSamples = sampleChrCounts %>% filter(ChrGainCount==1)
View(noSglCentroData %>% filter(SampleId %in% singleChrSamples$SampleId))
View(noSglCentroData %>% filter(SampleId %in% singleChrSamples$SampleId) %>% 
       group_by(CentroCNChg=pmin(round(CentroCNChg),4),
                LowerCN=pmin(pmax(round(pmin(CentroCnP,CentroCnQ)),1),4)) %>% count %>% spread(CentroCNChg,n))

View(centroGainSglData %>% filter(HasCentroGain&SampleId %in% singleChrSamples$SampleId) %>% 
       group_by(HasSGLs,
                CentroCNChg=pmin(round(CentroCNChg),4),
                LowerCN=pmin(pmax(round(pmin(CentroCnP,CentroCnQ)),1),4)) %>% count %>% spread(CentroCNChg,n))




noSglCentroChrData = merge(noSglCentroChrData,sampleChrCounts,by='SampleId',all.x=T)
noSglCentroChrData = noSglCentroData %>% group_by(SampleId,CentroCNChg=pmax(pmin(round(CentroCNChg),3),1)) %>% summarise(ChrCount=n())
View(noSglCentroChrData)

View(noSglCentroChrData %>% group_by(SampleId) %>% summarise(CnChgBuckets=n(),ChrCount=sum(ChrCount)) %>% group_by(CnChgBuckets,ChrCount) %>% count %>% spread(CnChgBuckets,n))
#View(noSglCentroChrData %>% filter(ChrGainCount==1) %>% group_by(CentroCNChg) %>% count)
View(noSglCentroChrData %>% group_by(ChrGainCount,CentroCNChg) %>% count %>% spread(CentroCNChg,n))


View(noSglCentroData %>% group_by(SampleId,CentroCNChg=pmax(pmin(round(CentroCNChg),3),1)) %>% summarise(ChrCount=n()) %>% 
       group_by(ChrCount,CentroCNChg) %>% summarise(Count=n()))



######
## SGLs without repeat mask set, possibly mappable in HG38
sglUnmapped = svData %>% filter(Type=='SGL'&RepeatType==''&RepeatClass==''&stri_length(InsertSeq)>=25)
sglUnmapped = sglUnmapped %>% filter(!grepl('TTTTTTTTTTTTTTTT',InsertSeq)&!grepl('AAAAAAAAAAAAAAAA',InsertSeq)
                                     &!grepl('CCCCCCCCCCCCCCCC',InsertSeq)&!grepl('GGGGGGGGGGGGGGGG',InsertSeq))

nrow(sglUnmapped)
View(sglUnmapped %>% group_by(ChrStart,ArmStart) %>% count %>% spread(ArmStart,n))

View(sglUnmapped %>% select(SampleId,Id,InsertSeq))

write.csv(sglUnmapped[1:2000,] %>% mutate(Id=paste(SampleId,Id,sep='_')) %>% select(Id,InsertSeq),'~/logs/sgl_unmap_ins_seq_01.fa',row.names = F, quote = F)
write.csv(sglUnmapped[2001:4000,] %>% mutate(Id=paste(SampleId,Id,sep='_')) %>% select(Id,InsertSeq),'~/logs/sgl_unmap_ins_seq_02.fa',row.names = F, quote = F)
write.csv(sglUnmapped[4001:6000,] %>% mutate(Id=paste(SampleId,Id,sep='_')) %>% select(Id,InsertSeq),'~/logs/sgl_unmap_ins_seq_03.fa',row.names = F, quote = F)
write.csv(sglUnmapped[6001:8000,] %>% mutate(Id=paste(SampleId,Id,sep='_')) %>% select(Id,InsertSeq),'~/logs/sgl_unmap_ins_seq_04.fa',row.names = F, quote = F)
write.csv(sglUnmapped[8001:nrow(sglUnmapped),] %>% mutate(Id=paste(SampleId,Id,sep='_')) %>% select(Id,InsertSeq),'~/logs/sgl_unmap_ins_seq_05.fa',row.names = F, quote = F)


sglUnmapData = unmapBlatResults %>% arrange(QName,-Score) %>% group_by(QName) %>% 
  summarise(ChrEnd=first(TName),Strand=first(Strand),ObePosStart=first(TStart),ObePosEnd=first(TEnd),
            TopScore=first(Score),NextScore=first(SecondScore),NextChrScore=max(ifelse(Score==TopScore,0,Score)),
            InsSeqLength=first(QSize),Match=first(Match),Mismatch=first(Mismatch)) %>% ungroup()

View(sglUnmapData)

sglUnmapData = sglUnmapData %>% mutate(HighConf=TopScore>900&(TopScore-NextChrScore>25))
sglUnmapData = sglUnmapData %>% separate(QName,c('SampleId','Id'),sep='_')
sglUnmapData$ObeChr = stri_replace_all_fixed(sglUnmapData$ObeChr,'chr','')

sglUnmapData = sglUnmapData %>% filter(SampleId %in% sampleCancerTypes$SampleId)

# merge with SGL SV data and centromere change data
View(sglUnmapData)
sglUnmapSvData = merge(sglUnmapData,
                       sglUnmapped %>% select(SampleId,Id,ChrStart,PosStart,OrientStart,ArmStart,ClusterId,ClusterCount,ResolvedType,
                                              Ploidy,CNChgStart,CNStart,GeneStart),
                        by=c('SampleId','Id'),all.x=T)

View(sglUnmapSvData %>% filter(HighConf&ChrStart==UnmapChr) %>%group_by(OrientStart,UnmapStrand) %>% count %>% spread(UnmapStrand,n))

# inferred SV types
sglUnmapSvData = sglUnmapSvData %>% mutate(OrientEnd=ifelse((OrientStart==1&Strand=='+')|(OrientStart==-1&Strand=='-'),-1,1),
                                             PosEnd=ifelse(OrientEnd==-1,ObePosStart,ObePosEnd),
                                             SglType=ifelse(ChrStart!=ObeChr,'BND',
                                                            ifelse(OrientStart==OrientEnd,'INV',
                                                                   ifelse(OrientStart==1&OrientEnd==-1,'DEL','DUP'))),
                                             SglLength=ifelse(SglType=='BND',0,abs(PosStart-(ObePosStart+ObePosEnd)*0.5)))

names(sglUnmapSvData)[names(sglUnmapSvData) == 'ObeChr'] <- 'ChrEnd'
View(sglUnmapSvData)

View(sglUnmapSvData %>% group_by(ArmStart,OrientStart,OrientEnd) %>% count %>% spread(ArmStart,n))

View(sglUnmapSvData %>% group_by(SglType,LengthBucket=ifelse(SglType=='BND',0,2**round(log(SglLength,2)))) %>% count %>% spread(SglType,n))
View(sglUnmapSvData %>% group_by(SglType,HighConf) %>% count %>% spread(HighConf,n))

print(ggplot(data = sglUnmapSvData %>% filter(HighConf&SglType!='BND') %>% 
               group_by(SglType,LengthBucket=ifelse(SglType=='BND',0,2**round(log(SglLength,2)))) %>% count, 
             aes(x=LengthBucket, y=n))
      + geom_line()
      + scale_x_log10()
      + facet_wrap(~SglType)
      + labs(title = "Unmapped SGL Types and Lengths"))

# any SGLs ending in the heterochromatin region of chr 1 (125184587-143184587)
View(sglUnmapSvData %>% filter(ChrEnd==1&PosEnd>2e6) %>% group_by(Position=round(PosEnd,-6)) %>% count)

# &PosEnd>1.1e8&PosEnd<1.6e8
print(ggplot(sglUnmapSvData %>% filter(ChrEnd==1&PosEnd>4e6) %>% group_by(Position=round(PosEnd,-6),OrientEnd) %>% count, 
             aes(x=Position, y=n))
      + geom_bar(stat = "identity", colour = "black")
      + facet_wrap(~OrientEnd)
      + labs(title = "SGLs to chromosome 1"))


View(sglUnmapSvData)








## DEBUG



View(sampleChrData %>% group_by(SglCount) %>% count)
View(sampleChrData %>% filter(CentroChr==1))

View(sampleChrData %>% filter(CentroChr==1&QCentroGain&SglCount==1))
View(sampleChrData %>% filter(CentroChr==1&SglCount==1) %>% group_by(FirstCentroOrient,PCentroGain,QCentroGain) %>% count)
nrow(sampleChrData %>% filter(CentroChr==1&SglCount==1))
View(sampleChrData %>% filter(CentroChr==12&SglCount==1) %>% group_by(FirstCentroOrient,PCentroGain,QCentroGain) %>% count)
View(sampleChrData %>% filter(CentroChr==1) %>% group_by(NetOrient=ifelse(NetOrient>0,1,ifelse(NetOrient==0,0,-1)),PCentroGain,QCentroGain) %>% count)

View(sglCentroSvData %>% filter(HighConf&CentroChr==1) %>% group_by(OrientEnd,QCentroGain) %>% count)

View(sglCentroSvData %>% filter(CentroChr==1&ChrStart==19) %>% group_by(SampleId) %>% count)
View(sglCentroSvData %>% filter(CentroChr==1&ChrStart==5) %>% group_by(SampleId) %>% count)

View(sampleChrData %>% filter(SglCount==2))
View(sampleChrData %>% filter(SglCount==2) %>% group_by(NetOrient) %>% count)
View(sampleChrData %>% filter(SglCount==2) %>% group_by(SourceChrCount) %>% count)
View(sampleChrData %>% filter(CentroChr==1) %>% group_by(PCentroGain,QCentroGain) %>% count)

View(sglCentroData)
View(sglCentroSvData %>% group_by(CentroChr,OrientEnd) %>% count)
View(sglCentroData %>% group_by(CentroChr) %>% count)
View(sglCentroData %>% group_by(CentroChr,HighConf) %>% count %>% spread(HighConf,n))


# DEBUG and working

blatRawData = read.csv('~/dev/blat/output_processed.csv')


filteredBlatResults = load_blat_results('~/data/sv/sglCentroOutput_subset.csv')
View(filteredBlatResults)




blatRawData = read.csv('~/data/sv/sglCentroOutput_subset.csv')
colnames(blatRawData)
View(blatResults)

colnames(blatResults)

filteredBlatResults = blatRawData %>% filter(!grepl('chrUn',TName)&!grepl('_random',TName))
View(filteredBlatResults)

filteredBlatResults = filteredBlatResults %>% mutate(QSizeMult=9-floor(QSize/100),
                                                     Score=(1000-(QSizeMult*Mismatch+QGapCount+TGapCount))*pmin(Match/QSize,1))

View(filteredBlatResults %>% filter(QSize>200&QSize<300&Mismatch==1&(QGapCount<=2|TGapCount<=2)) %>% select(QName,QSize,QSizeMult,Score,Match,RepMatch,Mismatch,QGapCount,TGapCount))
View(filteredBlatResults %>% select(QName,QSize,QSizeMult,Score,Match,RepMatch,Mismatch,QGapCount,TGapCount))


View(filteredBlatResults %>% filter(QName=='CPCT02010003T_431') %>% select(QName,TName,QSize,QSizeMult,Score,Match,RepMatch,Mismatch,QGapCount,TGapCount))

# s = y(m + f(r)) - ym′ - gq - gt,
View(filteredBlatResults %>% select(QName,QSize,Score,Match,RepMatch,Mismatch,QGapCount,TGapCount))

blatResults = filteredBlatResults %>% group_by(QName,TName) %>% summarise(Records=n(),
                                                                          AvgScore=round(mean(Score)),
                                                                          MaxScore=round(max(Score)))

View(blatResults %>% group_by(QName,MaxScoreBucket=round(MaxScore,-2)) %>% count %>% spread(MaxScoreBucket,n))
View(blatResults %>% filter(QName=='CPCT02010122T_109'))


# related to genes
View(sglCentros %>% group_by(ChrStart,Position=round(PosStart,-6)) %>% count)

# by position
print(ggplot(data = sglCentros %>% group_by(ChrStart,Position=round(PosStart,-6)) %>% count, 
             aes(x=Position,y=n))
      + geom_line()
      + facet_wrap(~ChrStart))




View(sglCentros %>% group_by(PloidyBucket=2**round(log(Ploidy,2)),
                             ClusterCountBucket=2**round(log(ClusterCount,2))) 
     %>% count() %>% spread(ClusterCountBucket,n))







# Questions about centromeric SGLs
# - what sort of complex clusters are they in - shattering, BFB etc
# - are they close to foldbacks or other activity
# - are they close to the centromere

View(centroLengths)

sglCentros = sglCentros %>% mutate(CentroType=ifelse(RepeatClass=='Satellite/centr','Centro','PeriC'),
                                   SampleClusterId=paste(SampleId,ClusterId,sep='_'))

sglCentros = merge(sglCentros,centroLengths,by.x='ChrStart',by.y='Chromosome',all.x=T)
View(sglCentros)
sglCentros = sglCentros %>% mutate(CentroDistance=pmax(abs(PosStart-Centromere)-1.5e6,1),
                                   CentroDistPerc=ifelse(ArmStart=='P',CentroDistance/PArmLength,CentroDistance/QArmLength),
                                   CentroDistBucket=2**round(log(CentroDistance,2)),
                                   HasDB=DBLenStart>-1000)

sglCentros = merge(sglCentros,sglCentroClusters %>% select(SampleClusterId,BFB,Foldbacks,SglCentroCount,MaxPloidy),by='SampleClusterId',all.x=T)
View(sglCentros %>% select(SampleClusterId,CentroType,ChrStart,PosStart,OrientStart,Centromere,CentroDistance,SglCentroCount,Foldbacks,everything()))

sglCentros = merge(sglCentros,sampleCancerTypes,by='SampleId',all.x=T)

# types of clusters
View(sglCentros %>% filter(CentroType=='PeriC'&ResolvedType!='DUP_BE'&ResolvedType!='LOW_VAF') %>% group_by(ResolvedType) %>% count)

# where do peri-centromeric SGLs occur
View(sglCentros %>% filter(CentroType=='PeriC'&ResolvedType!='DUP_BE'&ResolvedType!='LOW_VAF') %>% 
       group_by(ResolvedType,PercToTelomere=round(CentroDistPerc,1)) %>% count %>% spread(ResolvedType,n))

View(sglCentros  %>% group_by(ResolvedType) %>% count)

View(sglCentros %>% filter(CentroType=='Centro'&ResolvedType!='DUP_BE'&ResolvedType!='LOW_VAF') %>% group_by(ResolvedType,CentroDistBucket) 
     %>% count %>% spread(ResolvedType,n))


View(sglCentros %>% group_by(SampleId) %>% summarise(SglCount=sum(CentroType=='Centro')))

sglGenes = sglCentros %>% mutate(Length=get_length(Type,PosStart,PosEnd),
                                             LengthMin=get_length_min(Type,Length),
                                             LengthMax=get_length_max(Type,Length,LengthMin),
                                             Type=ifelse(LengthMin==lengthLong,longDelDupInv,as.character(Type)))

sglGenes = getBreakendsInGenes(sglGenes)
View(sglGenes %>% group_by(GeneName) %>% count)


View(sglCentros %>% group_by(ResolvedType,InBFB) %>% count %>% spread(InBFB,n))

View(sglCentros %>% filter(ResolvedType=='COMPLEX') %>% group_by(SampleId,ClusterId,InBFB) %>% count %>% group_by(InBFB,n) %>% count %>% spread(InBFB,nn))

bfbClusters = clusters %>% filter(Foldbacks>2&ResolvedType=='COMPLEX')
bfbClusters = bfbClusters %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))
View(clusters %>% filter(ResolvedType=='COMPLEX') %>% group_by(BFB=(Foldbacks>2&ResolvedType=='COMPLEX')) %>% summarise(sum(ClusterCount)))


sglCentros = sglCentros %>% mutate(BoundsLOH=(OrientStart==1&MinorAPStartPost<0))






nrow(sglCentros)
sglCentroArms = sglCentros %>% group_by(SampleId,Chr=ChrStart,Arm=ArmStart) %>% summarise(SglCentroCount=n())
View(sglCentroArms)

fbArmData = foldbacks %>% group_by(SampleId,Chr,Arm) %>% summarise(FoldbackCount=sum(ifelse(IsChained,0.5,1)))
fbArmData = fbArmData %>% ungroup() %>% mutate(FoldbackCount=ifelse(FoldbackCount<1,1,round(FoldbackCount)))
View(fbArmData)

View(foldbacks %>% filter(SampleId=='CPCT02010386T'))

cmpSvStart = svData %>% filter(ResolvedType=='COMPLEX') %>% select(SampleId,ClusterId,Chr=ChrStart,Arm=ArmStart) 
cmpSvEnd = svData %>% filter(ResolvedType=='COMPLEX'&Type!='SGL'&Type!='INF') %>% select(SampleId,ClusterId,Chr=ChrEnd,Arm=ArmEnd) 
complexArmData = rbind(cmpSvStart,cmpSvEnd)
complexArmData = complexArmData %>% group_by(SampleId,Chr,Arm,ClusterId) %>% summarise(ComplexBreakendCount=n())
complexArmData = complexArmData %>% group_by(SampleId,Chr,Arm) %>% summarise(ComplexEventCount=n())
View(complexArmData)

unbalTransStart = svData %>% filter(ClusterCount==1&ResolvedType=='UNBAL_TRANS') %>% select(SampleId,ClusterId,Chr=ChrStart,Arm=ArmStart)
unbalTransEnd = svData %>% filter(ClusterCount==1&ResolvedType=='UNBAL_TRANS') %>% select(SampleId,ClusterId,Chr=ChrEnd,Arm=ArmEnd)
unbalTrans = rbind(unbalTransStart,unbalTransEnd)
unbalTrans = unbalTrans %>% group_by(SampleId,Chr,Arm) %>% summarise(UnbalTransCount=n())
View(unbalTrans)


nrow(cnChrData)
View(cnChrData)
cnChrSvData = merge(cnChrData,unbalTrans %>% filter(Arm=='P') %>% select(SampleId,Chromosome=Chr,UnbalTransP=UnbalTransCount),
                    by=c('SampleId','Chromosome'),all.x=T)
cnChrSvData = merge(cnChrSvData,unbalTrans %>% filter(Arm=='Q') %>% select(SampleId,Chromosome=Chr,UnbalTransQ=UnbalTransCount),
                    by=c('SampleId','Chromosome'),all.x=T)

cnChrSvData = merge(cnChrSvData,complexArmData %>% filter(Arm=='P') %>% select(SampleId,Chromosome=Chr,CmpEventP=ComplexEventCount),
                    by=c('SampleId','Chromosome'),all.x=T)
cnChrSvData = merge(cnChrSvData,complexArmData %>% filter(Arm=='Q') %>% select(SampleId,Chromosome=Chr,CmpEventQ=ComplexEventCount),
                    by=c('SampleId','Chromosome'),all.x=T)

cnChrSvData = merge(cnChrSvData,fbArmData %>% filter(Arm=='P') %>% select(SampleId,Chromosome=Chr,FoldbacksP=FoldbackCount),
                    by=c('SampleId','Chromosome'),all.x=T)
cnChrSvData = merge(cnChrSvData,fbArmData %>% filter(Arm=='Q') %>% select(SampleId,Chromosome=Chr,FoldbacksQ=FoldbackCount),
                    by=c('SampleId','Chromosome'),all.x=T)

cnChrSvData = merge(cnChrSvData,sglCentroArms %>% filter(Arm=='P') %>% select(SampleId,Chromosome=Chr,SglCentroP=SglCentroCount),
                    by=c('SampleId','Chromosome'),all.x=T)
cnChrSvData = merge(cnChrSvData,sglCentroArms %>% filter(Arm=='Q') %>% select(SampleId,Chromosome=Chr,SglCentroQ=SglCentroCount),
                    by=c('SampleId','Chromosome'),all.x=T)

View(cnChrSvData %>% filter(UnbalTransQ>0))
View(cnChrSvData %>% filter(CmpEventQ>0))
View(cnChrSvData %>% filter(FoldbacksP>0))
View(cnChrSvData %>% filter(SglCentroP>0))
nrow(cnChrSvData)
cnChrSvData[is.na(cnChrSvData)] = 0
colnames(cnChrSvData)

write.csv(cnChrSvData,'~/data/sv/chr_arm_cn_sv_data.csv', row.names = F,quote = F)

sglFbCentroSummary = cnChrSvData %>% filter(!Chromosome %in% c(13,14,15,21,22,'X','Y')) %>% group_by(SampleId) %>% 
  summarise(CentroGain=sum(PCentroGain|QCentroGain),
            SglCentro=sum(SglCentroQ>0|SglCentroP>0),
            FBArms=sum(FoldbacksQ>0|FoldbacksP>0))

View(sglFbCentroSummary)


# gains in centromere vs number of centromeric SGLs per arm
print(ggplot(data = sglFbCentroSummary %>% filter(CentroGain>0|SglCentro>0|FBArms>0), aes(x=CentroGain,y=SglCentro))
      + geom_point(position="jitter")
      + geom_smooth(,method=lm,se=FALSE, fullrange=F)
      + labs(x='Arms with Centromeric Copy Number Change', y='Arms with Centromeric SGLs',title = "Chromosomal Arms per Sample with Centromeric Copy Number Change vs Centromeric SGLs"))

# centromeric SGLs vs foldback arms
print(ggplot(data = sglFbCentroSummary %>% filter(CentroGain>0|SglCentro>0|FBArms>0), aes(x=FBArms,y=SglCentro))
      + geom_point(position="jitter")
      + geom_smooth(,method=lm,se=FALSE, fullrange=F)
      + labs(x='Arms with Foldbacks', y='Arms with Centromeric SGLs', title = "Chromosomal Arms per Sample with Foldbacks vs Centromeric SGLs"))


print(ggplot(data = sglFbCentroSummary %>% filter(CentroGain>0|SglCentro>0|FBArms>0), aes(x=CentroGain))
      + geom_point(aes(y=SglCentro,colour='SglCentro'),position="jitter")
      + geom_point(aes(y=FBArms,colour='FBArms'),position="jitter")
      + geom_smooth(aes(y=SglCentro,colour='SglCentro'),method=lm,se=FALSE, fullrange=F)
      + geom_smooth(aes(y=FBArms,colour='FBArms'),method=lm,se=FALSE, fullrange=F)
      + labs(title = "Centromere Gain vs SGL-Centros & FB Arms"))


print(ggplot(data = sglFbCentroSummary %>% filter(CentroGain>0|SglCentro>0|FBArms>0), aes(x=CentroGain,y=SglCentro))
      + geom_hex()
      + labs(title = "Centromere Gain vs SGL-Centros & FB Arms"))

print(ggplot(data = sglFbCentroSummary %>% filter(CentroGain>0|SglCentro>0|FBArms>0), aes(x=CentroGain))
      + geom_hex(aes(y=SglCentro,colour='SglCentro'))
      + geom_hex(aes(y=FBArms,colour='FBArms'))
      + geom_smooth(aes(y=SglCentro,colour='SglCentro'),method=lm,se=FALSE, fullrange=F)
      + geom_smooth(aes(y=FBArms,colour='FBArms'),method=lm,se=FALSE, fullrange=F)
      + labs(title = "Centromere Gain vs SGL-Centros & FB Arms"))

View(cnChrSvData %>% filter(!Chromosome %in% c(13,14,15,21,22,'X','Y')) %>% group_by(PCentroGain,QCentroGain,HasSglCentro=SglCentroQ>0|SglCentroP>0) %>% count() %>% spread(HasSglCentro,n))
View(cnChrSvData %>% group_by(PCentroGain|QCentroGain,HasFB=FoldbacksQ>0|FoldbacksP>0,HasSglCentro=SglCentroQ+SglCentroP>0) %>% count() %>% spread(HasSglCentro,n))

View(cnChrSvData %>% group_by(HasFB=FoldbacksQ>0|FoldbacksP>0,HasSglCentro=SglCentroQ+SglCentroP>0) %>% count() %>% spread(HasSglCentro,n))
View(cnChrSvData %>% group_by(PCentroGain|QCentroGain,HasSglCentro=SglCentroQ+SglCentroP>0) %>% count() %>% spread(HasSglCentro,n))

View(cnChrSvData %>% filter(!Chromosome %in% c(13,14,15,21,22,'X','Y')) %>% group_by(PCentroGain|QCentroGain,HasFB=FoldbacksQ>0|FoldbacksP>0,HasCMP=CmpEventQ>0|CmpEventP>0,HasSglCentro=SglCentroQ+SglCentroP) %>% count() %>% spread(HasSglCentro,n))
View(cnChrSvData %>% filter(!Chromosome %in% c(13,14,15,21,22,'X','Y')) %>% group_by(PCentroGain|QCentroGain,HasFB=FoldbacksQ>0|FoldbacksP>0,HasCMP=CmpEventQ>0|CmpEventP>0,HasSglCentro=SglCentroQ+SglCentroP==1) %>% count() %>% spread(HasSglCentro,n))




#####
## Telomeric SGLs

sampleCancerAndSubTypes = load_cancer_types('~/data/hpc_sample_cancer_types.csv',T,10)

sglTelos = svData %>% filter(Type=='SGL'&(RepeatClass=='Simple_repeat'&(RepeatType=='(TTAGGG)n'|RepeatType=='(CCCTAA)n')))
nrow(sglTelos)

sglTelos = sglTelos %>% mutate(HasGGGTTA=grepl('GGGTTA',InsertSeq),
                               HasCCCTAA=grepl('CCCTAA',InsertSeq),
                               TowardsTelo=(HasGGGTTA&OrientStart==-1)|(HasCCCTAA&OrientStart==1),
                               InsertSeqLength=stri_length(InsertSeq),
                               RepeatCount=pmax(str_count(InsertSeq,'GGGTTA'),str_count(InsertSeq,'CCCTAA')))

View(sglTelos %>% group_by(HasGGGTTA,HasCCCTAA,OrientStart,TowardsTelo) %>% count())

sglPatternCounts = sglTelos %>% mutate(GGGTTA=str_count(InsertSeq,'GGGTTA'),
                         CCCTAA=str_count(InsertSeq,'CCCTAA'),
                         GGGGTT=str_count(InsertSeq,'GGGGTT'),
                         GGGTCA=str_count(InsertSeq,'GGGTCA'),
                         GGGTGA=str_count(InsertSeq,'GGGTGA'),
                         GGGGTTC=str_count(InsertSeq,'GGGGTTC'),
                         CCCTGA=str_count(InsertSeq,'CCCTGA'),
                         CCCCGA=str_count(InsertSeq,'CCCCGA'),
                         CCCTCA=str_count(InsertSeq,'CCCTCA'),
                         CCCAGC=str_count(InsertSeq,'CCCAGC'),
                         CCCACA=str_count(InsertSeq,'CCCACA'),
                         CCCCAA=str_count(InsertSeq,'CCCCAA'),
                         InsertSeqLength=stri_length(InsertSeq),
                         Total=GGGTTA+CCCTAA+GGGGTT+CCCTGA+GGGTCA+GGGGTTC+CCCCGA+CCCTCA+GGGTGA+CCCAGC+CCCACA+CCCCAA,
                         Unaccounted=InsertSeqLength-Total*6)

sglPatternCounts2 = sglPatternCounts %>% select(SampleId,Id,TeloSglCount,TowardsTelo,GGGTTA,CCCTAA,
                                                GGGGTT,GGGTCA,GGGTGA,GGGGTTC,CCCTGA,CCCCGA,CCCTCA,CCCAGC,CCCACA,CCCCAA)

colCount=ncol(sglPatternCounts2)
sglPatternCounts2 = sglPatternCounts2 %>% gather('Pattern','Count',5:colCount)
View(sglPatternCounts2)

print(ggplot(sglPatternCounts2 %>% group_by(Pattern,TeloSglCount,Direction=ifelse(TowardsTelo,'Faces Telomere','Faces Centromere')) %>% summarise(Count=sum(Count)), 
             aes(x=TeloSglCount, y=Count, fill=Pattern))
      + geom_bar(stat = "identity", colour = "black")
      + facet_wrap(~Direction)
      + labs(x='SGLs per Sample', y='Count per Pattern',title='# of Telomeric SGLs by Insert Sequence Pattern by Directionality'))

sglPatternCounts3 = sglPatternCounts %>% select(SampleId,Id,InsertSeqLength,TowardsTelo,GGGTTA,CCCTAA,
                                                GGGGTT,GGGTCA,GGGTGA,GGGGTTC,CCCTGA,CCCCGA,CCCTCA,CCCAGC,CCCACA,CCCCAA)

colCount=ncol(sglPatternCounts3)
sglPatternCounts3 = sglPatternCounts3 %>% gather('Pattern','Count',5:colCount)
View(sglPatternCounts3)

sglPatternCounts3 = sglPatternCounts3 %>% mutate(InsertSeq=10*round(InsertSeqLength/10))
sglPatternCounts3Totals = sglPatternCounts3 %>% group_by(InsertSeq,TowardsTelo) %>% summarise(PatternTotal=sum(Count))
sglPatternCounts3 = merge(sglPatternCounts3,sglPatternCounts3Totals,by=c('InsertSeq','TowardsTelo'),all.x=T)
View(sglPatternCounts3)

print(ggplot(sglPatternCounts3 %>% group_by(Pattern,InsertSeq,TowardsTelo) %>% summarise(Count=sum(Count)), 
             aes(x=InsertSeq, y=Count, fill=Pattern))
      + geom_bar(stat = "identity", colour = "black")
      + facet_wrap(~TowardsTelo)
      + labs(x='Insert Sequence Length', y='Count per Pattern'))

print(ggplot(sglPatternCounts3 %>% group_by(Pattern,InsertSeq=pmin(InsertSeq,150),TowardsTelo) %>% summarise(Percent=sum(Count)/first(PatternTotal)), 
             aes(x=InsertSeq, y=Percent, fill=Pattern))
      + geom_bar(stat = "identity", colour = "black")
      + facet_wrap(~TowardsTelo)
      + labs(x='Insert Sequence Length', y='% per Pattern'))


# %>% 
#       select(SampleId,Id,GGGTTA,CCCTAA,GGGGTT,GGGGTTC,GGGTCA,GGGTGA,CCCTGA,CCCCGA,CCCTCA,CCCAGC,CCCACA,CCCCAA,
#              Total,InsertSeqLength,Unaccounted,InsertSeq)) # %>% filter(GGGTTA==0&CCCTAA==0)

possibleSgls = svData %>% filter(Type=='SGL') %>% 
  mutate(Identified=(RepeatClass=='Simple_repeat'&(RepeatType=='(TTAGGG)n'|RepeatType=='(CCCTAA)n')),
    GGGTTA=str_count(InsertSeq,'GGGTTA'),
         CCCTAA=str_count(InsertSeq,'CCCTAA'),
          GGGGTT=str_count(InsertSeq,'GGGGTT'),
          GGGTCA=str_count(InsertSeq,'GGGTCA'),
          GGGTGA=str_count(InsertSeq,'GGGTGA'),
          GGGTTA=str_count(InsertSeq,'GGGTTA'),
          GGGGTTC=str_count(InsertSeq,'GGGGTTC'),
          CCCTGA=str_count(InsertSeq,'CCCTGA'),
          CCCCGA=str_count(InsertSeq,'CCCCGA'),
          CCCTCA=str_count(InsertSeq,'CCCTCA'),
          CCCAGC=str_count(InsertSeq,'CCCAGC'),
          CCCACA=str_count(InsertSeq,'CCCACA'),
          CCCTAA=str_count(InsertSeq,'CCCTAA'),
          CCCCAA=str_count(InsertSeq,'CCCCAA'),
          InsertSeqLength=stri_length(InsertSeq),
          Total=GGGTTA+CCCTAA+GGGGTT+CCCTGA+GGGTCA+GGGGTTC+CCCCGA+CCCTCA+GGGTGA+CCCAGC+CCCACA+CCCTAA+GGGTTA+CCCCAA) %>% 
  filter(Total>=4)

View(possibleSgls %>% group_by(Identified) %>% count)
View(possibleSgls %>% group_by(Identified,PatternCount=5*round(Total/5)) %>% count %>% spread(Identified,n))
View(possibleSgls %>% filter(!Identified&Total<10) %>% select(Total,InsertSeq))
View(possibleSgls %>% filter(!Identified&Total<10) %>% select(Total,InsertSeq))


sglTelos = sglTelos %>% filter(SampleId %in% sampleCancerAndSubTypes$SampleId)
sglTelos = merge(sglTelos,sampleCancerAndSubTypes,by='SampleId',all.x=T)

colnames(sglTelos)
print(ggplot(data = sglTelos %>% filter(!TowardsTelo==T) %>% group_by(InsertSeqLength=pmin(InsertSeqLength,250),TeloSglCount) %>% count, aes(x=InsertSeqLength,y=n))
      + geom_line()
      + facet_wrap(~TeloSglCount))

View(sglTelos %>% group_by(InsertSeqLength=pmin(InsertSeqLength,250),TeloSglCount,Direction=ifelse(TowardsTelo,'TowardsTelo','TowardsCentro')) 
     %>% count %>% spread(Direction,n,fill=0))

print(ggplot(data = sglTelos %>% group_by(InsertSeqLength=round(pmin(InsertSeqLength,250),-1),
                                          TeloSglCount,Direction=ifelse(TowardsTelo,'TowardsTelo','TowardsCentro')) 
             %>% count %>% spread(Direction,n,fill=0),aes(x=InsertSeqLength))
      + geom_line(aes(y=TowardsTelo,color='Towards Telomere'))
      + geom_line(aes(y=TowardsCentro,color='Towards Centromere'))
      + labs(x='Insert Sequennce',y='# of SGLs',color='Telomeric Sequence Direction',
             title='Insert Sequence Length Distribution of Telomeric SGLs by Facing Telomere vs Centromere')
      + facet_wrap(~TeloSglCount))

write.csv(sglTelos %>% select(SampleId,Id,TowardsTelo,InsertSeq),'~/logs/sgl_insert_seqs.csv',row.names = F, quote = F)


View(sglTelos %>% select(TowardsTelo,InsertSeqLength,RepeatCount,InsertSeq,everything()))
View(sglTelos %>% filter(InsertSeqLength>50&InsertSeqLength<80&TowardsTelo) %>% select(TowardsTelo,InsertSeq,InsertSeqLength,everything()))
View(sglTelos %>% group_by(TowardsTelo,InsertSeq=5*round(InsertSeqLength/5)) %>% count %>% spread(TowardsTelo,n))

View(sglTelos %>% group_by(TowardsTelo) %>% summarise(Count=n(),
                                                      AvgInsSeq=mean(InsertSeqLength),
                                                      MedInsSeq=median(InsertSeqLength)))

# non-SGL with telomeric inserts
nonSglTelos = svData %>% filter(Type!='SGL'&(RepeatClass=='Simple_repeat'&(RepeatType=='(TTAGGG)n'|RepeatType=='(CCCTAA)n')))
nrow(nonSglTelos)

nonSglTelos = nonSglTelos %>% mutate(HasGGGTTA=grepl('GGGTTA',InsertSeq),
                                     HasCCCTAA=grepl('CCCTAA',InsertSeq),
                                     InsertSeqLen=stri_length(InsertSeq))

View(nonSglTelos %>% group_by(HasGGGTTA,HasCCCTAA) %>% count())
View(nonSglTelos %>% group_by(Type,InsertSeqLen) %>% count())
View(nonSglTelos %>% group_by(HasGGGTTA,HasCCCTAA) %>% count())

View(sglTelos %>% group_by(HasGGGTTA,HasCCCTAA) %>% count())
View(sglTelos %>% group_by(HasGGGTTA,HasCCCTAA,OrientStart,TowardsTelo) %>% count())
View(sglTelos %>% group_by(CancerType,SampleId,TowardsTelo) %>% count() %>% spread(TowardsTelo,n))
View(sglTelos %>% filter(!(ChrStart %in% c(13,14,15,21,22))) %>% group_by(HasGGGTTA,HasCCCTAA,OrientStart,TowardsTelo,ArmStart) %>% count())
View(sglTelos %>% group_by(HasGGGTTA,HasCCCTAA,ChrStart) %>% count() %>% spread(ChrStart,n))
View(sglTelos %>% group_by(TowardsTelo,ChrStart,ArmStart) %>% count() %>% spread(TowardsTelo,n))
View(sglTelos %>% select(HasGGGTTA,HasCCCTAA,OrientStart,InsertSeq,RepeatType,everything()))


# Analysis
# - generally not clustered with other variants
# - >97% of SVs with telomeric  insert sequences are SGL type breakends
# - A small number of samples are highly enriched in telomeric single breakends,  in particular 84 samples have >5 observations.
# - Greatly enriched in Sarcoma, CNS,Skin and NET 
# - include subtypes for Skin,CNS and Sarcoma]
# - 42% of these 84 samples have ATRX drivers compared to 2% in rest of cohort
# - Test if enrichment is correlated with any other driver per cancer type??
# - Show orientation - how to summarise?
# - Telomeric inserts are evenly spread across chromosomes and are not typically clustered together.
# - XX% cause LOH CHARLES
# - Telomeric singles appear to cause genuine copy number change and are not generally concentrated 

View(sglTelos %>% mutate(ClusterType=ifelse(ResolvedType=='COMPLEX',
                                            ifelse(ClusterCount==1,'1',ifelse(ClusterCount<=3,'2-3','4+')))) 
     %>% group_by(ClusterType) %>% count)

print(ggplot(sglTelos %>% filter(ResolvedType!='DUP_BE'&ResolvedType!='LOW_VAF') %>% mutate(ClusterType=ifelse(ResolvedType=='COMPLEX',
                                                    ifelse(ClusterCount<=3,'COMPLEX <= 3 SVs','COMLPEX 5+ SVs'),as.character(ResolvedType))) 
                                 %>% group_by(ClusterType) %>% count, 
             aes(x=reorder(ClusterType,-n), y=n))
      + geom_bar(stat = "identity", colour = "black")
      + labs(x='Cluster Type', y='# Telomeric SGLs',title='# Telomeric SGLs by Cluster Type and Size'))


teloSampleData = sglTelos %>% group_by(SampleId,CancerType,CancerSubtype) %>% summarise(TotalTeloSgls=n(),
                                                                                        TeloFacing=sum(TowardsTelo),
                                                                                        CentroFacing=sum(!TowardsTelo))

View(teloSampleData)
nrow(teloSampleData %>% filter(TotalTeloSgls>5))

teloSampleData = teloSampleData %>% mutate(TeloSglCount=ifelse(TotalTeloSgls<=1,'1',ifelse(TotalTeloSgls<=3,'2-3',ifelse(TotalTeloSgls<=7,'4-7','8+'))))

# certain cancer types and cancer subtypes have an excess of telomeric SGLs
cancerTypeSummary = teloSampleData %>% group_by(CancerType,CancerSubtype,TeloSglCount) %>% count %>% spread(TeloSglCount,n)
cancerTypeSummary[is.na(cancerTypeSummary)] = 0
cancerTypeSummary = cancerTypeSummary %>% mutate(Total=`1`+`2-3`+`4-7`+`8+`)
View(cancerTypeSummary)

# TS are spread evenly across chromosomes
sampleChrData = sglTelos %>% group_by(SampleId,ChrStart) %>% summarise(SglCount=n()) %>% group_by(SampleId) %>% summarise(Chromosomes=n(),
                                                                                                                      TotalTeloSGLs=sum(SglCount))
View(sampleChrData)

print(ggplot(data = sampleChrData, aes(x=TotalTeloSGLs, y=Chromosomes))
      + geom_point()
      + labs(title = "# of Telomeric SGLs per Sample vs # of Chromosomes with Telomeric SGLs"))



altTertDrivers = read.csv('~/logs/atrx_samples.csv')
View(altTertDrivers)
atrxSamples = altTertDrivers %>% filter(Gene=='ATRX')
daxxSamples = altTertDrivers %>% filter(Gene=='DAXX')
View(atrxSamples)
View(daxxSamples)
sglTelos = sglTelos %>% mutate(HasAtrxDriver=SampleId %in% atrxSamples$SampleId,
                               HasDaxxDriver=SampleId %in% daxxSamples$SampleId)

sampleDriverData = sglTelos %>% group_by(SampleId,HasAtrxDriver) %>% summarise(SglCount=n())

View(teloSampleData %>% group_by(TeloSglCount,HasAtrxDriver) %>% count %>% spread(HasAtrxDriver,n))

print(ggplot(data = teloSampleData, aes(x=TeloSglCount))
      + geom_bar()
      + facet_wrap(~HasAtrxDriver)
      + labs(title = "# of Telomeric SGLs per Sample by ATRX driver present"))

nrow(teloSampleData)
View(teloSampleData)

# test enrichment by all possible driver genes
driverGenes = read.csv('~/logs/driver_catalog.csv')
relevantDriverGenes = driverGenes %>% filter(SampleId %in% teloSampleData$SampleId)
nrow(relevantDriverGenes)

View(driverGenes %>% group_by(Gene) %>% count())
enrichedSamples = teloSampleData %>% filter(TotalTeloSgls>=6)
nrow(enrichedSamples)

View(driverGenes %>% filter(SampleId %in% enrichedSamples$SampleId) %>% group_by(Gene) %>% count %>% filter(n>=5))


teloSglSampleCounts = merge(sampleCancerAndSubTypes,teloSampleData %>% ungroup() %>% select(SampleId,TotalTeloSgls,TeloSglCount),by='SampleId',all.x=T)
teloSglSampleCounts[is.na(teloSglSampleCounts)] = 0 
View(teloSglSampleCounts)
teloSglSampleCounts = teloSglSampleCounts %>% mutate(HasAtrxDriver=SampleId %in% atrxSamples$SampleId)

atrxDriverTotals = teloSglSampleCounts %>% group_by(HasAtrxDriver) %>% summarise(AtrxDriverTotal=n())
View(atrxDriverTotals)
teloSglSampleCounts = merge(teloSglSampleCounts,atrxDriverTotals,by='HasAtrxDriver',all.x=T)
teloSglAtrxSummary = teloSglSampleCounts %>% group_by(TeloSglCount,HasAtrxDriver) %>% summarise(Count=n(),
                                                                                                Percent=round(n()/first(AtrxDriverTotal),4))
print(ggplot(data = teloSglAtrxSummary, aes(x=ifelse(HasAtrxDriver,'ATRX Driver','No ATRX Driver'),y=Percent,fill=TeloSglCount))
      + geom_bar(stat = "identity", colour = "black")
      + scale_fill_manual(values = tcColours)
      + labs(x='', y='% of Samples', title = "# of Telomeric SGLs per Sample by ATRX driver present"))

print(ggplot(data = teloSglSampleCounts, aes(x=TeloSglCount))
      + geom_bar()
      + facet_wrap(~HasAtrxDriver)
      + scale_y_log10()
      + labs(x='Count of Telomeric SGLs', y='Samples', title = "# of Telomeric SGLs per Sample by ATRX driver present"))


teloSglCancerSummary = teloSglSampleCounts %>% group_by(CancerType,CancerSampleCount,TeloSglCount) %>% count() # %>% spread(TeloSglCount,n)
teloSglCancerSummary = teloSglCancerSummary %>% mutate(Percent=round(n/CancerSampleCount,4),
                                                       CancerTypeLabel=sprintf("%s (%d)",CancerType,CancerSampleCount))
View(teloSglCancerSummary)

print(ggplot(teloSglCancerSummary, aes(x=reorder(CancerTypeLabel,-CancerSampleCount), y=Percent, fill=TeloSglCount))
      + geom_bar(stat = "identity", colour = "black")
      + scale_fill_manual(values = tcColours)
      + labs(x='Cancer Type', y='Percent of Samples')
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10)))


teloSglCancerSubtypeSummary = teloSglSampleCounts %>% filter(CancerType=='Bone/Soft tissue') %>% 
  group_by(CancerSubtype,CancerSubtypeSampleCount,TeloSglCount) %>% count() # %>% spread(TeloSglCount,n)

teloSglCancerSubtypeSummary = teloSglCancerSubtypeSummary %>% mutate(Percent=round(n/CancerSubtypeSampleCount,4),
                                                                     CancerSubtypeLabel=sprintf("%s (%d)",CancerSubtype,CancerSubtypeSampleCount))
teloSglCancerSubtypeSummary = teloSglCancerSubtypeSummary %>% filter(CancerSubtype!='')
View(teloSglCancerSubtypeSummary)

print(ggplot(teloSglCancerSubtypeSummary, aes(x=reorder(CancerSubtypeLabel,-CancerSubtypeSampleCount), y=Percent, fill=TeloSglCount))
      + geom_bar(stat = "identity", colour = "black")
      + scale_fill_manual(values = tcColours)
      + labs(x='Cancer Subtype', y='Percent of Samples')
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10)))
      
            
calc_fisher_et(sampleCancerTypes,
               driverGenes %>% filter(Gene=='ATRX'),
               enrichedSamples,
               'ATRX','TeloSgl')
  
applicableGenes = relevantDriverGenes %>% group_by(Gene) %>% count

for(gene in applicableGenes$Gene)
{
  prob = calc_fisher_et(sampleCancerTypes, driverGenes %>% filter(Gene==gene), enrichedSamples, gene,'TeloSgl',F)
  
  if(prob<0.001)
  {
    print(sprintf("Gene(%s) prob=%f",gene,prob))
  }
}



teloSampleData = teloSampleData %>% mutate(HasAtrxDriver=SampleId %in% atrxSamples$SampleId)
View(teloSampleData)

print(ggplot(data = sampleDriverData, aes(x=SampleId, y=SglCount))
      + geom_point()
      + facet_wrap(~HasAtrxDriver)
      + labs(title = "# of Telomeric SGLs per Sample vs # of Chromosomes"))


#####
## facing telomere or centromere
View(sglTelos %>% group_by(TowardsTelo,ArmStart,OrientStart) %>% count() %>% spread(TowardsTelo,n))
View(sglTelos %>% group_by(TowardsTelo,ChrStart,ArmStart) %>% count() %>% spread(TowardsTelo,n))
View(sglTelos %>% group_by(TowardsTelo,ResolvedType) %>% count() %>% spread(TowardsTelo,n))
View(sglTelos %>% group_by(TowardsTelo,ClusterSize=2**round(log(ClusterCount,2))) %>% count() %>% spread(TowardsTelo,n))


sglTelos = sglTelos %>% mutate(HasAtrxDriver=SampleId %in% atrxSamples$SampleId)
sglTelos = merge(sglTelos,centroLengths,by.x='ChrStart',by.y='Chromosome',all.x=T)
sglTelos = sglTelos %>% mutate(TeloDistance=ifelse(ArmStart=='P',PosStart,Length-PosStart),
                               TeloDistBucket=2**round(log(TeloDistance,2)))
# View(centroLengths)
View(sglTelos %>% group_by(TowardsTelo,TeloDistBucket) %>% count() %>% spread(TowardsTelo,n))

sglTelos = sglTelos %>% mutate(SglFacesTelo=(ArmStart=='P')==(OrientStart==1))

View(sglTelos %>% group_by(TowardsTelo,SglFacesTelo) %>% count() %>% spread(TowardsTelo,n))
View(sglTelos %>% group_by(TowardsTelo,TeloDistBucket,SglFacesTelo) %>% count() %>% spread(TowardsTelo,n))
View(sglTelos %>% group_by(TowardsTelo,TeloDistBucket,ResolvedType) %>% count() %>% spread(TowardsTelo,n))

View(sglTelos)
View(sglTelos %>% group_by(SampleId,CancerType,HasAtrxDriver,TowardsTelo) %>% count() %>% spread(TowardsTelo,n))
View(sglTelos %>% group_by(SampleId,CancerType,HasDaxxDriver,TowardsTelo) %>% count() %>% spread(TowardsTelo,n))

View(sglTelos %>% group_by(SampleId,CancerType,ChrStart,TowardsTelo) %>% count() %>% spread(ChrStart,n))

View(sglTelos %>% group_by(SampleId,TowardsTelo) %>% count() %>% spread(TowardsTelo,n,fill=0) %>% mutate(total=pmin(10,`TRUE`+`FALSE`)) %>% 
       group_by(total) %>% summarise(sum(`TRUE`),sum(`FALSE`)))


## Summary figures for telomeric SGLs
colnames(sglTelos)

sglTelos = merge(sglTelos,teloSampleData %>% select(SampleId,TeloSglCount,CancerType,CancerSubtype),by='SampleId',all.x=T)

# facing telomere or not by enrichment 
View(sglTelos %>% group_by(TeloSglCount,TowardsTelo) %>% count %>% spread(TowardsTelo,n))
View(sglTelos %>% group_by(TeloSglCount,TowardsTelo,HasDriver=(HasAtrxDriver|HasDaxxDriver)) %>% count %>% spread(TowardsTelo,n))

View(sglTelos %>% group_by(TeloSglCount,TowardsTelo,HasDriver=(HasAtrxDriver|HasDaxxDriver)) %>% summarise(SglCount=n(),Samples=n_distinct(SampleId)))
View(sglTelos %>% group_by(TeloSglCount,HasDriver=(HasAtrxDriver|HasDaxxDriver)) %>% summarise(SglCount=n(),Samples=n_distinct(SampleId)))

print(ggplot(data = sglTelos %>% group_by(TeloSglCount,Direction=ifelse(TowardsTelo,'Faces Telomere','Faces Centromere')) %>% count, 
             aes(x=TeloSglCount,y=n,fill=Direction))
      + geom_bar(stat = "identity", colour = "black")
      + scale_fill_manual(values = tcColours)
      + labs(x='', y='Telomeric SGL count', title = "Telomeric SGLs: # Facing T or C by # per sample"))

# distance to telomere - makes no difference
View(sglTelos %>% group_by(TeloDistBucket,TowardsTelo) %>% count %>% spread(TowardsTelo,n) %>% mutate(FacesTeloPercent=`TRUE`/(`TRUE`+`FALSE`)))


# linking to LOH and/or drivers


# looking into samples with enrichment but not too many SVs
View(teloSampleData %>% filter(TotalTeloSgls>=10&SvCount<200))

nrow(sglTelos)



## telomeric SGLs forming an LOH or exhibiting loss out to the telomere
View(svData %>% filter(SampleId=='CPCT02010035T'&(ChrStart==6|ChrEnd==6)))
svDataWithTeloSgls = svData %>% filter(SampleId %in% teloSampleData$SampleId)
nrow(svDataWithTeloSgls)

pArmSvs = svDataWithTeloSgls %>% filter(ArmStart=='P'|(Type=='BND'&ArmEnd=='P')) %>% select(SampleId,Id,Type,ChrStart,ChrEnd,ArmStart,ArmEnd,PosStart,PosEnd,OrientStart,OrientEnd) %>% 
  mutate(Chr=ifelse(Type!='BND'|ArmStart=='P',as.character(ChrStart),as.character(ChrEnd)),Arm='P',
         Position=ifelse(Type!='BND'|ArmStart=='P',PosStart,PosEnd),
         Orient=ifelse(Type!='BND'|ArmStart=='P',OrientStart,OrientEnd)) %>%
  arrange(SampleId,Chr,Position)

#View(pArmSvs %>% filter(SampleId=='CPCT02060181T'))
#View(pArmSvs %>% filter(Type=='BND'&ChrStart!=Chr))
outerPSv = pArmSvs %>% group_by(SampleId,Chr,Arm) %>% summarise(Id=first(Id),Type=first(Type),Position=first(Position),Orient=first(Orient)) %>% ungroup()
outerPSgls = outerPSv %>% filter(Type=='SGL'&Orient==-1)
View(outerPSgls)
View(outerPSv %>% filter(SampleId=='CPCT02060181T'))

qArmSvs = svDataWithTeloSgls %>% filter(ArmEnd=='Q'|((Type=='BND'|Type=='SGL'|Type=='INF')&ArmStart=='Q')) %>% select(SampleId,Id,Type,ChrStart,ChrEnd,ArmStart,ArmEnd,PosStart,PosEnd,OrientStart,OrientEnd) %>% 
  mutate(Chr=ifelse((Type!='BND'&Type!='SGL'&Type!='INF')|ArmEnd=='Q',as.character(ChrEnd),as.character(ChrStart)),Arm='Q',
         Position=ifelse((Type!='BND'&Type!='SGL'&Type!='INF')|ArmEnd=='Q',PosEnd,PosStart),
         Orient=ifelse((Type!='BND'&Type!='SGL'&Type!='INF')|ArmEnd=='Q',OrientEnd,OrientStart)) %>%
  arrange(SampleId,Chr,-Position)

View(qArmSvs %>% filter(SampleId=='CPCT02060181T'))
View(qArmSvs %>% filter(Type=='BND'&ChrStart!=Chr))

outerQSv = qArmSvs %>% group_by(SampleId,Chr,Arm) %>% summarise(Id=first(Id),Type=first(Type),Position=first(Position),Orient=first(Orient)) %>% ungroup()
View(outerQSv)
outerQSgls = outerQSv %>% filter(Type=='SGL'&Orient==1)
View(outerQSgls)

sglTelosPlus = rbind(merge(sglTelos %>% filter(ArmStart=='P'),outerPSgls %>% select(SampleId,Id,Position),by=c('SampleId','Id'),all.x=T),
                     merge(sglTelos %>% filter(ArmStart=='Q'),outerQSgls %>% select(SampleId,Id,Position),by=c('SampleId','Id'),all.x=T))

View(sglTelosPlus %>% filter(!is.na(Position)))
View(sglTelosPlus %>% filter(!is.na(Position)) %>% group_by(SingleCluster=ClusterCount==1,TowardsTelo) %>% count)
View(sglTelosPlus %>% filter(!is.na(Position)&ClusterCount==1) %>% group_by(TeloDistBucket,TowardsTelo) %>% count %>% spread(TowardsTelo,n))


lohEvents = read.csv('~/data/sv/CN_LOH_EVENTS.csv')
teloLohs = lohEvents %>% filter((SegStart=='TELOMERE'&SegEnd=='SGL')|(SegStart=='SGL'&SegEnd=='TELOMERE'))
View(teloLohs)

centroLohs = centroLohs %>% mutate(LengthBucket=2**round(log(Length,2)),
                                   Arm=ifelse(SegStart=='CENTROMERE','Q','P'),
                                   ChrArm=paste(Chromosome,Arm,sep='_'),
                                   WholeArm=(SegStart=='TELOMERE'|SegEnd=='TELOMERE'))




# telomeric SGLs per sample and total sample SV counts
nrow(svData)
sampleSvCounts = svData %>% group_by(SampleId) %>% summarise(SvCount=n())
View(sampleSvCounts)

teloSampleData = merge(teloSampleData,sampleSvCounts,by='SampleId',all.x=T)

View(teloSampleData %>% filter(SvCount<100))

lowSvSglTeloSamples = teloSampleData %>% filter(SvCount<200&TeloSglCount==1)
View(lowSvSglTeloSamples)
View(lowSvSglTeloSamples %>% group_by(CancerType,CancerSubtype) %>% count)

lowSstsSgls = sglTelos %>% filter(SampleId %in% lowSvSglTeloSamples$SampleId)
View(lowSstsSgls %>% select(SampleId,Id,TeloDistance,TowardsTelo,Ploidy,ClusterCount,ResolvedType,ChainCount,DBLenStart,NearestLen,everything()))


View(lowSstsSgls %>% group_by(TeloDistBucket,TowardsTelo) %>% count() %>% spread(TowardsTelo,n))

View(teloSglSampleCounts %>% filter(CancerSubtype=='Leiomyosarcoma'))
leisSamples = teloSglSampleCounts %>% filter(CancerSubtype=='Leiomyosarcoma'&TeloSglCount=='8+')
View(leisSamples)

leisSgls = sglTelos %>% filter(SampleId %in% leisSamples$SampleId)
nrow(leisSgls)
View(leisSgls)

View(leisSgls %>% filter(HasAtrxDriver) %>% group_by(SampleId,CancerType,CancerSubtype) %>% count() %>% group_by(CancerType,CancerSubtype) %>% count)


View(leisSgls %>% group_by(SampleId,ChrStart) %>% count())
View(leisSgls %>% group_by(SampleId,ChrStart,ClusterId,ClusterCount) %>% count())
# View(leisSgls %>% group_by(SampleId,ClusterId,ClusterCount) %>% count() %>% group_by(ClusterCount) %>% summarise(Count=n(),


# about half are unclustered, 15% in cluster 2s, 5% in cluster 3s

# cluster 1s
View(leisSgls %>% filter(ClusterCount==1) %>% select(SampleId,ClusterId,ChrStart,TeloDistance,TowardsTelo,everything()))
View(leisSgls %>% filter(ClusterCount==1) %>% group_by(HasAtrxDriver,ChrStart,TeloDistance,TowardsTelo) %>% count())
View(sglTelos %>% filter(ClusterCount==1) %>% group_by(TowardsTelo,Enriched=SampleId %in% leisSamples$SampleId) %>% count() %>% spread(TowardsTelo,n))

# most SGLs close to the telomere face the centromere, whereas >10M becomes more even (2:1)
View(leisSgls %>% filter(ClusterCount==1) %>% group_by(TeloDistBucket,TowardsTelo) %>% count() %>% spread(TowardsTelo,n))
View(sglTelos %>% filter(ClusterCount==1) %>% group_by(TeloDistBucket,TowardsTelo) %>% count() %>% spread(TowardsTelo,n))

View(leisSgls %>% filter(ClusterCount==1) %>% group_by(HasAtrxDriver,TeloDistBucket,TowardsTelo) %>% count() %>% spread(TowardsTelo,n))
View(sglTelos %>% filter(ClusterCount==1) %>% group_by(HasAtrxDriver,TeloDistBucket,TowardsTelo) %>% count() %>% spread(TowardsTelo,n))

# no enrichment by chromosome
View(leisSgls %>% filter(ClusterCount==1) %>% group_by(ChrStart,TowardsTelo) %>% count() %>% spread(TowardsTelo,n))

# cluster 2s

# 20 are pairs of SGLs with a possible telomereic insertion or vice versa and there are 70 such clusters in the full telomeric SGL dataset
View(leisSgls %>% filter(ClusterDesc=='SGL=2') %>% group_by(SampleId,ClusterId) %>% summarise(SvCount=n(),
                                                                                              SvProximity=abs(first(PosStart)-last(PosStart)),
                                                                                              TeloDistance=min(TeloDistance),
                                                                                              TowardsTeloCount=sum(TowardsTelo)) %>% filter(SvCount==2))

View(leisSgls %>% filter(ClusterCount==2) %>% group_by(ClusterDesc,TowardsTelo) %>% count())

# larger clusters - 26% have a SGL facing the telomere, more than the 16% of cluster 1s which do
View(leisSgls %>% group_by(SampleId,ClusterId,ClusterSize=ifelse(ClusterCount<=3,ClusterCount,4)) %>% 
       summarise(TeloDistBucket=mean(TeloDistBucket),TowardsTelo=sum(TowardsTelo)/n()) %>%
       group_by(ClusterSize) %>% summarise(Clusters=n(),
                                           TeloDistBucket=mean(TeloDistBucket),
                                           TowardsTelo=mean(TowardsTelo)))


View(clusters)

sglTelos = sglTelos %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))

sglTeloClusters = merge(sglTelos,clusters %>% select(SampleClusterId,ChainCount,Consistency,UnchainedSVs),by='SampleClusterId',all.x=T)
View(sglTeloClusters)

View(sglTeloClusters %>% filter(ChainCount.y==1&UnchainedSVs==0&ClusterCount>=3) %>% select(SampleId,ClusterId,TowardsTelo,TeloDistance,everything()))


View(sglTelos %>% filter(SampleId=='CPCT02020691T'))






# Samples with 1 or few telomeric CN gain vs SGLs going to telomeres




## DEBUG



View(sglCentros %>% filter(ClusterCount>3&CentroType=='Centro') %>% group_by(BFB,CentroDistBucket,ClusterCount=2**round(log(ClusterCount,2))) 
     %>% count %>% spread(ClusterCount,n))

View(sglCentros %>% filter(ClusterCount>3) %>% group_by(ChrStart,round(CentroDistPerc,1)) %>% count %>% spread(ChrStart,n))
View(sglCentros %>% filter(ClusterCount>3) %>% group_by(BFB,round(CentroDistPerc,1)) %>% count %>% spread(BFB,n))
View(sglCentros %>% group_by(BFB,SglCentroCount=2**round(log(SglCentroCount,2)),round(CentroDistPerc,1)) %>% count %>% spread(SglCentroCount,n))

View(sglCentros %>% filter(ClusterCount>3&CentroType=='Centro') %>% group_by(BFB,CentroDistBucket) %>% count %>% spread(BFB,n))
View(sglCentros %>% filter(ClusterCount>3&CentroType=='Centro') %>% group_by(HasDB,CentroDistBucket) %>% count %>% spread(HasDB,n))
View(sglCentros %>% filter(ClusterCount>3&CentroType=='Centro') %>% group_by(SglCentroCount=2**round(log(SglCentroCount,2)),BFB,HasDB) %>% count %>% spread(SglCentroCount,n))


print(min(14, log2(3e9)/2 - 1))

# bfbClusters = clusters %>% filter(Foldbacks>2&ResolvedType=='COMPLEX')
clusters = clusters %>% mutate(SampleClusterId=paste(SampleId,ClusterId,sep='_'))

sglCentroClusters = clusters %>% filter(SampleClusterId %in% sglCentros$SampleClusterId)
nrow(sglCentroClusters)

View(sglCentroClusters %>% group_by(ResolvedType,ClusterSize=2**round(log(ClusterCount,2))) %>% count() %>% spread(ResolvedType,n))

scClusterCounts = sglCentros %>% group_by(SampleClusterId) %>% summarise(SglCentroCount=n())
sglCentroClusters = merge(sglCentroClusters,scClusterCounts,by='SampleClusterId',all.x=T)
sglCentroClusters = sglCentroClusters %>% mutate(HasFoldbacks=Foldbacks>0,
                                                 BFB=HasFoldbacks&ResolvedType=='COMPLEX')

View(sglCentroClusters %>% select(SampleClusterId,ClusterCount,ResolvedType,ClusterDesc,SglCentroCount,everything()))

View(sglCentroClusters %>% filter(ClusterCount>=3) %>% mutate(SglCentroPerc=round(SglCentroCount/ClusterCount,1),
                                                              HasFoldbacks=Foldbacks>0) %>% 
       select(SampleClusterId,ClusterCount,SglCentroPerc,HasFoldbacks,ClusterDesc,SglCentroCount,everything()))

View(sglCentroClusters %>% filter(ClusterCount>=3) %>% mutate(SglCentroPerc=round(SglCentroCount/ClusterCount,1),
                                                              HasFoldbacks=Foldbacks>0) %>% 
       group_by(SglCentroPerc,HasFoldbacks,ClusterSize=2**round(log(ClusterCount,2))) %>% count %>% spread(SglCentroPerc,n))


lohEvents = read.csv('~/data/sv/CN_LOH_EVENTS.csv')
centroLohs = lohEvents %>% filter(SegStart=='CENTROMERE'|SegEnd=='CENTROMERE')
View(centroLohs)

centroLohs = centroLohs %>% mutate(LengthBucket=2**round(log(Length,2)),
                                   Arm=ifelse(SegStart=='CENTROMERE','Q','P'),
                                   ChrArm=paste(Chromosome,Arm,sep='_'),
                                   WholeArm=(SegStart=='TELOMERE'|SegEnd=='TELOMERE'))

View(centroLohs %>% group_by(Chromosome,Arm,WholeArm) %>% count() %>% spread(WholeArm,n))
View(centroLohs %>% group_by(Chromosome,Arm,LengthBucket) %>% count())

plot_length_facetted(centroLohs,'!WholeArm','LengthBucket,ChrArm','LengthBucket','ChrArm','LOH Length by ChrArm',F)


## Gather foldback data

fbStart = svData %>% filter(FoldbackLenStart>=0) %>% select(SampleId,Id,ClusterId,Chr=ChrStart,Arm=ArmStart,FoldbackLength=FoldbackLenStart,OtherId=FoldbackLnkStart)
fbEnd = svData %>% filter(FoldbackLenEnd>=0) %>% select(SampleId,Id,ClusterId,Chr=ChrEnd,Arm=ArmEnd,FoldbackLength=FoldbackLenEnd,OtherId=FoldbackLnkEnd)
foldbacks = rbind(fbStart,fbEnd)
View(foldbacks)
foldbacks = foldbacks %>% mutate(IsChained=(OtherId!=Id),
                                 SingleBreakend=(OtherId==Id&FoldbackLength==0),
                                 FoldbackId=ifelse(Id<OtherId,Id,OtherId))

# unique foldback data
View(foldbacks %>% group_by(SampleId,ClusterId,FoldbackId) %>% summarise(Chr=first(Chr),Arm=first(Arm),FoldbackLength=first(FoldbackLength)))

fbArmData = foldbacks %>% group_by(SampleId,Chr,Arm) %>% summarise(FoldbackCount=sum(ifelse(IsChained,0.5,1)))
fbArmData = fbArmData %>% ungroup() %>% mutate(FoldbackCount=ifelse(FoldbackCount<1,1,round(FoldbackCount)))

foldbackSummary = foldbacks %>% group_by(SampleId,Chr,Arm) %>% summarise(FoldbackCount=n()/2)
View(foldbackSummary)



sglCentroFoldbackArms = merge(sglCentroArms,foldbackSummary,by=c('SampleId','Chr','Arm'),all=T)
sglCentroFoldbackArms = sglCentroFoldbackArms %>% mutate(FoldbackCount=ifelse(is.na(FoldbackCount),0,FoldbackCount),
                                                         SglCentroCount=ifelse(is.na(SglCentroCount),0,SglCentroCount),
                                                         HasFoldbacks=FoldbackCount>0,
                                                         HasSglCentro=SglCentroCount>0)

chr1Gains = teloCentroCn %>% filter(Chromosome==1&QCentroGain)
View(chr1Gains)
View(chr1Gains %>% group_by(CentroGainGroup2=ifelse(CentroGain<3.5,0.2*round(CentroGain/0.2),3.5)) %>% count())
View(teloCentroCn %>% filter(CentroGain>0.5) %>%
       group_by(ChrArm,CentroGainGroup2=ifelse(CentroGain<3.5,0.2*round(CentroGain/0.2),3.5)) %>% count() %>% spread(ChrArm,n))


View(sglCentroFoldbackArms)

totalChr1Arms = nrow(svData %>% filter(SampleId %in% chr1Gains$SampleId) %>% group_by(SampleId,ChrStart,ArmStart) %>% count())
print(totalChr1Arms)

chr1SglCentroFbArms = sglCentroFoldbackArms %>% filter(SampleId %in% chr1Gains$SampleId)
View(chr1SglCentroFbArms)

View(chr1SglCentroFbArms %>% group_by(HasFoldbacks,HasSglCentro) %>% count())

expected = nrow(chr1SglCentroFbArms %>% filter(HasFoldbacks)) / totalChr1Arms * nrow(chr1SglCentroFbArms %>% filter(HasSglCentro))
print(expected)

print(ggplot(data = chr1SglCentroFbArms, aes(x=SglCentroCount,y=FoldbackCount))
      + geom_point(position='jitter')
      + labs(title = "Foldback Count vs SGL-centro Count with Chr 1 Gain"))



View(sglCentroFoldbackArms %>% group_by(HasFoldbacks,HasSglCentro) %>% count())
totalNone = nrow(sglCentroFoldbackArms %>% filter(HasFoldbacks|HasSglCentro))
print(totalNone)

total = nrow(svData %>% group_by(SampleId,ChrStart,ArmStart) %>% count())
print(total)

withFb = nrow(sglCentroFoldbackArms %>% filter(HasFoldbacks))
withSglCentro = nrow(sglCentroFoldbackArms %>% filter(HasSglCentro))
withFbWithSglCentro = 5331
withFbNoSglCentro = 14598
noFbWithSglCentro = 8875
noFbNoSglCentro = total-withFbWithSglCentro-withFbNoSglCentro-noFbWithSglCentro
print(noFbNoSglCentro)
expected = withFb/total * withSglCentro
print(expected)
fishMatrix = rbind(c(withFbWithSglCentro,noFbWithSglCentro), c(withFbNoSglCentro,noFbNoSglCentro))

fetProb = fisher.test(fishMatrix, alternative="less")$p.value
print(fetProb)

print(ggplot(data = sglCentroFoldbackArms, aes(x=SglCentroCount,y=FoldbackCount))
      + geom_point()
      + labs(title = "Foldback Count vs SGL-centromere Count"))


teloCentroSampleSummary2 = teloCentroCn %>% filter(PCentroGain|QCentroGain) %>% 
  group_by(CancerType,SampleId,Chromosome,CentroGainGroup,PCentroGain,QCentroGain) %>% count()

View(teloCentroSampleSummary2)
teloCentroSampleSummary2 = teloCentroSampleSummary2 %>% mutate(PCentroGainGrp=ifelse(PCentroGain,CentroGainGroup,-CentroGainGroup),
                                                               QCentroGainGrp=ifelse(QCentroGain,CentroGainGroup,-CentroGainGroup))

View(teloCentroSampleSummary2)

foldbackCentroData1 = merge(foldbackSummary %>% filter(Arm=='P'),
                            teloCentroSampleSummary2 %>% ungroup() %>% 
                              mutate(Chr=Chromosome,Arm='P',CentroGain=PCentroGainGrp) %>% select(CancerType,SampleId,Chr,Arm,CentroGain),
                            by=c('SampleId','Chr','Arm'),all=T)

foldbackCentroData2 = merge(foldbackSummary %>% filter(Arm=='Q'),
                            teloCentroSampleSummary2 %>% ungroup() %>% 
                              mutate(Chr=Chromosome,Arm='Q',CentroGain=QCentroGainGrp) %>% select(CancerType,SampleId,Chr,Arm,CentroGain),
                            by=c('SampleId','Chr','Arm'),all=T)

foldbackCentroData = rbind(foldbackCentroData1,foldbackCentroData2)
foldbackCentroData = foldbackCentroData %>% mutate(FoldbackCount=ifelse(is.na(FoldbackCount),0,FoldbackCount),
                                                   CentroGain=ifelse(is.na(CentroGain),0,CentroGain),
                                                   CentroGainGrp=as.character(CentroGain))

# foldbackCentroData = foldbackCentroData %>% mutate(CentroGainGrp=as.character(CentroGain))
View(foldbackCentroData)
View(foldbackCentroData %>% group_by(CentroGain) %>% summarise(FoldbackCount=sum(FoldbackCount)))

ggplot(foldbackCentroData %>% group_by(CentroGain) %>% summarise(FoldbackCount=sum(FoldbackCount)), aes(CentroGain,FoldbackCount)) +
  geom_violin(scale="area",fill="#6baed6")

ggplot(foldbackCentroData, aes(CentroGainGrp,FoldbackCount)) +
  geom_violin(scale="area",fill="#6baed6")

print(ggplot(data = combinedSummary %>% filter(BfbClusters>0), aes(x=SglArmCount,y=ArmGainCount))
      + geom_point()
      + labs(title = "ArmGainCount vs SGL-centromere Count"))

print(ggplot(data = combinedSummary, aes(x=SglCentroCount,y=ArmGainCount))
      + geom_point()
      + labs(title = "ArmGainCount vs SGL-centromere Count"))


bfbClusters = clusters %>% filter(Foldbacks>0&ResolvedType=='COMPLEX')

bfbClusters = bfbClusters %>% mutate(SingleChain=FullyChained=='true'&ChainCount==1,
                                     Resolved=SingleChain&Consistency==0,
                                     HasSgls=(SglCount+InfCount>0))

bfbClusters = merge(bfbClusters,sglCentroClusters,by=c('SampleId','ClusterId'),all.x=T)

View(bfbClusters %>% group_by(SingleChain,NoSgls,Resolved) %>% count())
View(bfbClusters %>% group_by(SingleChain,HasSgls,Resolved,HasSglCentro=!is.na(SglCentroCount)) %>% count())
View(bfbClusters)
