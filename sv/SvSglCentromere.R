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




