
#####
## Telomeric SGLs

sampleCancerAndSubTypes = load_cancer_types('~/data/hpc_sample_cancer_types.csv',T,10)

#View(egSvData %>% filter(Type=='SGL'))
# sglTelos = egSvData %>% filter(Type=='SGL'&(RepeatClass=='Simple_repeat'&(RepeatType=='(TTAGGG)n'|RepeatType=='(CCCTAA)n')))
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
