

## Create sample cancer type data
sampleCancerTypes = read.csv('~/data/cup/sample_cancer_types_20200727.csv')
nrow(sampleCancerTypes) # 4609
View(sampleCancerTypes)

## write cohort sampleIds
write.csv(sampleCancerTypes %>% select(SampleId),'~/data/cup/sample_ids_4610.csv',row.names = F,quote = F)

## De-deup for reference file creation
sampleDedupSummary = sampleCancerTypes %>% group_by(HmfPatientId) %>% summarise(Samples=n(),SampleId=first(SampleId))
nrow(sampleDedupSummary) # 4180
View(sampleDedupSummary)

sampleCancerTypesDD = sampleCancerTypes %>% filter(SampleId %in% sampleDedupSummary$SampleId)
View(sampleCancerTypesDD)

majorCancerTypes = c('Bone/Soft tissue','Neuroendocrine','Biliary','Uterus','Head and neck','Breast','Lung','Colon/Rectum','Prostate','Skin','Urinary tract','Kidney','Ovary',
                     'Esophagus','Pancreas','Nervous system','Liver','Stomach','Mesothelioma','Lymphoid','Anus','Thyroid','Unknown')


sampleCancerTypes = sampleCancerTypes %>% mutate(MajorCancerType=ifelse(CancerType %in% majorCancerTypes,as.character(CancerType),'Other'))
View(sampleCancerTypes)
View(sampleCancerTypes %>% group_by(MajorCancerType) %>% count)

cancerSampleCounts = sampleCancerTypes %>% group_by(MajorCancerType) %>% summarise(SampleCount=n())
View(cancerSampleCounts)

write.csv(sampleCancerTypes %>% select(SampleId,CancerType=MajorCancerType,CancerSubtype),'~/data/cup/cup_sample_data_4610.csv',row.names = F,quote = F)

# strip bucket name from sample counts
snvSampleCounts = read.csv('~/data/cup/snv_sample_counts_4610.csv')
ncol(snvSampleCounts)
snvSampleCounts = snvSampleCounts %>% select(-BucketName,-pancreatic.or.biliary.tract)
write.csv(snvSampleCounts,'~/data/cup/snv_matrix_data_4610.csv',row.names = F,quote = F)
snvSampleCounts = read.csv('~/data/cup/snv_sample_counts_4610.csv')

# adding a specific sample
specSampleCounts = read.csv('~/data/cup/snv_test_sample_counts.csv')

snvSampleCountsExtra = cbind(snvSampleCounts,specSampleCounts %>% select(WIDE01010879T))
write.csv(snvSampleCountsExtra,'~/data/cup/snv_matrix_data_4610_extras.csv',row.names = F,quote = F)

## SNV Cosine Similarity
snvCss = read.csv('~/data/cup/css_data.csv')
View(snvCss)

snvCss = snvCss %>% mutate(CancerSubtype=ifelse(is.na(CancerSubtype),'Other',as.character(CancerSubtype)))

## sunmary stats
snvCss = snvCss %>% arrange(SampleId,-MatchedCancerCss)

sampleSummary = snvCss %>% group_by(SampleId,CancerType,CancerSubtype,SnvCount) %>% 
  summarise(TopMatch=first(MatchedCancerType),
            TopMatch=first(MatchedCancerType),
            TopCss=first(MatchedCancerCss),
            TotalMatches=sum(MatchedCancerType!='NONE'))

sampleSummary = sampleSummary %>% mutate(TopIsSame=as.character(CancerType)==as.character(TopMatch))
View(sampleSummary)

View(sampleSummary %>% group_by(CancerType,CssGroup=ifelse(TopCss>0.2,2**round(log(TopCss,2)),0.2)) %>% count %>% spread(CssGroup,n,fill=0))

View(sampleSummary %>% filter(SnvCount<=5e3) %>% group_by(CancerType,CssGroup=ifelse(TopCss>0.2,2**round(log(TopCss,2)),0.2)) %>% count %>% spread(CssGroup,n,fill=0))

ctStats = sampleSummary %>% filter(CancerType!='Other') %>%
  mutate(SnvCountBucket=ifelse(SnvCount<5e3,'1_MB_<5K',ifelse(SnvCount<2e4,'2_MB_5-20K','3_MB_>20K'))) %>%
  group_by(CancerType,TopIsSame,SnvCountBucket) %>% count

View(ctStats %>% spread(SnvCountBucket,n,fill=0) %>% mutate(Total=`1_MB_<5K`+`2_MB_5-20K`+`3_MB_>20K`))



# rm(snvCssAll)
snvCssAll=snvCss
nrow(snvCss)

# snvCss = snvCss %>% filter(CSS>=0.95)

snvCss = merge(snvCss,sampleCancerTypes %>% select(SampleId,CancerType=MajorCancerType),by='SampleId',all.x=T)
snvCss = merge(snvCss,sampleCancerTypes %>% select(OtherSampleId=SampleId,OtherCancerType=MajorCancerType),by='OtherSampleId',all.x=T)
snvCss = snvCss %>% arrange(-CSS)
snvCss = snvCss %>% filter(!is.na(CancerType)&!is.na(OtherCancerType))
nrow(snvCss)
View(head(snvCss,1000) %>% mutate(Same=CancerType==OtherCancerType))

snvCss = snvCss %>% select(-CancerType1)
View(snvCss %>% group_by(SampleId,CancerType) %>% count)

snvCss = snvCss %>% mutate(CssBucket=ifelse(CSS>=0.99,'Css_99_100',ifelse(CSS>=0.98,'Css_98_99',ifelse(CSS>=0.97,'Css_97_98',ifelse(CSS>=0.96,'Css_96_97','Css_95_96')))))

snvSampleData = read.csv('~/data/cup/known_ct_sample_ids.csv')
snvSampleData = merge(snvSampleData,sampleCancerTypes %>% select(SampleId,MajorCancerType,CancerSubtype),by='SampleId',all.x=T)
snvSampleData = snvSampleData %>% mutate(MajorCancerType=ifelse(is.na(MajorCancerType),'Other',as.character()))
View(snvSampleData %>% filter(is.na(MajorCancerType)))
write.csv(snvSampleData %>% select(SampleId,CancerType=MajorCancerType,CancerSubtype),'~/data/cup/cup_sample_data_3700.csv',row.names = F,quote = F)

View(snvSampleData)
View(snvCss)
View(snvCss %>% filter(is.na(CancerType)))
View(snvCss %>% filter(is.na(OtherCancerType)))

sampleCssBuckets = snvCss %>% group_by(SampleId,CssBucket) %>% summarise(PanCancerCount=n())
View(sampleCssBuckets)

snvCssSummary = snvCss %>% group_by(SampleId,CancerType,SampleTotal,OtherCancerType,CssBucket) %>% count %>% spread(OtherCancerType,n,fill=0)

snvCssSummary = merge(snvCssSummary,sampleCssBuckets,by=c('SampleId','CssBucket'),all.x=T)
View(snvCssSummary %>% select(SampleId,CancerType,SampleTotal,CssBucket,PanCancerCount,everything()))

View(snvCssSummary %>% select(SampleId,CancerType,SampleTotal,CssBucket,PanCancerCount,everything()) %>% arrange(SampleId,-PanCancerCount))

sampleId='CPCT02020619T'
sampleCss = snvCss %>% filter(SampleId1==sampleId|SampleId2==sampleId)
View(sampleCss %>% arrange(-CSS))

# percentile summary
sampleCssSummary = sampleCss %>% 
  filter(CSS>0.9) %>%
  mutate(SampleId=sampleId,OtherSampleId=ifelse(SampleId1==sampleId,as.character(SampleId2),as.character(SampleId1)),
         CancerType=ifelse(SampleId1==sampleId,as.character(CancerType1),as.character(CancerType2)),
         OtherCancerType=ifelse(SampleId1==sampleId,as.character(CancerType2),as.character(CancerType1)),
         CssBucket=ifelse(CSS>=0.995,'Css_0.995',ifelse(CSS>=0.99,'Css_0.99',ifelse(CSS>=0.98,'Css_0.98',
                          ifelse(CSS>=0.97,'Css_0.97',ifelse(CSS>=0.96,'Css_0.96',ifelse(CSS>=0.95,'Css_0.95','Css_0.90'))))))) %>%
  select(SampleId,CancerType,OtherSampleId,OtherCancerType,CSS,CssBucket)

View(sampleCssSummary)
View(sampleCssSummary %>% group_by(CancerType,OtherCancerType,CssBucket) %>% count %>% spread(CssBucket,n,fill=0))

## possible mis-classifications
View(sampleCss %>% filter(CSS>0.99&Sample))


## Drivers
allDrivers = read.csv('~/data/cup/drivers_4600.csv')
allDrivers = allDrivers %>% filter(SampleId %in% sampleCancerTypes$SampleId)
View(allDrivers %>% group_by(SampleId) %>% count)

write.csv(allDrivers %>% filter(Driverlikelihood>=0.8) %>% select(SampleId,Gene,Driver),'~/data/cup/sample_driver_data.csv',row.names = F,quote = F)


allDrivers = merge(allDrivers,sampleCancerTypes %>% select(SampleId,CancerType=MajorCancerType,CancerSubtype),by='SampleId',all.x=T)
View(allDrivers)

driverCtGenes = allDrivers %>% filter(CancerType!='Unknown') %>% group_by(CancerType,Gene,Driver) %>% summarise(SampleCount=round(sum(Driverlikelihood))) %>% filter(SampleCount>0)
driverCtGenes = merge(driverCtGenes,cancerSampleCounts %>% select(CancerType=MajorCancerType,CancerSampleCount=SampleCount),by='CancerType',all.x=T)
driverCtGenes = driverCtGenes %>% mutate(SamplePerc=round(SampleCount/CancerSampleCount,6))
View(driverCtGenes)

driverCtGenes2 = allDrivers %>% filter(CancerType!='Unknown') %>% group_by(CancerType,SampleId,Gene) %>% summarise(DriverCount=sum(Driverlikelihood))
driverCtGenes2 = driverCtGenes2 %>% mutate(DriverCount=pmin(DriverCount,1))
driverCtGenes2 = driverCtGenes2 %>% group_by(CancerType,Gene) %>% summarise(SampleCount=sum(DriverCount)) %>% filter(SampleCount>0)
View(driverCtGenes2)
driverCtGenes2 = merge(driverCtGenes2,cancerSampleCounts %>% select(CancerType=MajorCancerType,CancerSampleCount=SampleCount),by='CancerType',all.x=T)
driverCtGenes2 = driverCtGenes2 %>% mutate(SamplePerc=round(SampleCount/CancerSampleCount,6))
View(driverCtGenes2)
driverCtGenes2 = driverCtGenes2 %>% mutate(Driver='ALL')

View(driverCtGenes2 %>% group_by(CancerType) %>% summarise(sum(SampleCount>0)))
View(driverCtGenes2 %>% group_by(Gene) %>% count)
View(driverCtGenes2 %>% group_by(CancerType,Gene) %>% summarise(SampleCount=n()) %>% group_by(CancerType) %>% summarise(GeneCount=n()))

driverSigs = rbind(driverCtGenes,driverCtGenes2)
View(driverSigs)

minCohortPrev = 0.0001;

View(driverSigs %>% filter(SamplePerc>=0.001) %>% filter(Driver=='ALL') %>% group_by(Gene) %>% count)

write.csv(driverSigs %>% filter(SamplePerc>=0.0001) %>% select(CancerType,Gene,Driver,SamplePerc),'~/data/cup/driver_prevalence.csv',row.names = F,quote = F)

write.csv(driverCtGenes2 %>% filter(SamplePerc>0.001) %>% 
          filter(Gene %in% c('TERT','CDKN2A','BRAF','NRAS','PTEN','TP53','NF1','FAT4','ARID2','APC','MAP2K1','RB1','KMT2B',
                             'DPYD','B2M','RAC1','PPP6C','IDH1','FAT1','PBRM1','PTPRD','SETD2','BAP1','CDH10','PIK3CA')) %>%
            select(CancerType,Gene,Driver,SamplePerc),'~/data/cup/driver_prev_test.csv',row.names = F,quote = F)

driverSigs2 = driverSigs %>% select(Gene,Driver,CancerType,SamplePerc) %>% spread(CancerType,SamplePerc,fill=0)

View(driverSigs2 %>% select(Gene,Driver,Breast) %>% filter(Driver=='ALL'&Breast>0.01))
View(driverSigs2 %>% select(Gene,Driver,Breast,Prostate) %>% filter(Driver=='ALL'&(Breast>0.01|Prostate>0.01)))

View(driverSigs2 %>% select(Gene,Driver,Skin,Lung,Breast,Prostate) %>% filter(Driver=='ALL'&(Skin>0.01|Lung>0.01|Breast>0.01|Prostate>0.01)))

driverSigs3 = driverSigs %>% mutate(SamplePerc=pmax(SamplePerc,0.0001)) %>% select(Gene,Driver,CancerType,SamplePerc) %>% spread(CancerType,SamplePerc,fill=0.0001)

View(driverSigs3 %>% filter(Driver=='ALL'))


write.csv(driverSigs2,'~/data/cup/driver_prev_sigs.csv',row.names = F,quote = F)

View(allDrivers %>% group_by(Category,Driver,Likelihood=0.2*round(Driverlikelihood/0.2)) %>% count)

View(allDrivers %>% group_by(Category,Driver,Gene) %>% summarise(SampleCount=n()) %>% group_by(Gene) %>% summarise(Types=n()))
View(allDrivers %>% group_by(Category,Driver,Gene) %>% summarise(SampleCount=n()))


## SNV Signatures
leastSqFit = read.csv('~/logs/snv_ut_sample_contribs.csv')
View(leastSqFit[,1:10])

expMaxFit = read.csv('~/logs/snv_ut_emax_sample_contribs.csv')
View(expMaxFit[,1:10])

cosmicSigs = read.csv('~/data/sigs/snv_cosmic_sigs.csv')
View(cosmicSigs)

cosmicSubset = cosmicSigs %>% select(Signature.2,Signature.4,Signature.6,Signature.7,Signature.13,Signature.17)
colnames(cosmicSubset) = c('Sig2','Sig4','Sig6','Sig7','Sig13','Sig17')
write.csv(cosmicSubset,'~/data/cup/cosmic_sig_subset.csv',row.names = F,quote = F)
