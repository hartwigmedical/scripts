


## Driver results
driverResults = read.csv('~/data/cup/CUP.DRIVERS.csv')
View(driverResults %>% filter(MatchedCancerType!='NONE'))

driverResults = driverResults %>% filter(MatchedCancerType!='NONE')
driverResults = merge(driverResults,sampleCancerTypes %>% select(SampleId,CancerType=MajorCancerType),by='SampleId',all.x=T)

driverMatch = driverResults %>% arrange(SampleId,-Probability) %>% group_by(SampleId,CancerType) %>% 
  summarise(Matches=n(),
            TopCancerType=first(MatchedCancerType),
            TopProb=first(Probability),
            ProbSum=sum(Probability),
            DriverGenes=first(DriverGenes))

driverMatch = driverMatch %>% mutate(Matched=as.character(CancerType)==as.character(TopCancerType))

View(driverMatch %>% filter(CancerType!='Unknown') %>% group_by(TopProb=round(TopProb,1),Matched,Matches=pmin(Matches,10)) %>% count %>% spread(Matches,n,fill=0))

View(driverMatch %>% filter(CancerType=='Unknown') %>% group_by(TopProb=round(TopProb,1),Matches=pmin(Matches,10)) %>% count %>% spread(Matches,n,fill=0))

View(driverResults %>% arrange(SampleId,-Probability) %>% group_by(SampleId,CancerType) %>% 
       summarise(Matches=n(),
                 TopCancerType=first(MatchedCancerType),
                 TopProb=first(Probability),
                 ProbSum=sum(Probability),
                 DriverGenes=first(DriverGenes)))


View(allDrivers %>% filter(Gene %in% c('TERT','MEN1','DAXX','BRAF')) %>% group_by(SampleId,CancerType) %>% summarise(Drivers=n()))

View(allDrivers)





## Signatures
# adding a specific sample
specSampleCounts = read.csv('~/data/cup/snv_test_sample_counts.csv')
write.csv(specSampleCounts %>% select(WIDE01010879T),'~/data/cup/snv_test_sample_matrix_data.csv',row.names = F,quote = F)
View(specSampleCounts)

signatures = as.matrix(read.csv(sigsFile), stringsAsFactors=F)

# manual fit in R
sigFitResult = fit_to_signatures(specSampleCounts %>% select(WIDE01010879T), as.matrix(cosmicSigs, stringsAsFactors=F))
sigFitContributions = as.data.frame(sigFitResult$contribution)
sampleTotal = sum(specSampleCounts$WIDE01010879T)
sigNames = rownames(sigFitContributions)
sigFitContributions = cbind(sigNames,sigFitContributions)
colnames(sigFitContributions) = c('Signature','SigContribution')
sigFitContributions = sigFitContributions %>% mutate(SigPercent=round(SigContribution/sampleTotal,4),
                                                     SigContribution = round(SigContribution,4))

View(sigFitContributions)

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