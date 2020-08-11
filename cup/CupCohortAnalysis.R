
## Sample Set
## write cohort sampleIds
sampleCancerTypes = read.csv('~/data/cup/sample_cancer_types_20200727.csv')
nrow(sampleCancerTypes) # 4609
View(sampleCancerTypes)

write.csv(sampleCancerTypes %>% select(SampleId),'~/data/cup/sample_ids_4610.csv',row.names = F,quote = F)


## Sample Results
cupSampleId = 'CPCT02010440T' # Skin
cupSampleId = 'CPCT02050268T' # Prostate
cupSampleId = 'CPCT02020492T' # Uterus with HPV
cupSampleResults = read.csv(paste('~/data/cup/samples/',cupSampleId,'.cup.data.csv',sep=''))
# colnames(cupSampleResults )
View(cupSampleResults)

View(cupSampleResults %>% group_by(SampleId,CancerType,Category,DataType,Value,RefCancerType) %>% summarise(CtValue=first(RefValue)) 
     %>% spread(RefCancerType,CtValue,fill=0) %>% select(SampleId,CancerType,Category,DataType,Value,Uterus,everything()))

View(cupSampleResults %>% filter(Category=='CLASSIFIER') %>% arrange(SampleId,DataType,-RefValue) %>%
       group_by(SampleId,CancerType,Category,DataType,Value) %>% summarise(TopRefCancerType=first(RefCancerType),
                                                                           TopRefValue=first(RefValue)))


cupCohortResults = read.csv('~/data/cup/samples/CUP.SAMPLE_DATA.csv')
nrow(cupCohortResults)
View(cupCohortResults)

cupClassData = cupCohortResults %>% filter(Category=='CLASSIFIER')
nrow(cupClassData)
View(cupClassData)

cupPredictions = cupClassData %>% arrange(SampleId,DataType,-RefCancerTypeValue) %>%
  group_by(SampleId,CancerType,DataType) %>% summarise(TopRefCancerType=first(RefCancerType),
                                                       TopRefValue=first(RefCancerTypeValue))

cupPredictions = cupPredictions %>% mutate(PredictionBucket=round(TopRefValue/0.1)*0.1,
                                           IsCorrect=as.character(CancerType)==as.character(TopRefCancerType))

View(cupPredictions %>% group_by(DataType,PredictionBucket) %>% summarise(Total=n(),
                                                                          Correct=sum(IsCorrect),
                                                                          CorrectPerc=round(Correct/Total,3)))

cupPredictionsResult = cupPredictions %>% group_by(DataType,CancerType,PredictionBucket) %>% 
  summarise(Total=n(),
            Correct=sum(IsCorrect),
            Incorrect=sum(IsCorrect==F),
            CorrectPerc=round(Correct/Total,3))

View(cupPredictionsResult)

print(ggplot(cupPredictions %>% group_by(DataType,PredictionBucket) %>% 
               summarise(Total=n(),Correct=sum(IsCorrect),CorrectPerc=round(Correct/Total,3)), aes(x=PredictionBucket,y=CorrectPerc,fill=DataType))
      + geom_bar(stat="identity",colour="black",position='dodge'))

cupPredictionsResult2 = cupPredictions %>% group_by(DataType,CancerType,PredictionBucket,Result=ifelse(IsCorrect,'Correct','Incorrect')) %>% count
View(cupPredictionsResult2)

print(ggplot(cupPredictionsResult2 %>% filter(CancerType=='Liver'), aes(x=PredictionBucket,y=n,fill=Result))
      + geom_bar(stat="identity",colour="black",position='dodge')
      + facet_wrap(~DataType))

print(ggplot(cupPredictions %>% group_by(DataType,PredictionBucket,Result=ifelse(IsCorrect,'Correct','Incorrect')) %>% count, aes(x=PredictionBucket,y=n,fill=Result))
      + geom_bar(stat="identity",colour="black",position='dodge')
      + facet_wrap(~DataType))

print(ggplot(cupPredictionsResult, aes(x=PredictionBucket,y=CorrectPerc,fill=DataType))
      + geom_bar(stat="identity",colour="black",position='dodge')
      + labs(x='Prediction Rate', y='Success Rate')
      + facet_wrap(~CancerType))



#View(cupCohortResults %>% group_by(SampleId,CancerType,Category,DataType,Value,RefCancerType) %>% summarise(CtValue=first(RefCancerTypeValue)) 
#     %>% spread(RefCancerType,CtValue,fill=0))

prostateSamples = cupCohortResults %>% filter(CancerType=='Prostate')
nrow(prostateSamples)
View(prostateSamples %>% filter(DataType=='FUSION'))

cupSampleResults = prostateSamples %>% filter(SampleId=='CPCT02070192T')
nrow(cupSampleResults)

View(cupSampleResults %>%
       group_by(SampleId,CancerType,Category,DataType,Value,RefCancerType) %>% 
       summarise(CtValue=first(RefCancerTypeValue)) %>% 
       spread(RefCancerType,CtValue,fill=0) %>% 
       select(SampleId,CancerType,Category,DataType,Value,everything()))



View(prostateSamples %>%  filter(Category=='CLASSIFIER') %>% 
       arrange(SampleId,DataType,-RefCancerTypeValue) %>%
       group_by(SampleId,CancerType,Category,DataType,Value) %>% 
       summarise(TopRefCancerType=first(RefCancerType),TopRefValue=first(RefCancerTypeValue)))

View(prostateSamples %>% filter(Category=='CLASSIFIER'&)
       group_by(SampleId,CancerType,Category,DataType,Value,RefCancerType) %>% summarise(CtValue=first(RefCancerTypeValue)) 
     %>% spread(RefCancerType,CtValue,fill=0))

## Driver data
allDrivers = read.csv('~/data/cup/drivers_4600.csv')

allDrivers = allDrivers %>% filter(SampleId %in% sampleCancerTypes$SampleId)
View(allDrivers %>% group_by(SampleId) %>% count)

write.csv(allDrivers %>% filter(Driverlikelihood>=0.8) %>% select(SampleId,Gene,Driver),'~/data/cup/sample_driver_data.csv',row.names = F,quote = F)
View(allDrivers)

View(refVirusData)

# write with fusions as well
sampleDriversFusions = rbind(allDrivers %>% filter(Driverlikelihood>=0.8) %>% mutate(Type="DRIVER") %>% select(SampleId,Gene,Type),
                             fusions %>% mutate(Type='FUSION') %>% select(SampleId,Gene=Fusion,Type),
                             refVirusData %>% filter(VirusType!='OTHER') %>% mutate(Type='VIRUS') %>% select(SampleId,Gene=VirusType,Type))

View(sampleDriversFusions)
write.csv(sampleDriversFusions,'~/data/cup/sample_driver_fusion_data.csv',row.names = F,quote = F)


refDriverCombined = rbind(refDriverGenes,refDriverSubtypes)
View(refDriverCombined)

write.csv(driverSigs %>% filter(SamplePerc>=0.0001) %>% select(CancerType,Gene,Driver,SamplePerc),'~/data/cup/driver_prevalence.csv',row.names = F,quote = F)

refDriverPrev = refDriverCombined %>% filter(Driver=='ALL') %>% select(Gene,Driver,CancerType,SamplePerc) %>% spread(CancerType,SamplePerc,fill=0)
View(refDriverPrev)
View(driverSigs2 %>% select(Gene,Driver,Breast) %>% filter(Driver=='ALL'&Breast>0.01))
View(driverSigs2 %>% select(Gene,Driver,Breast) %>% filter(Driver=='ALL'))
View(driverSigs2 %>% select(Gene,Driver,Breast,Prostate) %>% filter(Driver=='ALL'&(Breast>0.01|Prostate>0.01)))

View(driverSigs2 %>% select(Gene,Driver,Skin,Lung,Breast,Prostate) %>% filter(Driver=='ALL'&(Skin>0.01|Lung>0.01|Breast>0.01|Prostate>0.01)))

View(driverSigs2 %>% filter(Driver=='ALL'&Gene %in% c('DAXX','MEN1','DMD','OR4A5','AKT1'))) # AKT1;DAXX;DMD;MEN1;OR4A5

driverSigs3 = driverSigs %>% mutate(SamplePerc=pmax(SamplePerc,0.0001)) %>% select(Gene,Driver,CancerType,SamplePerc) %>% spread(CancerType,SamplePerc,fill=0.0001)

View(driverSigs3 %>% filter(Driver=='ALL'))

View(driverSigs2)

write.csv(driverSigs2,'~/data/cup/driver_prev_sigs.csv',row.names = F,quote = F)

View(allDrivers %>% group_by(Category,Driver,Likelihood=0.2*round(Driverlikelihood/0.2)) %>% count)

View(allDrivers %>% group_by(Category,Driver,Gene) %>% summarise(SampleCount=n()) %>% group_by(Gene) %>% summarise(Types=n()))
View(allDrivers %>% group_by(Category,Driver,Gene) %>% summarise(SampleCount=n()))


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





# fit to subset of elevated sigs and sig 1
View(cosmicSigs)
cosmicSubset = cosmicSigs %>% select(Signature.1,Signature.2,Signature.4,Signature.6,Signature.7,Signature.10,Signature.11,Signature.13,Signature.17)
colnames(cosmicSubset) = c('Sig1','Sig2','Sig4','Sig6','Sig7','Sig10','Sig11','Sig13','Sig17')

View(cosmicSubset)
write.csv(cosmicSubset,'~/data/sigs/cosmic_subset.csv',row.names = F,quote = F)

sampleSubset = head(refSnvCounts %>% filter(SnvCount>2e4&SnvCount<=4e4),10)
View(sampleSubset)

subsetSnvCounts = snvSampleCounts2 %>% filter(SampleId %in% sampleSubset$SampleId)
View(subsetSnvCounts)
subsetSnvMatrixData = subsetSnvCounts %>% spread(SampleId,SnvCount,fill=0)
subsetSnvMatrixData = subsetSnvMatrixData %>% select(-BucketName)
View(subsetSnvMatrixData)
write.csv(subsetSnvMatrixData,'~/data/sigs/subset_10_matrix_data.csv',row.names = F,quote = F)


ssFitResults = read.csv('~/logs/SIG_FIT_RESULTS.csv')
View(ssFitResults)

ssBucketFit = read.csv('~/logs/SIG_FIT_RESIDUALS.csv')
ssBucketFit = ssBucketFit %>% mutate(FitDiffPerc=round(abs(FitCount-Count)/Count,4),
                                     SigOverfitPerc=ifelse(SigName!='NONE',round(SigFitCount/Count,4),0))
View(ssBucketFit)
View(ssBucketFit %>% filter(SigName!='NONE'))

View(ssBucketFit %>% filter(SigName!='NONE') %>% group_by(SigName,Bucket) %>% summarise(Count=n(),
                                                                                        SigBucketTotal=sum(SigFitCount),
                                                                                        Max=max(SigOverfitPerc),
                                                                                        Median=median(SigOverfitPerc)))



sigFitResult2 = fit_to_signatures(specSampleCounts %>% select(WIDE01010879T), as.matrix(cosmicSubset, stringsAsFactors=F))
sigFitContributions2 = as.data.frame(sigFitResult2$contribution)
sigFitContributions2 = cbind(rownames(sigFitContributions2),sigFitContributions2)
colnames(sigFitContributions2) = c('Signature','SigContribution')
sigFitContributions2 = sigFitContributions2 %>% mutate(SigPercent=round(SigContribution/sampleTotal,4),
                                                     SigContribution = round(SigContribution,4))

View(sigFitContributions2)


snvSampleCountsExtra = cbind(snvSampleCounts,specSampleCounts %>% select(WIDE01010879T))
write.csv(snvSampleCountsExtra,'~/data/cup/snv_matrix_data_4610_extras.csv',row.names = F,quote = F)

# other fit methods
leastSqFit = read.csv('~/logs/snv_ut_sample_contribs.csv')
View(leastSqFit[,1:10])

expMaxFit = read.csv('~/logs/snv_ut_emax_sample_contribs.csv')
View(expMaxFit[,1:10])



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

## Background Sigs
skinSamples = refSampleData %>% filter(CancerType=='Skin')

write.csv(refSampleData %>% filter(SampleId %in% skinSamples$SampleId) %>% select(SampleId,CancerType,CancerSubtype),
          '~/data/sigs/skin_sample_sig_data.csv',row.names = F,quote = F)

skinSampleCounts = snvSampleCounts2 %>% filter(SampleId %in% skinSamples$SampleId) %>% spread(SampleId,SnvCount,fill=0)
View(snvSampleCounts2)
View(skinSampleCounts[,1:10])
skinSampleCounts = skinSampleCounts %>% select(-BucketName)

write.csv(skinSampleCounts,'~/data/sigs/skin_sample_matrix_data.csv',row.names = F,quote = F)

bgSigs = read.csv('~/logs/snv_bg_ba_ct_bg_sigs.csv')
View(bgSigs)


## SNV Sig Fitting Methods
cosmicSigSubset = cosmicSigs %>% select(Signature.1,Signature.2,Signature.4,Signature.6,Signature.7,Signature.13,Signature.17)
colnames(cosmicSigSubset) = c('Sig1','Sig2','Sig4','Sig6','Sig7','Sig13','Sig17')
write.csv(cosmicSigSubset,'~/data/sigs/cosmic_sigs_subset.csv',row.names = F,quote = F)




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