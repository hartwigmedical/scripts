
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

cupShortDups = cupCohortResults %>% filter(DataType=='SIMPLE_DUP_32B_200B')
View(cupShortDups)

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
sampleDriversFusions = rbind(allDrivers %>% filter(Driverlikelihood>=0.5) %>% mutate(Type="DRIVER") %>% select(SampleId,Gene,Type,Likelihood=Driverlikelihood),
                             fusions %>% mutate(Type='FUSION',Likelihood=1) %>% select(SampleId,Gene=Fusion,Type,Likelihood),
                             refVirusData %>% filter(VirusType!='OTHER') %>% mutate(Type='VIRUS',Likelihood=1) %>% select(SampleId,Gene=VirusType,Type,Likelihood))

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
ssFitResults = merge(ssFitResults,refSampleData %>% select(SampleId,CancerType),by='SampleId',all.x=T)
View(ssFitResults)

View(ssFitResults %>% select(-Percent,-FitMethod) %>% spread(Signature,Allocation,fill=0))
View(ssFitResults %>% filter(SampleId=='CPCT02010440T') %>% select(-Percent) %>% spread(Signature,Allocation,fill=0))

sigMethodCompare = ssFitResults %>% 
  filter(SigName!='EXCESS'&SigName!='UNALLOC') %>% 
  filter(Percent>0.05) %>%
  select(SampleId,FitMethod,SigName,Percent) %>% spread(FitMethod,Percent,fill=0)

View(sigMethodCompare)

print(ggplot(sigMethodCompare,aes(x=LEAST_SQUARES,y=SIG_OPTIMISER))
      + geom_point()
      + facet_wrap(~SigName))
             

ssBucketFit = read.csv('~/logs/sig_fit_bucket_data.csv')
ssBucketFit = read.csv('~/logs/sig_fit_residuals.subset.csv')

nrow(ssBucketFit)
ssBucketFit = ssBucketFit %>% mutate(FitDiffPerc=round(abs(FitCount-Count)/Count,4),
                                     SigOverfitPerc=ifelse(SigName!='NONE',round(SigFitCount/Count,4),0),
                                     Prob=ifelse(FitCount>Count,1-ppois(FitCount-1,Count,T),1))
View(ssBucketFit)
View(ssBucketFit %>% filter(SampleId=='CPCT02010440T'))
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




#####
## SNV position buckets
snvPosFreq = read.csv('~/data/cup/snv_position_freq.csv')
nrow(snvPosFreq)
View(head(snvPosFreq,100))
colnames(snvPosFreq)
snvPosFreq = snvPosFreq %>% filter(SampleId %in% refSampleData$SampleId)

View(snvPosFreq %>% group_by(Count=pmin(Count,20)) %>% count)


snvPosFreq2 = snvPosFreq %>% group_by(Chromosome,Position) %>% 
  summarise(Samples=n(),Total=sum(Count),Median=median(Count),Max=max(Count))
View(snvPosFreq2)

snvPosFreq3 = snvPosFreq %>% spread(SampleId,Count,fill=0) 
snvPosBuckets = snvPosFreq3 %>% select(Chromosome,Position)
nrow(snvPosBuckets)
write.csv(snvPosFreq3 %>% select(-Chromosome,-Position),'~/data/cup/ref/cup_ref_snv_pos_matrix_data.csv',row.names = F,quote = F)
View(snvPosFreq3[,1:10])

View(refSnvCounts)

snvCounts = snvPosFreq %>% group_by(SampleId) %>% summarise(SnvCount=sum(Count))
refSnvCounts = snvCounts %>% filter(SampleId %in% refSampleData$SampleId)
View(snvCounts)

snvPosFreq4 = merge(snvPosFreq,refSnvCounts,by='SampleId',all.x=T)
snvPosFreq4 = snvPosFreq4 %>% mutate(ExpCount=round(SnvCount/3e3,3))
snvPosFreq5 = snvPosFreq4 %>% filter(Count>ExpCount)
snvPosFreq5 = snvPosFreq5 %>% mutate(Prob=ppois(Count,ExpCount,F))

View(snvPosFreq5 %>% filter(SampleId=='CPCT02010003T'|SampleId=='CPCT02050386T'))

View(snvPosFreq5 %>% filter(Prob<0.001) %>% group_by(SampleId) %>% summarise(Total=sum(Count)))
snvPosFreq6 = snvPosFreq5 %>% filter(Prob<0.001&Chromosome!='X'&Chromosome!='Y') %>% select(Chromosome,Position,SampleId,Count) %>% spread(SampleId,Count,fill=0)
View(snvPosFreq6[,1:10])
ncol(snvPosFreq6)
nrow(snvPosFreq6)
write.csv(snvPosFreq6 %>% select(-Chromosome,-Position),'~/data/cup/ref/cup_ref_snv_pos_matrix_data_elev.csv',row.names = F,quote = F)


# per-cancer background counts capped at a max
maxSnvCount=2e4
maxSnvCount=5e4
snvPosFreq = merge(snvPosFreq,refSampleData %>% select(SampleId,CancerType),by='SampleId',all.x=T)
snvPosFreq = merge(snvPosFreq,refSnvCounts %>% select(SampleId,SnvCount),by='SampleId',all.x=T)

snvPosCancerRefData = snvPosFreq %>% 
  mutate(Count=ifelse(SnvCount>maxSnvCount,maxSnvCount/SnvCount*Count,Count)) %>%
  group_by(CancerType,Chromosome,Position) %>% summarise(CtCount=sum(Count))

View(snvPosCancerRefData)

snvPosCancerRefSigs = snvPosCancerRefData %>% select(Chromosome,Position,CancerType,CtCount) %>% spread(CancerType,CtCount,fill=0)
View(snvPosCancerRefSigs)
nrow(snvPosCancerRefSigs)
write.csv(snvPosCancerRefSigs %>% ungroup() %>% select(-Chromosome,-Position),'~/data/cup/ref/cup_ref_snv_pos_ct_data.csv',row.names = F,quote = F)
snvPosCancerRefSigs = read.csv('~/data/cup/ref/cup_ref_snv_pos_ct_data.csv')
snvPosCancerRefSigs = cbind(snvPosBuckets,snvPosCancerRefSigs)
View(snvPosCancerRefSigs)

# 10MB and 2MB positions
snvPosFreq10M = snvPosFreq %>% mutate(Position=1e7*round(Position/1e7)) %>% group_by(CancerType,SnvCount,SampleId,Chromosome,Position) %>% summarise(Count=sum(Count))
snvPosFreq10M = snvPosFreq %>% mutate(Position=2e6*round(Position/2e6)) %>% group_by(CancerType,SnvCount,SampleId,Chromosome,Position) %>% summarise(Count=sum(Count))
View(snvPosFreq10M)
View(head(snvPosFreq10M,100))
nrow(snvPosFreq10M)
nrow(snvPosFreq)

snvPosFreqCounts10M = snvPosFreq10M %>% ungroup() %>% select(SampleId,Chromosome,Position,Count) %>% spread(SampleId,Count,fill=0) 
View(snvPosFreqCounts10M[,1:10])
snvPosBuckets10M = snvPosFreqCounts10M %>% select(Chromosome,Position)
nrow(snvPosBuckets10M)
write.csv(snvPosFreqCounts10M %>% select(-Chromosome,-Position),'~/data/cup/ref/cup_ref_snv_pos_matrix_data_2m.csv',row.names = F,quote = F)


snvPosCancerRefData10M = snvPosFreq10M %>% 
  mutate(Count=ifelse(SnvCount>maxSnvCount,maxSnvCount/SnvCount*Count,Count)) %>%
  group_by(CancerType,Chromosome,Position) %>% summarise(CtCount=sum(Count))

View(snvPosCancerRefData10M)

snvPosCancerRefSigs10M = snvPosCancerRefData10M %>% select(Chromosome,Position,CancerType,CtCount) %>% spread(CancerType,CtCount,fill=0)
View(snvPosCancerRefSigs10M)
nrow(snvPosCancerRefSigs10M)
write.csv(snvPosCancerRefSigs10M %>% ungroup() %>% select(-Chromosome,-Position),'~/data/cup/ref/cup_ref_snv_pos_ct_data_2m.csv',row.names = F,quote = F)


snvPosCancerRefData = snvPosCancerRefSigs %>% gather('CancerType','Count',3:ncol(snvPosCancerRefSigs))
View(snvPosCancerRefData)

snvPosCancerRefData = merge(snvPosCancerRefData,snvPosCancerRefData %>% group_by(CancerType) %>% summarise(CancerTotal=sum(Count)),by='CancerType',all.x=T)
snvPosCancerRefData = snvPosCancerRefData %>% mutate(PosPercent=Count/CancerTotal)

snvPosCancerRefSigPerc = snvPosCancerRefData %>% select(Chromosome,Position,CancerType,PosPercent) %>% spread(CancerType,PosPercent,fill=0)
View(snvPosCancerRefSigPerc)

print(ggplot(snvPosCancerRefSigPerc,aes(x=Breast,y=Colon.Rectum))
      + geom_point())

print(ggplot(snvPosCancerRefSigPerc,aes(x=Lung,y=Neuroendocrine))
      + geom_point())

print(ggplot(snvPosCancerRefSigPerc,aes(x=Esophagus,y=Stomach))
      + geom_point())


cssVsRef = read.csv('~/logs/sig_css.ct_ref.csv')
cssVsRef = merge(cssVsRef,refSampleData %>% select(SampleId,CancerType),by='SampleId',all.x=T)
cssVsRef = cssVsRef %>% mutate(Matched=as.character(CancerType)==as.character(RefName),CssBucket=round(CSS,2))
View(cssVsRef %>% group_by(Matched,CssBucket) %>% count %>% spread(Matched,n,fill=0))

# remove Other in both a sample and as a ref
cssVsRef = cssVsRef %>% filter(CancerType!='Other'&RefName!='Other')

cssVsRefResults = cssVsRef %>% arrange(SampleId,-CSS) %>% group_by(SampleId,CancerType) %>% 
  summarise(TopCancerType=first(RefName),
            TopCSS=first(CSS),NextCSS=nth(CSS,2),
            Total=n())

cssVsRef10M = read.csv('~/logs/sig_css.ct_ref_10m.csv')
cssVsRef10M = read.csv('~/logs/sig_css.ct_ref_2m.csv')
cssVsRef10M = merge(cssVsRef10M,refSampleData %>% select(SampleId,CancerType),by='SampleId',all.x=T)
cssVsRef10M = cssVsRef10M %>% mutate(Matched=as.character(CancerType)==as.character(RefName),CssBucket=round(CSS,2))
cssVsRef10M = cssVsRef10M %>% filter(CancerType!='Other'&RefName!='Other')

cssVsRefResults10M = cssVsRef10M %>% arrange(SampleId,-CSS) %>% group_by(SampleId,CancerType) %>% 
  summarise(TopCancerType=first(RefName),
            TopCSS=first(CSS),NextCSS=nth(CSS,2),
            Total=n())

cssVsRefResults = cssVsRefResults %>% mutate(Matched=as.character(CancerType)==as.character(TopCancerType))
cssVsRefResults = merge(cssVsRefResults,refSnvCounts,by='SampleId',all.x=T)
cssVsRefResults = merge(cssVsRefResults,refSampleData %>% select(SampleId,CancerSubtype),by='SampleId',all.x=T)

cssVsRefResults20K = cssVsRefResults
cssVsRefResults50K = cssVsRefResults

View(rbind(cssVsRefResults20K %>% group_by(CancerType) %>% summarise(Sample=n(),Matched=sum(Matched),MatchedPerc=sum(Matched)/n()) %>% mutate(Pos='20K'),
           cssVsRefResults50K %>% group_by(CancerType) %>% summarise(Sample=n(),Matched=sum(Matched),MatchedPerc=sum(Matched)/n()) %>% mutate(Pos='50K')) %>%
     select(CancerType,Sample,MatchedPerc,Pos) %>% spread(Pos,MatchedPerc))

View(cssVsRefResults)
View(cssVsRefResults %>% group_by(CancerType) %>% summarise(Sample=n(),Matched=sum(Matched),MatchedPerc=sum(Matched)/n()))
View(cssVsRefResults %>% group_by(CancerType,CancerSubtype) %>% summarise(Sample=n(),Matched=sum(Matched),MatchedPerc=sum(Matched)/n()) %>% filter(Sample>=5))
View(cssVsRefResults %>% group_by(CancerType,MutLoad=ifelse(SnvCount<0.5e5,'Low','High')) %>% summarise(Sample=n(),Matched=sum(Matched),MatchedPerc=sum(Matched)/n()))
View(cssVsRefResults %>% group_by(CancerType,MutLoad=ifelse(SnvCount<0.5e5,'Low','High')) %>% summarise(MatchedPerc=sum(Matched)/n()) %>% spread(MutLoad,MatchedPerc,fill=0))


cssVsRefResults10M = cssVsRefResults10M %>% mutate(Matched=as.character(CancerType)==as.character(TopCancerType))
cssVsRefResults10M = merge(cssVsRefResults10M,refSnvCounts,by='SampleId',all.x=T)
cssVsRefResults10M = merge(cssVsRefResults10M,refSampleData %>% select(SampleId,CancerSubtype),by='SampleId',all.x=T)
View(cssVsRefResults10M %>% group_by(CancerType) %>% summarise(Sample=n(),Matched=sum(Matched),MatchedPerc=sum(Matched)/n()))

View(rbind(cssVsRefResults %>% group_by(CancerType) %>% summarise(Sample=n(),Matched=sum(Matched),MatchedPerc=sum(Matched)/n()) %>% mutate(Pos='1M'),
           cssVsRefResults10M %>% group_by(CancerType) %>% summarise(Sample=n(),Matched=sum(Matched),MatchedPerc=sum(Matched)/n()) %>% mutate(Pos='2M')))


# by next score
cssVsRefResults = cssVsRefResults %>% mutate(NextCssDiff=round(TopCSS-NextCSS,3))
View(cssVsRefResults)
View(cssVsRefResults %>% group_by(CancerType,NextCssDiff=round(TopCSS-NextCSS,2)) %>% summarise(MatchedPerc=sum(Matched)/n()) %>% spread(NextCssDiff,MatchedPerc,fill=0))

View(cssVsRefResults %>% group_by(NextCssDiff,MutLoad=ifelse(SnvCount<0.5e5,'Low','High')) %>% summarise(Samples=n(),Matched=sum(Matched),MatchedPerc=sum(Matched)/n()))


View(cssVsRefResults %>% filter(Matched==F) %>% group_by(CancerType,TopCancerType) %>% count)

cssVsRefSpread = merge(cssVsRef %>% select(-CssBucket,-Matched),cssVsRefResults %>% select(SampleId,CancerSubtype,SnvCount,TopCSS,TopCancerType,NextCSS,Matched),by='SampleId',all.x=T)
View(cssVsRefSpread)
write.csv(cssVsRefSpread,'~/logs/snv_pos_css_results.csv',row.names = F,quote = F)
View(cssVsRefSpread %>% spread(RefName,CSS,fill=0))

# snvVsCancerTypeCss = read.csv('~/logs/SIG_CSS_RESULTS.csv')

snvVsCancerTypeCss = merge(snvVsCancerTypeCss,refSampleData %>% select(SampleId1=SampleId,CancerType1=CancerType),by='SampleId1',all.x=T)
snvVsCancerTypeCss = merge(snvVsCancerTypeCss,refSnvCounts %>% select(SampleId1=SampleId,SnvCount1=SnvCount),by='SampleId1',all.x=T)
View(snvVsCancerTypeCss)

snvVsCancerTypeCss = snvVsCancerTypeCss %>% mutate(Matched=as.character(CancerType1)==as.character(SampleId2),CssBucket=round(CSS,2))
View(snvVsCancerTypeCss %>% group_by(Matched,CssBucket) %>% count %>% spread(Matched,n,fill=0))

## pairwise comparison using elevated buckets
cssPairs = read.csv('~/logs/sig_css.pair_elev.csv')

cssPairs = merge(cssPairs,refSampleData %>% select(SampleId1=SampleId,CancerType1=CancerType),by='SampleId1',all.x=T)
cssPairs = merge(cssPairs,refSampleData %>% select(SampleId2=SampleId,CancerType2=CancerType),by='SampleId2',all.x=T)
cssPairs = merge(cssPairs,refSnvCounts %>% select(SampleId1=SampleId,SnvCount1=SnvCount),by='SampleId1',all.x=T)
cssPairs = merge(cssPairs,refSnvCounts %>% select(SampleId2=SampleId,SnvCount2=SnvCount),by='SampleId2',all.x=T)

cssPairs = cssPairs %>% mutate(Matched=as.character(CancerType1)==as.character(CancerType2),CssBucket=round(CSS,2))
View(cssPairs %>% group_by(Matched,CssBucket) %>% count %>% spread(Matched,n,fill=0))

cssSamples = rbind(cssPairs %>% select(SampleId=SampleId1,CancerType=CancerType1,OtherCancerType=CancerType2,CSS,Matched),
                  cssPairs %>% select(SampleId=SampleId2,CancerType=CancerType2,OtherCancerType=CancerType1,CSS,Matched))

View(refCancerSampleCounts)
cssSamples = merge(cssSamples,refCancerSampleCounts %>% select(CancerType,CancerTotal=SampleCount),by='CancerType',all.x=T)
nrow(cssSamples)

cssSamples = cssSamples %>% mutate(CssWeight=2^(-100*(1-CSS)),
                                   WeightedCss=CSS*CssWeight/sqrt(CancerTotal))

# double cssWeight = pow(2, -100 * (1 - css));
# weightedCss = css * cssWeight * sqrt(cancerTypeCount);

cssSampleResults = cssSamples %>% group_by(SampleId,CancerType,OtherCancerType) %>% 
  summarise(OtherSamples=n(),WeightedCssTotal=sum(WeightedCss))

View(cssSampleResults)

cssPairResults = cssSampleResults %>% arrange(SampleId,-WeightedCssTotal) %>% group_by(SampleId,CancerType) %>% 
  summarise(TopCancerType=first(OtherCancerType),
            OtherSamples=sum(OtherSamples))

cssPairResults = cssPairResults %>% mutate(Matched=as.character(CancerType)==as.character(TopCancerType))
cssPairResults = merge(cssPairResults,refSnvCounts,by='SampleId',all.x=T)
# cssPairResults = merge(cssPairResults,refSampleData %>% select(SampleId,CancerSubtype),by='SampleId',all.x=T)

View(cssPairResults)

View(cssPairResults %>% group_by(CancerType) %>% summarise(Sample=n(),Matched=sum(Matched),MatchedPerc=sum(Matched)/n()))
View(cssPairResults %>% group_by(OtherSamples=2**round(log(OtherSamples,2))) %>% summarise(Sample=n(),Matched=sum(Matched),MatchedPerc=sum(Matched)/n()))
View(cssPairResults %>% group_by(CancerType,CancerSubtype) %>% summarise(Sample=n(),Matched=sum(Matched),MatchedPerc=sum(Matched)/n()) %>% filter(Sample>=5))
View(cssPairResults %>% group_by(CancerType,MutLoad=ifelse(SnvCount<0.5e5,'Low','High')) %>% summarise(Sample=n(),Matched=sum(Matched),MatchedPerc=sum(Matched)/n()))



write.csv(snvPosFreq %>% filter(SampleId %in% topCss$SampleId1 | SampleId %in% topCss$SampleId2),'~/logs/snv_positions_data.csv',row.names = F,quote = F)

modMutLoadSamples = head(refSnvCounts %>% filter(SnvCount>1e4&SnvCount<2e4),50)
nrow(modMutLoadSamples)
write.csv(snvPosFreq %>% filter(SampleId %in% modMutLoadSamples$SampleId),'~/logs/snv_positions_10-20K_MB.csv',row.names = F,quote = F)



igRegions = snvPosFreq %>% filter((Chromosome==14&Position>=106032614&Position<=107288051)
                                  |(Chromosome==2&Position>=89890568&Position<=90274235)
                                  |(Chromosome==22&Position>=22380474&Position<=23265085))
  

igRegions = merge(igRegions,refSampleData,by='SampleId',all.x=T)
View(igRegions)

View(igRegions %>% group_by(CancerType,SampleId) %>% summarise(SnvCount=sum(Count)) %>%
       group_by(CancerType) %>% summarise(Samples=n(),Total=sum(SnvCount),Median=median(SnvCount),Max=max(SnvCount)))

# rm(snvCssAll)
snvCssAll=snvCss
nrow(snvCss)



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