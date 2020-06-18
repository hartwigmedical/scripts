
#####
## Urinary Tract and AHR

utSampleData = read.csv('~/data/urinary_tract_sample_data.csv')

utSampleData = utSampleData %>% mutate(AHRMutationType=ifelse(CanonicalHgvsProteinImpact=='p.Gln383His','SNV',
                                                              ifelse(FusionName=='AHR_AHR','FUSION',ifelse(RelativeCN>=3,'AMP','NONE'))),
                                       HasAHRMutation=(AHRMutationType!='NONE'))

View(utSampleData)
View(utSampleData %>% group_by(AHRMutationType,HasAHRMutation) %>% count)

write.csv(utSampleData %>% select(SampleId), '~/logs/ut_sample_ids.csv',row.names = F,quote = F)

View(utSampleData %>% group_by(SampleId) %>% summarise(MutType=first(AHRMutationType)))
View(utSampleData %>% group_by(SampleId) %>% summarise(MutType=first(AHRMutationType)) %>% group_by(MutType) %>% count)


utGenePanelExp = read.csv('~/data/rna/cohort/isofox_ut.sample_gene_perc_data.csv')
View(utGenePanelExp)

utAhrExp = utGenePanelExp %>% filter(GeneName=='AHR')
nrow(utAhrExp)
View(utAhrExp %>% arrange(TPM))
View(utAhrExp %>% group_by(SampleId) %>% count)
utAhrExp = merge(utAhrExp,utSampleData %>% select(SampleId,AHRMutationType),by='SampleId',all.x=T)
utAhrExp = utAhrExp %>% mutate(AHRMutationType=ifelse(is.na(AHRMutationType),'NONE',as.character(AHRMutationType)))

View(utAhrExp %>% group_by(AHRMutationType) %>% 
       summarise(Samples=n(),MinTPM=min(TPM),MedianTPM=median(TPM),MaxTPM=max(TPM),
                 MinPerc=min(CancerPercentile),MedianPerc=median(CancerPercentile),MaxPerc=max(CancerPercentile)))


write.csv(utAhrExp %>% arrange(TPM) %>% select(-GeneId), '~/data/rna/cohort/ut_ahr_tpm_vs_cohort.csv',row.names = F, quote = F)


## co-occurrence of gene expression with AHR mutation (Fusion and AMP)

utGenePanelExp = merge(utGenePanelExp,utSampleData %>% select(SampleId,AHRMutationType),by='SampleId',all.x=T)
utGenePanelExp = utGenePanelExp %>% mutate(AHRMutationType=ifelse(is.na(AHRMutationType),'NONE',as.character(AHRMutationType)))

# define high relative expression as >= 90th percentile, low as <= 10th percentile
View(utGenePanelExp)

# testing
utGenePanelExp = utGenePanelExp %>% mutate(PercBand=ifelse(CancerPercentile<=10,'Low',ifelse(CancerPercentile>=90,'High','Middle')))

View(genesList)
geneName='HRAS'
geneExpData = utGenePanelExp %>% filter(GeneName==geneName)
nrow(geneExpData)
mutationType='AMP'
samplesWithMutation = geneExpData %>% filter(AHRMutationType==mutationType)
band='High'
samplesInBand = geneExpData %>% filter(PercBand==band)
testName=sprintf("%s-%s-%s",geneName,mutationType,band)
utResults = calc_fisher_et(geneExpData,samplesWithMutation,samplesInBand,"Mutation","Band",log=T,returnDF=T,testLabel=testName)
View(utResults)

utResults = data.frame()

mutationTypes = c('AMP','FUSION')
bands = c('High','Low')
genesList = unique(utGenePanelExp$GeneName)

for(mutationType in mutationTypes)
{
  for(band in bands)
  {
    for(geneName in genesList)
    {
      geneExpData = utGenePanelExp %>% filter(GeneName==geneName)
      samplesWithMutation = geneExpData %>% filter(AHRMutationType==mutationType)
      samplesInBand = geneExpData %>% filter(PercBand==band)
      testName=sprintf("%s-%s-%s",geneName,mutationType,band)
      geneResult = calc_fisher_et(geneExpData,samplesWithMutation,samplesInBand,"Mutation","Band",log=F,returnDF=T,testLabel=testName)
      
      if(nrow(utResults) == 0)
        utResults = geneResult
      else
        utResults = rbind(utResults,geneResult)
    }
  }
}

View(utResults)



e




# driver catalog
utDrivers = read.csv('~/data/ut_ahr_driver_catalog.csv')

View(utDrivers %>% group_by(Gene) %>% count)


View(utDrivers)

View(utGenesList)


geneName = 'TERT'
geneSamples = utDrivers %>% filter(Gene==geneName)
ahrSamples = utSampleData %>% filter(HasAHRMutation)
calc_fisher_et(utSampleData,geneSamples,ahrSamples,geneName,'AHR_Mutation',T,F,testLabel='AHR Mutation vs DriverGene')


utGenesList = unique(utDrivers$Gene)
# utGenesList = c('TERT','TP53','EGFR')

#tmp = calc_fisher_et(utSampleData,utDrivers %>% filter(Gene=='TERT'),ahrSamples,'Gene','AHR_Mutation',F,T,testLabel='AHR_vs_DriverGene')
#View(tmp)

ahrSamples = utSampleData %>% filter(HasAHRMutation)
utCocResults = data.frame(matrix(ncol = 10, nrow = 0))
colnames(utCocResults) = colnames(tmp)

for(geneName in utGenesList)
{
  geneSamples = utDrivers %>% filter(Gene==geneName)
  
  geneResults = calc_fisher_et(utSampleData,geneSamples,ahrSamples,'Gene','AHR_Mutation',T,T,testLabel=geneName)
  
  utCocResults = rbind(utCocResults,geneResults)
}

View(utCocResults)



#####
## Signatures

snvMatrixData = read.csv('~/data/sigs/snv_ut_db_sample_counts.csv')
View(snvMatrixData[,1:10])
ncol(snvMatrixData)
nrow(snvMatrixData)

# remove BucketName if exists
snvBuckets = snvMatrixData$BucketName
colnames(snvBuckets) = c('Bucket')
View(snvBuckets)
snvMatrixData = within(snvMatrixData,rm(BucketName))
snvMatrixData = as.matrix(snvMatrixData)
write.csv(snvMatrixData,'~/data/sigs/snv_ut_db_sample_counts.csv',row.names = F,quote = F)

#tmp = read.csv('~/logs/snv_db_sample_counts.csv')
#View(tmp)
#snvBuckets = tmp %>% mutate(Bucket=BucketName) %>% select(Bucket)
#write.csv(snvBuckets,'~/data/sigs/snv_buckets.csv',row.names = F,quote = F)
snvBuckets = read.csv('~/data/sigs/snv_buckets.csv')

snvSummaryCounts = matrix_to_sample_counts(snvMatrixData,snvBuckets)
View(snvSummaryCounts)

utSampleCancerTypes = utSampleData %>% mutate(CancerType=ifelse(HasAHRMutation,'UT_AHR','UT_NoAHR')) %>% select(SampleId,CancerType) 
View(utSampleCancerTypes)
  
cosmicSigs = read.csv('~/data/sigs/snv_cosmic_sigs.csv')
cosmicSigs = as.matrix(cosmicSigs)
contributions = read.csv('~/logs/snv_ut_sample_contribs.csv')
contributions = as.matrix(contributions)
View(contributions[,1:10])

sigNamesNamed = get_signame_list(30,T)
View(sigNamesNamed)

# use MP package least squares fit instead of sig_analyser (or call calc_lsqnonneg)
fitResult = fit_to_signatures(snvMatrixData, cosmicSigs)
mpContributions = fitResult$contribution
View(fitResult$contribution)

View(mpContributions[,1:10])

# 

evaluate_signature_fit('SNV', 'UT_mut_pat', cosmicSigs, mpContributions, snvMatrixData, snvSummaryCounts, utSampleCancerTypes, snvBuckets, sigNamesNamed, 
                         plotByCancerType=T, viewResults=T, bgSigCount=0, printAllPlots=T)

evaluate_signature_fit('SNV', 'UT_siga_lq', cosmicSigs, contributions, snvMatrixData, snvSummaryCounts, utSampleCancerTypes, snvBuckets, sigNamesNamed, 
                       plotByCancerType=T, viewResults=T, bgSigCount=0, printAllPlots=F)

# using bucket analyser fit
baContribs = read.csv('~/logs/snv_ut_ba_sample_contribs.csv')
baContribs = as.matrix(baContribs)
View(baContribs[,1:10])

evaluate_signature_fit('SNV', 'UT_siga_ba', cosmicSigs, baContribs, snvMatrixData, snvSummaryCounts, utSampleCancerTypes, snvBuckets, sigNamesNamed, 
                       plotByCancerType=T, viewResults=F, bgSigCount=0, printAllPlots=F)



sampleTotals = snvSummaryCounts %>% group_by(SampleId) %>% summarise(SampleCount=sum(Count))
signatures = cosmicSigs
contribution = contributions
sigCount = ncol(cosmicSigs)
sigNamesUnamed = get_signame_list(sigCount, F)
colnames(signatures) = sigNamesUnamed
View(sigNamesUnamed)

sigNamesCombined = as.data.frame(cbind(sigNamesUnamed, sigNamesNamed))
colnames(sigNamesCombined) <- c("Signature", "SigName")

# 1. Bucket Evaluation
sigBucketData = get_bucket_data(signatures, contribution, bucketNames)
sigBucketData = merge(sigBucketData,bucketNamesIndexed,by="Bucket",all.x=T)
sigBucketData = merge(sigBucketData,sigNamesCombined,by="Signature",all.x=T)

sigBucketStats = get_sig_bucket_stats(sigBucketData)

# Signature Discovery, by looking at relative contribution of buckets

# Top Bucket Counts per Sample - for now get all buckets
# sampleBucketTopN = get_top_buckets_by_sample(summaryCounts, origSampleCounts, sampleCancerTypes, 0)
sampleBucketData = get_sample_bucket_data(matrixData, origSampleCounts, bucketNames)

sampleSigData = get_sig_data(signatures, contribution, sigNamesNamed)
sampleSigData = append_residuals(contribution, signatures, snvMatrixData, snvBuckets, sampleSigData)

# get cancer type and SV Count
sampleSigData = merge(sampleSigData,utSampleCancerTypes,by='SampleId',all.x=T)
sampleSigData = merge(sampleSigData,sampleTotals,by='SampleId',all.x=T)

sampleSigData = merge(sampleSigData,sigNamesCombined, by="SigName",all.x=T)
View(sampleSigData)

sampleSigData = sampleSigData %>% mutate(AHRStatus=(CancerType=='UT_AHR'))
write.csv(sampleSigData %>% mutate(SigCount=round(Count,1),
                                   SigPercent=round(SigPercent,3),
                                   Signature=SigName) %>% select(SampleId,AHRStatus,SampleCount,Signature,SigCount,SigPercent),
          '~/logs/urinary_tract_ahr_signature_data.csv',row.names = F,quote=F)





