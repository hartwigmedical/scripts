


#####
## RNA Samples

rnaSampleIds_2131 = read.csv('~/data/rna/samples/rna_samples_2131.txt')
nrow(rnaSampleIds_2131)

rnaSampleData = read.csv('~/data/rna/samples/rna_sample_data.csv')
nrow(rnaSampleData) # 1613 at Apr 2020
View(rnaSampleData)
View(rnaSampleData %>% group_by(CancerType) %>% count)
colnames(rnaSampleData)
write.csv(rnaSampleData,'~/data/rna/samples/rna_sample_data.csv',row.names = F,quote = F)


rnaSamples20200503 = read.csv('~/data/rna/samples/rna_samples_20200503.csv')
nrow(rnaSamples20200503) # 1845
View(rnaSamples20200503) # 1845
write.csv(rnaSamples20200503,'~/data/rna/samples/rna_samples_20200503_1845.txt',row.names = F,quote = F)

# new samples up to 1845
newSamples = rnaSamples20200503 %>% filter(!(SampleId %in% rnaSampleData$SampleId))
View(newSamples)
write.csv(newSamples,'~/data/rna/samples/gcp_batch_star_232.txt',row.names = F,quote = F)

newSampleCTs = read.csv('~/logs/rna_232_sample_ct.csv')
View(newSampleCTs)
newSampleCTs = newSampleCTs %>% mutate(ReadLength=151,
                                       FastqDir='Unknown')
rnaSampleData1845 = rbind(rnaSampleData,newSampleCTs)
nrow(rnaSampleData1845)
View(rnaSampleData1845)
write.csv(rnaSampleData1845,'~/data/rna/samples/rna_sample_data_1845.csv',row.names = F,quote = F)
write.csv(rnaSampleData1845 %>% select(SampleId),'~/data/rna/samples/rna_sample_ids_1845.txt',row.names = F,quote = F)
write.csv(rnaSampleData1845 %>% select(SampleId,ReadLength),'~/data/rna/samples/gcp_batch_isofox_1845.csv',row.names = F,quote = F)

View(rnaSampleData1845 %>% group_by(CancerType) %>% count)

write.csv(rnaSampleData %>% select(SampleId),'~/data/rna/samples/rna_all_samples.txt',row.names=F,quote=F)


# new samples up to 2131 - 1845 + 286
nrow(sampleDedup2)
newSamples20200516 = sampleDedup2 %>% filter(!(SampleId %in% rnaSampleData1845$SampleId))
View(newSamples20200516)
write.csv(newSamples20200516 %>% select(SampleId),'~/data/rna/samples/rna_sample_new_20200516.txt',row.names = F,quote = F)

# cancer type not set correctly
rnaSampleData2131 = rbind(rnaSampleData1845,newSamples20200516 %>% select(SampleId) %>% mutate(ReadLength=151,
                                                                                               FastqDir='Unknown',
                                                                                               CancerType='Unknown'))
View(rnaSampleData2131)
View(rnaSampleData2131 %>% group_by(ReadLength) %>% count)
write.csv(rnaSampleData2131 %>% select(SampleId,ReadLength),'~/data/rna/samples/gcp_batch_isofox_2131.csv',row.names = F,quote = F)

View(rnaSampleData2131 %>% filter(!(SampleId %in% rnaSampleData$SampleId)))
write.csv(rnaSampleData2131 %>% filter(!(SampleId %in% rnaSampleData$SampleId)) %>% select(SampleId,ReadLength),
          '~/data/rna/samples/gcp_batch_isofox_extra_500.csv',row.names = F,quote = F)



## Cohort analysis
fusionCohort = read.csv('~/logs/isofox_fusion_cohort.csv')
View(fusionCohort)

sampleDedup = read.csv('~/data/rna/samples/rna_sample_dedup.csv')
View(sampleDedup)
nrow(sampleDedup)
sampleDedupSummary = sampleDedup %>% group_by(PatientId) %>% summarise(Samples=n(),SampleId=first(SampleId))
View(sampleDedupSummary)

View(sampleDedup %>% group_by(PatientId) %>% summarise(Samples=n(),FirstSample=first(SampleId),SecondSample=last(SampleId)) %>% filter(Samples>=2))

sampleDedupSummary = sampleDedup %>% group_by(PatientId) %>% summarise(Samples=n(),SampleId=first(SampleId))

sampleDedup2 = read.csv('~/data/rna/samples/rna_sample_dedup2.csv')
View(sampleDedup2)
View(sampleDedup2 %>% group_by(HmfPatientId) %>% summarise(Samples=n()))
nrow(sampleDedup2)
sampleDedupSummary2 = sampleDedup2 %>% group_by(HmfPatientId) %>% summarise(Samples=n(),SampleId=first(SampleId))
nrow(sampleDedupSummary2) #2030

View(sampleDedup2 %>% group_by(HmfPatientId) %>% summarise(Samples=n(),FirstSample=first(SampleId),SecondSample=last(SampleId)) %>% filter(Samples>=2))

write.csv(rnaSampleData1845 %>% filter(SampleId %in% sampleDedupSummary2$SampleId) %>% select(SampleId,CancerType),
          '~/data/rna/samples/dl_rna_sample_cts_dedup.csv',row.names = F,quote = F)

# new samples as of 2020-05-16
write.csv(sampleDedup %>% filter(!(SampleId %in% rnaSampleData1845$SampleId)) %>% select(SampleId),
          '~/data/rna/samples/gcp_batch_star_20200516.csv',row.names = F,quote = F)


write.csv(sampleDedup2 %>% select(SampleId),'~/data/rna/samples/rna_samples_2131.txt',row.names = F,quote = F)

write.csv(rnaSampleData2131 %>% filter(SampleId %in% sampleDedupSummary2$SampleId) %>% select(SampleId,CancerType),
          '~/data/rna/samples/dl_rna_sample_dedup_2030.csv',row.names = F,quote = F)

write.csv(rnaSampleData2131 %>% select(SampleId,CancerType),
          '~/data/rna/samples/dl_rna_samples_2131.csv',row.names = F,quote = F)


# fastqFile='/Users/charlesshale/data/rna/rna_sample_fastqs_20200318.txt'

# write all samples for data-loading, eg for summary data
write.csv(rnaSampleData %>% select(SampleId,MajorCancerType),'~/data/rna/samples/dl_all_samples.csv', row.names = F,quote = F)


## Cancer Cohorts for Expression Analysis
View(rnaSampleData %>% group_by(CancerType) %>% count %>% arrange(-n))
cancerSampleCutoff=15
majorCancerTypes = rnaSampleData %>% group_by(CancerType) %>% count %>% filter(n>=cancerSampleCutoff)
View(majorCancerTypes)

rnaSampleData = rnaSampleData %>% mutate(MajorCancerType=ifelse(CancerType %in% majorCancerTypes$CancerType,as.character(CancerType),'Other'))
View(rnaSampleData %>% group_by(MajorCancerType) %>% count)

dataDir='~/data/rna/samples'
for(cancerType in unique(rnaSampleData$MajorCancerType))
{

  fileName = stri_replace_all_fixed(cancerType,'/','')
  fileName = stri_replace_all_fixed(fileName,' ','')
  fileName = stri_trans_tolower(fileName)

  print(sprintf("%s: %s", cancerType, fileName))
  
  write.csv(rnaSampleData %>% filter(MajorCancerType==cancerType) %>% select(SampleId,MajorCancerType),
          paste(sprintf("%s/dl_%s_samples.csv",dataDir,fileName)), row.names = F,quote = F)
}



## Cohort Summary Stats
starStats = read.csv('~/data/rna/cohort/star_stats.csv')
View(starStats)
colnames(starStats)
View(starStats %>% group_by(SampleId) %>% count)
View(rnaSampleData %>% group_by(SampleId) %>% count)

# both complete
View(rnaSampleData %>% filter(!(SampleId %in% starStats$SampleId)))
View(starStats %>% filter(!(SampleId %in% rnaSampleData$SampleId)))

View(starStats %>% filter(is.na(TooShortReads)|is.na(TooManyLociReads)))

cohortSummary = read.csv('~/data/rna/samples/isofox_summary.csv')
View(cohortSummary)

# validation
View(cohortSummary %>% filter(DuplicateFragments>TotalFragments))

cohortCombinedStats = merge(cohortSummary,starStats,by='SampleId',all.x=T)
View(cohortCombinedStats)
write.csv(cohortCombinedStats,'~/data/rna/cohort/star_isofox_summary.csv',row.names = F,quote = F)

View(ensemblGeneData)
View(ensemblTransExonData %>% filter(GeneId=='ENSG00000170242'))

# has been deduped
View(rnaSampleData %>% group_by(SampleId) %>% count)


# re-run STAR for the samples which may have included their TII fastqs??
View(rnaSampleData %>% filter(grepl('TII',SampleId)|grepl('TIII',SampleId)))
write.csv(rnaSampleData %>% filter(grepl('TII',SampleId)|grepl('TIII',SampleId)) %>% select(SampleId),'~/data/rna/samples/gcp_star_rerun_t2.txt',row.names = F, quote = F)
primeSamples = read.csv('~/data/rna/samples/gcp_star_rerun_t2.txt')
View(primeSamples)
View(primeSamples %>% filter(SampleId %in% rnaSampleData$SampleId))
write.csv(primeSamples %>% filter(SampleId %in% rnaSampleData$SampleId),'~/data/rna/samples/gcp_batch_star_rerun_t2.txt',row.names = F, quote = F)
write.csv(rnaSampleData %>% filter(SampleId %in% primeSamples$SampleId) %>% select(SampleId,ReadLength),'~/data/rna/samples/gcp_batch_isofox_rerun_t2.txt',row.names = F, quote = F)

# UT 138 - first set
write.csv(rnaSampleData %>% filter(CancerType=='Urinary tract'&ReadLength==76) %>% select(SampleId,CancerType),
          '~/data/rna/dl_ut_138_samples.csv', row.names = F,quote = F)

View(utSamples)



#############

library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)

qcStats = cohortCombinedStats
qcStats = read.csv('~/hmf/analyses/RNA/star_isofox_summary.csv') 

qcStats = qcStats %>% mutate(DuplicatePerc=DuplicateFragments/TotalFragments,UniqueFragments=TotalFragments-DuplicateFragments) %>% 
  mutate(ReadLength=as.factor(ReadLength), TooManyLociReads=TooManyLociReads/100,MultipleLoci=MultipleLoci/100,TooShortReads=TooShortReads/100,
         MapFailPerc=TooManyLociReads+TooShortReads,InitialFragments=TotalFragments/(1-MapFailPerc))

# 151 base samples have higher fragment counts, but a higher percentage of duplicates, and higher number of mapping fails
plot_grid(ggplot(qcStats, aes(TotalFragments,UniqueFragments,colour = ReadLength)) + geom_point(),
          ggplot(qcStats, aes(InitialFragments,UniqueFragments,colour = ReadLength)) + geom_point() + scale_x_log10()+ scale_y_log10(),
          ggplot(qcStats, aes(DuplicatePerc,colour = ReadLength)) + stat_ecdf(pad=FALSE,geom='step'),
          ggplot(qcStats, aes(MapFailPerc,colour = ReadLength)) + stat_ecdf(pad=FALSE,geom='step'))

# 10 samples have very low numbers of unique fragments both due to map fail and low number of initial fragments
View(qcStats %>% filter(UniqueFragments<2e6))

# 151 base reads have longer median fragment lengths, but still an average of ~210 bases meaning most bases are duplicated in paired fragments.   The 95th% needs investigation.  A significantly higher percentage of 151 base reads are
plot_grid(ggplot(qcStats, aes(FragLength5th,colour = ReadLength)) + stat_ecdf(pad=FALSE,geom='step'),
          ggplot(qcStats, aes(FragLength50th,colour = ReadLength)) + stat_ecdf(pad=FALSE,geom='step'),
          ggplot(qcStats, aes(FragLength95th,colour = ReadLength)) + stat_ecdf(pad=FALSE,geom='step'),
          ggplot(qcStats, aes(TooShortReads,colour = ReadLength)) + stat_ecdf(pad=FALSE,geom='step'))

# 151 base samples have more spliced and chimeric fragments, less unspliced and unspliced fragments
plot_grid(ggplot(qcStats, aes(ChimericFragmentPerc,colour = ReadLength)) + stat_ecdf(pad=FALSE,geom='step'),
          ggplot(qcStats, aes(SplicedFragmentPerc,colour = ReadLength)) + stat_ecdf(pad=FALSE,geom='step'),
          ggplot(qcStats, aes(UnsplicedFragmentPerc,colour = ReadLength)) + stat_ecdf(pad=FALSE,geom='step'),
          ggplot(qcStats, aes(AltFragmentPerc,colour = ReadLength)) + stat_ecdf(pad=FALSE,geom='step'))

#  5 of the 76 base samples have extremely high rates of chimeric fragments and 4 have very high rates of Alt fragments
View(qcStats %>% filter(ChimericFragmentPerc>0.2|AltFragmentPerc>0.1))

plot_grid(ggplot(qcStats, aes(EnrichedGenePercent,colour = ReadLength)) + stat_ecdf(pad=FALSE,geom='step'),
          ggplot(qcStats, aes(MedianGCRatio,colour = ReadLength)) + stat_ecdf(pad=FALSE,geom='step'))

#  11 of the 76 base samples have extremely high rates of enrichment of the 6 highly expresse genes
View(qcStats %>% filter(EnrichedGenePercent>0.5))

# More 76 base reads map to too many loci (set to >3 in STAR, likely due  to shorter length) but more 151 base reads map to multiple loci.
plot_grid(ggplot(qcStats, aes(TooManyLociReads,colour = ReadLength)) + stat_ecdf(pad=FALSE,geom='step'),
          ggplot(qcStats, aes(MultipleLoci,colour = ReadLength)) + stat_ecdf(pad=FALSE,geom='step'),
          ggplot(qcStats, aes(MultipleLoci,TooManyLociReads,colour = ReadLength)) + geom_point())

#High enriched gene percent associated with high GC, long 95th%, short 50th%, short 5th% 
plot_grid(ggplot(qcStats, aes(EnrichedGenePercent,MedianGCRatio,colour = ReadLength)) + geom_point(),
          ggplot(qcStats, aes(EnrichedGenePercent,FragLength5th,colour = ReadLength)) + geom_point(),
          ggplot(qcStats, aes(EnrichedGenePercent,FragLength50th,colour = ReadLength)) + geom_point(),
          ggplot(qcStats, aes(EnrichedGenePercent,FragLength95th,colour = ReadLength)) + geom_point())


# Long 95th% is associated with short 5th and 50th% size band often high Total Fragments
plot_grid(ggplot(qcStats, aes(FragLength95th,TotalFragments,colour = ReadLength)) + geom_point(),
          ggplot(qcStats, aes(FragLength95th,FragLength5th,colour = ReadLength)) + geom_point(),
          ggplot(qcStats, aes(FragLength95th,FragLength50th,colour = ReadLength)) + geom_point(),
          ggplot(qcStats, aes(FragLength95th,TooShortReads,colour = ReadLength)) + geom_point())

# Samples with few total fragments highly associated with high rates of multiple loci
plot_grid(ggplot(qcStats, aes(TotalFragments,MultipleLoci,colour = ReadLength)) + geom_point())


View(cohortCombinedStats %>% filter(FragLength95th==3000) %>% arrange(TotalFragments))


