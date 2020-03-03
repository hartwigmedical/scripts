
load_cancer_types<-function(sourceFile,keepSubtypes=T,minSampleCount=10)
{
  sampleCancerTypes = read.csv(sourceFile)
  
  if(!keepSubtypes)
  {
    sampleCancerTypes = sampleCancerTypes %>% select(SampleId,CancerType)
  }
  
  if(minSampleCount>0)
  {
    minorTypes = sampleCancerTypes %>% group_by(CancerType) %>% count %>% filter(n<minSampleCount)
    sampleCancerTypes = sampleCancerTypes %>% mutate(CancerType=ifelse(CancerType %in% minorTypes$CancerType,'Other',as.character(CancerType)))
  }
  
  return (sampleCancerTypes)
}

# load current cancer type assignments for 3784 samples
sampleCancerTypes = load_cancer_types('~/data/hpc_sample_cancer_types.csv',F,10)
nrow(sampleCancerTypes)
View(sampleCancerTypes)
View(sampleCancerTypes %>% group_by(CancerType) %>% count)


#####
## HPC (May 2019)
load('~/data/hmf_cohort_may_2019.RData')
hpcCancerTypes = cohort %>% select(SampleId=sampleId,CancerType=cancerType)
View(hpcCancerTypes)
# View(hpcCancerTypes %>% filter(is.na(CancerType)|CancerType=='NULL'))
#View(highestPurityCohort)
View(hpcCancerTypes %>% group_by(CancerType) %>% count)


load('~/data/hpc_driver_paper.RData')
View(highestPurityCohort)
dpHpcData = highestPurityCohort %>% select(SampleId=sampleId,DpCancerType=cancerType,DpCancerSubtype=cancerSubtype)
nrow(dpHpcData)
View(dpHpcData %>% group_by(DpCancerType) %>% count)

dbCancerTypes = read.csv('~/data/db_sample_cancer_types.csv')
nrow(dbCancerTypes)
dbCancerTypes = dbCancerTypes %>% filter(SampleId %in% hpcCancerTypes$SampleId)
View(dbCancerTypes)


mergedCancerTypes = merge(hpcCancerTypes,dpHpcData,by='SampleId',all.x=T)
mergedCancerTypes = merge(mergedCancerTypes,dbCancerTypes %>% 
                            select(SampleId,DbCancerType=CancerType,DbCancerSubtype=CancerSubtype),by='SampleId',all.x=T)

mergedCancerTypes = mergedCancerTypes %>% mutate(DpCancerType=ifelse(is.na(DpCancerType),'',DpCancerType),
                                                 DpCancerSubtype=ifelse(is.na(DpCancerSubtype),'',DpCancerSubtype),
                                                 DbCancerType=ifelse(DbCancerType=='NULL','',as.character(DbCancerType)),
                                                 DbCancerSubtype=ifelse(DbCancerSubtype=='NULL','',as.character(DbCancerSubtype)))

mergedCancerTypes = mergedCancerTypes %>% mutate(DpCancerSubtype=ifelse(DpCancerSubtype=='Colorectal','Colon/Rectum',DpCancerSubtype),
                                                 DpCancerSubtype=ifelse(DpCancerSubtype=='Small intestine','Small Intestinal',DpCancerSubtype),
                                                 DpCancerSubtype=ifelse(DpCancerSubtype=='Renal Cell','Renal cell',DpCancerSubtype),
                                                 DbCancerType=ifelse(DbCancerType=='Nervous system','CNS',DbCancerType))

mergedCancerTypes = mergedCancerTypes %>% mutate(FinalCancerType=ifelse(DbCancerType!='',DbCancerType,
                                                                        ifelse(DpCancerType!='',DpCancerType,CancerType)),
                                                 FinalCancerSubtype=ifelse(DbCancerSubtype!='',DbCancerSubtype,
                                                                           ifelse(DpCancerSubtype!='',DpCancerSubtype,'')))

# final name corrections
mergedCancerTypes = mergedCancerTypes %>% mutate(FinalCancerType=ifelse(is.na(FinalCancerType),'Unknown',FinalCancerType),
                                                 FinalCancerType=ifelse(FinalCancerType=='CUP','Unknown',FinalCancerType),
                                                 FinalCancerType=ifelse(FinalCancerType=='Bone/soft tissue','Bone/Soft tissue',FinalCancerType))

View(mergedCancerTypes %>% group_by(FinalCancerType) %>% count)
nrow(mergedCancerTypes)

write.csv(mergedCancerTypes %>% select(SampleId,CancerType=FinalCancerType,CancerSubtype=FinalCancerSubtype),
          '~/data/hpc_sample_cancer_types.csv',row.names = F,quote = F)


sampleCancerTypes = load_cancer_types('~/data/hpc_sample_cancer_types.csv',F,10)
View(sampleCancerTypes)
View(sampleCancerTypes %>% group_by(CancerType) %>% count)

set_cancer_sample_counts<-function(sampleCancerTypes)
{
  cstSampleCounts = sampleCancerTypes %>% group_by(CancerType,CancerSubtype) %>% summarise(CancerSubtypeSampleCount=n())
  ctSampleCounts = sampleCancerTypes %>% group_by(CancerType) %>% summarise(CancerSampleCount=n())
  cstSampleCounts = merge(cstSampleCounts,ctSampleCounts,by='CancerType',all.x=T)

  return (cstSampleCounts)
}

cancerSampleCounts = set_cancer_sample_counts(sampleCancerAndSubTypes)
View(cancerSampleCounts)


View(mergedCancerTypes %>% filter((DbCancerType!=''&toupper(CancerType)!=toupper(DbCancerType))))
View(mergedCancerTypes %>% filter((DbCancerType!=''&DbCancerType!='Unknown'&toupper(CancerType)!=toupper(DbCancerType))))

View(mergedCancerTypes)
View(mergedCancerTypes %>% filter(CancerType!=FinalCancerType))
View(mergedCancerTypes %>% filter(DbCancerType!=FinalCancerType))

View(mergedCancerTypes %>% filter((DbCancerType!=''&DbCancerType!='Unknown'&CancerType!=DbCancerType)
                                  |(DpCancerType!=''&DpCancerType!='Other'&CancerType!=DpCancerType)
                                  |(DpCancerSubtype!=''&DbCancerSubtype!=''&DpCancerSubtype!=DbCancerSubtype)))


View(mergedCancerTypes %>% filter((DbCancerType!=''&DbCancerType!='Unknown'&toupper(CancerType)!=toupper(DbCancerType))))

View(mergedCancerTypes %>% filter((DbCancerType!=''&DbCancerType!='Unknown'
                                   &DpCancerType!='Other'&DpCancerType!='CUP'&DpCancerSubtype!=''&toupper(DpCancerSubtype)!=toupper(DbCancerSubtype))))

mergedCancerTypes = mergedCancerTypes %>% mutate(Matched=!is.na(CancerType)&!is.na(DpCancerType)&CancerType==DpCancerType,
                                                 NoDpMatch=!is.na(CancerType)&is.na(DpCancerType),
                                                 DpReplacesNull=CancerType=='Unknown'&!is.na(DpCancerType)&DpCancerType!='Unknown')

View(mergedCancerTypes %>% group_by(Matched,NoDpMatch,DpReplacesNull) %>%)

View(mergedCancerTypes %>% filter(!Matched&!NoDpMatch))