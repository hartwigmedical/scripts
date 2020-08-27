
#####
## CUP Reference Data Creation

## SANPLE DATA

## Create sample cancer type data
sampleCancerTypes = read.csv('~/data/cup/sample_cancer_types_20200727.csv')
# sampleCancerTypes = read.csv('~/data/cup/sample_cancer_types_20200820.csv')
nrow(sampleCancerTypes) # 4609
View(sampleCancerTypes)

# TO-DO - case insensitive (eg Bone/soft tissue)
sampleCompare = merge(sampleCancerTypes,sampleCancerTypes2 %>% select(SampleId,PrimaryTumorLocation),by='SampleId',all.x=T)
View(sampleCompare %>% filter(as.character(CancerType)!=as.character(PrimaryTumorLocation)))

sampleCompare = merge(refSampleData,sampleCancerTypes2 %>% select(SampleId,PrimaryTumorLocation),by='SampleId',all.x=T)
View(sampleCompare %>% filter(as.character(OrigCancerType)!=as.character(PrimaryTumorLocation)) %>% group_by(OrigCancerType,PrimaryTumorLocation) %>% count)

sampleCancerTypes = sampleCancerTypes %>% filter(!is.na(HmfPatientId)&HmfPatientId!='NULL')
nrow(sampleCancerTypes) # 4481

# correct for unset cancer types
# View(sampleCancerTypes %>% filter(is.na(CancerType)|CancerType==''|CancerType=='NULL'))
sampleCancerTypes = sampleCancerTypes %>% mutate(CancerType=ifelse(CancerType=='NULL','Unknown',as.character(CancerType)))
View(sampleCancerTypes %>% group_by(CancerType) %>% count)
View(sampleCancerTypes %>% filter(CancerType=='Unknown') %>% group_by(CancerType) %>% count)


# link to set of major cancer types
majorCancerTypes = c('Anus','Bone/Soft tissue','Neuroendocrine','Biliary','Uterus','Head and neck','Breast','Lung','Colon/Rectum','Prostate','Skin','Urinary tract','Kidney','Ovary',
                     'Esophagus','Pancreas','Nervous system','Liver','Stomach','Mesothelioma','Lymphoid','Thyroid')

length(majorCancerTypes) # 22

sampleCancerTypes = sampleCancerTypes %>% mutate(OrigCancerType=CancerType,CancerType=ifelse(CancerType %in% majorCancerTypes,as.character(CancerType),'Other'))
View(sampleCancerTypes)

# Form composite cancer types using some sub-types
# Gastrointestinal stromal tumor (GIST) -> GIST

# sub-type inclusion
cancersWithSubtypes = c('Biliary','Bone/Soft tissue','Head and neck','Neuroendocrine','Skin','Uterus')

includedSubtypes = c('Cholangiocarcinoma','Gastrointestinal stromal tumor (GIST)','Leiomyosarcoma','Salivary Gland',
                     'NET: Pancreatic','NEC: Pancreatic','NET: Small Intestinal','Melanoma','Cervical')

View(sampleCancerTypes %>% group_by(CancerType,CancerSubtype) %>% count)

sampleCancerTypes = sampleCancerTypes %>% 
  mutate(OrigCancerType=CancerType,
         CancerType=ifelse(CancerType %in% majorCancerTypes,as.character(CancerType),'Other'),
         CancerKnownSubtype=ifelse(CancerType %in% cancersWithSubtypes,
                        ifelse(CancerSubtype %in% includedSubtypes,ifelse(grepl('GIST',CancerSubtype),'GIST',
                                                                   ifelse(grepl('Pancreatic',CancerSubtype),'Pancreatic',
                                                                   ifelse(grepl('Small Intestinal',CancerSubtype),'Small Intestinal',as.character(CancerSubtype)))),'Other'),''))

View(sampleCancerTypes %>% group_by(CancerType,CancerKnownSubtype,CancerSubtype) %>% count)
View(sampleCancerTypes %>% group_by(CancerType,CancerKnownSubtype) %>% count)


refSampleData = sampleCancerTypes %>% filter(CancerType!='Unknown'&CancerType!='Unknown '&CancerType!='Other')
nrow(refSampleData) # 4311
View(refSampleData)

incorrectSamples = c('CPCT02070117T','CPCT02080240T','CPCT02340030T','CPCT02230013T','CPCT02060124T','CPCT02360014T','CPCT02030520T','CPCT02290050T','CPCT02050181T','CPCT02030243T',
                     'WIDE01010622T','CPCT02130154T','CPCT02340067T','CPCT02010240T','CPCT02040037T','CPCT02050018T','CPCT02060214T','CPCT02220023T','CPCT02030500T',
                     'CPCT02060121T','CPCT02070422T','CPCT02040329T','CPCT02010866T','CPCT02050349T','CPCT02080048T','CPCT02040165T','CPCT02020497T','CPCT02190039T')

nrow(refSampleData %>% filter(SampleId %in% incorrectSamples)) # 28

## De-deup 
refSampleDataDD = refSampleData %>% arrange(HmfPatientId,SampleId) %>% group_by(HmfPatientId) %>% 
  summarise(Samples=n(),SampleId=first(SampleId),IncorrectCount=sum(SampleId %in% incorrectSamples))

nrow(refSampleDataDD) # 4028
View(refSampleDataDD)
View(refSampleDataDD %>% filter(SampleId %in% incorrectSamples | SecondSampleId %in% incorrectSamples))

refSampleDataDD = refSampleDataDD %>% filter(IncorrectCount==0)
nrow(refSampleDataDD) # 4001, was 4150 when cancer type 'Other' was included
View(refSampleDataDD)

refSampleData = refSampleData %>% filter(SampleId %in% refSampleDataDD$SampleId)
View(refSampleData)
nrow(refSampleData) # 4001
View(refSampleData %>% group_by(CancerType,CancerKnownSubtype) %>% count)

refSampleData = refSampleData %>% mutate(DisplayType=ifelse(CancerKnownSubtype!='',paste(CancerType,CancerKnownSubtype,sep=': '),as.character(CancerType)))
View(refSampleData %>% group_by(DisplayType,CancerType,CancerKnownSubtype) %>% count)

#refSampleDataRaw = refSampleData
refSampleData = refSampleData %>% select(SampleId,MajorCancerType=CancerType,CancerType=DisplayType,CancerKnownSubtype,CancerSubtype,OrigCancerType)

# write.csv(refSampleData,'~/data/cup/ref/cup_ref_sample_data_4150.csv',row.names = F,quote = F)
write.csv(refSampleData,'~/data/cup/ref/cup_ref_sample_data_4001.csv',row.names = F,quote = F)

refSampleData = read.csv('~/data/cup/ref/cup_ref_sample_data_4001.csv')
View(refSampleData)

refCancerSampleCounts = refSampleData %>% group_by(CancerType) %>% summarise(SampleCount=n())
View(refCancerSampleCounts) # 30 distinct types
write.csv(refCancerSampleCounts,'~/data/cup/ref/cup_ref_cancer_sample_counts.csv',row.names = F,quote = F)


## SNV Sample Counts

# strip bucket name from sample counts
snvSampleCounts = read.csv('~/data/cup/snv_sample_counts_4610.csv')
ncol(snvSampleCounts)

snvBuckets = snvSampleCounts %>% select(BucketName)
rowIndex = data.frame(as.numeric(as.character(rownames(snvBuckets))))
colnames(rowIndex) <- c("BucketId")
snvBuckets = cbind(rowIndex,snvBuckets)
snvBuckets = snvBuckets %>% mutate(BucketId=BucketId-1)
View(snvBuckets)

snvSampleCounts = snvSampleCounts %>% select(-pancreatic.or.biliary.tract)
snvSampleCounts2 = snvSampleCounts %>% gather('SampleId','SnvCount',2:ncol(snvSampleCounts))
View(snvSampleCounts2)

snvSampleCounts2 = snvSampleCounts2 %>% filter(SampleId %in% refSampleData$SampleId)
refSnvCounts = snvSampleCounts2 %>% group_by(SampleId) %>% summarise(SnvCount=sum(SnvCount))
nrow(refSnvCounts) # 4001
View(refSnvCounts)
refSnvSampleCounts = snvSampleCounts2 %>% spread(SampleId,SnvCount,fill=0)
refSnvSampleCounts = refSnvSampleCounts %>% select(-BucketName)
ncol(refSnvSampleCounts)
View(refSnvSampleCounts[,1:10])

write.csv(refSnvSampleCounts,'~/data/cup/ref/cup_ref_snv_counts.csv',row.names = F,quote = F)

# refSnvCounts = merge(refSnvCounts,refSampleData %>% select(SampleId,CancerType),by='SampleId',all.x=T)

View(refSnvCounts %>% arrange(CancerType,SnvCount) %>% group_by(CancerType) %>% summarise(Samples=n(),Min=min(SnvCount),Max=max(SnvCount),
                                                         Median=median(SnvCount),Avg=mean(SnvCount),
                                                         Rank5=nth(SnvCount,5),
                                                         Rank10=nth(SnvCount,10),
                                                         Rank25=nth(SnvCount,25),
                                                         Rank50=nth(SnvCount,50),
                                                         Rank100=nth(SnvCount,100)))
View(refSnvCounts)

cosmicSigs = read.csv('~/data/cup/build_ref/snv_cosmic_sigs.csv')
View(cosmicSigs)

# Signature Fit Percentages by Cancer Type
library(MutationalPatterns)

refSigFit = fit_to_signatures(refSnvSampleCounts, as.matrix(cosmicSigs, stringsAsFactors=F))
refFitContributions = as.data.frame(refSigFit$contribution)
View(rownames(refFitContributions))
refFitContributions = cbind(rownames(refFitContributions),refFitContributions)
# colnames(refFitContributions) = c('Signature','SigContribution')
View(refFitContributions[,1:10])
refFitContributions2 = refFitContributions %>% gather('SampleId','SigContrib',2:ncol(refFitContributions))
refFitContributions2 = merge(refFitContributions2,refSnvCounts %>% select(SampleId,SnvCount),by='SampleId',all.x=T)
View(refFitContributions2)
colnames(refFitContributions2) = c('SampleId','SigName','SigContrib','SampleTotal')
View(refFitContributions2 %>% group_by(SampleId) %>% count)

refFitContributions2 = refFitContributions2 %>% mutate(SigPercent=pmin(round(SigContrib/SampleTotal,4),1))
View(refFitContributions2)
nrow(refFitContributions2 %>% group_by(SampleId) %>% count) # still 4001..

View(refFitContributions2 %>% group_by(SigName) %>% summarise(Max=max(SigContrib),
                                                              Median=median(SigContrib),
                                                              Mean=mean(SigContrib),
                                                              Total=sum(SigContrib),
                                                              Above10KAnd20Pct=sum(SigContrib>1e4&SigPercent>=0.2),
                                                              Above25KAnd35Pct=sum(SigContrib>2.5e4&SigPercent>=0.35),
                                                              Above50KAnd50Pct=sum(SigContrib>5e4&SigPercent>=0.5)))

sigsIncluded = c('Sig1','Sig2','Sig4','Sig6','Sig7','Sig10','Sig11','Sig13','Sig17')

# write for Cuppa to turn into percentiles
write.csv(refFitContributions2 %>% filter(SigName %in% sigsIncluded) %>% 
            select(SampleId,SigName,SigContrib,SigPercent),'~/data/cup/build_ref/cup_ref_sig_contribs.csv',row.names = F,quote=F)


refSigContribPercentiles = read.csv('~/data/cup/ref/cup_ref_sig_percentiles.csv')
View(refSigContribPercentiles)


## Sample Purity Data
refPurityData = read.csv('~/data/cup/purity_data.csv')
nrow(refPurityData) # 4641
View(refPurityData)

refPurityData = refPurityData %>% filter(SampleId %in% refSampleData$SampleId)
nrow(refPurityData) # 4001

## add in Chord scores
chordScores = read.csv('~/data/cup/chord_scores.csv')
chordScores = chordScores %>% filter(SampleId %in% refSampleData$SampleId)
chordScores = merge(chordScores,refSampleData %>% select(SampleId,CancerType),by='SampleId',all.x=T)
nrow(chordScores)
write.csv(chordScores,'~/data/cup/ref/cup_ref_chord_scores.csv',row.names = F,quote = F)

View(refPurityData)
refPurityData = merge(refPurityData,chordScores %>% select(SampleId,ChordHrd=HRD),by='SampleId',all.x=T)

# write for Cuppa to turn into percentiles
write.csv(refPurityData,'~/data/cup/build_ref/cup_ref_purity_data.csv',row.names = F,quote = F)

refDataPercentiles = read.csv('~/data/cup/ref/cup_ref_sample_trait_percentiles.csv')
View(refDataPercentiles)

refPurityData = merge(refPurityData,refSampleData %>% select(SampleId,CancerType),by='SampleId',all.x=T)
View(refPurityData)
View(refPurityData %>% group_by(Gender) %>% count)

refPurityCt = refPurityData %>% group_by(CancerType) %>% 
  summarise(IsFemale=sum(Gender=='FEMALE'),
            WGD=sum(WholeGenomeDuplication==1))

refPurityCt = merge(refPurityCt,refCancerSampleCounts,by='CancerType',all.x=T)
refPurityCt = refPurityCt %>% mutate(GenderFemalePerc=round(IsFemale/SampleCount,3),
                                     WGDPerc=round(WGD/SampleCount,3))

View(refPurityCt)

write.csv(refPurityCt,'~/data/cup/ref/cup_ref_sample_trait_rates.csv',row.names=F,quote=F)

## Drivers
allDrivers = read.csv('~/data/cup/drivers_4600.csv')
nrow(allDrivers)
View(allDrivers)

refDriverData = allDrivers %>% filter(SampleId %in% refSampleData$SampleId)
nrow(refDriverData %>% group_by(SampleId) %>% count) # 3967


refDriverData = merge(refDriverData,refSampleData %>% select(SampleId,CancerType),by='SampleId',all.x=T)
View(refDriverData)

#refDriverSubtypes = refDriverData %>% group_by(CancerType,Gene,Driver) %>% summarise(SampleCount=round(sum(Driverlikelihood))) %>% filter(SampleCount>0)
#refDriverSubtypes = merge(refDriverSubtypes,cancerSampleCounts %>% select(CancerType,CancerSampleCount=SampleCount),by='CancerType',all.x=T)
#refDriverSubtypes = refDriverSubtypes %>% mutate(SamplePerc=round(SampleCount/CancerSampleCount,6))
#View(refDriverSubtypes)

refDriverGenes = refDriverData %>% group_by(CancerType,SampleId,Gene) %>% summarise(DriverCount=sum(Driverlikelihood))
refDriverGenes = refDriverGenes %>% mutate(DriverCount=pmin(DriverCount,1))
refDriverGenes = refDriverGenes %>% group_by(CancerType,Gene) %>% summarise(SampleCount=sum(DriverCount)) %>% filter(SampleCount>0)
refDriverGenes = merge(refDriverGenes,refCancerSampleCounts %>% select(CancerType,CancerSampleCount=SampleCount),by='CancerType',all.x=T)
refDriverGenes = refDriverGenes %>% mutate(SamplePerc=round(SampleCount/CancerSampleCount,6))
refDriverGenes = refDriverGenes %>% mutate(Driver='ALL')
View(refDriverGenes)

#write.csv(refDriverGenes %>% filter(Driver=='ALL'&SamplePerc>=0.0001) %>% select(CancerType,Gene,Driver,SamplePerc),
#          '~/data/cup/ref/cup_ref_driver_prev_data.csv',row.names = F,quote = F)


## Fusions
fusions = read.csv('~/data/sv/cohort/LNX_FUSIONS.csv')
fusions = fusions %>% filter(Reportable=='true')
fusions = fusions %>% mutate(Fusion=paste(GeneNameUp,GeneNameDown,sep='_'))
refFusions = fusions %>% filter(SampleId %in% refSampleData$SampleId)
View(refFusions)

refFusionData = merge(refFusions %>% select(SampleId,Fusion),refSampleData %>% select(SampleId,CancerType),by='SampleId',all.x=T)

View(refFusionData)
View(refFusionData %>% group_by(SampleId) %>% count)

refFusionPrev = refFusionData %>% group_by(CancerType,Fusion) %>% summarise(SampleCount=n())
refFusionPrev = merge(refFusionPrev,refCancerSampleCounts %>% select(CancerType,CancerSampleCount=SampleCount),by='CancerType',all.x=T)
refFusionPrev = refFusionPrev %>% mutate(SamplePerc=round(SampleCount/CancerSampleCount,6))
View(refFusionPrev)

## Virus Data
virusList = read.csv('~/data/viral_host_ref.csv')
#View(virusList)

virusList = virusList %>% mutate(MajorType=ifelse(grepl('Alphapapillo',virus_name)|grepl('papillo',virus_name),'HPV',
                                                  ifelse(grepl('HBV',virus_name)|grepl('Hepatitis B',virus_name),'HBV',
                                                         ifelse(grepl('Merkel',virus_name),'MERKEL',
                                                                ifelse(grepl('adenovirus',virus_name),'AAV',
                                                                       ifelse(grepl('herpes',virus_name),'HERPES','OTHER'))))))
virusList = virusList %>% select(VirusName=virus_name,VirusType=MajorType)

refVirusData = read.csv('~/data/cup/cup_virus_data.csv')
View(refVirusData)

refVirusData = refVirusData %>% filter(SampleId %in% refSampleData$SampleId)
refVirusData = merge(refVirusData,virusList,by='VirusName',all.x=T)

# remove duplicates
refVirusData = refVirusData %>% group_by(SampleId,VirusType) %>% count %>% ungroup() %>% select(-n)
refVirusData = merge(refVirusData,refSampleData %>% select(SampleId,CancerType),by='SampleId',all.x=T)
View(refVirusData)

refVirusPrev = refVirusData %>% filter(VirusType!='OTHER') %>% group_by(CancerType,VirusType) %>% summarise(SampleCount=n())
refVirusPrev = merge(refVirusPrev,refCancerSampleCounts %>% select(CancerType,CancerSampleCount=SampleCount),by='CancerType',all.x=T)
refVirusPrev = refVirusPrev %>% mutate(SamplePerc=round(SampleCount/CancerSampleCount,6))
View(refVirusPrev)


# write and treat drivers, fusions and viral insertions the same
refFeatures = rbind(refDriverGenes %>% mutate(Type='DRIVER') %>% select(CancerType,Gene,Type,SamplePerc),
                    refFusionPrev %>% mutate(Type='FUSION',Gene=Fusion) %>% select(CancerType,Gene,Type,SamplePerc),
                    refVirusPrev %>% mutate(Type='VIRUS',Gene=VirusType) %>% select(CancerType,Gene,Type,SamplePerc))
  
View(refFeatures)
write.csv(refFeatures %>% filter(SamplePerc>=0.0001),'~/data/cup/ref/cup_ref_dri_fus_vir_prev_data.csv',row.names = F,quote = F)

refFeatures = read.csv('~/data/cup/ref/cup_ref_dri_fus_vir_prev_data.csv')
View(refFeatures)

# number of features per sample 
refFeaturesPerSample = refDriverData %>% group_by(CancerType,SampleId) %>% summarise(DriverCount=sum(Driverlikelihood))
refFeaturesPerSample = refFeaturesPerSample %>% group_by(CancerType) %>% summarise(AvgFeatures=round(mean(DriverCount),2))
View(refFeaturesPerSample)

panCancerFeaturesPerSample = refDriverData %>% group_by(SampleId) %>% summarise(DriverCount=sum(Driverlikelihood)) %>% summarise(AvgFeatures=round(mean(DriverCount),2))
View(panCancerFeaturesPerSample)

write.csv(rbind(refFeaturesPerSample,panCancerFeaturesPerSample %>% mutate(CancerType='ALL')),'~/data/cup/ref/cup_ref_feature_per_sample.csv',row.names = F,quote = F)


## SVs
svData = read.csv('~/data/sv/cohort/LNX_SVS.csv')
nrow(svData)
View(svData)

# types
# LINE SV
# SIMPLE_DEL_20kb-1Mb
# SIMPLE_DUP_32b-200b
# SIMPLE_DUP_100kb-5mb
# Max SV in 1 event
# Telomeric SGL breakend
refSvData = svData %>% filter(SampleId %in% refSampleData$SampleId) %>% 
  mutate(Length=ifelse(Type=='DEL'|Type=='DUP',PosEnd-PosStart,0)) %>%
  group_by(SampleId) %>%
  summarise(LINE=sum(ResolvedType=='LINE'),
            SIMPLE_DEL_20KB_1MB=sum(Type=='DEL'&Length>=2e4&Length<=1e6),
            SIMPLE_DUP_32B_200B=sum(Type=='DUP'&Length>=32&Length<=200),
            SIMPLE_DUP_100KB_5MB=sum(Type=='DUP'&Length>=1e5&Length<=5e6),
            MAX_COMPLEX_SIZE=max(ifelse(ResolvedType=='COMPLEX',ClusterCount,0)),
            TELOMERIC_SGL=sum(Type=='SGL'&(RepeatClass=='Simple_repeat'&(RepeatType=='(TTAGGG)n'|RepeatType=='(CCCTAA)n'))))

nrow(refSvData)
View(refSvData)

# write for cuppa to generate percentiles
write.csv(refSvData,'~/data/cup/build_ref/cup_ref_sv_data.csv',row.names = F,quote=F)
refSvData = read.csv('~/data/cup/build_ref/cup_ref_sv_data.csv')
refSvData = refSvData %>% filter(SampleId %in% refSampleData$SampleId)

refSvPercentiles = read.csv('~/data/cup/ref/cup_ref_sv_percentiles.csv')
View(refSvPercentiles)
write.csv(refSvPercentiles %>% filter(CancerType!='Other'),'~/data/cup/ref/cup_ref_sv_percentiles.csv',row.names = F,quote = F)


# HPV, HBV, Merkel cell, AAV, Herpes
View(virusList)
View(virusList %>% group_by(MajorType) %>% count)
write.csv(virusList,'~/data/cup/virus_types.csv',row.names = F,quote = F)


