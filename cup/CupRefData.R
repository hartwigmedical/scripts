
#####
## CUP Reference Data Creation

## SANPLE DATA

## Create sample cancer type data
sampleCancerTypes = read.csv('~/data/cup/sample_cancer_types_20200727.csv')
nrow(sampleCancerTypes) # 4609
View(sampleCancerTypes)

sampleCancerTypes = sampleCancerTypes %>% filter(!is.na(HmfPatientId)&HmfPatientId!='NULL')
nrow(sampleCancerTypes) # 4481

# correct for unset cancer types
View(sampleCancerTypes %>% filter(is.na(CancerType)|CancerType==''|CancerType=='NULL'))
sampleCancerTypes = sampleCancerTypes %>% mutate(CancerType=ifelse(CancerType=='NULL','Unknown',as.character(CancerType)))
View(sampleCancerTypes %>% group_by(CancerType) %>% count)
View(sampleCancerTypes %>% filter(CancerType=='Unknown') %>% group_by(CancerType) %>% count)


# link to set of major cancer types
majorCancerTypes = c('Bone/Soft tissue','Neuroendocrine','Biliary','Uterus','Head and neck','Breast','Lung','Colon/Rectum','Prostate','Skin','Urinary tract','Kidney','Ovary',
                     'Esophagus','Pancreas','Nervous system','Liver','Stomach','Mesothelioma','Lymphoid','Thyroid')


sampleCancerTypes = sampleCancerTypes %>% mutate(OrigCancerType=CancerType,CancerType=ifelse(CancerType %in% majorCancerTypes,as.character(CancerType),'Other'))
View(sampleCancerTypes)

refSampleData = sampleCancerTypes %>% filter(CancerType!='Unknown'&CancerType!='Unknown ')
nrow(refSampleData) # 4481
View(refSampleData)

## De-deup 
refSampleDataDD = refSampleData %>% group_by(HmfPatientId) %>% summarise(Samples=n(),SampleId=first(SampleId))
nrow(refSampleDataDD) # 4178
View(refSampleDataDD)

refSampleData = refSampleData %>% filter(SampleId %in% refSampleDataDD$SampleId)
View(refSampleData)
nrow(refSampleData) # 4178

write.csv(refSampleData,'~/data/cup/ref/cup_ref_sample_data_4178.csv',row.names = F,quote = F)
refSampleData = read.csv('~/data/cup/ref/cup_ref_sample_data_4178.csv')
View(refSampleData)

# CPCT02010123T

refCancerSampleCounts = refSampleData %>% group_by(CancerType) %>% summarise(SampleCount=n())
View(refCancerSampleCounts) # 22 distinct types
write.csv(refCancerSampleCounts,'~/data/cup/ref/cup_ref_cancer_sample_counts.csv',row.names = F,quote = F)


## SNV Sample Counts

# strip bucket name from sample counts
snvSampleCounts = read.csv('~/data/cup/snv_sample_counts_4610.csv')
ncol(snvSampleCounts)
snvSampleCounts = snvSampleCounts %>% select(-pancreatic.or.biliary.tract)
snvSampleCounts2 = snvSampleCounts %>% gather('SampleId','SnvCount',2:ncol(snvSampleCounts))
View(snvSampleCounts2)

snvSampleCounts2 = snvSampleCounts2 %>% filter(SampleId %in% refSampleData$SampleId)
refSnvCounts = snvSampleCounts2 %>% group_by(SampleId) %>% summarise(SnvCount=sum(SnvCount))
nrow(refSnvCounts)
View(refSnvCounts)
refSnvSampleCounts = snvSampleCounts2 %>% spread(SampleId,SnvCount,fill=0)
refSnvSampleCounts = refSnvSampleCounts %>% select(-BucketName)
ncol(refSnvSampleCounts)
View(refSnvSampleCounts)

write.csv(refSnvSampleCounts,'~/data/cup/ref/snv_ref_matrix_data_4178.csv',row.names = F,quote = F)

refSnvCounts = merge(refSnvCounts,refSampleData %>% select(SampleId,CancerType),by='SampleId',all.x=T)

View(refSnvCounts %>% arrange(CancerType,SnvCount) %>% group_by(CancerType) %>% summarise(Samples=n(),Min=min(SnvCount),Max=max(SnvCount),
                                                         Median=median(SnvCount),Avg=mean(SnvCount),
                                                         Rank5=nth(SnvCount,5),
                                                         Rank10=nth(SnvCount,10),
                                                         Rank25=nth(SnvCount,25),
                                                         Rank50=nth(SnvCount,50),
                                                         Rank100=nth(SnvCount,100)))
View(refSnvCounts)
View(cosmicSigs)

# Signature Fit Percentages by Cancer Type
refSigFit = fit_to_signatures(refSnvSampleCounts, as.matrix(cosmicSigs, stringsAsFactors=F))
refFitContributions = as.data.frame(refSigFit$contribution)
refFitContributions = cbind(rownames(refFitContributions),refFitContributions)
colnames(refFitContributions) = c('Signature','SigContribution')
View(refFitContributions[,1:10])
refFitContributions2 = refFitContributions %>% gather('SampleId','SigContrib',2:ncol(refFitContributions))
refFitContributions2 = merge(refFitContributions2,refSnvCounts,by='SampleId',all.x=T)
View(refFitContributions2)
colnames(refFitContributions2) = c('SampleId','Signature','SigContrib','SampleTotal')

refFitContributions2 = refFitContributions2 %>% mutate(SigPercent=pmin(round(SigContrib/SampleTotal,4),1))
refFitContributions2 = refFitContributions2 %>% mutate(SigName=stri_replace_all_fixed(Signature,'Signature.','Sig'))

View(refFitContributions2 %>% group_by(SampleId) %>% count)

View(refFitContributions2 %>% group_by(SigName) %>% summarise(Max=max(SigContrib),
                                                              Median=median(SigContrib),
                                                              Mean=mean(SigContrib),
                                                              Total=sum(SigContrib),
                                                              Above10KAnd20Pct=sum(SigContrib>1e4&SigPercent>=0.2),
                                                              Above25KAnd35Pct=sum(SigContrib>2.5e4&SigPercent>=0.35),
                                                              Above50KAnd50Pct=sum(SigContrib>5e4&SigPercent>=0.5)))

sigsIgnored = c('Sig20','Sig21','Sig22','Sig23','Sig24','Sig25','Sig26','Sig27','Sig28','Sig29','Sig30','Sig14','Sig15')
write.csv(refFitContributions2 %>% filter(!(SigName %in% sigsIgnored)) %>% 
            select(SampleId,SigName,SigContrib,SigPercent),'~/data/cup/ref/cup_ref_sig_contribs.csv',row.names = F,quote=F)


refSigContribPercentiles = read.csv('~/data/cup/ref/cup_ref_sig_percentiles.csv')
View(refSigContribPercentiles)


## Sample Purity Data
refPurityData = read.csv('~/data/cup/purity_data.csv')
nrow(refPurityData) # 4641
View(refPurityData)

refPurityData = refPurityData %>% filter(SampleId %in% refSampleData$SampleId)
nrow(refPurityData) # 4178
refPurityData = merge(refPurityData,refSnvCounts,by='SampleId',all.x=T)

## add in Chord scores
chordScores = read.csv('~/data/cup/chord_scores.csv')
chordScores = chordScores %>% filter(SampleId %in% refSampleData$SampleId)
chordScores = merge(chordScores,refSampleData %>% select(SampleId,CancerType),by='SampleId',all.x=T)
nrow(chordScores)
write.csv(chordScores,'~/data/cup/ref/cup_ref_chord_scores.csv',row.names = F,quote = F)

View(refPurityData)
refPurityData = merge(refPurityData,chordScores %>% select(SampleId,ChordHrd=HRD),by='SampleId',all.x=T)

# write for Cuppa to turn into percentiles
write.csv(refPurityData,'~/data/cup/ref/ref_cup_purity_data.csv',row.names = F,quote = F)

refDataPercentiles = read.csv('~/data/cup/ref/cup_ref_sample_trait_percentiles.csv')
View(refDataPercentiles)

refPurityData = merge(refPurityData,refSampleData %>% select(SampleId,CancerType),by='SampleId',all.x=T)
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

refDriverData = allDrivers %>% filter(SampleId %in% refSampleData$SampleId)
nrow(refDriverData %>% group_by(SampleId) %>% count) # 4141


refDriverData = merge(refDriverData,refSampleData %>% select(SampleId,CancerType),by='SampleId',all.x=T)
View(refDriverData)

refDriverSubtypes = refDriverData %>% group_by(CancerType,Gene,Driver) %>% summarise(SampleCount=round(sum(Driverlikelihood))) %>% filter(SampleCount>0)
refDriverSubtypes = merge(refDriverSubtypes,cancerSampleCounts %>% select(CancerType,CancerSampleCount=SampleCount),by='CancerType',all.x=T)
refDriverSubtypes = refDriverSubtypes %>% mutate(SamplePerc=round(SampleCount/CancerSampleCount,6))
View(refDriverSubtypes)

refDriverGenes = refDriverData %>% group_by(CancerType,SampleId,Gene) %>% summarise(DriverCount=sum(Driverlikelihood))
refDriverGenes = refDriverGenes %>% mutate(DriverCount=pmin(DriverCount,1))
refDriverGenes = refDriverGenes %>% group_by(CancerType,Gene) %>% summarise(SampleCount=sum(DriverCount)) %>% filter(SampleCount>0)
refDriverGenes = merge(refDriverGenes,refCancerSampleCounts %>% select(CancerType,CancerSampleCount=SampleCount),by='CancerType',all.x=T)
refDriverGenes = refDriverGenes %>% mutate(SamplePerc=round(SampleCount/CancerSampleCount,6))
refDriverGenes = refDriverGenes %>% mutate(Driver='ALL')
View(refDriverGenes)

# for now only write driver gene, not the type (ie AMP, DEL etc)
write.csv(refDriverGenes %>% filter(Driver=='ALL'&SamplePerc>=0.0001) %>% select(CancerType,Gene,Driver,SamplePerc),
          '~/data/cup/ref/cup_ref_driver_prev_data.csv',row.names = F,quote = F)


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

#virusList = read.csv('~/data/viral_host_ref.csv')
#View(virusList)

virusList = virusList %>% mutate(MajorType=ifelse(grepl('Alphapapillo',virus_name)|grepl('papillo',virus_name),'HPV',
                                                  ifelse(grepl('HBV',virus_name)|grepl('Hepatitis B',virus_name),'HBV',
                                                         ifelse(grepl('Merkel',virus_name),'MERKEL',
                                                                ifelse(grepl('adenovirus',virus_name),'AAV',
                                                                       ifelse(grepl('herpes',virus_name),'HERPES','OTHER'))))))
virusList = virusList %>% select(VirusName=virus_name,VirusType=MajorType)

# write and treat drivers, fusions and viral insertions the same
refFusionDrivers = rbind(refDriverGenes %>% mutate(Type='DRIVER') %>% select(CancerType,Gene,Type,SamplePerc),
                         refFusionPrev %>% mutate(Type='FUSION',Gene=Fusion) %>% select(CancerType,Gene,Type,SamplePerc),
                         refVirusPrev %>% mutate(Type='VIRUS',Gene=VirusType) %>% select(CancerType,Gene,Type,SamplePerc))
  
View(refFusionDrivers)
write.csv(refFusionDrivers %>% filter(SamplePerc>=0.0001),'~/data/cup/ref/cup_ref_driver_fusion_prev_data.csv',row.names = F,quote = F)




## SVs
svData = read.csv('~/data/sv/cohort/LNX_SVS.csv')
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
write.csv(refSvData,'~/data/cup/ref/cup_ref_sv_data.csv',row.names = F,quote=F)

refSvPercentiles = read.csv('~/data/cup/ref/cup_ref_sv_percentiles.csv')
View(refSvPercentiles)

# check short DUPs
View(merge(refSvData %>% filter(SIMPLE_DUP_32B_200B>20),refSampleData,by='SampleId',all.x=T))


# HPV, HBV, Merkel cell, AAV, Herpes
View(virusList)
View(virusList %>% group_by(MajorType) %>% count)
write.csv(virusList,'~/data/cup/virus_types.csv',row.names = F,quote = F)



######
## Specific Sample Data Prep

# Purity data
ssPurityData = read.csv('~/data/cup/samples/cup_purity_raw_test1.csv')
View(ssPurityData)
ssPurityData = merge(ssPurityData,ssSnvCounts,by='SampleId',all.x=T)
write.csv(ssPurityData,'~/data/cup/samples/cup_purity_data_test1.csv',row.names = F,quote = F)

write.csv(ssPurityData %>% mutate(CancerType='Unknown') %>% select(SampleId,CancerType),'~/data/cup/samples/cup_sample_data_test1.csv',row.names = F,quote=F)

# SVs
ssSvData = read.csv('~/data/cup/samples/sv_data_raw_test1.csv')
View(ssSvData)

ssSvData = ssSvData %>% mutate(Length=ifelse(Type=='DUP'|Type=='DEL',as.numeric(as.character(PosEnd))-as.numeric(as.character(PosStart)),0))

ssSvSampleData = ssSvData %>% group_by(SampleId) %>%
  summarise(LINE=sum((LEStart!='NONE'&LEStart!='None')|(LEEnd!='NONE'&LEEnd!='None')),
            SIMPLE_DEL_20KB_1MB=sum(Type=='DEL'&Length>=2e4&Length<=1e6),
            SIMPLE_DUP_32B_200B=sum(Type=='DUP'&Length>=32&Length<=200),
            SIMPLE_DUP_100KB_5MB=sum(Type=='DUP'&Length>=1e5&Length<=5e6),
            MAX_EVENT_SIZE=max(ClusterCount),
            TELOMERIC_SGL=sum(Type=='SGL'&(RepeatClass=='Simple_repeat'&(RepeatType=='(TTAGGG)n'|RepeatType=='(CCCTAA)n'))))

View(ssSvSampleData)

write.csv(ssSvSampleData,'~/data/cup/samples/cup_sv_data_test1.csv',row.names = F,quote=F)


# SVN counts and sigs
ssSnvMatrixData = read.csv('~/data/cup/samples/snv_test_sample_counts.csv')

write.csv(ssSnvMatrixData %>% select(-BucketName),'~/data/cup/samples/cup_snv_matrix_data_test1.csv',row.names = F,quote=F)

View(ssSnvCounts)

ssSnvCounts = ssSnvMatrixData %>% gather('SampleId','Count',2:ncol(ssSnvMatrixData))
ssSnvCounts = ssSnvCounts %>% group_by(SampleId) %>% summarise(SnvCount=sum(Count))

View(ssSnvCounts)

ssSigFit = fit_to_signatures(ssSnvCounts %>% select(-BucketName), as.matrix(cosmicSigs, stringsAsFactors=F))
ssFitContributions = as.data.frame(ssSigFit$contribution)
ssFitContributions = cbind(rownames(ssFitContributions),ssFitContributions)
View(ssFitContributions)
ssFitContributions2 = ssFitContributions %>% gather('SampleId','SigContrib',2:ncol(ssFitContributions))
ssFitContributions2 = merge(ssFitContributions2,ssSnvCounts,by='SampleId',all.x=T)
colnames(ssFitContributions2) = c('SampleId','SigName','SigContrib','SampleTotal')
ssFitContributions2 = ssFitContributions2 %>% mutate(SigName=stri_replace_all_fixed(SigName,'Signature.','Sig'))
View(ssFitContributions2)

write.csv(ssFitContributions2 %>% filter(!(SigName %in% sigsIgnored)),'~/data/cup/samples/cup_sig_contribs_test1.csv',row.names = F,quote=F)


# SampleId,SigName,SigContrib,SigPercent


# Drivers and Fusions
ssDriverFusions = read.csv('~/data/cup/samples/.csv')





