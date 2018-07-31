# Multiple Biopsy data creation

load(file="~/data/multipleBiopsySomaticsWithScope.RData")
View(head(multipleBiopsySomaticsWithScope,100))
nrow(multipleBiopsySomaticsWithScope %>% filter(type=='INDEL'))
nrow(multipleBiopsySomaticsWithScope %>% filter(type=='SNP'))
nrow(multipleBiopsySomaticsWithScope %>% filter(type=='MNP'))

mnvMBData = multipleBiopsySomaticsWithScope %>% filter(type=='MNP')
# View(mnvMBData)
save(mnvMBData, file="~/data/mnvMultBiopData.RData")

indelMBData = multipleBiopsySomaticsWithScope %>% filter(type=='INDEL')
# View(indelMBData)
save(indelMBData, file="~/data/indelMultBiopData.RData")

snvMBData = multipleBiopsySomaticsWithScope %>% filter(type=='SNP')
# View(snvMBData)
save(snvMBData, file="~/data/snvMultBiopData.RData")

rm(multipleBiopsySomaticsWithScope)



# Multiple Biopsy Evaluation

# Inputs:
# variants annotated with PatientId (to link SampleIds) and Scope (Shared, Sample1, Sample2 etc)

load(file="~/data/r_data/mnvMultBiopData.RData")
nrow(mnvMBData)
View(mnvMBData)

# conform to existing MNV data
mnvMBData = mnvMBData %>% select(sampleId, patientId, ref, alt, clonality, scope)
colnames(mnvMBData) = c("SampleId", "PatientId", "Ref", "Alt", "Clonality", "Scope")

mnvMBData = create_mnv_fields(mnvMBData)
mnvMBData$Bucket = mnvMBData$MutationAdj

View(mnvMBData %>% filter(PatientId=="CPCT02020611"))

# mnvBuckets = mnvMBData %>% group_by(Bucket) %>% count()
# View(mnvBuckets)

# dbSampleCancerTypes = read.csv("~/data/patient_cancertypes.csv")
dr022Samples = read.csv("~/data/DR-022_metadata.tsv", sep='\t')
View(dr022Samples)

dr022SampleData = dr022Samples %>% select(patientId,sampleId,primaryTumorLocation)
colnames(dr022SampleData) = c("PatientId","SampleId","CancerType")
View(dr022SampleData)

# check that only patients with 2 samples are included
mbPatientSamples = mnvMBData %>% group_by(PatientId,SampleId) %>% summarise(Count=n())
mbPatientSampleCounts = mbPatientSamples %>% group_by(PatientId) %>% summarise(SampleCount=n(),VariantCount=sum(Count))
View(mbPatientSampleCounts)
mbPatientTwoInstanceIds = mbPatientSampleCounts %>% filter(SampleCount==2)

dr022PatientCancerTypes = dr022SampleData %>% select(PatientId,CancerType) %>% group_by(PatientId) %>% summarise(CancerType=first(CancerType))

mbPatientTwoInstanceIds = merge(mbPatientTwoInstanceIds, dr022PatientCancerTypes, by.x="PatientId", by.y="PatientId", all.x=T)
View(mbPatientTwoInstanceIds)
mbPatientCancerTypes = mbPatientTwoInstanceIds %>% select(PatientId,CancerType)

# now get cancer type for these samples using patientId (not sampleId) to guarantee it is the same

mnvMBData = mnvMBData %>% filter(PatientId %in% mbPatientTwoInstanceIds$PatientId)
View(mnvMBData)
View(head(mnvMBData,100))
nrow(mnvMBData)

# the evaluate function expects a Count field corresponding to bucket
mnvMBData$Count = 1

mnvBuckets = mnvMBData %>% group_by(Bucket) %>% count()
View(mnvBuckets)

# need to ensure that at every MNV bucket (which match the sigs) is represented in at least 1 sample
mnvDefinedBuckets = read.csv("~/data/r_data/mnv_buckets.csv")
View(mnvDefinedBuckets)


viewData=F
plotByCancerType=F
writeToPDF=T
evaluate_nmf_multiple_biopsies("MNV", "multi_biop", mnvMBData, mnvSignatures, mnvDefinedBuckets, mbPatientCancerTypes,
                               viewData, plotByCancerType, writeToPDF)




# INDELs

load(file="~/data/r_data/indelMultBiopData.RData")
nrow(indelMBData)
View(head(indelMBData,100))

# conform to existing INDEL data
indelMBData = indelMBData %>% select(sampleId, patientId, ref, alt, clonality, scope, microhomology, repeatSequence, repeatCount)
colnames(indelMBData) = c("SampleId", "PatientId", "Ref", "Alt", "Clonality", "Scope", "Microhomology", "RepeatSequence", "RepeatCount")
View(head(indelMBData,100))

# View(indelMBData %>% filter(grepl(',',alt)))

indelMBData$AltLength = apply(indelMBData[,c('Alt')], 1, function(x) calcSplitAltLength(x[1], ','))
indelMBData = create_indel_fields(indelMBData)
View(head(indelMBData,100))

indelMBSampleCounts = indelMBData %>% group_by(PatientId,SampleId,Scope,SubType,LengthGrp,MH,RepeatCountHigh) %>% summarise(Count=n())

# make a combined bucket name: DEL vs INS, LengthGroup, RepeatCount high and for DELs whether has MH or not
indelMBSampleCounts = create_indel_bucket(indelMBSampleCounts)
View(indelMBSampleCounts)


# check that only patients with 2 samples are included
imbPatientSamples = indelMBData %>% group_by(PatientId,SampleId) %>% summarise(Count=n())
imbPatientSampleCounts = imbPatientSamples %>% group_by(PatientId) %>% summarise(SampleCount=n(),VariantCount=sum(Count))
View(imbPatientSampleCounts)
imbPatientTwoInstanceIds = imbPatientSampleCounts %>% filter(SampleCount==2)

dr022PatientCancerTypes = dr022SampleData %>% select(PatientId,CancerType) %>% group_by(PatientId) %>% summarise(CancerType=first(CancerType))

imbPatientTwoInstanceIds = merge(imbPatientTwoInstanceIds, dr022PatientCancerTypes, by.x="PatientId", by.y="PatientId", all.x=T)
View(imbPatientTwoInstanceIds)
imbPatientCancerTypes = imbPatientTwoInstanceIds %>% filter(!is.na(CancerType)) %>% select(PatientId,CancerType)
View(imbPatientCancerTypes)

# now get cancer type for these samples using patientId (not sampleId) to guarantee it is the same

indelMBData = indelMBData %>% filter(PatientId %in% imbPatientCancerTypes$PatientId)
View(indelMBData)
nrow(indelMBData)


# need to ensure that at every MNV bucket (which match the sigs) is represented in at least 1 sample
indelDefinedBuckets = read.csv("~/data/r_data/indel_buckets.csv")
View(indelDefinedBuckets)
View(indelBuckets)
View(indelSignatures)

viewData=T
plotByCancerType=F
writeToPDF=T
evaluate_nmf_multiple_biopsies("INDEL", "multi_biop", indelMBSampleCounts, indelSignatures, indelSigNamesStr, indelDefinedBuckets, imbPatientCancerTypes,
                               viewData, plotByCancerType, writeToPDF)



# SNVs
load(file="~/data/r_data/snvMultBiopData.RData")
nrow(snvMBData)
View(head(snvMBData,100))

# conform to existing INDEL data
snvMBData = snvMBData %>% select(sampleId, patientId, ref, alt, clonality, scope, trinucleotideContext)
colnames(snvMBData) = c("SampleId", "PatientId", "Ref", "Alt", "Clonality", "Scope", "Context")
View(head(snvMBData,100))

snvMBData = snvMBData %>% filter(stri_length(Alt)==1)

nrow(snvMBData %>% filter(stri_length(Alt)>1))

snvMBSampleCounts = snvMBData %>% group_by(SampleId,PatientId,Scope,Type,Context) %>% summarise(Count=n())
nrow(snvMBSampleCounts)

# convert to standard type and context
snvMBSampleCounts$StdType = apply(snvMBSampleCounts[,c("Type")], 1, function(x) standard_mutation(x[1]))
snvMBSampleCounts$StdContext = apply(snvMBSampleCounts[,c("Type","StdType","Context")], 1, function(x) standard_context(x[1],x[2],x[3]))

# emptyBuckets = create_empty_signature()
# colnames(emptyBuckets) = c("Type", "Context")
# snvMBSampleCounts = merge(emptyBuckets, snvMBSampleCounts, all=T)
# snvMBSampleCounts[is.na(snvMBSampleCounts)] <- 0
# View(snvMBSampleCounts)
snvMBSampleCounts$Bucket = paste(snvMBSampleCounts$StdType, snvMBSampleCounts$StdContext, sep='_')


# are all buckets represented?
snvBucketData = snvMBSampleCounts %>% group_by(Bucket) %>% summarise(Count=sum(Count))
View(snvBucketData)
# snvMBSampleCounts = merge(emptySigs, snvMBSampleCounts, all=T)
# snvMBSampleCounts[is.na(snvMBSampleCounts)] <- 0

# check that only patients with 2 samples are included
snvPatientSamples = snvMBSampleCounts %>% group_by(PatientId,SampleId) %>% summarise(Count=n())
snvPatientSampleCounts = snvPatientSamples %>% group_by(PatientId) %>% summarise(SampleCount=n(),VariantCount=sum(Count))
View(snvPatientSampleCounts)
snvPatientTwoInstanceIds = snvPatientSampleCounts %>% filter(SampleCount==2)

snvPatientTwoInstanceIds = merge(snvPatientTwoInstanceIds, dr022PatientCancerTypes, by.x="PatientId", by.y="PatientId", all.x=T)
View(snvPatientTwoInstanceIds)
snvPatientCancerTypes = snvPatientTwoInstanceIds %>% filter(!is.na(CancerType)) %>% select(PatientId,CancerType)
View(snvPatientCancerTypes)

# now get cancer type for these samples using patientId (not sampleId) to guarantee it is the same

snvMBSampleCounts = snvMBSampleCounts %>% filter(PatientId %in% snvPatientCancerTypes$PatientId)
nrow(snvMBSampleCounts)

# snvBucketNames = create_empty_signature()
# snvBucketNames = unite(snvBucketNames, "Bucket", type, context, sep='_')
View(snvBucketNames)

View(snvBuckets)
View(snvSignatures)
View(snvSigNamesStr)

viewData=T
plotByCancerType=F
writeToPDF=T
evaluate_nmf_multiple_biopsies("SNV", "multi_biop", snvMBSampleCounts, snvSignatures, snvSigNamesStr, snvBucketNames, snvPatientCancerTypes,
                               viewData, plotByCancerType, writeToPDF)





