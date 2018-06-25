library(purple);
library(devtools);
library(RMySQL)
library(data.table)
library(IRanges)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(stringi)
library(NMF)
library(MutationalPatterns)
library(gridExtra)
library(ggpubr)


detach("package:svnmf", unload=TRUE);
library(svnmf)

# DATA PREPARATION
svData = read.csv('~/logs/CLUSTER_V23.csv')
nrow(svData)

# filter out multiple biopsy (approximately)
# svData = svData %>% filter(!grepl("DRUP", SampleId)&!grepl("TIII", SampleId)&!grepl("TII", SampleId))
# svData = svData %>% filter(PONCount<2) # should already be taken out

# restrict to samples in high-purity and DR022 set
load("~/data/highestPurityCohortSummary.RData")
dr022Samples = read.csv("~/data/DR-022_metadata.tsv", sep='\t')
hpcSamples = highestPurityCohortSummary %>% filter(sampleId %in% dr022Samples$sampleId)
# hpcSamples = highestPurityCohortSummary
nrow(hpcSamples)

svData = svData %>% filter(SampleId %in% hpcSamples$sampleId)
nrow(svData)
n_distinct(svData$SampleId)


# patientCancerTypes = read.csv('~/data/patient_cancertypes.csv')
# View(patientCancerTypes)
# svData = (merge(svData, patientCancerTypes, by.x="SampleId", by.y="SampleId", all.x=TRUE))
# svData$CancerType = ifelse(is.na(svData$CancerType), 'N/A', paste(svData$CancerType, sep="")) # set 'N/A' for unknowns
#
# cancerTypes = svData %>% group_by(CancerType,SampleId) %>% summarise(Count=n())
# cancerTypes = cancerTypes %>% group_by(CancerType) %>% summarise(SampleCount=n(), SvCount=sum(Count))
# View(cancerTypes)

cancerTypes = hpcSamples %>% group_by(cancerType) %>% count()
cancerTypes = setNames(cancerTypes, c("CancerType", "Count"))
sampleCancerTypes = hpcSamples %>% select(sampleId, cancerType)
sampleCancerTypes = setNames(sampleCancerTypes, c("SampleId", "CancerType"))


# prepare additional fields and counts for NMF

svData$Length = ifelse(svData$Type=='BND'|svData$Type=='INS'|svData$ArmEnd!=svData$ArmStart, -1, svData$PosEnd-svData$PosStart)

svData$ArmExpStart = ifelse(svData$ArmExpStart>0,svData$ArmExpStart,0.1)
svData$ArmExpEnd = ifelse(svData$ArmExpEnd>0,svData$ArmExpEnd,0.1)
svData$StressedPPStart = round(1 - ppois(svData$ArmCountStart - 1, svData$ArmExpStart),4)
svData$StressedPPEnd = round(1 - ppois(svData$ArmCountEnd - 1, svData$ArmExpEnd),4)
svData$IsStressed = ifelse((svData$ArmCountStart >= 10 & svData$StressedPPStart <= 0.001)|(svData$ArmCountEnd >= 10 & svData$StressedPPEnd <= 0.001), 1, 0)
# nrow(svData %>% filter(IsStressed==0))

svData$IsLINE = ifelse(svData$LEStart!='false'|svData$LEEnd!='false',1,0)

# filter out invalid SVs (eg cross-chromosomal-arm SVs)
svData = svData %>% filter(svData$Type=='BND'|Length>=0)
nrow(svData)

svData$IsTI = ifelse(svData$LnkTypeStart=='TI'|svData$LnkTypeEnd=='TI',1,0)
svData$IsDB = ifelse(svData$LnkTypeStart=='DB'|svData$LnkTypeEnd=='DB',1,0)

svData$ClusterNone=ifelse(svData$ClusterCount==1,1,0)
svData$ClusterSmall=ifelse(svData$ClusterCount>1&svData$ClusterCount<=3,1,0)
svData$ClusterLarge=ifelse(svData$ClusterCount>3,1,0)

# group all samples into relevant counts - for now lump all stressed-arm variants together
svSampleCounts = (svData %>% group_by(SampleId)
             %>% summarise(LINE=sum(IsLINE==1),
                           Stressed=sum(IsLINE==0&IsStressed==1),
                           Del_LT10K_CN=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length<=1e4&ClusterNone==1),
                           Del_LT10K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length<=1e4&ClusterSmall==1&IsDB==1),
                           Del_LT10K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length<=1e4&ClusterSmall==1&IsTI==1),
                           Del_LT10K_CL=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length<=1e4&ClusterLarge==1),
                           Del_10Kto100K_CN=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e4&Length<=1e5&ClusterNone==1),
                           Del_10Kto100K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e4&Length<=1e5&ClusterSmall==1&IsDB==1),
                           Del_10Kto100K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e4&Length<=1e5&ClusterSmall==1&IsTI==1),
                           Del_10Kto100K_CL=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e4&Length<=1e5&ClusterLarge==1),
                           Del_100Kto500K_CN=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e5&Length<=5e5&ClusterNone==1),
                           Del_100Kto500K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e5&Length<=5e5&ClusterSmall==1&IsDB==1),
                           Del_100Kto500K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e5&Length<=5e5&ClusterSmall==1&IsTI==1),
                           Del_100Kto500K_CL=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>1e5&Length<=5e5&ClusterLarge==1),
                           Del_500Kto5M_CN=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e5&Length<=5e6&ClusterNone==1),
                           Del_500Kto5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e5&Length<=5e6&ClusterSmall==1&IsDB==1),
                           Del_500Kto5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e5&Length<=5e6&ClusterSmall==1&IsTI==1),
                           Del_500Kto5M_CL=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e5&Length<=5e6&ClusterLarge==1),
                           Del_GT5M_CN=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e6&ClusterNone==1),
                           Del_GT5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e6&ClusterSmall==1&IsDB==1),
                           Del_GT5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e6&ClusterSmall==1&IsTI==1),
                           Del_GT5M_CL=sum(IsLINE==0&IsStressed==0&Type=='DEL'&Length>5e6&ClusterLarge==1),
                           Dup_LT10K_CN=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length<=1e4&ClusterNone==1),
                           Dup_LT10K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length<=1e4&ClusterSmall==1&IsDB==1),
                           Dup_LT10K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length<=1e4&ClusterSmall==1&IsTI==1),
                           Dup_LT10K_CL=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length<=1e4&ClusterLarge==1),
                           Dup_10Kto100K_CN=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e4&Length<=1e5&ClusterNone==1),
                           Dup_10Kto100K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e4&Length<=1e5&ClusterSmall==1&IsDB==1),
                           Dup_10Kto100K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e4&Length<=1e5&ClusterSmall==1&IsTI==1),
                           Dup_10Kto100K_CL=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e4&Length<=1e5&ClusterLarge==1),
                           Dup_100Kto500K_CN=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e5&Length<=5e5&ClusterNone==1),
                           Dup_100Kto500K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e5&Length<=5e5&ClusterSmall==1&IsDB==1),
                           Dup_100Kto500K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e5&Length<=5e5&ClusterSmall==1&IsTI==1),
                           Dup_100Kto500K_CL=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>1e5&Length<=5e5&ClusterLarge==1),
                           Dup_500Kto5M_CN=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e5&Length<=5e6&ClusterNone==1),
                           Dup_500Kto5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e5&Length<=5e6&ClusterSmall==1&IsDB==1),
                           Dup_500Kto5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e5&Length<=5e6&ClusterSmall==1&IsTI==1),
                           Dup_500Kto5M_CL=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e5&Length<=5e6&ClusterLarge==1),
                           Dup_GT5M_CN=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e6&ClusterNone==1),
                           Dup_GT5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e6&ClusterSmall==1&IsDB==1),
                           Dup_GT5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e6&ClusterSmall==1&IsTI==1),
                           Dup_GT5M_CL=sum(IsLINE==0&IsStressed==0&Type=='DUP'&Length>5e6&ClusterLarge==1),
                           Inv_LT10K_CN=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length<=1e4&ClusterNone==1),
                           Inv_LT10K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length<=1e4&ClusterSmall==1&IsDB==1),
                           Inv_LT10K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length<=1e4&ClusterSmall==1&IsTI==1),
                           Inv_LT10K_CL=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length<=1e4&ClusterLarge==1),
                           Inv_10Kto100K_CN=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e4&Length<=1e5&ClusterNone==1),
                           Inv_10Kto100K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e4&Length<=1e5&ClusterSmall==1&IsDB==1),
                           Inv_10Kto100K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e4&Length<=1e5&ClusterSmall==1&IsTI==1),
                           Inv_10Kto100K_CL=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e4&Length<=1e5&ClusterLarge==1),
                           Inv_100Kto500K_CN=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e5&Length<=5e5&ClusterNone==1),
                           Inv_100Kto500K_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e5&Length<=5e5&ClusterSmall==1&IsDB==1),
                           Inv_100Kto500K_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e5&Length<=5e5&ClusterSmall==1&IsTI==1),
                           Inv_100Kto500K_CL=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>1e5&Length<=5e5&ClusterLarge==1),
                           Inv_500Kto5M_CN=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e5&Length<=5e6&ClusterNone==1),
                           Inv_500Kto5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e5&Length<=5e6&ClusterSmall==1&IsDB==1),
                           Inv_500Kto5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e5&Length<=5e6&ClusterSmall==1&IsTI==1),
                           Inv_500Kto5M_CL=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e5&Length<=5e6&ClusterLarge==1),
                           Inv_GT5M_CN=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e6&ClusterNone==1),
                           Inv_GT5M_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e6&ClusterSmall==1&IsDB==1),
                           Inv_GT5M_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e6&ClusterSmall==1&IsTI==1),
                           Inv_GT5M_CL=sum(IsLINE==0&IsStressed==0&Type=='INV'&Length>5e6&ClusterLarge==1),
                           Bnd_CN=sum(IsLINE==0&IsStressed==0&Type=='BND'&ClusterNone==1),
                           Bnd_CS_DB=sum(IsLINE==0&IsStressed==0&Type=='BND'&ClusterSmall==1&IsDB==1),
                           Bnd_CS_TI=sum(IsLINE==0&IsStressed==0&Type=='BND'&ClusterSmall==1&IsTI==1),
                           Bnd_CL=sum(IsLINE==0&IsStressed==0&Type=='BND'&ClusterLarge==1)))

nrow(svSampleCounts)

# inputSigCount[[runNumber]] = 12
# inputSampleCounts[[runNumber]] = svSummary
# inputSigNames[[runNumber]] = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
# inputSigNamesNamed[[runNumber]] = c("01_LongDEL_INV", "02_Stressed", "03_LINE", "04_MidDEL", "05_BND_CN", "06_LongDUP", "07_BND_CL", "08_ShortINV", "09_ShortDEL", "10_DBs", "11_ShortDUP", "12_INV_CL")

svSigCount = 11
svSigNamesNum = get_signame_list(11, F)
svSigNamesStr = get_signame_list(11, T)
# sigNamesNamed = c("01_LongDEL_INV", "02_Stressed", "03_LINE", "04_MidDEL", "05_BND_CN", "06_LongDUP", "07_BND_CL", "08_", "09_ShortDEL", "10_DBs", "11_ShortDUP", "12_INV_CL")
# sigNamesNamed = c("01_LINE", "02_BND_CL", "03_Stressed", "04_ShortINV", "05_LongDEL_INV", "06_DBs", "07_ShortDEL", "08_LongDUP", "09_ShortDUP", "10_BND_CN", "11_MidDEL", "12_INV_CL")
# sigNamesNamed = c("01_LINE", "02_BND_CL", "03_Stressed", "04_ShortINV", "05_LongDEL_INV", "06_DBs", "07_ShortDEL", "08_LongDUP", "09_ShortDUP", "10_BND_CN", "11_MidDEL")
svSigNamesNamed = c("01_BND_CS", "02_Med_DUP", "03_Short_INV", "04_Short_DEL", "05_BND_CL", "06_LINE", "07_Long_DUP", "08_Mid_DEL", "09_Short_DUP", "10_BND_CN", "11_Stressed")

svMatrixData = svnmf::convert_summary_counts_to_nmf(svSampleCounts)
View(svMatrixData[,1:10])
svBucketNames = data.frame(svMatrixData$CountField)
colnames(svBucketNames) <- c("Bucket")
View(svBucketNames)
svMatrixData = within(svMatrixData, rm(CountField))
nrow(svMatrixData)
ncol(svMatrixData)

write.csv(svMatrixData, file="~/logs/r_output/sv_nmf_counts_dr22.csv", row.names=F, quote=F)

View(nmfMatrixData)

#nmfEstimate <- nmf(svMatrixData, rank=6:15, method="brunet", nrun=4, seed=123456)
#plot(nmfEstimate)

# generate the actual NMF results
svNmfResult <- nmf(svMatrixData, rank=svSigCount, method="brunet", nrun=5, seed=123456, .opt='vp6')
# save(svNnmfResult, file="~/logs/r_output/nmfResult_all_ploidyGT05.RData")

load(file="~/data/svNmfResult_sig11_dr22.RData")

evaluate_nmf_run("SV", "sig11_DR22", svSigCount, svNmfResult, svMatrixData, svSampleCountsGrouped,
                 sampleCancerTypes, svBucketNames, svSigNamesNum, svSigNamesNamed, FALSE, FALSE, TRUE)


# save sig date for driver-gene correlations
svSignatures = NMF::basis(svNmfResult)
svContribution = NMF::coef(svNmfResult)
svSampleNames = colnames(svContribution)

svSampleSigData = get_sig_data(svSignatures, svContribution, svSigNamesNamed, svSampleNames)
View(svSampleSigData)
write.csv(svSampleSigData, "~/logs/r_output/svSampleSigData.csv", row.names=F, quote=F)


# multiple biopsy data
load(file="~/data/multipleBiopsyStructuralVariantsWithScope.RData")
View(multipleBiopsyStructuralVariantsWithScope)






# RUN 1: Non-stressed SVs (LINE and Stressed removed)

runNumber = 1

nsData = svData %>% filter(IsLINE==0&IsStressed==0)
nrow(nsData)

nsData$IsDB=ifelse(nsData$NearestDBLen>-1&(nsData$NearestDBLen<nsData$NearestTILen|nsData$NearestTILen<30),1,0)
nsData$IsTI=ifelse(nsData$NearestTILen>=30&nsData$IsDB==0,1,0)
nsData$ClusterNone=ifelse(nsData$ClusterCount==1,1,0)
nsData$ClusterSmall=ifelse(nsData$ClusterCount>1&nsData$ClusterCount<=3,1,0)
nsData$ClusterLarge=ifelse(nsData$ClusterCount>3,1,0)

# optionally filter by cancer type
cancerTypes = nsData %>% group_by(CancerType,SampleId) %>% summarise(Count=n())
cancerTypes = cancerTypes %>% group_by(CancerType) %>% summarise(Count=n(), SvCount=sum(Count))
View(cancerTypes)
# nsData = nsData %>% filter(CancerType=='Melanoma')

# group all samples into relevant counts - for now lump all stressed-arm variants together
nsSummary = (nsData %>% group_by(SampleId)
                 %>% summarise(Del_LT10K_CN=sum(Type=='DEL'&Length<=1e4&ClusterNone==1),
                               Del_LT10K_CS_DB=sum(Type=='DEL'&Length<=1e4&ClusterSmall==1&IsDB==1),
                               Del_LT10K_CS_TI=sum(Type=='DEL'&Length<=1e4&ClusterSmall==1&IsTI==1),
                               Del_LT10K_CL=sum(Type=='DEL'&Length<=1e4&ClusterLarge==1),
                               Del_10Kto100K_CN=sum(Type=='DEL'&Length>1e4&Length<=1e5&ClusterNone==1),
                               Del_10Kto100K_CS_DB=sum(Type=='DEL'&Length>1e4&Length<=1e5&ClusterSmall==1&IsDB==1),
                               Del_10Kto100K_CS_TI=sum(Type=='DEL'&Length>1e4&Length<=1e5&ClusterSmall==1&IsTI==1),
                               Del_10Kto100K_CL=sum(Type=='DEL'&Length>1e4&Length<=1e5&ClusterLarge==1),
                               Del_100Kto500K_CN=sum(Type=='DEL'&Length>1e5&Length<=5e5&ClusterNone==1),
                               Del_100Kto500K_CS_DB=sum(Type=='DEL'&Length>1e5&Length<=5e5&ClusterSmall==1&IsDB==1),
                               Del_100Kto500K_CS_TI=sum(Type=='DEL'&Length>1e5&Length<=5e5&ClusterSmall==1&IsTI==1),
                               Del_100Kto500K_CL=sum(Type=='DEL'&Length>1e5&Length<=5e5&ClusterLarge==1),
                               Del_500Kto5M_CN=sum(Type=='DEL'&Length>5e5&Length<=5e6&ClusterNone==1),
                               Del_500Kto5M_CS_DB=sum(Type=='DEL'&Length>5e5&Length<=5e6&ClusterSmall==1&IsDB==1),
                               Del_500Kto5M_CS_TI=sum(Type=='DEL'&Length>5e5&Length<=5e6&ClusterSmall==1&IsTI==1),
                               Del_500Kto5M_CL=sum(Type=='DEL'&Length>5e5&Length<=5e6&ClusterLarge==1),
                               Del_GT5M_CN=sum(Type=='DEL'&Length>5e6&ClusterNone==1),
                               Del_GT5M_CS_DB=sum(Type=='DEL'&Length>5e6&ClusterSmall==1&IsDB==1),
                               Del_GT5M_CS_TI=sum(Type=='DEL'&Length>5e6&ClusterSmall==1&IsTI==1),
                               Del_GT5M_CL=sum(Type=='DEL'&Length>5e6&ClusterLarge==1),
                               Dup_LT10K_CN=sum(Type=='DUP'&Length<=1e4&ClusterNone==1),
                               Dup_LT10K_CS_DB=sum(Type=='DUP'&Length<=1e4&ClusterSmall==1&IsDB==1),
                               Dup_LT10K_CS_TI=sum(Type=='DUP'&Length<=1e4&ClusterSmall==1&IsTI==1),
                               Dup_LT10K_CL=sum(Type=='DUP'&Length<=1e4&ClusterLarge==1),
                               Dup_10Kto100K_CN=sum(Type=='DUP'&Length>1e4&Length<=1e5&ClusterNone==1),
                               Dup_10Kto100K_CS_DB=sum(Type=='DUP'&Length>1e4&Length<=1e5&ClusterSmall==1&IsDB==1),
                               Dup_10Kto100K_CS_TI=sum(Type=='DUP'&Length>1e4&Length<=1e5&ClusterSmall==1&IsTI==1),
                               Dup_10Kto100K_CL=sum(Type=='DUP'&Length>1e4&Length<=1e5&ClusterLarge==1),
                               Dup_100Kto500K_CN=sum(Type=='DUP'&Length>1e5&Length<=5e5&ClusterNone==1),
                               Dup_100Kto500K_CS_DB=sum(Type=='DUP'&Length>1e5&Length<=5e5&ClusterSmall==1&IsDB==1),
                               Dup_100Kto500K_CS_TI=sum(Type=='DUP'&Length>1e5&Length<=5e5&ClusterSmall==1&IsTI==1),
                               Dup_100Kto500K_CL=sum(Type=='DUP'&Length>1e5&Length<=5e5&ClusterLarge==1),
                               Dup_500Kto5M_CN=sum(Type=='DUP'&Length>5e5&Length<=5e6&ClusterNone==1),
                               Dup_500Kto5M_CS_DB=sum(Type=='DUP'&Length>5e5&Length<=5e6&ClusterSmall==1&IsDB==1),
                               Dup_500Kto5M_CS_TI=sum(Type=='DUP'&Length>5e5&Length<=5e6&ClusterSmall==1&IsTI==1),
                               Dup_500Kto5M_CL=sum(Type=='DUP'&Length>5e5&Length<=5e6&ClusterLarge==1),
                               Dup_GT5M_CN=sum(Type=='DUP'&Length>5e6&ClusterNone==1),
                               Dup_GT5M_CS_DB=sum(Type=='DUP'&Length>5e6&ClusterSmall==1&IsDB==1),
                               Dup_GT5M_CS_TI=sum(Type=='DUP'&Length>5e6&ClusterSmall==1&IsTI==1),
                               Dup_GT5M_CL=sum(Type=='DUP'&Length>5e6&ClusterLarge==1),
                               Inv_LT10K_CN=sum(Type=='INV'&Length<=1e4&ClusterNone==1),
                               Inv_LT10K_CS_DB=sum(Type=='INV'&Length<=1e4&ClusterSmall==1&IsDB==1),
                               Inv_LT10K_CS_TI=sum(Type=='INV'&Length<=1e4&ClusterSmall==1&IsTI==1),
                               Inv_LT10K_CL=sum(Type=='INV'&Length<=1e4&ClusterLarge==1),
                               Inv_10Kto100K_CN=sum(Type=='INV'&Length>1e4&Length<=1e5&ClusterNone==1),
                               Inv_10Kto100K_CS_DB=sum(Type=='INV'&Length>1e4&Length<=1e5&ClusterSmall==1&IsDB==1),
                               Inv_10Kto100K_CS_TI=sum(Type=='INV'&Length>1e4&Length<=1e5&ClusterSmall==1&IsTI==1),
                               Inv_10Kto100K_CL=sum(Type=='INV'&Length>1e4&Length<=1e5&ClusterLarge==1),
                               Inv_100Kto500K_CN=sum(Type=='INV'&Length>1e5&Length<=5e5&ClusterNone==1),
                               Inv_100Kto500K_CS_DB=sum(Type=='INV'&Length>1e5&Length<=5e5&ClusterSmall==1&IsDB==1),
                               Inv_100Kto500K_CS_TI=sum(Type=='INV'&Length>1e5&Length<=5e5&ClusterSmall==1&IsTI==1),
                               Inv_100Kto500K_CL=sum(Type=='INV'&Length>1e5&Length<=5e5&ClusterLarge==1),
                               Inv_500Kto5M_CN=sum(Type=='INV'&Length>5e5&Length<=5e6&ClusterNone==1),
                               Inv_500Kto5M_CS_DB=sum(Type=='INV'&Length>5e5&Length<=5e6&ClusterSmall==1&IsDB==1),
                               Inv_500Kto5M_CS_TI=sum(Type=='INV'&Length>5e5&Length<=5e6&ClusterSmall==1&IsTI==1),
                               Inv_500Kto5M_CL=sum(Type=='INV'&Length>5e5&Length<=5e6&ClusterLarge==1),
                               Inv_GT5M_CN=sum(Type=='INV'&Length>5e6&ClusterNone==1),
                               Inv_GT5M_CS_DB=sum(Type=='INV'&Length>5e6&ClusterSmall==1&IsDB==1),
                               Inv_GT5M_CS_TI=sum(Type=='INV'&Length>5e6&ClusterSmall==1&IsTI==1),
                               Inv_GT5M_CL=sum(Type=='INV'&Length>5e6&ClusterLarge==1),
                               Bnd_CN=sum(Type=='BND'&ClusterNone==1),
                               Bnd_CS_DB=sum(Type=='BND'&ClusterSmall==1&IsDB==1),
                               Bnd_CS_TI=sum(Type=='BND'&ClusterSmall==1&IsTI==1),
                               Bnd_CL=sum(Type=='BND'&ClusterLarge==1)))

inputSigCount[[runNumber]] = 10
inputSampleCounts[[runNumber]] = nsSummary
inputSigNames[[runNumber]] = c("Sig1", "Sig2", "Sig3", "Sig4", "Sig5", "Sig6", "Sig7", "Sig8", "Sig9", "Sig10")
# sig10Names <- c('BndNC', 'MidDUP', 'LongDUP', 'BNDClust', 'RES2:1:1', 'BndDB', 'RS3-ShortDUP', 'ShortDEL', "ShortINV", 'LongDEL-FS')

# View(inputSampleCounts[[runNumber]])
# View(inputSigCount[[runNumber]])


# RUN 2: Include all stressed and LINE SVs

runNumber = 1

nrow(svData)


# extract the results
print(runNumber)
signatures = NMF::basis(nmfResult)
contribution = NMF::coef(nmfResult)
View(signatures)
View(contribution)
sampleNames = colnames(contribution)
View(sigNames)
View(sigNamesNamed)
View(sampleNames)


# begin evaluation routines
# to fix
View(sigCount)
# svnmf::evaluate_run(runNumber, sigCount, nmfMatrixData, sampleCounts, signatures, contribution, bucketNames, sigNames)

# Bucket Evaluation
sigBucketData = svnmf::get_bucket_data(signatures, contribution, bucketNames)
View(sigBucketData)
sigBucketStats = svnmf::get_sig_bucket_stats(sigBucketData)
View(sigBucketStats)

# Signature Discovery, by looking at relative contribution of buckets

# report top contributing buckets to aid with signature naming
sigNamesNamed = c("01_CL", "02_LINE", "03_LongDUP", "04_Stressed", "05_DBs", "06_ShortDUP", "07_BND_CN", "08_LongDEL_INV", "09_ShortDEL", "10_MidDUP", "11_MidDEL")

sigBucketTopN = svnmf::get_top_buckets(sigBucketData, sigNamesUnamed, sigNamesNamed)
View(sigBucketTopN)

# key bucket stats
bucketSummaryData = svnmf::get_bucket_stats(sigBucketData)
View(bucketSummaryData)

# least contributing 10 buckets
leastContribBuckets = svnmf::get_least_contrib_buckets(bucketSummaryData)
View(leastContribBuckets)

# compare to signature graphs


# 2 Signature Evaluation

# Insert meaningful signature names here if required

# optionally name signatues for subsequent output
sigNames = sigNamesNamed

rm(sampleSigData)
sampleSigData = svnmf::get_sig_data(signatures, contribution, sigNames, sampleNames)
View(sampleSigData)
View(sampleNames)


# key stats per signature
sigStats = svnmf::get_sig_stats(sampleSigData)
View(sigStats)

# run again, this time bucketing samples into mutational load and cancer types

# get cancer type and SV Count
sampleSigData = (merge(sampleSigData, patientCancerTypes,by.x="SampleId",by.y="SampleId",all.x=TRUE))
sampleSigData$CancerType = ifelse(is.na(sampleSigData$CancerType), 'N/A', paste(sampleSigData$CancerType, sep=""))

sampleSvCounts = sampleSigData %>% group_by(SampleId) %>% summarise(SampleSvCount=sum(SvCount))
sampleSigData = (merge(sampleSigData, sampleSvCounts, by.x="SampleId",by.y="SampleId",all.x=TRUE))
View(sampleSigData)


detach("package:svnmf", unload=TRUE);
library(svnmf)

# DATA OUTPUT TO PDF
View(runNumber)
fileId = "all"
pdf(file=paste("~/logs/r_output/svnmf_", fileId, ".pdf", sep = ""), height = 14, width = 20)

par(mar=c(1,1,1,1))

# 1. NMF estimate data
# plot(nmfEstimate)

# 2. bucket data
title = textGrob("Bucket Summary Data & Top-N Buckets", gp=gpar(fontface="bold", fontsize=16))
grid.arrange(tableGrob(head(bucketSummaryData, 40), rows=NULL),
             tableGrob(sigBucketTopN, rows=NULL),
             ncol = 2, newpage = TRUE, top=title)


bucketSummaryPlot = svnmf::get_bucket_summary_plot(bucketSummaryData)
# View(bucketSummaryData)
grid.arrange(bucketSummaryPlot, ncol = 1, nrow = 1, newpage = TRUE)


title = textGrob("Signature-Bucket Stats & Least Important Buckets", gp=gpar(fontface="bold", fontsize=16))
grid.arrange(tableGrob(sigBucketStats, rows=NULL),
             tableGrob(leastContribBuckets, rows = NULL),
             nrow = 2, newpage = TRUE, top=title)


# 3. default signature-bucket plot
sigBucketsPlot = svnmf::get_bucket_signatures_plot(bucketNames, signatures, sigNames)
print(sigBucketsPlot)
grid.arrange(sigBucketsPlot, ncol = 1, nrow = 1, newpage = TRUE)

# 4. Signature data
title = textGrob("Signature Stats", gp=gpar(fontface="bold", fontsize=16))
grid.arrange(tableGrob(sigStats, rows = NULL), ncol = 1, nrow = 1, top=title, newpage = TRUE)

# 5. Top 50 samples by signature, but include all other signatures as well
plot_top_n_samples_by_sig(sampleSigData, sigNames)

# 6. Sigs with Samples by cancer type
svnmf::plot_sig_samples(sampleSigData, "") # all samples

# detach("package:svnmf", unload=TRUE);
# library(svnmf)

for(cancerType in cancerTypes$CancerType)
{
  if(!is.na(cancerType))
  {
    svnmf::plot_sig_samples(sampleSigData, cancerType)
  }
}

dev.off()



# sample & signature correlation
View(sampleSigData)

modMutLoad = sampleSigData %>% filter(SampleSvCount <= 500)
View(modMutLoad)

sampleSig0 = modMutLoad %>% select(SampleId, SigName, SigPercent)
View(sampleSig0)

sampleSig1 = sampleSig0 %>% spread(SigName, SigPercent)
View(sampleSig1)

# add cancer type back in
sampleSigPerCols = (merge(sampleSig1, patientCancerTypes, by.x="SampleId", by.y="SampleId", all.x=TRUE))
View(sampleSigPerCols)
write.csv(sampleSigPerCols, "~/logs/r_output/sample_sig_percents.csv")

View(cancerTypes)
View(sigNamesNamed)

for(ct in cancerTypes$CancerType) {

  canSamSigPerCols = sampleSigPerCols %>% filter(CancerType==ct)

  if(nrow(canSamSigPerCols) == 0)
  {
    continue
  }

  maxCorrelation = 0
  maxCombo = ""
  minCorrelation = 0
  minCombo = 0

  for(i in 1:length(sigNames))
  {
    for(j in 1:length(sigNames))
    {
      if(j > i)
      {
        pairCols = data.frame(canSamSigPerCols[,i+1], canSamSigPerCols[,j+1])
        colnames(pairCols) <- c('C1', 'C2')
        correlation = round(cor(pairCols$C1, pairCols$C2),3)

        if(!is.na(correlation))
        {
          # print(paste(i, ", ", j, ", correlation ", correlation, sep=""))

          if(correlation > maxCorrelation) {
            maxCorrelation = correlation
            maxCombo = paste(colnames(canSamSigPerCols[i+1])," and ", colnames(canSamSigPerCols[j+1]), sep="")
            # print(paste(ct, " has new maxCorrelation ", maxCorrelation, " with ", maxCombo, sep=""))
          }

          if(correlation < minCorrelation) {
            minCorrelation = correlation
            minCombo = paste(colnames(canSamSigPerCols[i+1])," and ", colnames(canSamSigPerCols[j+1]), sep="")
            # print(paste(ct, " has new minCorrelation ", minCorrelation, " with ", minCombo, sep=""))
          }
        }
      }
    }
  }

  if(maxCorrelation > 0.5)
  {
    print(paste(ct, " has maxCorrelation ", maxCorrelation, " with ", maxCombo, sep=""))
  }

  if(minCorrelation < -0.5)
  {
    print(paste(ct, " has minCorrelation ", minCorrelation, " with ", minCombo, sep=""))
  }
}




