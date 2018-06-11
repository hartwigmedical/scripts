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
svData = read.csv('~/logs/CLUSTER_V16.csv')

# filter out multiple biopsy (approximately)
svData = svData %>% filter(!grepl("DRUP", SampleId)&!grepl("TIII", SampleId)&!grepl("TII", SampleId))
svData = svData %>% filter(PONCount<2) # should already be taken out


# take cancer type from a CSV file
patientCancerTypes = read.csv('~/data/patient_cancertypes.csv')
View(patientCancerTypes)
svData = (merge(svData, patientCancerTypes, by.x="SampleId", by.y="SampleId", all.x=TRUE))
svData$CancerType = ifelse(is.na(svData$CancerType), 'N/A', paste(svData$CancerType, sep="")) # set 'N/A' for unknowns

cancerTypes = svData %>% group_by(CancerType,SampleId) %>% summarise(Count=n())
cancerTypes = cancerTypes %>% group_by(CancerType) %>% summarise(SampleCount=n(), SvCount=sum(Count))
View(cancerTypes)


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



# PREPARE INPUT COUNTS

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

svData$IsDB=ifelse(svData$NearestDBLen>-1&(svData$NearestDBLen<svData$NearestTILen|svData$NearestTILen<30),1,0)
svData$IsTI=ifelse(svData$NearestTILen>=30&svData$IsDB==0,1,0)
svData$ClusterNone=ifelse(svData$ClusterCount==1,1,0)
svData$ClusterSmall=ifelse(svData$ClusterCount>1&svData$ClusterCount<=3,1,0)
svData$ClusterLarge=ifelse(svData$ClusterCount>3,1,0)

# run again but with ploidy >= 0.5
ploidyCutoffData = svData %>% filter(Ploidy >= 0.5)
nrow(ploidyCutoffData)
svFullData = svData
nrow(svFullData)
svData = ploidyCutoffData

# group all samples into relevant counts - for now lump all stressed-arm variants together
svSummary = (svData %>% group_by(SampleId)
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

# inputSigCount[[runNumber]] = 12
# inputSampleCounts[[runNumber]] = svSummary
# inputSigNames[[runNumber]] = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
# inputSigNamesNamed[[runNumber]] = c("01_LongDEL_INV", "02_Stressed", "03_LINE", "04_MidDEL", "05_BND_CN", "06_LongDUP", "07_BND_CL", "08_ShortINV", "09_ShortDEL", "10_DBs", "11_ShortDUP", "12_INV_CL")

sigCount = 11
sampleCounts = svSummary
sigNamesUnamed = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
sigNamesUnamed = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
sigNames = sigNamesUnamed
# sigNamesNamed = c("01_LongDEL_INV", "02_Stressed", "03_LINE", "04_MidDEL", "05_BND_CN", "06_LongDUP", "07_BND_CL", "08_", "09_ShortDEL", "10_DBs", "11_ShortDUP", "12_INV_CL")
# sigNamesNamed = c("01_LINE", "02_BND_CL", "03_Stressed", "04_ShortINV", "05_LongDEL_INV", "06_DBs", "07_ShortDEL", "08_LongDUP", "09_ShortDUP", "10_BND_CN", "11_MidDEL", "12_INV_CL")
sigNamesNamed = c("01_LINE", "02_BND_CL", "03_Stressed", "04_ShortINV", "05_LongDEL_INV", "06_DBs", "07_ShortDEL", "08_LongDUP", "09_ShortDUP", "10_BND_CN", "11_MidDEL")


# BUCKET AND SIGNATURE ANALYSIS

# all following code ought to be run-agnostic

#runNumber = 2
#View(runNumber)

#sampleCounts = inputSampleCounts[[runNumber]]
#sigCount = inputSigCount[[runNumber]]
sampleIds = sampleCounts[,1]
print(sigCount)
View(sampleCounts)
View(sampleIds)
View(sampleNames)

nmfMatrixData = svnmf::convert_summary_counts_to_nmf(sampleCounts)
bucketNames = svnmf::get_bucket_names(nmfMatrixData)
View(bucketNames)
#inputBucketNames[[runNumber]] = bucketNames
nmfMatrixData = svnmf::remove_bucket_names(nmfMatrixData)

View(nmfMatrixData)

#nmfEstimate <- nmf(nmfMatrixData, rank=6:15, method="brunet", nrun=4, seed=123456)
#plot(nmfEstimate)


# generate the actual NMF results
nmfResult <- nmf(nmfMatrixData, rank=sigCount, method="brunet", nrun=5, seed=123456, .opt='vp6')
save(nmfResult, file="~/logs/r_output/nmfResult_all_ploidyGT05.RData")


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

# experimenting
View(cancerTypes)


svnmf::plot_sig_samples(sampleSigData, "Breast") # all samples

cancerType = ""
cancerSigData = sampleSigData

if(cancerType != "") {
  cancerSigData = cancerSigData %>% filter(CancerType==cancerType)
}

cancerSampleSigData = cancerSigData %>% arrange(-SampleSvCount, SampleId) %>% select('SampleId', 'SigName', 'SvCount')

# only plot 50 samples at a time
numSamples = n_distinct(cancerSampleSigData$SampleId)
numSigs = n_distinct(cancerSampleSigData$SigName)
samplesPerPlot = 50
rowsPerPlot = samplesPerPlot * numSigs
plotCount = ceiling(numSamples/samplesPerPlot)

View(plotCount)
View(numSamples)
View(numSigs)

sigCancerPlots = list()
plotIndex = 1

for (n in 1:plotCount) {
  rowStart = ((n-1) * rowsPerPlot + 1)
  rowEnd = min((n * rowsPerPlot), numSamples * numSigs)

  if(cancerType == "")
  {
    title = "Sig SV Counts by Sample"
  }
  else
  {
    title = paste("Sig SV Counts by Sample for ", cancerType, sep="")
  }

  sampleSigPlot <- (ggplot(cancerSampleSigData[rowStart:rowEnd,], aes(x = reorder(SampleId, -SvCount), y = SvCount, fill = SigName))
                    + geom_bar(stat = "identity", colour = "black")
                    + labs(x = "", y = "SV Count by Sample")
                    + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                    + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                    + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))

  if(n == 1)
  {
    sampleSigPlot <- sampleSigPlot + ggtitle(title)
  }
  if(n < plotCount)
  {
    sampleSigPlot <- sampleSigPlot + theme(legend.position="none")
  }

  sigCancerPlots[[plotIndex]] <- sampleSigPlot

  if(plotIndex >=4)
  {
    print(paste("Printing 4 current plots with n=", n, ", plotIndex=", plotIndex, sep=""))
    multiplot(plotlist = sigCancerPlots, cols = 2)
    sigCancerPlots = list()
    plotIndex = 1
  }
  else
  {
    plotIndex = plotIndex + 1
  }

  if(n >= 8)
  {
    break
  }
}

multiplot(plotlist = sigCancerPlots, cols = 2)



myCOLORS = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
             "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a",
             "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
             "#dea185","#a0729d","#8a392f")

sampleSigPlot <- (ggplot(cancerSampleSigData[rowStart:rowEnd,], aes(x = reorder(SampleId, -SvCount), y = SvCount, fill = SigName))
                  + geom_bar(stat = "identity", colour = "black")
                  + labs(x = "", y = "SV Count by Sample")
                  + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                  + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                  + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))


# other NMF output

# estimate functions
nmfEstimate <- nmf(nmfMatrixData, rank=6:12, method="brunet", nrun=4, seed=123456)
plot(nmfEstimate)

# The most common approach is to choose the smallest rank for which cophenetic correlation coefficient starts decreasing.
# Another approach is to choose the rank for which the plot of the residual sum of squares (RSS) between the input matrix
# and its estimate shows an inflection point.”


contribsWithSigNames = contribution
rownames(contribsWithSigNames) <- sigNames
View(contribsWithSigNames)
plot_contribution_heatmap(contribsWithSigNames, cluster_samples=T)


cos_sim_samples_signatures = cos_sim_matrix(nmfMatrixData, cancer_signatures)

# Plot heatmap with specified signature order
plot_cosine_heatmap(nmfMatrixData, col_order = cosmic_order, cluster_rows = TRUE)


# Compare the reconstructed mutational profile with the original mutational profile:
# doesn't work since assumes SNV buckets
plot_compare_profiles(nmfMatrixData[,1], nmfResult$reconstructed[,1], profile_names = c("Original", "Reconstructed"), condensed = TRUE)





# prepare PDF for output
pdf(file=paste("~/logs/r_output/svnmf_", runNumber, ".pdf", sep = ""), height = 14, width = 20)

grid.arrange(ggparagraph(text = "Bucket Stats", size = 16, color = "black"),
             ggparagraph(text = "Top-3 Buckets", size = 16, color = "black"),
             tableGrob(sigBucketStats, rows=NULL),
             tableGrob(sigBucketTop3, rows=NULL),
             ncol = 2, nrow = 2, heights = c(0.1, 1), newpage = TRUE)


grid.arrange(tableGrob(head(bucketSummaryData, 50), rows=NULL), newpage = TRUE)

dev.off()

ggarrange(ggtexttable(sigBucketStats, rows = NULL),
          ggtexttable(sigBucketTop3, rows = NULL),
          ggtexttable(bucketSummaryData, rows = NULL),
          ggtexttable(leastContribBuckets, rows = NULL),
          ncol = 1, nrow = 4, heights = c(1, 1, 1, 1))

ggarrange(sigBucketsPlot, ncol = 1, nrow = 1)

ggarrange(ggtexttable(sigStats, rows = NULL), ncol = 1, nrow = 1)

dev.off()


someText = "Here is an explanation"
someText <- ggparagraph(text = someText, face = "italic", size = 11, color = "black")
print(someText)

sigDataTable <- ggtexttable(sigStats, rows = NULL, theme = ttheme("mOrange"))


# Arrange the plots on the same page

ggarrange(sigBucketsPlot, sigDataTable, someText,
          ncol = 1, nrow = 3,
          heights = c(1, 0.3, 0.1))

ggarrange(sigDataTable, someText,
          ncol = 1, nrow = 2,
          heights = c(0.8, 0.2))

ggarrange(sigBucketsPlot, sigDataTable, someText,
          ncol = 1, nrow = 3,
          heights = c(2, 0.5, 0.1))


# grid.table(sigStats, rows = NULL)

dev.off()







