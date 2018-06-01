library(purple);
library(RMySQL)
library(data.table)
library(IRanges)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(stringi)
library(devtools)
library(NMF)
library(MutationalPatterns)

detach("package:svnmf", unload=TRUE);
library(svnmf)

####################
### CODE STARTS HERE
####################

# SV data file
svData = read.csv('~/logs/CLUSTER_V16.csv')

# filter out multiple biopsy (approximately)
svData = svData %>% filter(!grepl("DRUP", SampleId)&!grepl("TIII", SampleId)&!grepl("TII", SampleId))

svData$IsLINE = ifelse(svData$LEStart!='false'|svData$LEEnd!='false',1,0)
svData$IsFS = ifelse(svData$FSStart!='false'|svData$FSEnd!='false',1,0)
svData$Length = ifelse(svData$Type=='BND'|svData$Type=='INS'|svData$ArmEnd!=svData$ArmStart, -1, svData$PosEnd-svData$PosStart)
svData$IsStressed = ifelse((svData$ArmCountStart>=svData$ArmExpStart*2.5&svData$ArmCountStart>=10)|(svData$ArmCountEnd>=svData$ArmExpEnd*2.5&svData$ArmCountEnd>=10),1,0)


svData$ClusterNone=ifelse(svData$ClusterCount==1,1,0)
svData$ClusterSmall=ifelse(svData$ClusterCount>1&svData$ClusterCount<=3,1,0)
svData$ClusterLarge=ifelse(svData$ClusterCount>3,1,0)
# View(svData)

# first group SVs by distinct chromosomal arm
svArmData = svData
svArmData$Chr = svArmData$ChrStart
svArmData$Arm = svArmData$ArmStart
svArmData$ArmCount = svArmData$ArmCountStart
svArmData$ArmExpected = svArmData$ArmExpStart
svArmData$Position = svArmData$PosStart

svArmBndData = svData %>% filter(Type=='BND')
svArmBndData$Chr = svArmBndData$ChrEnd
svArmBndData$Arm = svArmBndData$ArmEnd
svArmBndData$ArmCount = svArmBndData$ArmCountEnd
svArmBndData$ArmExpected = svArmBndData$ArmExpEnd
svArmBndData$Position = svArmBndData$PosEnd

# merge rows prior to arm grouping
combinedArmData = rbind(svArmData, svArmBndData)

# determine IsStressed using poisson distribution using expected vs actual SV counts per arm
combinedArmData$ArmExpected = ifelse(combinedArmData$ArmExpected>0,combinedArmData$ArmExpected,0.1)
combinedArmData$StressedPoissonProb = round(1 - ppois(combinedArmData$ArmCount - 1, combinedArmData$ArmExpected),4)


## STRESSED ARM ANALYSIS

combinedArmData$StressedRate = round(combinedArmData$ArmCount/combinedArmData$ArmExpected,1)
combinedArmData$StressedPPA = round(get_poisson_prob(combinedArmData$ArmCount, combinedArmData$ArmExpected),4)

# look at SV Arm Count vs SV Stressed Rate across cohort
armStressedCounts =  (combinedArmData %>% group_by(SampleId, Chr, Arm)
                      %>% summarise(SvCount=n(),
                                    StressedRate=min(StressedRate),
                                    StressedPP=min(StressedPoissonProb),
                                    StressedPPA=min(StressedPPA))
                      %>% arrange(SampleId, Chr, Arm))

View(armStressedCounts)

armStressedCounts$SvCountBucket = (ifelse(armStressedCounts$SvCount <= 20, armStressedCounts$SvCount,
                                          ifelse(armStressedCounts$SvCount <= 50, round(armStressedCounts$SvCount/5)*5,
                                                 ifelse(armStressedCounts$SvCount <= 100, round(armStressedCounts$SvCount/10)*10, round(armStressedCounts$SvCount/100)*100)
                                          )))

armStressedCounts$StressedRateBucket = (ifelse(armStressedCounts$StressedRate < 1, 0,
                                               ifelse(armStressedCounts$StressedRate <= 3.0, round(armStressedCounts$StressedRate/0.2)*0.2, 4)
))

View(armStressedCounts)

nrow(armStressedCounts %>% filter(SvCountBucket >= 10 & StressedPP <= 0.001))


armStressBucketed = (armStressedCounts %>% filter(SvCountBucket >= 6 & StressedRateBucket >= 1) %>% group_by(SvCountBucket, StressedRateBucket)
                     %>% summarise(ArmCount=n(),
                                   SvCount=sum(SvCount),
                                   AvgPoisson=round(sum(StressedPP)/n(),4),
                                   AvgPoissonA=round(sum(StressedPPA)/n(),4))
                     %>% arrange(SvCountBucket, StressedRateBucket))


View(armStressBucketed)

write.csv(armStressBucketed, "~/logs/r_output/armStressBucketed.csv")


stressedArmsByPP = (armStressedCounts %>% filter(SvCountBucket >= 10 & StressedPP <= 0.001) %>% group_by(SvCountBucket, StressedRateBucket)
                    %>% summarise(ArmCount=n(),
                                  SvCount=sum(SvCount))
                    %>% arrange(SvCountBucket, StressedRateBucket))

View(stressedArmsByPP)

sum(stressedArmsByPP$ArmCount)
sum(stressedArmsByPP$SvCount)

write.csv(stressedArmsByPP, "~/logs/r_output/stressedArmsByPP.csv")

# stressed arms by sample
stressedArmsBySample = (armStressedCounts %>% filter(SvCountBucket >= 10 & StressedPP <= 0.001) %>% group_by(SampleId)
                        %>% summarise(ArmCount=n(),
                                      SvCount=sum(SvCount))
                        %>% arrange(-ArmCount))

View(stressedArmsBySample)

write.csv(stressedArmsBySample, "~/logs/r_output/stressedArmsBySample.csv")

stressedBySampleAndArm = (stressedArmsBySample %>% group_by(ArmCount)
                          %>% summarise(SampleCount=n(),
                                        SvCount=sum(SvCount))
                          %>% arrange(-ArmCount))

View(stressedBySampleAndArm)

write.csv(stressedBySampleAndArm, "~/logs/r_output/stressedBySampleAndArm.csv")

# old method
combinedArmData$IsStressedOld = ifelse(combinedArmData$ArmCount>=10 & combinedArmData$ArmCount >= 2.5 * combinedArmData$ArmExpected,1,0)
nrow(combinedArmData %>% filter(IsStressedOld==1))

oldGroupedArms = combinedArmData %>% filter(IsStressedOld==1) %>% group_by(SampleId,Chr,Arm) %>% summarise(SampleArms=n())
nrow(oldGroupedArms)


## NMF STRESSED ARM ANALYSIS

combinedArmData$IsStressed = ifelse(combinedArmData$ArmCount >= 10 & combinedArmData$StressedPoissonProb <= 0.001, 1, 0)
saData = combinedArmData %>% filter(IsStressed==1)


# run if skipped the stressed arm logic above
# svData$ArmExpStart = ifelse(svData$ArmExpStart>0,svData$ArmExpStart,0.1)
# svData$ArmExpEnd = ifelse(svData$ArmExpEnd>0,svData$ArmExpEnd,0.1)
# svData$SPPStart = round(1 - ppois(svData$ArmCountStart - 1, svData$ArmExpStart),4)
# svData$SPPEnd = round(1 - ppois(svData$ArmCountEnd - 1, svData$ArmExpEnd),4)
# svData$IsStressed = ifelse((svData$ArmCountStart>=10&svData$SPPStart<=0.001)|(svData$ArmCountEnd>=10&svData$SPPEnd<=0.001),1,0)
# saData = svData %>% filter(IsStressed==1)

nrow(saData)



# prepare extra fields
saData$Length = ifelse(saData$Type=='BND'|saData$Type=='INS'|saData$ArmEnd!=saData$ArmStart, -1, saData$PosEnd-saData$PosStart)
saData$LengthBucket = ifelse(saData$Length <= 0, 0, 10**(round(log(saData$Length,10),0)))
saData$LengthBucket = ifelse(saData$LengthBucket > 1e6,1e6,saData$LengthBucket)
saData$LengthBucket = ifelse(saData$LengthBucket <= 1e3,1e3,saData$LengthBucket)
saData$IsDB=ifelse(saData$NearestDBLen>-1&(saData$NearestDBLen<saData$NearestTILen|saData$NearestTILen<30)&saData$NearestDBLen<=1e2,1,0)
saData$IsTI=ifelse(saData$NearestTILen>=30&saData$IsDB==0,1,0)
saData$ClusterSize = ifelse(saData$ClusterCount==1, 'CN', ifelse(saData$ClusterCount < 4, 'CS', 'CL'))

saData$LinkType=ifelse(saData$IsDB==1&saData$ClusterSize!='None','DB',ifelse(saData$IsTI==1&saData$ClusterSize!='None','TI','NA'))

saData$CN = (saData$AdjCNStart + saData$AdjCNEnd)*0.5
# saData$CN = ifelse(saData$CN <= 0.8, 1.0, ifelse(saData$CN > 64,64,saData$CN))
# saData$CNBucket = 2**(round(log(saData$CN,2),0))
saData$CNBucket = ifelse(saData$CN <= 0.8,1, ifelse(saData$CN < 2.4,2, ifelse(saData$CN <=8,8,16)))
View(saData %>% group_by(CNBucket) %>% summarise(SvCount=n()))

saByLength = (saData %>% group_by(LengthBucket) %>% summarise(SvCount=n()))
View(saByLength)

# optionally remove standard SVs from stressed arms
saData$IsStdSig = ifelse(saData$IsLINE==1|(saData$Type=='DEL'&saData$Length<=1e4)|(saData$Type=='DUP'&saData$Length<=1e5)
                            |(saData$Type=='DEL'&saData$Length>1e4&saData$Length<=5e5)|(saData$Type=='DUP'&saData$Length>1e5&saData$Length<=5e6)
                            |(saData$Type=='INV'&saData$Length<=1e4),1,0)

saData$ArmCountNI = saData$ArmCount - (saData$ArmExpected * 0.7) # based on observed rate of SI
saData$StressedPPNI = round(1 - ppois(saData$ArmCountNI - 1, saData$ArmExpected),4)
saData$IsStressedNI = ifelse(saData$ArmCountNI >= 10 & saData$StressedPPNI <= 0.001, 1, 0)

saNIData = saData %>% filter(IsStressedNI==1&IsStdSig==0)
nrow(saNIData)
saData = saNIData

# prepare arm stats
saSampleCounts = (saData %>% group_by(SampleId, Chr, Arm, LengthBucket, CNBucket, ClusterSize, LinkType)
            %>% summarise(SvCount=n(),
                          DelCount=sum(Type=='DEL'),
                          DupCount=sum(Type=='DUP'),
                          InvCount=sum(Type=='INV'),
                          BndCount=sum(Type=='BND')))

View(saSampleCounts)

# group without samples
saSummaryCounts = (saSampleCounts %>% group_by(LengthBucket, CNBucket, ClusterSize, LinkType)
                   %>% summarise(ArmsAffected=n(),
                                 SvCount=sum(SvCount),
                                 # SvAvgCount=round(sum(SvCount)/n(),1),
                                 DelCount=sum(DelCount),
                                 DelPerc=round(sum(DelCount)/sum(SvCount),3),
                                 DuplCount=sum(DupCount),
                                 DupPerc=round(sum(DupCount)/sum(SvCount),3),
                                 InvCount=sum(InvCount),
                                 InvPerc=round(sum(InvCount)/sum(SvCount),3),
                                 BndCount=sum(BndCount),
                                 BndPerc=round(sum(BndCount)/sum(SvCount),3),
                                 DBCount=sum(DBCount),
                                 TICount=sum(TICount))
                   %>% arrange(-SvCount))

View(saSummaryCounts)

sum(saSummaryCounts$ArmCount)
sum(saSummaryCounts$SvCount)


# NMF prep
saNmfSampleArmCounts = saData %>% group_by(SampleId, Chr, Arm, LengthBucket, CNBucket, ClusterSize, LinkType, Type) %>% summarise(Count=n())
View(saNmfSampleArmCounts)

# collapse arms
saNmfSampleCounts = saNmfSampleArmCounts %>% group_by(SampleId, LengthBucket, CNBucket, ClusterSize, LinkType, Type) %>% summarise(Count=sum(Count))
View(saNmfSampleCounts)

# convert groupings to concatenated buckets
saNmfBucketed = unite(saNmfSampleCounts, "Bucket", LengthBucket, CNBucket, ClusterSize, LinkType, Type, sep="_")
View(saNmfBucketed)
nrow(saNmfBucketed)
sum(saNmfBucketed$Count)

# bucket analysis
saNmfBucketed$CountBucket = 2**round(log(saNmfBucketed$Count,2),0)
saCountBucketed = saNmfBucketed %>% group_by(CountBucket) %>% summarise(Items=n(), SvCount=sum(Count))
View(saCountBucketed)

bucketCounts = saNmfBucketed %>% group_by(Bucket) %>% summarise(SampleCount=n(), SvCount=sum(Count))
View(bucketCounts)
nrow(bucketCounts)
write.csv(bucketCounts, "~/logs/r_output/bucketCounts.csv")

includedBuckets = bucketCounts %>% filter(SampleCount >= 80)
View(includedBuckets)
nrow(includedBuckets)
includedBucketNames = includedBuckets$Bucket

saNmfBucketedLimited = saNmfBucketed %>% filter(Bucket %in% includedBucketNames)
nrow(saNmfBucketedLimited)
sum(saNmfBucketedLimited$Count)


bndBucketCounts = saNmfBucketed %>% filter(grepl("BND", Bucket)) %>% group_by(Bucket) %>% summarise(SampleCount=n(), SvCount=sum(Count))
View(bndBucketCounts)

saNmfBucketed = saNmfBucketedLimited
nrow(saNmfBucketed)

# adjust length bucket names
saNmfBucketed$Bucket = stri_replace_all_fixed(saNmfBucketed$Bucket, "1000_", "1K_")
saNmfBucketed$Bucket = stri_replace_all_fixed(saNmfBucketed$Bucket, "10000_", "10K_")
saNmfBucketed$Bucket = stri_replace_all_fixed(saNmfBucketed$Bucket, "1e+05_", "100K_")
saNmfBucketed$Bucket = stri_replace_all_fixed(saNmfBucketed$Bucket, "1e+06_", "1M_")
saNmfBucketed$Bucket = stri_replace_all_fixed(saNmfBucketed$Bucket, "1e+07_", "10M_")
View(saNmfBucketed)

saNmfBucketed = within(saNmfBucketed, rm(CountBucket))

saNmfTemp1 = saNmfBucketed %>% spread(SampleId, Count)
View(saNmfTemp1)
nrow(saNmfTemp1)

# saNmfTemp2 = apply(saNmfTemp1, 2, function(x) ifelse(is.na(x), 0, x))

saNmfTemp2 = saNmfTemp1

# try this instead:
sampleResult[is.na(sampleResult)] <- 0

# for(i in 1:nrow(saNmfTemp2))
# {
#   for(j in 2:ncol(saNmfTemp2))
#   {
#     if(is.na(saNmfTemp2[i,j]))
#     {
#       saNmfTemp2[i,j] = 0
#     }
#   }
# }
# nrow(saNmfTemp2)

# View(saNmfTemp2)
# write.csv(saNmfTemp2, "~/logs/r_output/saNMF2.csv")

# extract bucket names and then remove them
saNmfTemp3 = as.data.frame(saNmfTemp2)
bucketNames = saNmfTemp3$Bucket
View(bucketNames)
saNmfMatrixData = within(saNmfTemp3, rm(Bucket))
nrow(saNmfMatrixData)
#View(nmfMatrixData)


saNmfEstimate <- nmf(saNmfMatrixData, rank=6:12, method="brunet", nrun=4, seed=123456, .opt='vp6')
plot(saNmfEstimate)
save(saNmfEstimate, file="~/logs/r_output/nmfEstimate_run3_stressedArms.RData")

# generate the actual NMF results
sigCount = 10
saNmfResult <- nmf(saNmfMatrixData, rank=sigCount, method="brunet", nrun=5, seed=123456, .opt='vp6')
save(saNmfResult, file="~/logs/r_output/nmfResult_run3_stressedArms.RData")


# NMF EVALUATION

detach("package:svnmf", unload=TRUE);
library(svnmf)


sampleIds = colnames(saNmfMatrixData)
sampleNames = sampleIds
View(sampleIds)
View(bucketNames)
signatures = NMF::basis(saNmfResult)
contribution = NMF::coef(saNmfResult)
# sampleNames = colnames(contribution)

sigNames = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
sigNamesNamed[[runNumber]] = c("01_LongDEL_INV", "02_Stressed", "03_LINE", "04_MidDEL", "05_BND_CN", "06_LongDUP", "07_BND_CL", "08_ShortINV", "09_ShortDEL", "10_DBs", "11_ShortDUP", "12_INV_CL")
View(sigNames)

# Bucket Evaluation
sigBucketData = svnmf::get_bucket_data(signatures, contribution, bucketNames)
View(sigBucketData)
sigBucketStats = svnmf::get_sig_bucket_stats(sigBucketData)
View(sigBucketStats)

# Signature Discovery, by looking at relative contribution of buckets

# report on top 3 contributing buckets only
# sigBucketTop3 = svnmf::get_top_buckets(sigBucketData, 3)
# View(sigBucketTop3)
sigBucketTopN = svnmf::get_top_buckets(sigBucketData, 5, 0.9)
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
sigNamesNamed = c("LongDEL_INV", "Stressed", "LINE", "MidDEL", "BND_CN", "LongDUP", "BND_CL", "ShortINV", "ShortDEL", "BND_DB", "ShortDUP", "INV_CL")
View(sigNamesNamed)
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
# sampleSigData = within(sampleSigData, rm(SampleIdStripped))
sampleSigData$CancerType = ifelse(is.na(sampleSigData$CancerType), 'N/A', paste(sampleSigData$CancerType, sep=""))

sampleSvCounts = sampleSigData %>% group_by(SampleId) %>% summarise(SampleSvCount=sum(SvCount))
sampleSigData = (merge(sampleSigData, sampleSvCounts, by.x="SampleId",by.y="SampleId",all.x=TRUE))
View(sampleSigData)

# cancerTypes = sampleSigData %>% group_by(CancerType) %>% summarise(Count=n())
# View(cancerTypes)


detach("package:svnmf", unload=TRUE);
library(svnmf)


# DATA OUTPUT TO PDF
View(runNumber)
pdf(file=paste("~/logs/r_output/svnmf_", runNumber, ".pdf", sep = ""), height = 14, width = 20)

par(mar=c(1,1,1,1))

# 1. NMF estimate data
# plot(nmfEstimate)

# 2. bucket data
title = textGrob("Bucket Summary Data & Top-N Buckets", gp=gpar(fontface="bold", fontsize=16))
grid.arrange(tableGrob(head(bucketSummaryData, 40), rows=NULL),
             tableGrob(sigBucketTopN, rows=NULL),
             ncol = 2, newpage = TRUE, top=title)


bucketSummaryPlot = svnmf::get_bucket_summary_plot(bucketSummaryData)
print(bucketSummaryPlot)
grid.arrange(bucketSummaryPlot, ncol = 1, nrow = 1, newpage = TRUE)


title = textGrob("Signature-Bucket Stats & Least Important Buckets", gp=gpar(fontface="bold", fontsize=16))
grid.arrange(tableGrob(sigBucketStats, rows=NULL),
             tableGrob(leastContribBuckets, rows = NULL),
             nrow = 2, newpage = TRUE, top=title)


# 3. default signature-bucket plot
sigBucketsPlot = svnmf::get_bucket_signatures_plot(bucketNames, signatures)
print(sigBucketsPlot)
grid.arrange(sigBucketsPlot, ncol = 1, nrow = 1, newpage = TRUE)

# 4. Signature data
title = textGrob("Signature Stats", gp=gpar(fontface="bold", fontsize=16))
grid.arrange(tableGrob(sigStats, rows = NULL), ncol = 1, nrow = 1, top=title, newpage = TRUE)

# 5. Sigs with Samples by cancer type
svnmf::plot_sig_samples(sampleSigData, "") # all samples

for(cancerType in cancerTypes$CancerType)
{
  if(!is.na(cancerType))
  {
    svnmf::plot_sig_samples(sampleSigData, cancerType)
  }
}

dev.off()


# Testing if stressed arms are still stressed when standard independent signature (SI) SVs are removed

indSvData = combinedArmData

indArmData = (indSvData %>% group_by(SampleId, Chr, Arm, IsStressed)
                 %>% summarise(SvCount=n(),
                               ShortDel=sum(IsLINE==0&Type=='DEL'&Length<=1e4),
                               ShortDup=sum(IsLINE==0&Type=='DUP'&Length<=1e5),
                               MidDel=sum(IsLINE==0&Type=='DEL'&Length>1e4&Length<=5e5),
                               LongDup=sum(IsLINE==0&Type=='DUP'&Length>1e5&Length<=5e6),
                               ShortInv=sum(IsLINE==0&Type=='INV'&Length<=1e4),
                               LINE=sum(IsLINE==1),
                               BndNC=sum(IsLINE==0&Type=='BND'&ClusterNone==1))
                 %>% arrange(-SvCount))

View(indArmData)

indSummaryData = (indArmData %>% group_by(IsStressed)
                 %>% summarise(ArmCount=n(),
                               SvCount=sum(SvCount),
                               Sv_Avg=round(sum(SvCount)/n(),1),
                               ShortDel=sum(ShortDel),
                               ShortDel_Avg=round(sum(ShortDel)/n(),1),
                               ShortDel_Perc=round(sum(ShortDel)/SvCount,3),
                               ShortDup=sum(ShortDup),
                               ShortDup_Avg=round(sum(ShortDup)/n(),1),
                               ShortDup_Perc=round(sum(ShortDup)/SvCount,3),
                               MidDel=sum(MidDel),
                               MidDel_Avg=round(sum(MidDel)/n(),1),
                               MidDel_Perc=round(sum(MidDel)/SvCount,3),
                               LongDup=sum(LongDup),
                               LongDup_Avg=round(sum(LongDup)/n(),1),
                               LongDup_Perc=round(sum(LongDup)/SvCount,3),
                               ShortInv=sum(ShortInv),
                               ShortInv_Avg=round(sum(ShortInv)/n(),1),
                               ShortInv_Perc=round(sum(ShortInv)/SvCount,3),
                               LINE=sum(LINE),
                               LINE_Avg=round(sum(LINE)/n(),1),
                               LINE_Perc=round(sum(LINE)/SvCount,3),
                               BndNC=sum(BndNC),
                               BndNC_Avg=round(sum(BndNC)/n(),1),
                               BndNC_Perc=round(sum(BndNC)/SvCount,3))
                 %>% arrange(-SvCount))

View(indSummaryData)

# now assess changes to stress-rates after removing these standard independent signature SVs
indSvData$IsStdSig = ifelse(indSvData$IsLINE==1|(indSvData$Type=='DEL'&indSvData$Length<=1e4)|(indSvData$Type=='DUP'&indSvData$Length<=1e5)
                            |(indSvData$Type=='DEL'&indSvData$Length>1e4&indSvData$Length<=5e5)|(indSvData$Type=='DUP'&indSvData$Length>1e5&indSvData$Length<=5e6)
                            |(indSvData$Type=='INV'&indSvData$Length<=1e4)|(indSvData$Type=='BND'&indSvData$ClusterNone==1),1,0)

indAssArmData = (indSvData %>% group_by(SampleId, Chr, Arm, IsStressed)
              %>% summarise(SvCount=n(),
                            StdSig=sum(IsStdSig==1),
                            NonBndStdSig=sum(IsStdSig==1&Type!='BND'),
                            BndStdSig=sum(IsStdSig==1&Type=='BND'),
                            NonStdSig=sum(IsStdSig==0),
                            ArmExpected=min(ArmExpected),
                            ArmCount=min(ArmCount))
              %>% arrange(-SvCount))

View(indAssArmData)

indAssSummaryData = (indAssArmData %>% group_by(IsStressed)
                  %>% summarise(ArmCount=n(),
                                SvCount=sum(SvCount),
                                Sv_Avg=round(sum(SvCount)/n(),1),
                                StdSig=sum(StdSig),
                                StdSig_Avg=round(sum(StdSig)/n(),1),
                                StdSig_Perc=round(sum(StdSig)/SvCount,3),
                                NonStdSig=sum(NonStdSig),
                                NonStdSig_Avg=round(sum(NonStdSig)/n(),1),
                                NonStdSig_Perc=round(sum(NonStdSig)/SvCount,3))
                  %>% arrange(-SvCount))

View(indAssSummaryData)

# Non-Independent (NI) - ie removing the effects of SIs SVs
indAssArmData$ArmCountNI = indAssArmData$ArmCount - (indAssArmData$ArmExpected * 0.7) # based on observed rate of SI
indAssArmData$StressedPPNI = round(1 - ppois(indAssArmData$ArmCountNI - 1, indAssArmData$ArmExpected),4)
indAssArmData$StressedRateNI = round(indAssArmData$ArmCountNI/indAssArmData$ArmExpected,1)
indAssArmData$IsStressedNI = ifelse(indAssArmData$ArmCountNI >= 10 & indAssArmData$StressedPPNI <= 0.001, 1, 0)

View(indAssArmData %>% filter(IsStressedNI==0&IsStressed==1))

# now look at how the stressed arm rate has changed
indNIStressedComparison = (indAssArmData %>% group_by(IsStressed,IsStressedNI)
                     %>% summarise(ArmCount=n(),
                                   SvCount=sum(SvCount),
                                   StdSigCount=sum(StdSig),
                                   NonStdSigCount=sum(NonStdSig),
                                   StressedPPNI_Avg=round(sum(StressedPPNI)/n(),3),
                                   StressedRate_Avg=round(sum(StressedRateNI)/n(),3))
                     %>% arrange(-SvCount))

View(indNIStressedComparison)


# look at SV Arm Count vs SV Stressed Rate across cohort
armStressedCounts =  (combinedArmData %>% group_by(SampleId, Chr, Arm)
                      %>% summarise(SvCount=n(),
                                    StressedRate=min(StressedRate),
                                    StressedPP=min(StressedPoissonProb),
                                    StressedPPA=min(StressedPPA))
                      %>% arrange(SampleId, Chr, Arm))

View(armStressedCounts)


