library(VGAM)   # zeta function


# Cosine Similarity
cosine_sim<-function(vec1, vec2)
{
  cosineSim = (vec1 %*% vec2) / (sqrt(vec1 %*% vec1)*sqrt(vec2 %*% vec2))
  return (cosineSim)
}

calc_cosine_sims<-function(dataMatrix, cssCutoff, itemName, logMatches=F)
{
  # compare all entries in the columns
  numItems = ncol(dataMatrix)
  itemNames = colnames(dataMatrix)

  cssResults = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(cssResults) <- c(paste(itemName,'1',sep=''), paste(itemName,'2',sep=''), "CSS")

  for(i in 1:numItems)
  {
    data1 = dataMatrix[,i]

    for(j in i+1:numItems)
    {
      if(j>numItems)
        break

      data2 = dataMatrix[,j]

      css = cosine_sim(data1,data2)

      if(css >= cssCutoff)
      {
        if(logMatches)
        {
          print(paste(i, "=", itemNames[i], ", count=", sum(data1), ", ", j, "=", itemNames[j], ", count=", sum(data2), ", CSS=", round(css,4), sep=''))
        }

        rowIndex = nrow(cssResults)+1
        cssResults[rowIndex,1] = itemNames[i]
        cssResults[rowIndex,2] = itemNames[j]
        cssResults[rowIndex,3] = round(css,6)
      }
    }

    if(i %% 100 == 0)
    {
      print(paste("processed ", i, sep=''))
    }
  }

  return (cssResults)
}


# bucket analyser output
baSigs = as.matrix(read.csv(file="~/dev/nmf/logs/snv_ba_sigs.csv", stringsAsFactors=F))
baContribs = as.matrix(read.csv(file="~/dev/nmf/logs/snv_ba_contribs.csv", stringsAsFactors=F))
baSigCount = nrow(baContribs)
baSigsNamesStr = get_signame_list(baSigCount, T)
baSigsNamesNum = get_signame_list(baSigCount, F)
colnames(baSigs) = baSigsNamesNum

View(baSigs)
View(baContribs[,1:10])

sampleSigData = get_sig_data(baSigs, baContribs, baSigsNamesStr, colnames(baContribs))
sampleSigData = sampleSigData %>% filter(!is.na(SigPercent))
sigStats = get_sig_stats(sampleSigData)
View(sigStats)
View(sampleSigData)

View(sampleCancerTypes)
sampleSigData = merge(sampleSigData, sampleCancerTypes,by.x="SampleId",by.y="SampleId",all.x=TRUE)

View(sampleSigData %>% filter(SigPercent==1))

View(sampleSigData %>% filter(SigPercent==1) %>% group_by(CancerType,SigName) %>% count())


snvSampleCounts = read.csv(file="~/data/r_data/snv_nmf_sample_counts.csv")
View(snvSampleCounts)
nrow(snvSampleCounts)
ncol(snvSampleCounts)
View(snvSampleCounts %>% group_by(SampleId) %>% count())

sampleCancerTypes = highestPurityCohortSummary %>% select(sampleId,cancerType)
colnames(sampleCancerTypes) = c("SampleId", "CancerType")
View(sampleCancerTypes)

snvSampleCounts = merge(snvSampleCounts,sampleCancerTypes,all.X=T)
View(snvSampleCounts)

cancerTypesList = unique(sampleCancerTypes$CancerType)
View(cancerTypesList)

nrow(snvSampleCounts %>% filter(Count==0))

View(snvBucketNames)

percInterval = 0.05
percValues = seq(0, 1, percInterval)

View(percValues)

bucketStatsCols = c("CancerType", "Bucket", "Samples", "Min", "Max", "Mean", "Median")

for(i in 1:length(percValues))
{
  colIndex = length(bucketStatsCols)+1
  bucketStatsCols[colIndex] = paste("P_", percValues[i]*100, sep='')
}
View(bucketStatsCols)

bucketStats = data.frame(matrix(ncol=7+length(percValues), nrow = 0))
colnames(bucketStats) = bucketStatsCols
View(bucketStats)

ctTmp = c("Skin", "Prostate", "Breast", "Lung", "Colon/Rectum", "Kidney")
View(cancerTypesList)
cancerTypesList = ctTmp

for(cancerType in cancerTypesList)
{
  scData = snvSampleCounts %>% filter(CancerType==cancerType)

  print(paste("CancerType=", cancerType, ", sampeCount=", nrow(scData), sep=''))

  # now gather stats for each bucket
  for(bucketName in snvBucketNames$Bucket)
  {
    bucketCounts = scData %>% filter(Bucket==bucketName) %>% arrange(Count)

    sampleCount = nrow(bucketCounts)

    rowIndex = nrow(bucketStats)+1
    bucketStats[rowIndex,1] = cancerType
    bucketStats[rowIndex,2] = bucketName
    bucketStats[rowIndex,3] = sampleCount
    bucketStats[rowIndex,4] = min(bucketCounts$Count)
    bucketStats[rowIndex,5] = max(bucketCounts$Count)
    bucketStats[rowIndex,6] = round(mean(bucketCounts$Count),1)
    bucketStats[rowIndex,7] = round(median(bucketCounts$Count),1)

    rowsPerInterval = floor(percInterval * sampleCount)
    for(i in 1:length(percValues))
    {
      startRow = (i-1)*rowsPerInterval+1
      endRow = min(startRow + rowsPerInterval-1, sampleCount)
      # print(paste(i, ": start=", startRow, ", endRow=", endRow, sep=''))
      bcSubset = slice(bucketCounts, startRow:endRow)
      bcMean = mean(bcSubset$Count)
      bucketStats[rowIndex,7+i] = round(bcMean,0)
    }
  }
}

View(bucketStats)

bucketStats$Total = bucketStats$Mean*bucketStats$Samples

write.csv(bucketStats, "~/logs/r_output/snvBucketPercentiles.csv", row.names=F, quote=F)

# find lowest rate of change buckets
bucketStatsLowCTs = (bucketStats %>% filter(CancerType!="Skin"&CancerType!="Lung") %>% group_by(Bucket)
                     %>% summarise(P_0=sum(P_0), P_5=sum(P_5), P_10=sum(P_10), P_15=sum(P_15), P_20=sum(P_20), P_25=sum(P_25), P_30=sum(P_30),
                                   P_35=sum(P_35), P_40=sum(P_40), P_45=sum(P_45), P_50=sum(P_50), P_55=sum(P_55), P_60=sum(P_60), P_65=sum(P_65), P_70=sum(P_70),
                                   P_75=sum(P_75), P_80=sum(P_80), P_85=sum(P_85), P_90=sum(P_90), P_95=sum(P_95), P_100=sum(P_100)))


# compute ratios
bucketStatsLowCTs$R_5_0 = round(bucketStatsLowCTs$P_5 / bucketStatsLowCTs$P_0,3)
bucketStatsLowCTs$R_10_5 = round(bucketStatsLowCTs$P_10 / bucketStatsLowCTs$P_5,3)
bucketStatsLowCTs$R_15_10 = round(bucketStatsLowCTs$P_15 / bucketStatsLowCTs$P_10,3)
bucketStatsLowCTs$R_20_15 = round(bucketStatsLowCTs$P_20 / bucketStatsLowCTs$P_15,3)
bucketStatsLowCTs$R_25_20 = round(bucketStatsLowCTs$P_25 / bucketStatsLowCTs$P_20,3)
bucketStatsLowCTs$R_30_25 = round(bucketStatsLowCTs$P_30 / bucketStatsLowCTs$P_25,3)
bucketStatsLowCTs$R_35_30 = round(bucketStatsLowCTs$P_35 / bucketStatsLowCTs$P_30,3)
bucketStatsLowCTs$R_40_35 = round(bucketStatsLowCTs$P_40 / bucketStatsLowCTs$P_35,3)
bucketStatsLowCTs$R_45_40 = round(bucketStatsLowCTs$P_45 / bucketStatsLowCTs$P_40,3)
bucketStatsLowCTs$R_50_45 = round(bucketStatsLowCTs$P_50 / bucketStatsLowCTs$P_45,3)
bucketStatsLowCTs$R_55_50 = round(bucketStatsLowCTs$P_55 / bucketStatsLowCTs$P_50,3)
bucketStatsLowCTs$R_60_55 = round(bucketStatsLowCTs$P_60 / bucketStatsLowCTs$P_55,3)
bucketStatsLowCTs$R_65_60 = round(bucketStatsLowCTs$P_65 / bucketStatsLowCTs$P_60,3)
bucketStatsLowCTs$R_70_65 = round(bucketStatsLowCTs$P_70 / bucketStatsLowCTs$P_65,3)
bucketStatsLowCTs$R_75_70 = round(bucketStatsLowCTs$P_75 / bucketStatsLowCTs$P_70,3)
bucketStatsLowCTs$R_80_75 = round(bucketStatsLowCTs$P_80 / bucketStatsLowCTs$P_75,3)
bucketStatsLowCTs$R_85_80 = round(bucketStatsLowCTs$P_85 / bucketStatsLowCTs$P_80,3)
bucketStatsLowCTs$R_90_85 = round(bucketStatsLowCTs$P_90 / bucketStatsLowCTs$P_85,3)
bucketStatsLowCTs$R_95_90 = round(bucketStatsLowCTs$P_95 / bucketStatsLowCTs$P_90,3)
bucketStatsLowCTs$R_100_95 = round(bucketStatsLowCTs$P_100 / bucketStatsLowCTs$P_95,3)

bucketStatsLowCTs$AvgRatio = 0
colIndex = ncol(bucketStatsLowCTs)
colnames(bucketStatsLowCTs)

startCol = 25
endCol = 40
for(i in 1:nrow(bucketStatsLowCTs))
{
  ratioTotal = 0
  for(j in startCol:endCol)
  {
    ratioTotal = ratioTotal + bucketStatsLowCTs[i,j]
  }

  bucketStatsLowCTs[i,colIndex] = ratioTotal / (endCol - startCol + 1)
}


View(bucketStatsLowCTs)
View(cosmicSignatures[,1])
cosmicSigsNamed = cosmicSignatures
cosmicSigsNamed = round(cosmicSigsNamed,4)
cosmicSigsNamed = cbind(snvBucketNames,cosmicSigsNamed)
View(cosmicSigsNamed)


View(snvSampleCounts)

# extract specific cancer types for NMF

skinMatrixData2 = snvSampleCounts %>% filter(CancerType=="Skin") %>% select(SampleId,Bucket,Count) %>% spread(SampleId,Count)
skinMatrixData2[is.na(skinMatrixData2)] = 0
View(skinMatrixData2[,1:10])
ncol(skinMatrixData2)
write.csv(skinMatrixData2, "~/data/r_data/skin_matrix_data.csv", row.names=F, quote=F)
skinMatrixData = within(skinMatrixData2, rm(Bucket))
View(skinMatrixData)
ncol(skinMatrixData)
write.csv(skinMatrixData, "~/data/r_data/skin_matrix_data_all.csv", row.names=F, quote=F)

kidneyMatrixData = snvSampleCounts %>% filter(CancerType=="Kidney") %>% select(SampleId,Bucket,Count) %>% spread(SampleId,Count)
kidneyMatrixData[is.na(kidneyMatrixData)] = 0
View(kidneyMatrixData[,1:10])
ncol(kidneyMatrixData)
kidneyMatrixData = within(kidneyMatrixData, rm(Bucket))
write.csv(kidneyMatrixData, "~/data/r_data/kidney_matrix_data.csv", row.names=F, quote=F)

# MMR samples
load("~/data/r_data/highestPurityCohortSummary.RData")
View(highestPurityCohortSummary)
View(highestPurityCohortSummary %>% group_by(msiStatus) %>% count())
msiSamples = highestPurityCohortSummary %>% filter(msiStatus=="MSI")
View(msiSamples)
mmrSampleCounts = snvSampleCounts %>% filter(SampleId %in% msiSamples$sampleId) %>% select(SampleId,Bucket,Count)
mmrMatrixData = mmrSampleCounts %>% spread(SampleId,Count)
mmrMatrixData[is.na(mmrMatrixData)] = 0
View(mmrMatrixData)
ncol(mmrMatrixData)
mmrMatrixData = within(mmrMatrixData, rm(Bucket))
write.csv(mmrMatrixData, "~/data/r_data/mmr_matrix_data.csv", row.names=F, quote=F)


dr022Samples = read.csv("~/data/DR-022_metadata.tsv", sep='\t')
hpcSamples = highestPurityCohortSummary %>% filter(sampleId %in% dr022Samples$sampleId)
# hpcSamples = highestPurityCohortSummary




# just analyse similar samples - eg skin, high ML
skinSampleCounts = snvSampleCounts %>% filter(CancerType=="Skin") %>% arrange(-Count)
View(skinSampleCounts)
skinSampleCounts = merge(skinSampleCounts, origSampleCounts, all.x=T)
View(skinSampleCounts %>% group_by(SampleId) %>% summarise(Count=sum(Count)))

skinBucketCounts = snvSampleCounts %>% group_by(Bucket) %>% summarise(BucketTotal=sum(Count))
View(skinBucketCounts)
skinSampleCounts = merge(skinSampleCounts, skinBucketCounts, all.x=T)


skinSampleCountsHighML = skinSampleCounts %>% filter(SampleTotal>=20000,BucketTotal>=500000) # filtering small buckets as well
View(skinSampleCountsHighML)
View(skinSampleCountsHighML %>% group_by(Bucket) %>% summarise(BucketTotal=sum(Count)))

skinMatrixData = skinSampleCountsHighML %>% select(SampleId,Bucket,Count) %>% spread(SampleId,Count)
View(skinMatrixData)
skinMatrixData[is.na(skinMatrixData)] = 0
skinMatrixData = skinMatrixData %>% select(-Bucket)
write.csv(skinMatrixData, "~/data/r_data/skin_matrix_data.csv", row.names=F, quote=F)


# same again for mutational load on each sample
mutLoadStatsCols = c("CancerType", "Samples", "Mean", "Median")

for(i in 1:length(percValues))
{
  colIndex = length(mutLoadStatsCols)+1
  mutLoadStatsCols[colIndex] = paste("P_", percValues[i]*100, sep='')
}
View(mutLoadStatsCols)

mutLoadStats = data.frame(matrix(ncol=length(mutLoadStatsCols), nrow = 0))
colnames(mutLoadStats) = mutLoadStatsCols
View(mutLoadStats)

ctTmp = c("Skin", "Prostate", "Breast", "Lung", "Colon/Rectum", "Kidney")
View(cancerTypesList)
cancerTypesList = ctTmp

for(cancerType in cancerTypesList)
{
  scData = snvSampleCounts %>% filter(CancerType==cancerType) %>% group_by(SampleId) %>% summarise(Count=sum(Count)) %>% arrange(Count)

  print(paste("CancerType=", cancerType, ", sampleCount=", nrow(scData), sep=''))

  sampleCount = nrow(scData)

  rowIndex = nrow(mutLoadStats)+1
  mutLoadStats[rowIndex,1] = cancerType
  mutLoadStats[rowIndex,2] = sampleCount
  mutLoadStats[rowIndex,3] = round(mean(scData$Count),1)
  mutLoadStats[rowIndex,4] = round(median(scData$Count),1)

  rowsPerInterval = floor(percInterval * sampleCount)
  for(i in 1:length(percValues))
  {
    startRow = (i-1)*rowsPerInterval+1
    endRow = min(startRow + rowsPerInterval-1, sampleCount)
    # print(paste(i, ": start=", startRow, ", endRow=", endRow, sep=''))
    scSubset = slice(scData, startRow:endRow)
    scMean = mean(scSubset$Count)
    mutLoadStats[rowIndex,4+i] = round(scMean,0)
  }
}

View(mutLoadStats)





# plot bucket % per sample for each cancer type
origSampleCounts = snvSampleCounts %>% group_by(SampleId) %>% summarise(SampleTotal=sum(Count))
sampleBucketData = get_sample_bucket_data(snvMatrixData, origSampleCounts, snvBucketNames)
sampleBucketData = merge(sampleBucketData, sampleCancerTypes, all.x=T)
View(sampleBucketData)
View(sampleCancerTypes)
View(origSampleCounts)


cancerTypesList = unique(sampleCancerTypes$CancerType)


outputFile = "~/logs/r_output/pdfs/snv_bucket_perc_by_cancer_abs.pdf"

# DATA OUTPUT TO PDF
pdf(file=outputFile, height = 14, width = 20)

par(mar=c(1,1,1,1))

for(cancerType in cancerTypesList)
{
  cancerSBData = sampleBucketData %>% filter(CancerType==cancerType)
  plotTitle = paste("Bucket Counts for ", cancerType, sep='')
  plot_sample_bucket_contrib(cancerSBData, nrow(snvBucketNames), "SNV", 0, plotTitle, T)
}

dev.off()





snvMatrixData = read.csv(file="~/data/r_data/snv_nmf_matrix_data.csv")
nrow(snvMatrixData)
ncol(snvMatrixData)

# append bucket names back on, or just keep them handy
View(snvBucketNames)


# order each bucket's counts from lowest to highest to get idea of distribution
snvBucketCount = nrow(snvMatrixData)

nrow(snvSampleCounts %>% filter(Count>=10))
nrow(snvSampleCounts %>% filter(Count>=20))
nrow(snvSampleCounts %>% filter(Count>=50))
nrow(snvSampleCounts %>% filter(Count>=100))


# power distribution analysis
snvBucketCountSize = 100
snvSampleCounts$CountBucket = round(snvSampleCounts$Count/snvBucketCountSize)*snvBucketCountSize

View(snvSampleCounts %>% group_by(Bucket) %>% summarise(VarCount=sum(Count)))

sd(tmpBucketCounts$Count)

tmpBucketCounts = snvSampleCounts %>% filter(Bucket=="C>T_TCC") %>% arrange(Count)
View(tmpBucketCounts)
View(snvSampleCounts %>% filter(Bucket=="C>T_TCT") %>% group_by(CountBucket) %>% summarise(Count=n()))

View(tmpBucketCounts %>% group_by(round(Count/5,0)*5) %>% count())
View(tmpBucketCounts %>% group_by(round(Count/20,0)*20) %>% count())

tmpBucketCounts$CountRnd = round(tmpBucketCounts$Count/5,0)*5
bcDistrib = tmpBucketCounts %>% group_by(Count) %>% summarise(Freq=n()) %>% arrange(Freq)
bcDistribRnd = tmpBucketCounts %>% group_by(CountRnd) %>% summarise(Freq=n()) %>% arrange(Freq)
View(bcDistrib)

distPlot = (ggplot(data = bcDistribRnd %>% filter(CountRnd<1000), aes(x = CountRnd))
            + geom_line(aes(y=Freq, colour='Frequency'))
            # + scale_x_log10()
            # + coord_cartesian(xlim = c(-1e4, 1e6))
            # + scale_x_continuous()
            + ylab("Count group frequency"))

print(distPlot)


logDistrib = bcDistrib %>% select(Freq,Count)
logDistrib$Count = log(logDistrib$Count)
logDistrib$Freq = log(logDistrib$Freq)
View(logDistrib)
View(logDistrib %>% filter(Freq>0))



# evaluate simulated sample bucket counts
simSampleBucketData = as.matrix(read.csv("~/dev/nmf/logs/snv_sim_sample_counts.csv", stringsAsFactors=F))
simSampleBucketData = cbind(snvBucketNames,simSampleBucketData)
View(simSampleBucketData[,1:20])

simSampleBucketCounts = simSampleBucketData %>% gather("SampleId", "Count", -Bucket)
View(simSampleBucketCounts)

# look across just C->T buckets
CtoTBucketCounts = simSampleBucketCounts %>% filter(grepl("C>T",Bucket)) %>% group_by(SampleId) %>% summarise(Count=sum(Count)) %>% arrange(Count)
View(CtoTBucketCounts)
countRound = 10000
CtoTBucketCounts$CountRnd = round(CtoTBucketCounts$Count/countRound,0)*countRound
CtoTBucketCountsDist = CtoTBucketCounts %>% group_by(CountRnd) %>% summarise(Freq=n())
View(CtoTBucketCountsDist)
View(CtoTBucketCounts %>% filter(Bucket=="C>T_TCT") %>% group_by(CountBucket) %>% summarise(Count=n()))


distPlot = (ggplot(data = CtoTBucketCountsDist %>% filter(CountRnd<1e5), aes(x = CountRnd))
            + geom_line(aes(y=Freq, colour='Frequency'))
            # + scale_x_log10()
            # + coord_cartesian(xlim = c(-1e4, 1e6))
            # + scale_x_continuous()
            + ylab("Count group frequency"))

print(distPlot)





# log-likelihood comparison of sample counts



View(tmpKidney %>% select(CPCT02030356,CPCT02020580))
View(tmpKidney)
css = cosine_sim(vec1, vec2)
print(css)
colnames(tmpKidney)
tmp2 = cbind(tmpKidney[,9],tmpKidney[,12])
View(tmp2)



# compare CPCT02030356 (12) with CPCT02020580 (9) and CPCT02230015
tmpKidney = as.data.frame(snvSimMatrixData)
vec1 = tmpKidney$CPCT02030356
vec2 = tmpKidney$CPCT02020580
ratio = sum(vec2)/sum(vec1)
print(ratio)
vec3 = vec1
vec3 = vec3 * ratio
vec3 = round(vec3)
tmp3 = as.data.frame(cbind(vec3,vec2))
colnames(tmp3) = c("Sam1", "Sam2")
View(tmp3)

# now calc poisson probs for each pair and log them
poisLLDiff = sum(dpois(x=tmp3$Sam2, lambda=tmp3$Sam1, log=T))
print(poisLLDiff)
poisLLSame = sum(dpois(x=tmp3$Sam2, lambda=tmp3$Sam2, log=T))
print(poisLLSame)

prob = 1 - pchisq(2 * (poisLLSame - poisLLDiff), df=c(1,3))
prob = 1 - pchisq(2 * (poisLLDiff - poisLLSame), df=c(1,3))
print(prob)

prob = 1 - pchisq(1, df=c(1,3))
print(prob)

prob = 1 - pchisq(2 * (-400 - c(-200,-300)), df=c(1,3))
print(prob)


prob = 1 - pchisq(1, df=c(1,2))
print(prob)


if (outp > 1) {
  message("[4] Running dNdSloc...")

  selfun_loc = function(j) {
    y = as.numeric(genemuts[j,-1])
    x = RefCDS[[j]]

    # a. Neutral model: wmis==1, wnon==1, wspl==1
    mrfold = sum(y[1:4])/sum(y[5:8]) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
    ll0 = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,1,1,1),dim=c(4,numrates))), log=T)) # loglik null model

    # b. Missense model: wmis==1, free wnon, free wspl
    mrfold = max(1e-10, sum(y[c(1,2)])/sum(y[c(5,6)])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
    wfree = y[3:4]/y[7:8]/mrfold; wfree[y[3:4]==0] = 0
    llmis = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,1,wfree),dim=c(4,numrates))), log=T)) # loglik free wmis

    # c. free wmis, wnon and wspl
    mrfold = max(1e-10, y[1]/y[5]) # Correction factor of "t"
    w = y[2:4]/y[6:8]/mrfold; w[y[2:4]==0] = 0 # MLE of dN/dS based on the local rate (using syn muts as neutral)
    llall = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,w),dim=c(4,numrates))), log=T)) # loglik free wmis, wnon, wspl
    w[w>1e4] = 1e4

    p = 1-pchisq(2*(llall-c(llmis,ll0)),df=c(1,3))
    return(c(w,p))
  }






plfit(tmpBucketCounts$Count)

tmpX <- (1-runif(100))^(-1/(2.5-1))
View(tmpX)
plfit(tmpX)



# power dist fit method explained
x=rpareto(1000,10,2.5) # 1000 entries, 10 up to 143, randomly generated, smaller the s (dispersion value), the greater the range
method="limit"
value=c()
finite=FALSE
nowarn=FALSE
nosmall=FALSE


#init method value to NULL
vec <- c() ; sampl <- c() ; limit <- c(); fixed <- c()

#  estimate xmin and alpha in the discrete case
#

if( is.null(vec) )
{
  # create a sequence of values from 1.5 to 3.6 in 0.1 increments
  vec<-seq(1.5,3.5,.01)  # covers range of most practical scaling parameters
}
print(vec)

zvec <- zeta(vec) # apply zeta function, not sure what this does yet
print(zvec)

xmins <- sort(unique(x))
print(xmins)
xmins <- xmins[-length(xmins)]
print(xmins)
length(xmins)

if( !is.null(limit) ){
  limit <- round(limit)
  xmins <- xmins[xmins>=limit]
}
print(limit)

if( !is.null(fixed) ){
  xmins <- fixed
}

if( !is.null(sampl) ){
  xmins <- xmins[unique(round(seq(1,length(xmins),length.out=sampl)))]
}

if( is.null(xmins) || length(xmins) < 2){
  stop("(plfit) error: x must contain at least two unique values.")
}

if(length(which(xmins==0) > 0)){
  stop("(plfit) error: x must not contain the value 0.")
}

xmax <- max(x)
print(xmax)
dat <- matrix(0,nrow=length(xmins),ncol=2)
View(dat)
z <- x

# for( xm in 1:length(xmins) )
# {
xm = 1
  xmin <- xmins[xm] # xmins is the set of values ordered smallest to largest, so xmin is the first and smallest value
  # print(xmin)
  z <- z[z>=xmin] # z is the set of values above the current min value
  n <- length(z)
  # print(n)

  # mine:
  repM = rep(t(1:(xmin-1)),length(vec)) # creates a repeated array of values 1 through current xmin-1, to a length of the input vector
  View(repM)

  # mine:
  kronM = t(kronecker(t(array(1,xmin-1)),vec)) # a matrix of 1-9 (as per 1 -> xmin-1) by the 0.1 sequence above (1.5 to 3.5)
  View(kronM[,50:60])

  # estimate alpha via direct maximization of likelihood function
  #  vectorized version of numerical calculation
  # matlab: zdiff = sum( repmat((1:xmin-1)',1,length(vec)).^-repmat(vec,xmin-1,1) ,1);
  if(xmin==1){
    zdiff <- rep(0,length(vec))
  }else{
    zdiff <- apply(rep(t(1:(xmin-1)),length(vec))^-t(kronecker(t(array(1,xmin-1)),vec)),2,sum)
    print(zdiff)
  }
  # matlab: L = -vec.*sum(log(z)) - n.*log(zvec - zdiff);
  L <- -vec*sum(log(z)) - n*log(zvec - zdiff);
  print(L)
  I <- which.max(L)
  print(I)
  # compute KS statistic
  fit <- cumsum((((xmin:xmax)^-vec[I])) / (zvec[I] - sum((1:(xmin-1))^-vec[I])))
  cdi <- cumsum(hist(z,c(min(z)-1,(xmin+.5):xmax,max(z)+1),plot=FALSE)$counts/n)
  dat[xm,] <- c(max(abs( fit - cdi )),vec[I])

# repeat for all values xm in xmin


D     <- min(dat[,1])
I     <- which.min(dat[,1])
xmin  <- xmins[I]
n     <- sum(x>=xmin)
alpha <- dat[I,2]

print(D)
print(alpha)

if( finite ){
  alpha <- alpha*(n-1)/n+1/n # finite-size correction
}
if( n<50 && !finite && !nowarn){
  print("(plfit) Warning : finite-size bias may be present")
}

#  end discrete case






medianIndex = round(snvSampleCount*0.5,0)
topXBucketSize = 0.05
topXBucketCount = 1/topXBucketSize

bsColNames = c("Bucket", "Min", "Max")
for(j in 1:topXBucketCount)
{
  bsColNames[3+j] = j*topXBucketSize
}

View(bsColNames)
colnames(snvBucketSummary) = bsColNames
snvBucketSummary = data.frame(matrix(ncol = 3+topXBucketCount, nrow = 0))

for(i in 1:snvBucketCount)
{
  bucketData = snvSampleCounts[i,]
  bucketData = data.frame(t(bucketData))
  colnames(bucketData) = c("BucketCount")
  bucketData = bucketData %>% arrange(BucketCount)

  snvBucketSummary[i,1] = snvBucketNames[i,1]
  snvBucketSummary[i,2] = min(bucketData$BucketCount)
  snvBucketSummary[i,3] = max(bucketData$BucketCount)

  for(j in 1:topXBucketCount)
  {
    topIndex = j*topXBucketSize*snvSampleCount
    snvBucketSummary[i,3+j] = bucketData[topIndex,1]
  }
}

View(snvBucketSummary)


# then test bucket-pairing CSS

snvMatrixDataTrans = t(snvMatrixData)
colnames(snvMatrixDataTrans) = snvBucketNames$Bucket
snvBucketCss = calc_cosine_sims(snvMatrixDataTrans, 0.8, "Bucket", T)
View(snvBucketCss)

View(snvMatrixData[,1:10])

snvMatrixDataWithBuckets = cbind(snvBucketNames$Bucket,snvMatrixData)
View(snvMatrixDataWithBuckets[,1:10])

names(snvMatrixDataWithBuckets)[names(snvMatrixDataWithBuckets) == 'snvBucketNames$Bucket'] <- 'Bucket'

snvBucketSS = snvMatrixDataWithBuckets %>% filter(grepl("C>G",Bucket))
View(snvBucketSS[,1:10])

snvBucketSampleCtoG = gather(snvBucketSS, "SampleId", "Count", -Bucket)
View(snvBucketSampleCtoG)

snvBucketSampleCtoGSpec = snvBucketSampleCtoG %>% filter(Bucket=="C>G_TCT"|Bucket=="C>G_TCA")
View(snvBucketSampleCtoGSpec)
snvTmp1 = snvBucketSampleCtoGSpec %>% spread(Bucket,Count)
View(snvTmp1)

snvTmp1$Ratio = ifelse(snvTmp1$`C>G_TCT`>0,snvTmp1$`C>G_TCA` / snvTmp1$`C>G_TCT`,0)
snvTmp1$RatioBucket = round(snvTmp1$Ratio/0.01,0)*0.01

snvTmpBucketRatioStats = snvTmp1 %>% group_by(RatioBucket) %>% summarise(Count=n()) %>% arrange(RatioBucket)
View(snvTmpBucketRatioStats)

sam1 = snvTmp1$`C>G_TCT`
sam2 = snvTmp1$`C>G_TCA`
tmp1CSS = cosine_sim(sam1, sam2)
print(tmp1CSS)


View(snvSampleCounts)

snvBucketCountSize = 100
snvSampleCounts$CountBucket = round(snvSampleCounts$Count/snvBucketCountSize)*snvBucketCountSize
snvSampleCounts$CountBucket_20 = round(snvSampleCounts$Count/20)*20
snvSampleCounts$CountBucket_5 = round(snvSampleCounts$Count/5)*5

View(snvSampleCounts %>%group_by(Bucket) %>% summarise(Count=sum(Count)))

View(snvSampleCounts %>% filter(Bucket=="C>G_TCT") %>% group_by(CountBucket) %>% summarise(Count=n()))
View(snvSampleCounts %>% filter(Bucket=="C>G_TCT"&CountBucket_20<=10000) %>% group_by(CountBucket_20) %>% summarise(Count=n()))

View(snvSampleCounts %>% filter(Bucket=="C>T_TCT") %>% group_by(CountBucket) %>% summarise(Count=n()))
View(snvSampleCounts %>% filter(Bucket=="T>G_GTC") %>% group_by(CountBucket_2) %>% summarise(Count=n()))

View(snvSampleCounts %>% filter(CountBucket_5<=10000) %>% group_by(CountBucket_5) %>% summarise(Count=n()))

View(snvSampleCounts %>% filter(Bucket=="C>T_ACG"|Bucket=="C>T_GCG"|Bucket=="C>T_CCG"|Bucket=="C>T_TCG")
     %>% group_by(CountBucket_5) %>% summarise(Count=n()))

snvTmpBCs = (snvSampleCounts %>% filter(Bucket=="C>T_ACG"|Bucket=="C>T_GCG"|Bucket=="C>T_CCG"|Bucket=="C>T_TCG")
             %>% group_by(CountBucket_5) %>% summarise(Count=n()))

View(snvTmpBCs)
View(snvBucketNames)

i = 0
for(bucketName in snvBucketNames$Bucket)
{
  bucketCountData = snvSampleCounts %>% filter(Bucket==bucketName)
  bucketCounts = bucketCountData$Count

  print(paste("testing bucket=", bucketName, ", count=", sum(bucketCounts), sep=''))
  pdFit = plfit(bucketCounts)
  print(paste("result: bucket=", bucketName, ", alpha=", pdFit$alpha, ", xmin=", pdFit$xmin, ", D=", pdFit$D, sep=''))
  i = i + 1

  if(i > 10)
    break
}




TPL=function(par,x,y)
{
  C=par[1]
  beta=par[2]
  xo=par[3]

  est=(C*x^(-beta))*exp(-x/xo)
  sum((log(y) - log(est))^2)
}

f = optim(par=c(0.1,0.9,700),TPL,x=snvTmpBCs$CountBucket_5, y=snvTmpBCs$Count)
print(f)

x1=c(seq(0.1,1,0.1),seq(1,10000,1))
y1=(f$par[1]*x1^(-f$par[2]))*exp(-x1/f$par[3])

plot(x1,y1,log="xy",type="l",lwd=3,col="#CC6666",xlab="",ylab="",las=1,xaxt="n",yaxt="n",xlim=c(min(x),max(x)),ylim=c(min(y),max(y)))

dev.off()

hist(snvTmpBCs$Count, breaks=5)

dev.off()

h=hist(as.matrix(snvTmpBCs), breaks=5, plot=FALSE)
x=h$mids
y=h$density

plot(x, y, log="xy", pch=16, type="p", cex=2,col="steelblue3",xlab="",ylab="",las=1,xaxt="n",yaxt="n")


brea=seq(5,log(10000),10000)
brea=c(0, exp(brea))
h=hist(as.matrix(snvTmpBCs), breaks=brea, plot=FALSE)
x=h$mids
y=h$density

par(mar=c(5,6,2,2))
plot(x,y,log="xy",pch=16,type="p",cex=2,col="steelblue3",xlab="",ylab="",las=1,xaxt="n",yaxt="n"))
mtext(quote(X),1,line=3,cex=2)
mtext(quote(P(X)),2,line=4,cex=2)
eaxis(1,at=c(1,10,10^2,10^3,10^4),cex.axis=1.5)
eaxis(2,at=c(10^-9,10^-8,10^-7,10^-6,10^-5,10^-4,10^-3,10^-2,10^-1,10^0),cex.axis=1.5)

View(snvTmpBCs)

plot = pwrdist(snvTmpBCs2)

# age-related: C>T_ACG, C>T_GCG, C>T_CCG, C>T_TCG





pwrdist<-function(u,...)
{
  # u is vector of event counts
  fx <- table(u)
  i <- as.numeric(names(fx))
  y <- rep(0,max(i))
  y[i] <- fx
  m0 <- glm(y~log(1:max(i)),family=quasipoisson())
  print(summary(m0))
  sub <- paste("s=",round(m0$coef[2],2),"lambda=",sum(u),"/",length(u))
  plot(i,fx,log="xy",xlab="x",sub=sub,ylab="counts",...)
  grid()
  lines(1:max(i),(fitted(m0)),type="b")
  return(m0)
}

# check for bucket pairings excluding the top 10% of samples
View(snvSampleTotals)
snvBottomXPerc = tail(snvSampleTotals, nrow(snvSampleTotals)*0.95)
View(snvBottomXPerc)

snvLowestXPCounts = snvSampleCounts %>% filter(SampleId %in% snvBottomXPerc$SampleId)
snvLowestXPMatrixData = snvLowestXPCounts %>% spread(SampleId,Count)
snvLowestXPMatrixData[is.na(snvLowestXPMatrixData)] = 0
View(snvLowestXPMatrixData[,1:10]) # check a subset
ncol(snvLowestXPMatrixData)
snvLowestXPMatrixData = within(snvLowestXPMatrixData, rm(Bucket))
snvLowestXPMdt = t(snvLowestXPMatrixData)
View(snvLowestXPMdt[,1:10])
ncol(snvLowestXPMdt)

colnames(snvLowestXPMdt) = snvBucketNames$Bucket
snvLowestXPCss = calc_cosine_sims(snvLowestXPMdt, 0.8, "Bucket", T)
View(snvLowestXPCss)

# top10% of samples - should be very low due to lack of correlation and fewer samples
snvHighestXPerc = head(snvSampleTotals, nrow(snvSampleTotals)*0.10)
View(snvHighestXPerc)

snvHighestXPCounts = snvSampleCounts %>% filter(SampleId %in% snvHighestXPerc$SampleId)
snvHighestXPMatrixData = snvHighestXPCounts %>% spread(SampleId,Count)
snvHighestXPMatrixData[is.na(snvHighestXPMatrixData)] = 0
View(snvHighestXPMatrixData[,1:10]) # check a subset
ncol(snvHighestXPMatrixData)
snvHighestXPMatrixData = within(snvHighestXPMatrixData, rm(Bucket))
snvHighestXPMdt = t(snvHighestXPMatrixData)
View(snvHighestXPMdt[,1:10])
ncol(snvHighestXPMdt)

colnames(snvHighestXPMdt) = snvBucketNames$Bucket
snvHighestXPCss = calc_cosine_sims(snvHighestXPMdt, 0.8, "Bucket", T)
View(snvHighestXPCss)




# View(cssBucketResultGroups)
View(cssBucketResults)


# all samples, all CSS values - note this will be 2400^2/2 entries
snvAllCssResults = calc_cosine_sims(snvMatrixData, 0.5, "Sample", F)
View(snvAllCssResults)

snvHighMLCssResults = calc_cosine_sims(snvHighMatrixData, 0.5, "Sample", F)
View(snvHighMLCssResults)

snvHighMLCssResults = merge(snvHighMLCssResults, sampleCancerTypes, by.x="Sample1", by.y="SampleId", ,all.x=TRUE)
snvHighMLCssResults = merge(snvHighMLCssResults, sampleCancerTypes, by.x="Sample2", by.y="SampleId", ,all.x=TRUE)
snvHighMLCssResults = snvHighMLCssResults %>% arrange(-CSS)
colnames(snvHighMLCssResults) <- c("Sample1", "Sample2", "CSS", "CancerType1", "CancerType2")
write.csv(snvHighMLCssResults, "~/logs/r_output/snvHighMLCssResults.csv")

numSamples = ncol(snvHighMatrixData)
snvHighMLSampleNames = colnames(snvHighMatrixData)
View(snvHighMLSampleNames)

View(snvBucketNames)

cssRCount = 100
cssRInv = 1/cssRCount
cssResultGroups = data.frame(matrix(ncol = 2, nrow = cssRCount))
colnames(cssResultGroups) <- c("CSSBand", "Count")

for(i in 1:nrow(cssResultGroups))
{
  cssResultGroups[i,1] = i*cssRInv
  # cssResultGroups[i,2] = 0
}



View(cssResultGroups)


print(cosine_sim(snvHighMatrixData$CPCT02060041T, snvHighMatrixData$DRUP01230001T))

View(cosineSimResults)

View(snvHighMatrixData %>% select(CPCT02060041T,DRUP01230001T))
View(highMLSampleCancerTypes %>% filter(SampleId=="CPCT02060041T"|SampleId=="DRUP01230001T"))
View(sampleSigData %>% filter(SampleId=="CPCT02060041T"|SampleId=="DRUP01230001T"))

View(sampleSigData)

View(snvHighMatrixData)

# same again but for buckets

cssBucketResults = data.frame(matrix(ncol = 3, nrow = 0))
colnames(cssBucketResults) <- c("Bucket1", "Bucket2", "CSS")

cssRCount = 100
cssRInv = 1/cssRCount
cssBucketResultGroups = data.frame(matrix(ncol = 2, nrow = cssRCount))
colnames(cssBucketResultGroups) <- c("CSSBand", "Count")

for(i in 1:nrow(cssBucketResultGroups))
{
  cssBucketResultGroups[i,1] = i*cssRInv
  cssBucketResultGroups[i,2] = 0
}

numBuckets = nrow(snvBucketNames)
snvHighMatrixDataTrans = t(snvHighMatrixData)
# ncol(snvHighMatrixDataTrans)
# nrow(snvHighMatrixDataTrans)
# View(snvHighMatrixDataTrans[,1:10])
rownames(snvHighMatrixDataTrans) <- NULL

View(snvBucketNames)
snvBucketNamesVec = snvBucketNames$Bucket

for(i in 1:numBuckets)
{
  data1 = snvHighMatrixDataTrans[,i]

  for(j in i+1:numBuckets)
  {
    if(j>numBuckets)
      break

    data2 = snvHighMatrixDataTrans[,j]

    css = cosine_sim(data1,data2)

    # print(paste("s1=", snvHighMLSampleNames[i], ", count=", sum(sam1Data), ", s2=", snvHighMLSampleNames[j], ", count=", sum(sam2Data), sep=''))
    # print(paste(i, j, css, sep=', '))

    if(css >= 0.8)
    {
      print(paste("b1=", snvBucketNamesVec[i], ", count=", sum(data1), ", b2=", snvBucketNamesVec[j], ", count=", sum(data2), sep=''))
      print(paste(i, j, css, sep=', '))
      rowIndex = nrow(cssBucketResults)+1
      cssBucketResults[rowIndex,1] = snvBucketNamesVec[i]
      cssBucketResults[rowIndex,2] = snvBucketNamesVec[j]
      cssBucketResults[rowIndex,3] = round(css,6)
    }

    cssRounded = round(css/cssRInv)
    cssBucketResultGroups[cssRounded,2] = cssBucketResultGroups[cssRounded,2]+1
  }
}

View(cssBucketResultGroups)
View(cssBucketResults)


