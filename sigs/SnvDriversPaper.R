library(dplyr)
library(tidyr)
library(stringi)
library("pracma") # for ceil function

# plotting
library(grid)
library(gridExtra)
library(ggplot2)
library(cowplot)
# library(scater)

# Load sample counts in matrix form
snvMatrixData = read.csv('~/data/sigs/snv_prod_matrix_data_20190409.csv')
nrow(snvMatrixData)
ncol(snvMatrixData)
View(snvMatrixData[,1:10])
View(snvMatrixData[,1146:1147])

prodSamples = colnames(snvMatrixData)
View(prodSamples)

# re-create SNV buckets
snvBuckets = create_empty_signature()
snvBuckets = snvBuckets %>% mutate(Bucket=paste(type,context,sep='_'))
snvBuckets = snvBuckets %>% select(Bucket)
View(snvBuckets)

# create sample counts from matrix data
snvSampleCounts = matrix_to_sample_counts(snvMatrixData,snvBuckets)
View(snvSampleCounts)
View(snvSampleCounts %>% group_by(SampleId) %>% count()) # 4311 samples (at 09/04/2019), without QC

# limit to relevant samples in HPC
load('~/data/r_data/highestPurityCohortSummary.RData')
nrow(highestPurityCohortSummary)
View(highestPurityCohortSummary)

# take cancer types from DB or from HPC
sampleCancerTypes = read.csv('~/data/sigs/sample_cancer_types.csv')
sampleCancerTypes[is.na(sampleCancerTypes)] <- 'Unknown'

sampleCancerTypes = highestPurityCohortSummary %>% select(SampleId=sampleId,CancerType=cancerType)
View(sampleCancerTypes)
View(sampleCancerTypes %>% group_by(CancerType) %>% count())

snvDpSampleCounts = snvSampleCounts %>% filter(SampleId %in% highestPurityCohortSummary$sampleId)
View(snvDpSampleCounts %>% group_by(SampleId) %>% count()) 
snvDpSampleTotals = snvDpSampleCounts %>% group_by(SampleId) %>% count()
nrow(snvDpSampleTotals)

# no missing DP samples 
missingDpSamples = highestPurityCohortSummary %>% filter(!(sampleId %in% snvDpSampleTotals$SampleId))
View(missingDpSamples)
View(missingDpSamples$sampleId)

# write out DP sample matrix data
snvDpMatrixData = sample_counts_to_matrix(snvDpSampleCounts)
View(snvDpMatrixData[,1:10])
ncol(snvDpMatrixData)

write.csv(snvDpMatrixData, '~/data/sigs/snv_dp_matrix_data.csv', row.names = F, quote = F)
snvDpMatrixData = as.matrix(read.csv('~/data/sigs/snv_dp_matrix_data.csv'))

snvDpSampleCounts = matrix_to_sample_counts(snvDpMatrixData,snvBuckets)


# load in results from SigAnalyser NMF fit with PCAWG sigs
snvDpPcawgContribs = as.matrix(read.csv('~/data/sigs/logs/snv_pcawg_fit_nmf_contribs.csv', stringsAsFactors = F))
snvDpPcawgSigs = as.matrix(read.csv('~/data/sigs/logs/snv_pcawg_fit_nmf_sigs.csv', stringsAsFactors = F))

snvPcawgSigsNamesStr = get_pcawg_named_sigs()
length(snvPcawgSigsNamesStr)

View(snvDpPcawgSigs)
View(snvDpPcawgContribs[,1:20])


####
# Standard Signature Report
evaluate_nmf_data("SNV", "dp_pcawg_fit", snvDpPcawgSigs, snvDpPcawgContribs, snvDpMatrixData, snvDpSampleCounts,
                  sampleCancerTypes, snvBuckets, snvPcawgSigsNamesStr, F, F, 0, F)


# load in results from SigAnalyser NMF fit with COSMIC sigs
snvDpCosmicContribs = as.matrix(read.csv('~/data/sigs/logs/snv_cosmic_fit_nmf_contribs.csv', stringsAsFactors = F))
snvDpCosmicSigs = as.matrix(read.csv('~/data/sigs/logs/snv_cosmic_fit_nmf_sigs.csv', stringsAsFactors = F))

snvCosmicSigsNamesStr = get_signame_list(30,T)
length(snvCosmicSigsNamesStr)
View(snvCosmicSigsNamesStr)

View(snvDpCosmicSigs)
View(snvDpCosmicContribs[,1:20])

View(snvDpMatrixData[,1:10])

evaluate_nmf_data("SNV", "dp_cosmic_fit", cosmicSigs, snvDpCosmicContribs, snvDpMatrixData, snvDpSampleCounts,
                  sampleCancerTypes, snvBuckets, snvCosmicSigsNamesStr, T, F, 0, F)



evaluate_nmf_data("SNV", "dp_cosmic_nmf_fit", cosmicSigs, nmfFitContributions, snvDpMatrixData, snvDpSampleCounts,
                  sampleCancerTypes, snvBuckets, snvCosmicSigsNamesStr, T, F, 0, F)

#View(sampleCancerTypes)
#View(sampleCancerTypes %>% group_by(CancerType) %>% count())


####
# Signature prevalence by Cancer Type
snvDpSampleSigData = get_sig_summary(snvDpCosmicSigs, snvDpCosmicContribs, snvDpMatrixData, snvCosmicSigsNamesStr, snvBuckets)
snvDpSampleSigData = get_sig_summary(snvDpCosmicSigs, nmfFitContributions, snvDpMatrixData, snvCosmicSigsNamesStr, snvBuckets)
# snvDpSampleSigData = get_sig_summary(snvDpPcawgSigs, snvDpPcawgContribs, snvDpMatrixData, snvPcawgSigsNamesStr, snvBuckets)
View(snvDpSampleSigData)
# View(snvDpSampleSigData %>% filter(is.na(CancerType)))
# snvDpSampleSigData[is.na(snvDpSampleSigData)] <- "Unknown"

write.csv(snvDpSampleSigData, '~/data/sigs/dp/snvDpSampleSigData.csv', row.names = F, quote = F)
write.csv(nmfFitContributions, '~/data/sigs/dp/snvDpCosmicFitContributions.csv', row.names = F, quote = F)
write.csv(snvDpMatrixData, '~/data/sigs/dp/snvDpMatrixData.csv', row.names = F, quote = F)
write.csv(snvDpSampleCounts, '~/data/sigs/dp/snvDpSampleCounts.csv', row.names = F, quote = F)
write.csv(cosmicSigs, '~/data/sigs/dp/cosmicSigs.csv', row.names = F, quote = F)

# filter for samples with a sig contribution above 5%
sigPercentThreshold = 0.05
sigCountThreshold = 300
snvDpSampleSigData = snvDpSampleSigData %>% filter(SigPercent>=sigPercentThreshold&Count>=sigCountThreshold)

sampleTotals = snvDpSampleSigData %>% group_by(SampleId) %>% summarise(SampleTotal=round(sum(Count),0))
# View(sampleTotals)
snvDpSampleSigData = merge(snvDpSampleSigData,sampleTotals,by.x="SampleId",by.y="SampleId",all.x=T)

cancerSampleCounts = sampleCancerTypes %>% group_by(CancerType) %>% summarise(CancerSampleCount=n())
rowIndex = data.frame(as.numeric(as.character(rownames(cancerSampleCounts))))
colnames(rowIndex) <- c("CancerIndex")
cancerSampleCounts = cbind(rowIndex, cancerSampleCounts)
# View(cancerSampleCounts)

snvDpSampleSigData = merge(snvDpSampleSigData, sampleCancerTypes,by.x="SampleId",by.y="SampleId",all.x=T)

# View(snvDpSampleSigData)

cancerSampleSigSummary = (snvDpSampleSigData %>% filter(SigName!='Excess'&SigName!='Unalloc') %>% group_by(CancerType,SigName) 
                          %>% summarise(CancerSigSampleCount=sum(SigPercent>=sigPercentThreshold),
                                        CancerSigSampleMutLoad=round(sum(Count),0)))

cancerSampleSigSummary = merge(cancerSampleSigSummary,cancerSampleCounts,by='CancerType',all.x=T)
# View(cancerSampleSigSummary)


cancerTypesBySig = cancerSampleSigSummary %>% group_by(SigName) %>% summarise(SigCancerTypeCount=sum(CancerSigSampleCount>0))
rowIndex = data.frame(as.numeric(as.character(rownames(cancerTypesBySig))))
colnames(rowIndex) <- c("SigIndex")
cancerTypesBySig = cbind(rowIndex, cancerTypesBySig)
# View(cancerTypesBySig)

cancerSampleSigSummary = merge(cancerSampleSigSummary,cancerTypesBySig,by='SigName',all.x=T)
# View(cancerSampleSigSummary)

# work out median mutational load per sig and cancer type for samples with the signature
cancerSigMedianCounts = snvDpSampleSigData %>% filter(Count>0) %>% group_by(SigName,CancerType) %>% summarise(MedianSampleCount=round(median(Count),0))
# View(cancerSigMedianCounts)

cancerSampleSigSummary = merge(cancerSampleSigSummary,cancerSigMedianCounts,by=c('SigName','CancerType'),all.x=T)
# View(cancerSampleSigSummary)

mutLoadMin = 0.25
mutLoadMax = 16
cancerSampleSigSummary = cancerSampleSigSummary %>% mutate(SamplesWithSig=round(CancerSigSampleCount/CancerSampleCount,3),
                                                           MedianLoadPerMb=ifelse(MedianSampleCount>0,round(MedianSampleCount/3e3,3),0),
                                                           MedianLoadPerMb=ifelse(MedianLoadPerMb>mutLoadMax,mutLoadMax,
                                                                                  ifelse(MedianLoadPerMb<mutLoadMin,mutLoadMin,MedianLoadPerMb)))

View(cancerSampleSigSummary)
# Matching PCAWG's plot:
# Size of circle: Proportion of samples in this cancer type with this signature - use SamplesWithSig
# Colour of circle: Median mutations/Mb due to signature (among samples with the signature)  - use MedianLoadPerMb

# scale_fill_gradient(name = "Significance", trans = "log10", low = "#2171b5", high = "#bdd7e7", guide = "colourbar",limits=c(1e-14, 1e-2), breaks=c(1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2), labels=c(1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)) +
  
plotBreaks = c(0.25, 0.5, 1, 2, 4, 8, 16)

plot = (ggplot(cancerSampleSigSummary %>% filter(SamplesWithSig>0), aes(x=CancerType, y=reorder(SigName, -SigIndex))) 
          + geom_point(aes(color=MedianLoadPerMb, size=SamplesWithSig)) 
          + scale_color_gradient(name='MedianLoadPerMb', trans='log2', low="lightblue", high="midnightblue", guide="colourbar",limits=c(mutLoadMin,mutLoadMax), breaks=plotBreaks, labels=plotBreaks)
          + theme_gray()
          + theme(panel.grid.major = element_line(), panel.grid.minor = element_line())
          + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0,size=10))
          + scale_x_discrete(position = 'top')
          # + guides(color = guide_colourbar(barheight = 20, direction = "vertical", title.position="top", title.hjust = 0, title.vjust = 0.5, nbin = 50))
          + labs(subtitle="Signature Allocation by Cancer Types", y="Signature", x="CancerType"))

print(plot)



cancerSampleSigSummary = (cancerSampleSigSummary %>% mutate(MedianLoadBucket=round(0.5**round(log(MedianLoadPerMb,0.5)),2),
                                                            MedianLoadGroup=ifelse(MedianLoadBucket<0.075,0.05,
                                                                                   ifelse(MedianLoadBucket<0.017,0.1,
                                                                                   ifelse(MedianLoadBucket<0.37,0.25,
                                                                                   ifelse(MedianLoadBucket<0.75,0.5,
                                                                                   ifelse(MedianLoadBucket<1.75,1,
                                                                                   ifelse(MedianLoadBucket<3.75,2.5,
                                                                                   ifelse(MedianLoadBucket<7.5,5,
                                                                                   ifelse(MedianLoadBucket<17.5,10,25))))))))))
                                                            
cancerSampleSigSummary$MedianLoadGroup = as.factor(cancerSampleSigSummary$MedianLoadGroup)
View(cancerSampleSigSummary)
# View(cancerSampleSigSummary %>% group_by(SigName,CancerType,MedianLoadBucket) %>% summarise(SamplesWithSig=first(SamplesWithSig)) %>% spread(MedianLoadBucket,SamplesWithSig))

#mutLoadColours = c("cornsilk3", "hotpink", "darkorange", "seagreen3", "tomato3", "thistle2", "steelblue2", "darkgreen", "indianred", "honeydew2",
#              "turquoise3", "lightpink2", "goldenrod2", "darkslateblue", "yellowgreen", "wheat2", "violetred2", "ivory3", "coral1", "springgreen2")

mutLoadColours = c('wheat2','orange','darkorange','indianred','violet','purple','blue','darkblue')

plot = (ggplot(cancerSampleSigSummary %>% filter(SamplesWithSig>0), aes(x=CancerType, y=reorder(SigName, -SigIndex))) 
        + geom_point(aes(col=MedianLoadGroup, size=SamplesWithSig)) 
        + scale_colour_manual(values = mutLoadColours)
        # + geom_smooth(method="loess", se=F) 
        + theme_gray()
        + theme(panel.grid.major = element_line(), panel.grid.minor = element_line())
        + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0,size=10))
        + scale_x_discrete(position = 'top')
        +  labs(subtitle="Signature Allocation by Cancer Types", y="Signature", x="CancerType"))

print(plot)


####
# Residuals Report


residualsSummary = (snvDpSampleSigData %>% group_by(SampleId,CancerType) 
                    %>% summarise(SampleTotal=first(SampleTotal), 
                                  ResidualTotal=round(sum(ifelse(SigName=='Excess',-Count,ifelse(SigName=='Unalloc',Count,0))),0))
                    %>% mutate(ResidualPerc=round(ResidualTotal/SampleTotal,4)))

View(residualsSummary)

resPlot = (ggplot(residualsSummary, aes(CancerType, ResidualPerc))
           + theme(axis.text.x = element_text(angle = 90, hjust=1,size=7))
           + geom_boxplot(varwidth=T, fill="plum") +labs(title="Residuals % by Cancer Type", x="Cancer Type", y="Residuals %"))

print(resPlot)

# residuals percentage vs mutational load
tmp = snvDpSampleSigData %>% filter(SigName=='Excess'|SigName=='Unalloc') %>% group_by(SampleId) %>% summarise(ResidualsTotal=sum(abs(Count)))
tmp = merge(tmp,sampleTotals,by.x="SampleId",by.y="SampleId",all.x=T)
tmp$ResidualsTotal = round(tmp$ResidualsTotal/tmp$SampleTotal,3)

View(tmp)

print(ggplot(data=tmp, aes(x=SampleTotal, y=ResidualsTotal))
        + geom_point()
        + scale_x_log10())
        # + xlim(0,100))



# plot buckets for samples with residuals above 50%

worstSamples = c('CPCT02010772T','CPCT02010662T')
worstSamples = c('CPCT02010772T','CPCT02070053T','CPCT02070038T','CPCT02340067T','CPCT02010662T','CPCT02030415T')
worstSamplesData = snvDpSampleCounts %>% filter(SampleId %in% worstSamples)
worstSamplesMatixData = sample_counts_to_matrix(worstSamplesData)

plot_96_profile(worstSamplesMatixData)    
plot_96_profile(snvDpMatrixData[,c(1,7)])

View(snvDpMatrixData[,'CPCT02010772T'])
colnames()

plot_96_profile(snvDpMatrixData[,'CPCT02010772T'])

####
## Worst samples NMF sig discovery and report
worstSamples = residualsSummary %>% filter(ResidualPerc>0.3&SampleTotal>10000) %>% arrange(-ResidualTotal)
worst20Samples = head(worstSamples,20)
View(worst20Samples)

# which sigs are they using, however poorly?
View(snvDpSampleSigData %>% filter(SampleId %in% worst20Samples$SampleId & SigPercent>=0.05))

# obtain a limited data set of the worst samples for a de-novo NMF signature discovery
wsSampleCounts = snvDpSampleCounts %>% filter(SampleId %in% worst20Samples$SampleId)
wsSampleCancerTypes = wsSampleCounts %>% group_by(SampleId) %>% count()
wsSampleCancerTypes$CancerType = 'Unknown'
View(wsSampleCancerTypes)
write.csv(wsSampleCancerTypes %>% select(SampleId,CancerType), '~/data/sigs/ws_sample_ext_data.csv', row.names = F, quote = F)

wsSampleCounts = snvDpSampleCounts %>% filter(SampleId=='CPCT02010772T')
wsMatrixData = sample_counts_to_matrix(wsSampleCounts)
View(wsMatrixData)

write.csv(sampleCancerTypes, '~/data/sigs/dp_sample_cancer_types.csv', row.names = F, quote = F)
write.csv(wsMatrixData, '~/data/sigs/snv_ws_matrix_data.csv', row.names = F, quote = F)
write.csv(wsSampleCounts, '~/data/sigs/snv_ws_sample_counts.csv', row.names = F, quote = F)

baBackgroundAllocs = read.csv('~/dev/nmf/snv_ba_sam_calc_data_M12K.csv')
View(baBackgroundAllocs)

wsBackgroundAllocs = baBackgroundAllocs %>% filter(SampleName %in% worst20Samples$SampleId)
write.csv(wsBackgroundAllocs, '~/data/sigs/ws_ba_background_counts.csv', row.names = F, quote = F)


wsPcawgSigs = as.matrix(read.csv('~/data/sigs/logs/ws_pcawg_sigs.csv', stringsAsFactors = F))
wsDenovoSigs = as.matrix(read.csv('~/data/sigs/logs/snv_ws_ba_sigs.csv', stringsAsFactors = F))
colnames(wsDenovoSigs) = c('NEW01','NEW02','NEW03','NEW04','NEW05','NEW06','NEW07','NEW08','NEW09')

# remove small new sig
wsDenovoSigs = within(as.data.frame(wsDenovoSigs), rm('NEW09'))

View(wsDenovoSigs)
View(wsPcawgSigs)
wsCombinedSigs = cbind(wsPcawgSigs,wsDenovoSigs)
ncol(wsCombinedSigs)
write.csv(wsCombinedSigs, '~/data/sigs/logs/ws_combined_sigs.csv', row.names = F, quote = F)
wsCombinedSigs = as.matrix(wsCombinedSigs)


# load in results from SigAnalyser NMF fit with PCAWG + de-novo sigs
wsPcawgDenovoContribs = as.matrix(read.csv('~/data/sigs/logs/snv_ws_pcawg_denovo_fit_nmf_contribs.csv', stringsAsFactors = F))
wsSigsNamesStr = colnames(wsCombinedSigs)
print(wsSigsNamesStr)

evaluate_nmf_data("SNV", "dp_ws_pcawg_denovo_fit", wsCombinedSigs, wsPcawgDenovoContribs, wsMatrixData, wsSampleCounts,
                  sampleCancerTypes, snvBuckets, wsSigsNamesStr, F, T, 0, F)




worstSampleSigData = get_sig_summary(wsCombinedSigs, wsPcawgDenovoContribs, wsMatrixData, wsSigsNamesStr, snvBuckets)
worstSampleSigData$CancerType = "Unknown"

wsResSummary = (worstSampleSigData %>% group_by(SampleId,CancerType) 
                   %>% summarise(SampleTotal=sum(Count), 
                                 ResidualTotal=round(sum(ifelse(SigName=='Excess',-Count,ifelse(SigName=='Unalloc',Count,0))),0))
                   %>% mutate(ResidualPerc=round(ResidualTotal/SampleTotal,4)))

View(wsResSummary)

wsResPlot = (ggplot(wsResSummary %>% filter(SampleTotal>=1000), aes(CancerType, ResidualPerc))
                + theme(axis.text.x = element_text(angle = 90, hjust=1,size=7))
                + geom_boxplot(varwidth=T, fill="plum") +labs(title="Residuals % by Cancer Type", x="Cancer Type", y="Residuals %"))

print(wsResPlot)


# dev.off()

print(ggplot(data = residualsSummary, aes(ResidualTotal))
      + stat_ecdf()
      + scale_x_log10())

print(ggplot(data = residualsSummary, aes(ResidualPerc))
      + stat_ecdf()
      + facet_wrap(~CancerType))

resPlot = (ggplot(residualsSummary, aes(CancerType, ResidualPerc))
           + theme(axis.text.x = element_text(angle = 90, hjust=1,size=7))
           + geom_boxplot(varwidth=T, fill="plum") +labs(title="Residuals % by Cancer Type", x="Cancer Type", y="Residuals %"))

print(resPlot)



## Mutational Load by Cancer Type
snvDpSampleSummary = merge(snvDpSampleCounts, sampleCancerTypes,by="SampleId",all.x=T) %>% 
  group_by(SampleId,CancerType) %>% summarise(SampleTotal=sum(Count))

View(snvDpSampleSummary)

View(snvDpSampleSummary %>% group_by(CancerType) %>% summarise(Samples=n(),
                                                               CancerTotal=sum(SampleTotal),
                                                               Min=min(SampleTotal),
                                                               Max=max(SampleTotal),
                                                               Avg=round(sum(SampleTotal)/n(),0),
                                                               Med=round(median(SampleTotal),0),
                                                               StdDev=round(std(SampleTotal),1)))


print(ggplot(data = snvDpSampleSummary, aes(x=))
      + stat_ecdf()
      + facet_wrap(~CancerType))


View(sampleCancerTypes)


## DEBUG


# only skin and lung
View(sampleCancerTypes %>% group_by(CancerType) %>% count())

skinAndLungSamples = sampleCancerTypes %>% filter(CancerType=='Skin'|CancerType=='Lung')
nrow(skinAndLungSamples)
View(skinAndLungSamples %>% group_by(CancerType) %>% count())
snvSkinLungSampleCounts = snvDpSampleCounts %>% filter(SampleId %in% skinAndLungSamples$SampleId)
nrow(snvSkinLungSampleCounts %>% group_by(SampleId) %>% count())
rm(snvSkinLungSampleCounts)

# write out DP sample matrix data
snvSkinLungMatrixData = sample_counts_to_matrix(snvSkinLungSampleCounts)
View(snvSkinLungMatrixData[,1:10])
ncol(snvSkinLungMatrixData)
rm(snvSkinLungMatrixData)
write.csv(snvSkinLungMatrixData, '~/data/sigs/snv_dp_skin_lung_matrix_data.csv', row.names = F, quote = F)


snvSkinLungPcawgContribs = as.matrix(read.csv('~/data/sigs/logs/snv_pcawg_skin_lung_fit_nmf_contribs.csv', stringsAsFactors = F))
snvSkinLungPcawgSigs = as.matrix(read.csv('~/data/sigs/logs/snv_pcawg_skin_lung_fit_nmf_sigs.csv', stringsAsFactors = F))
rm(snvSkinLungPcawgContribs)
rm(snvSkinLungPcawgSigs)

View(snvSkinLungPcawgContribs[,1:20])
ncol(snvSkinLungPcawgContribs)

# snvSkinLungSampleCounts = matrix_to_sample_counts(snvSkinLungMatrixData, snvBuckets)

evaluate_nmf_data("SNV", "skin_lung_pcawg_fit", snvSkinLungPcawgSigs, snvSkinLungPcawgContribs, snvSkinLungMatrixData, snvSkinLungSampleCounts,
                  sampleCancerTypes, snvBuckets, snvPcawgSigsNamesStr, T, F)


specificSample = snvDpSampleCounts %>% filter(SampleId=='CPCT02010772T')
View(specificSample)
specificSampleMatrixData = sample_counts_to_matrix(specificSample)
View(specificSampleMatrixData)
write.csv(specificSampleMatrixData, '~/data/sigs/snv_dp_syd895_matrix_data.csv', row.names = F, quote = F)






  
# calculate and factor in residuals
residuals = calc_contrib_sample_residuals(contribution, signatures, matrixData, bucketNames)
View(residuals)

sampleSigData = append_residuals(contribution, signatures, matrixData, bucketNames, sampleSigData)
View(sampleSigData)


# get cancer type and SV Count
sampleSigData = merge(sampleSigData, sampleCancerTypes,by.x="SampleId",by.y="SampleId",all.x=T)
# sampleSigData$CancerType = ifelse(is.na(sampleSigData$CancerType), 'N/A', paste(sampleSigData$CancerType, sep=""))
sampleSigData = merge(sampleSigData, origSampleCounts, by.x="SampleId",by.y="SampleId",all.x=T)
sampleSigDataNoResiduals = sampleSigData %>% filter(SigName!="Unalloc"&SigName!="Excess")

snvDpSampleSigData = sampleSigData
View(snvDpSampleSigData)


cancerSampleSigSummary = (cancerSampleSigSummary %>% mutate(MedianLoadBucket=round(0.5**round(log(MedianLoadPerMb,0.5)),2),
                                                            MedianLoadGroup=ifelse(MedianLoadBucket<0.075,'#1 0.05',
                                                                                   ifelse(MedianLoadBucket<0.017,'#2 0.1',
                                                                                          ifelse(MedianLoadBucket<0.37,'#3 0.25',
                                                                                                 ifelse(MedianLoadBucket<0.75,'#4 0.5',
                                                                                                        ifelse(MedianLoadBucket<1.75,'#5 1',
                                                                                                               ifelse(MedianLoadBucket<3.75,'#6 2.5',
                                                                                                                      ifelse(MedianLoadBucket<7.5,'#7 5',
                                                                                                                             ifelse(MedianLoadBucket<17.5,'#8 10','#9 25'))))))))))



# DEBUG plotting colours

evaluate_nmf_data("SNV", "dp_cosmic_fit", cosmicSigs, snvDpCosmicContribs, snvDpMatrixData, snvDpSampleCounts,
                  sampleCancerTypes, snvBuckets, snvCosmicSigsNamesStr, T, F, 0, F)

sigCount = 30
sigColours = get_sig_colours(sigCount)
sigColours[sigCount+1] = "black"
sigColours[sigCount+2] = "grey30"
print(sigColours)

signatures = cosmicSigs
sigNamesUnamed = get_signame_list(sigCount, F)
colnames(signatures) = sigNamesUnamed

sampleSigData = get_sig_data(signatures, snvDpCosmicContribs, snvCosmicSigsNamesStr)
sampleSigData = append_residuals(snvDpCosmicContribs, signatures, snvDpMatrixData, snvBuckets, sampleSigData)

# get cancer type and SV Count
sampleSigData = merge(sampleSigData, sampleCancerTypes,by.x="SampleId",by.y="SampleId",all.x=T)
sampleSigData = merge(sampleSigData, snvDpSampleCounts %>% group_by(SampleId) %>% summarise(SampleCount=sum(Count)), by.x="SampleId",by.y="SampleId",all.x=T)

View(sampleSigData)

sigNamesCombined = as.data.frame(cbind(sigNamesUnamed, snvCosmicSigsNamesStr))
colnames(sigNamesCombined) <- c("Signature", "SigName")
View(sigNamesCombined)

sigNamesCombined = sigNamesCombined %>% add_row(Signature = sprintf("%d",sigCount+1), SigName = "Unalloc")
sigNamesCombined = sigNamesCombined %>% add_row(Signature = sprintf("%d",sigCount+2), SigName = "Excess")
sigNamesCombined$Signature = as.numeric(as.character(sigNamesCombined$Signature))
View(sigNamesCombined)


sampleSigData = merge(sampleSigData,sigNamesCombined, by="SigName",all.x=T)

plot_sig_samples(sampleSigData, "", sigColours, runType, 0, dropSigLegend) # all samples














