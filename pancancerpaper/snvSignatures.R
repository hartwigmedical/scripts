library(devtools)
library(purple); # for multiplot
library(data.table)
library(dplyr)
library(tidyr)
library(stringi)
library(MutationalPatterns)
library(ggplot2)
library(cowplot)
theme_set(theme_bw() + theme(axis.text=element_text(size=5),axis.title=element_text(size=7),legend.text = element_text(size=5)))

snvDpSampleSigData = read.csv('/Users/jon/hmf/RData/SIgnatures/snvDpSampleSigData.csv')
snvDpSampleCounts = read.csv('/Users/jon/hmf/RData/SIgnatures/snvDpSampleCounts.csv')
nmfFitContributions = read.csv('/Users/jon/hmf/RData/SIgnatures/snvDpCosmicFitContributions.csv')
snvDpMatrixData = read.csv('/Users/jon/hmf/RData/SIgnatures/snvDpMatrixData.csv')
cosmicSigs = read.csv('/Users/jon/hmf/RData/SIgnatures/cosmicSigs.csv')

load('~/hmf/RData/processed/highestPurityCohortSummary.RData')
sampleCancerTypes = highestPurityCohortSummary %>% select(SampleId=sampleId,CancerType=cancerType)

snvDpSampleSigData = merge(snvDpSampleSigData, sampleCancerTypes,by.x="SampleId",by.y="SampleId",all.x=T)

# merge in sample totals and cancer type
sampleTotals = snvDpSampleCounts %>% group_by(SampleId) %>% summarise(SampleTotal=round(sum(Count),0))
snvDpSampleSigData = merge(snvDpSampleSigData,sampleTotals,by.x="SampleId",by.y="SampleId",all.x=T)

########################################
# 1. Signature Prevalence by Cancer Type

# filter for samples with a sig contribution above 5% and 300 count
sigPercentThreshold = 0.05
sigCountThreshold = 500
highMLSampleSigData = snvDpSampleSigData %>% filter(SigPercent>=sigPercentThreshold&Count>=sigCountThreshold)

cancerSampleCounts = sampleCancerTypes %>% group_by(CancerType) %>% summarise(CancerSampleCount=n())
rowIndex = data.frame(as.numeric(as.character(rownames(cancerSampleCounts))))
colnames(rowIndex) <- c("CancerIndex")
cancerSampleCounts = cbind(rowIndex, cancerSampleCounts)


# prepare summary for plot

# filter out residuals categories
cancerSampleSigSummary = (highMLSampleSigData %>% filter(SigName!='Excess'&SigName!='Unalloc') %>% group_by(CancerType,SigName) 
                          %>% summarise(CancerSigSampleCount=sum(SigPercent>=sigPercentThreshold),
                                        CancerSigSampleMutLoad=round(sum(Count),0)))

cancerSampleSigSummary = merge(cancerSampleSigSummary,cancerSampleCounts,by='CancerType',all.x=T)

# work counts of sample by cancer and signature
cancerTypesBySig = cancerSampleSigSummary %>% group_by(SigName) %>% summarise(SigCancerTypeCount=sum(CancerSigSampleCount>0))
rowIndex = data.frame(as.numeric(as.character(rownames(cancerTypesBySig))))
colnames(rowIndex) <- c("SigIndex")
cancerTypesBySig = cbind(rowIndex, cancerTypesBySig)

cancerSampleSigSummary = merge(cancerSampleSigSummary,cancerTypesBySig,by='SigName',all.x=T)

# work out median mutational load per sig and cancer type for samples with the signature
cancerSigMedianCounts = highMLSampleSigData %>% filter(Count>0) %>% group_by(SigName,CancerType) %>% summarise(MedianSampleCount=round(median(Count),0))

cancerSampleSigSummary = merge(cancerSampleSigSummary,cancerSigMedianCounts,by=c('SigName','CancerType'),all.x=T)

# calculate counts per MB for each signature
mutLoadMin = 0.25
mutLoadMax = 16
cancerSampleSigSummary = cancerSampleSigSummary %>% 
  mutate(
    SamplesWithSig=round(CancerSigSampleCount/CancerSampleCount,3),
    MedianLoadPerMb=ifelse(MedianSampleCount>0,round(MedianSampleCount/3e3,3),0),
    MedianLoadPerMb=ifelse(MedianLoadPerMb>mutLoadMax,mutLoadMax,ifelse(MedianLoadPerMb<mutLoadMin,mutLoadMin,MedianLoadPerMb))) 

residualsSummary = (snvDpSampleSigData %>% group_by(SampleId,CancerType) 
                    %>% summarise(SampleTotal=first(SampleTotal), 
                                  ResidualTotal=round(sum(ifelse(SigName=='Excess',-Count,ifelse(SigName=='Unalloc',Count,0))),0))
                    %>% mutate(ResidualPerc=round(ResidualTotal/SampleTotal,4))) 

sampleIdMap = read.csv(file = "/Users/jon/hmf/secure/SampleIdMap.csv", stringsAsFactors = F)
worstSamples = c('HMF001562A','HMF002896A') # these are the 2 SYD985 samples
worstSamplesMatrixData = snvDpSampleCounts %>% left_join(sampleIdMap, by = c("SampleId" =  "sampleId")) %>% filter(hmfSampleId %in% worstSamples) %>% select(-SampleId) %>% spread(hmfSampleId,Count)
worstSamplesMatrixData[is.na(worstSamplesMatrixData)] <- 0
worstSamplesMatrixData = within(worstSamplesMatrixData, rm(Bucket))

colorGradients = c("#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#08519c","#08306b")

plotBreaks = c(0.5, 1, 2, 4, 8, 16)
p1 = ggplot(cancerSampleSigSummary %>% filter(SamplesWithSig>0), aes(x=CancerType, y=reorder(SigName, -SigIndex))) +
        geom_point(aes(color=MedianLoadPerMb, size=SamplesWithSig))  +
        scale_color_gradientn(name='Median SNV Load Per Mb', trans='log2', colors = colorGradients, guide="colourbar",limits=c(mutLoadMin,mutLoadMax), breaks=plotBreaks, labels=plotBreaks) +
        scale_size_continuous(name = "Proportion of Samples with Signature", range = c(0, 5)) +
        theme(
          legend.margin=margin(t = 3, b = 0, unit = "pt"), 
          plot.margin = margin(l = 5, r = 5, t = 10, b = 0, unit = "pt"),
          legend.title = element_text(size = 6),  legend.key.size = unit(0.4, "cm"),
          legend.background=element_blank(), 
          legend.key=element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, size = 7, hjust = 0.5, vjust = 0.5),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank()
              ) +
        #theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0,size=10)) +
        theme(legend.position = "top") +
        #scale_x_discrete(position = 'top') + 
        ylab("Signature") + 
        guides(color = guide_colourbar(title.vjust = 0.5, title.position = "left", nrow = 2, barheight = unit(4, "pt"), barwidth = unit(50, "pt")))  
    #labs(size="Proportion of Samples with Signature")
      

#####################
# 2. Residuals Report
p2 = ggplot(residualsSummary, aes(CancerType, ResidualPerc)) + 
  geom_boxplot(varwidth=T, fill="lightblue", size = 0.2, outlier.size = 0.1) +
  ylab("Residuals") +
  theme(
    axis.title.x = element_blank(), 
    axis.ticks = element_blank(), axis.text.x = element_blank(),
      panel.border = element_blank(), panel.grid.minor.x = element_blank(),  panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(),
      strip.background = element_blank(), strip.text = element_blank(), legend.title = element_blank(), 
      legend.margin=margin(t = 0,unit = "pt"), 
      plot.margin = margin(t = 0, b = 25, unit = "pt"),
      legend.background=element_blank(), legend.key=element_blank()) +
    scale_y_continuous(labels = percent) 




######################################################
# 3. Bucket plots for samples with residuals above 50%

# these sample IDs need to be replaced with the HMF IDs

p3 = plot_96_profile(worstSamplesMatrixData) + xlab("") +
  theme(
    strip.text.y = element_text(size = 5), axis.text.y = element_text(size = 5), axis.title.y = element_text(size = 7), strip.text.x = element_text(size = 7),
    panel.grid.minor = element_blank()
    ) + 
  scale_y_continuous(labels = percent, breaks = c(0, 0.1, 0.2)) + ylab("Contribution")

p4 = plot_grid(p1, p2, ncol = 1, align = "v", labels = "auto", rel_heights =  c(4, 1), label_size= 8)
p5 = plot_grid(p4, p3, labels = c("", "c"), rel_heights = c(4, 1.2), ncol = 1, label_size= 8)

ggplot2::ggsave("~/hmf/RPlot/Extended Figure 3.pdf", p5, width = 183, height = 217, units = "mm", dpi = 300)
ggplot2::ggsave("~/hmf/RPlot/Extended Figure 3.png", p5, width = 183, height = 217, units = "mm", dpi = 300)
ggplot2::ggsave("~/hmf/RPlot/Extended Figure 3.eps", p5, width = 183, height = 217, units = "mm", dpi = 300)

