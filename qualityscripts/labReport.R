# Parse and check inputs
args <- commandArgs(trailingOnly = TRUE)
print(args)

if(length(args) < 2)
{
  print("Requires argument 1=Outputdir and 2=ReportedDate (yyyy-mm-dd), which is the date from which samples should be highlighted")
  stop()
}

library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(stringi)
library(gtable)
library(DBI)
library(RMySQL)

## RETRIEVE DATA
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
dbPilot = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")

purity = dbGetQuery(dbProd, "select * from purity")
sample = dbGetQuery(dbPilot, "select * from sample")
metric = dbGetQuery(dbPilot, "select * from metric")
flagstat = dbGetQuery(dbPilot, "select * from flagstat")

dbDisconnect(dbProd)
dbDisconnect(dbPilot)
rm(dbProd,dbPilot)

## INPUT ARGS
outputDir= args[1]
reportedDateText= args[2]
#outputDir=
#reportedDateText='2021-05-17'

purityData = purity %>% inner_join(sample) %>% dplyr::filter(arrivalDate>='2020-01-01') %>% mutate(reportedDate>=reportedDateText) 
metricData = metric %>% inner_join(sample) %>% dplyr::filter(arrivalDate>='2020-01-01') %>% mutate(reportedDate>=reportedDateText) 
flagstatData = flagstat %>% inner_join(sample) %>% dplyr::filter(arrivalDate>='2020-01-01') %>% mutate(reportedDate>=reportedDateText) 

## PURITY 
print("Including purity")
score = ggplot(purityData, aes(x="x",y=score)) + geom_boxplot(col="black") + 
  scale_y_sqrt() + 
  geom_jitter(size=0.45, alpha=0.3,col="royalblue3") +
  geom_point(data=purityData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
  geom_text(data=purityData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5)

contamination = ggplot(purityData, aes(x="x",y=contamination)) + geom_boxplot(col="black") + 
  scale_y_sqrt() + 
  geom_jitter(size=0.45, alpha=0.3,col="royalblue3") +
  geom_point(data=purityData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
  geom_text(data=purityData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5) 

unsupported = ggplot(purityData, aes(x="x",y=unsupportedCopyNumberSegments)) + geom_boxplot(col="black") + 
  scale_y_sqrt() + 
  geom_jitter(size=0.45, alpha=0.3,col="royalblue3") +
  geom_point(data=purityData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
  geom_text(data=purityData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5) 

## METRIC REF EXCLUDED
print("Including metric - ref excluded reads")
refTotal = ggplot(metricData, aes(x="",y=refPctExcTotal)) + geom_boxplot(col="black") + 
  scale_y_sqrt() + 
  geom_jitter(size=0.45, alpha=0.3,col="royalblue3") +
  geom_point(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5)

refOverlap = ggplot(metricData, aes(x="",y=refPctExcOverlap)) + geom_boxplot(col="black") + 
  scale_y_sqrt() + 
  geom_jitter(size=0.45, alpha=0.3,col="royalblue3") +
  geom_point(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5) 

refCapped = ggplot(metricData, aes(x="",y=refPctExcCapped)) + geom_boxplot(col="black") + 
  scale_y_sqrt() + 
  geom_jitter(size=0.45, alpha=0.3,col="royalblue3") +
  geom_point(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5) 

refBaseQ = ggplot(metricData, aes(x="",y=refPctExcBaseQ)) + geom_boxplot(col="black") + 
  scale_y_sqrt() + 
  geom_jitter(size=0.45, alpha=0.3,col="royalblue3") +
  geom_point(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5) 

refUnpaired = ggplot(metricData, aes(x="",y=refPctExcUnpaired)) + geom_boxplot(col="black") + 
  scale_y_sqrt() + 
  geom_jitter(size=0.45, alpha=0.3,col="royalblue3") +
  geom_point(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5) 

refDupe = ggplot(metricData, aes(x="",y=refPctExcDupe)) + geom_boxplot(col="black") + 
  scale_y_sqrt() + 
  geom_jitter(size=0.45, alpha=0.3,col="royalblue3") +
  geom_point(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5) 

refMapQ = ggplot(metricData, aes(x="",y=refPctExcMapQ)) + geom_boxplot(col="black") + 
  scale_y_sqrt() + 
  geom_jitter(size=0.45, alpha=0.3,col="royalblue3") +
  geom_point(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5) 

## METRIC TUMOR EXCLUDED
print("Including metric - tumor excluded reads")
tumorTotal = ggplot(metricData, aes(x="",y=tumorPctExcTotal)) + geom_boxplot(col="black") + 
  scale_y_sqrt() + 
  geom_jitter(size=0.45, alpha=0.3,col="royalblue3") +
  geom_point(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5)

tumorOverlap = ggplot(metricData, aes(x="",y=tumorPctExcOverlap)) + geom_boxplot(col="black") + 
  scale_y_sqrt() + 
  geom_jitter(size=0.45, alpha=0.3,col="royalblue3") +
  geom_point(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5) 

tumorCapped = ggplot(metricData, aes(x="",y=tumorPctExcCapped)) + geom_boxplot(col="black") + 
  scale_y_sqrt() + 
  geom_jitter(size=0.45, alpha=0.3,col="royalblue3") +
  geom_point(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5) 

tumorBaseQ = ggplot(metricData, aes(x="",y=tumorPctExcBaseQ)) + geom_boxplot(col="black") + 
  scale_y_sqrt() + 
  geom_jitter(size=0.45, alpha=0.3,col="royalblue3") +
  geom_point(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5) 

tumorUnpaired = ggplot(metricData, aes(x="",y=tumorPctExcUnpaired)) + geom_boxplot(col="black") + 
  scale_y_sqrt() + 
  geom_jitter(size=0.45, alpha=0.3,col="royalblue3") +
  geom_point(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5) 

tumorDupe = ggplot(metricData, aes(x="",y=tumorPctExcDupe)) + geom_boxplot(col="black") + 
  scale_y_sqrt() + 
  geom_jitter(size=0.45, alpha=0.3,col="royalblue3") +
  geom_point(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5) 

tumorMapQ = ggplot(metricData, aes(x="",y=tumorPctExcMapQ)) + geom_boxplot(col="black") + 
  scale_y_sqrt() + 
  geom_jitter(size=0.45, alpha=0.3,col="royalblue3") +
  geom_point(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5) 


##METRIC COVERAGES
print("Including metric - coverages")
refMean = ggplot(metricData, aes(x=refMeanCoverage,color=reportedDate>=reportedDateText)) + geom_histogram(binwidth=0.5) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(y=2, label=sampleId), hjust=-0.001, angle=45,size=2.5,position=position_jitter()) + 
  theme(legend.position = "none")
refSd = ggplot(metricData, aes(x=refSdCoverage,color=reportedDate>=reportedDateText)) + geom_histogram(binwidth=0.1) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(y=2, label=sampleId), hjust=-0.001, angle=45,size=2.5,position=position_jitter())+ 
  theme(legend.position = "none")
ref30x = ggplot(metricData, aes(x=refCoverage30xPercentage,color=reportedDate>=reportedDateText)) + geom_histogram(binwidth=0.005) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(y=2, label=sampleId), hjust=-0.001, angle=45,size=2.5,position=position_jitter())


tumorMean = ggplot(metricData, aes(x=tumorMeanCoverage,color=reportedDate>=reportedDateText)) + geom_histogram(binwidth=1) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(y=2, label=sampleId), hjust=-0.001, angle=45,size=2.5,position=position_jitter())+ 
  theme(legend.position = "none")
tumorSd = ggplot(metricData, aes(x=tumorSdCoverage,color=reportedDate>=reportedDateText)) + geom_histogram(binwidth=0.2) + 
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(y=2, label=sampleId), hjust=-0.001, angle=45,size=2.5,position=position_jitter())+ 
  theme(legend.position = "none")
tumor60x = ggplot(metricData, aes(x=tumorCoverage60xPercentage,color=reportedDate>=reportedDateText)) + geom_histogram(binwidth=0.002) +
  geom_text(data=metricData %>% dplyr::filter(reportedDate>=reportedDateText),aes(y=2, label=sampleId), hjust=-0.001, angle=45,size=2.5,position=position_jitter())



## FLAGSTAT [TODO]
#ggplot(flagstatData, aes(x="",y=tumorDuplicateProportion)) + geom_boxplot(col="black") + 
#  scale_y_sqrt() + 
#  geom_jitter(size=0.45, alpha=0.3) +
#  geom_point(data=flagstatData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
#  geom_text(data=flagstatData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5) 

#ggplot(flagstatData, aes(x="",y=tumorMappedProportion)) + geom_boxplot(col="black") + 
#  scale_y_sqrt() + 
#  geom_jitter(size=0.45, alpha=0.3) +
#  geom_point(data=flagstatData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.9) +
#  geom_text(data=flagstatData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5) 

#ggplot(flagstatData, aes(x="",y=tumorUniqueReadCount)) + geom_boxplot(col="black") + 
#  scale_y_sqrt() + 
#  geom_jitter(size=0.45, alpha=0.3) +
#  geom_point(data=flagstatData %>% dplyr::filter(reportedDate>=reportedDateText),color="RED",size=1.5,alpha=0.6) +
#  geom_text(data=flagstatData %>% dplyr::filter(reportedDate>=reportedDateText),aes(label=sampleId),hjust=-0.2, vjust=0.5, size=2.5) 

#ggplot(flagstatData, aes(x=tumorUniqueReadCount,color=reportedDate>=reportedDateText)) + geom_histogram(binwidth=5000000) +
#  geom_text(data=flagstatData %>% dplyr::filter(reportedDate>=reportedDateText),aes(y=2, label=sampleId), hjust=-0.001, angle=45,size=2.5,position=position_jitter())

#ggplot(flagstatData, aes(x=tumorSecondaryCount,color=reportedDate>=reportedDateText)) + geom_histogram(binwidth=50) +
#  geom_text(data=flagstatData %>% dplyr::filter(reportedDate>=reportedDateText),aes(y=2, label=sampleId), hjust=-0.001, angle=45,size=2.5,position=position_jitter())



## REPORT DETAILS
print("Start building lab report")

#titleHeight=5
#subtitleHeight=2
#purityplotsHeight=40
#plotHeights=c(titleHeight,subtitleHeight,purityplotsHeight)

outputFile = paste(outputDir,Sys.Date(),"_lab_report.pdf", sep='')  
print(paste("Writing output to file: ",outputFile, sep=''))

## BUILD PDF
pdf(file=outputFile,height=10,width=20)

title = textGrob(paste("Weekly lab report of",Sys.Date(),sep=' '))
subtitle = textGrob(paste("This report highlights data from samples from reported date",reportedDateText,sep=' '))

grid.arrange(plot_grid(title,subtitle,ncol=1))

purityPlots = arrangeGrob(score, contamination, unsupported, ncol=3,top=textGrob("Purity"))

metricRefExcPlots = arrangeGrob(refTotal, refOverlap, refCapped, refBaseQ, refUnpaired, refDupe, refMapQ, ncol=3,top=textGrob("Ref excluded reads"),
                                  layout_matrix = rbind(c(1, 6, 7),
                                                        c(1, 3, 2),
                                                        c(1, 4, 5)))

metricTumorExcPlots = arrangeGrob(tumorTotal, tumorOverlap, tumorCapped, tumorBaseQ, tumorUnpaired, tumorDupe, tumorMapQ, ncol=3,top=textGrob("Tumor excluded reads"),
                                    layout_matrix = rbind(c(1, 6, 7),
                                                          c(1, 3, 2),
                                                          c(1, 4, 5)))

metricCoveragePlots = arrangeGrob(refMean, refSd, ref30x, tumorMean, tumorSd, tumor60x, ncol=2,top=textGrob("Ref/tumor coverage details"),
                                    layout_matrix = rbind(c(1,1,2,2,3,3,3),
                                                          c(4,4,5,5,6,6,6))) 

grid.arrange(plot_grid(purityPlots,ncol=1),newpage=T)
grid.arrange(plot_grid(metricRefExcPlots,ncol=1),newpage=T)
grid.arrange(plot_grid(metricTumorExcPlots,ncol=1),newpage=T)
grid.arrange(plot_grid(metricCoveragePlots,ncol=1),newpage=T)

dev.off()
