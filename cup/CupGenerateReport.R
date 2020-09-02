library(dplyr)
library(grid)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(stringi)

theme_set(theme_bw())

args <- commandArgs(trailingOnly = TRUE)
print(args)

if(length(args) < 3)
{
  print("Requires arguments 1=sample counts, 2=signatures, 3=outputFile")
  stop()
}

sampleId = args[1]
inputFile = args[2]
outputDir = args[3]

# sampleId = 'CPCT02010440T'
# outputDir = '~/data/cup/reports/'    
# inputFile = paste('~/data/cup/samples/',sampleId,'.cup.data.csv',sep='')

cupSampleResults = read.csv(inputFile)

print(sprintf('sample(%s) loaded %d results',sampleId,nrow(cupSampleResults)))

cupPlotData = cupSampleResults %>% select(Category,DataType,Value,RefCancerType,RefValue)

cupPlotData = cupPlotData %>% mutate(RefValueLabel=sprintf('%.0f%%',RefValue*100),
                                     DataType=stri_replace_all_fixed(DataType,'_',' '))

# View(cupPlotData)

cupClassData = cupPlotData %>% filter(Category=='CLASSIFIER') %>% mutate(DataLabel=DataType)

cupGender = cupPlotData  %>% filter(DataType=='GENDER') %>% 
  mutate(DataLabel=sprintf('%s (%s)',DataType,Value),
         PrevColour=ifelse(RefValue==0,'high',ifelse(RefValue<=0.02,'low','norm')))

# ensure features are shown in alphabetical order
cupFeatures = cupPlotData %>% filter(Category=='FEATURE')
featureOrder = cupFeatures %>% group_by(DataType,Value) %>% count %>% arrange(DataType,Value) %>% ungroup()
rowIndex = data.frame(as.numeric(as.character(rownames(featureOrder))))
colnames(rowIndex) = c("FeatureIndex")
featureOrder = cbind(featureOrder,rowIndex)
cupFeatures = merge(cupFeatures,featureOrder %>% select(Value,FeatureIndex),by='Value',all.x=T)
cupFeatures = cupFeatures %>% mutate(DataLabel=Value)

cupOtherData = cupPlotData %>% filter(Category!='CLASSIFIER'&DataType!='GENDER'&Category!='FEATURE') %>%
  mutate(DataLabel=ifelse(Category=='SNV_SIG'|Category=='SV',sprintf('%s (%.0f)',DataType,as.numeric(as.character(Value))),
                          ifelse(DataType %in% c('PURITY','PLOIDY','MS INDELS TMB','CHORD HRD'),sprintf('%s (%.2f)',DataType,as.numeric(as.character(Value))),
                                 sprintf('%s (%s)',DataType,Value))),
         PercColour=ifelse(RefValue<(-2)|RefValue>2,'high',ifelse(RefValue<0|RefValue>1,'medium',ifelse(RefValue<=0.02|RefValue>=0.98,'low','norm'))),
         DataTypeOrder=ifelse(Category=='SAMPLE_TRAIT',0,ifelse(Category=='SV',1,2)))


## Generate Plots

font = 'sans'

gradColourMin='white'
gradColourMax='seagreen3'
prevColours = c('high'='indianred3','medium'='salmon2','low'='peachpuff','norm'='white')

summaryPlot = ggplot(cupClassData,aes(x=RefCancerType,y=DataLabel)) +
  geom_tile(aes(fill=RefValue),colour="grey",stat="identity",position="identity") +
  geom_text(aes(label=RefValueLabel),size=3) +
  scale_fill_gradient(low=gradColourMin,high=gradColourMax) +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.bottom=element_blank(),axis.ticks.x=element_blank(),
        axis.text.x.top=element_text(angle=90,hjust=0,size=10,face='bold',family=font),
        axis.ticks.y=element_blank(),
        legend.position='none') +
  theme(panel.background = element_blank(),panel.border = element_blank()) +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_text(size=10,face='bold',hjust=1,family=font)) +
  theme(legend.position='none') +
  labs(x='',y='',title='')

genderPlot = ggplot(cupGender,aes(x=RefCancerType,y=DataLabel)) +
  geom_tile(aes(fill=PrevColour),colour="grey",stat="identity",position="identity") +
  geom_text(aes(label=RefValueLabel),size=3) +
  scale_colour_manual(values=prevColours) +
  scale_fill_manual(values=prevColours,limits=names(prevColours)) +
  theme(panel.background = element_blank(),panel.border = element_blank()) +
  theme(axis.text.x.bottom=element_blank(),axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=10,face='bold',family=font)) +
  theme(legend.position='none') +
  labs(x='',y='',title='')

featurePlot = ggplot(cupFeatures,aes(x=RefCancerType,y=reorder(DataLabel,-FeatureIndex))) +
  geom_tile(aes(fill=RefValue),colour="grey",stat="identity",position="identity") + 
  geom_text(aes(label=RefValueLabel),size=3) +
  scale_fill_gradient(low=gradColourMin,high=gradColourMax) +
  theme(axis.text.x.bottom=element_blank(),axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position='none') +
  theme(panel.background = element_blank(),panel.border = element_blank()) +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_text(size=10,face='bold',hjust=1,family=font)) +
  theme(legend.position='none') +
  labs(x='',y='',title='FEATURE PREVELANCES')

svTraitsPlot = ggplot(cupOtherData %>% filter(Category=='SV'|Category=='SAMPLE_TRAIT'),aes(x=RefCancerType,y=reorder(DataLabel,-DataTypeOrder))) +
  geom_tile(aes(fill=PercColour),colour="grey",stat="identity",position="identity") +
  geom_text(aes(label=RefValueLabel),size=3) +
  scale_colour_manual(values=prevColours) +
  scale_fill_manual(values=prevColours,limits=names(prevColours)) +
  theme(axis.text.x.bottom=element_blank(),axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position='none') +
  theme(panel.background = element_blank(),panel.border = element_blank()) +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_text(size=10,face='bold',hjust=1,family=font)) +
  theme(legend.position='none') +
  labs(x='',y='',title='FEATURE PERCENTILES')

sigPlot = ggplot(cupOtherData %>% filter(Category=='SNV_SIG'),aes(x=RefCancerType,y=DataLabel)) +
  geom_tile(aes(fill=PercColour),colour="grey",stat="identity",position="identity") +
  geom_text(aes(label=RefValueLabel),size=3) +
  scale_colour_manual(values=prevColours) +
  scale_fill_manual(values=prevColours,limits=names(prevColours)) +
  theme(axis.text.x.bottom=element_blank(),axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position='none') +
  theme(panel.background = element_blank(),panel.border = element_blank()) +
  theme(axis.ticks.y=element_blank(),axis.text.y=element_text(size=10,face='bold',hjust=1,family=font)) +
  theme(legend.position='none') +
  labs(x='',y='',title='SNV SIGNATURES')

outputFile = paste(outputDir,sampleId,'_cup_report.pdf', sep='')
print(paste("writing output to file: ", outputFile, sep=''))

featureLimit = 15
featureCount = nrow(cupFeatures %>% group_by(Value) %>% count)

if(featureCount > featureLimit)
{
  separateFeaturePlot=T
  print(sprintf('features(%d) print separately',featureCount))
  plotHeights = c(10,195,45,110,90)
} else
{
  separateFeaturePlot=F
  featureHeight = 45 + (featureCount - 1) * 10
  print(sprintf('features(%d) featureHeight(%d)',featureCount,featureHeight))
  plotHeights = c(10,195,45,110,90,featureHeight)
}

pdf(file=outputFile,height=14,width=20)

par(mar=c(1,1,1,1))

title = textGrob(paste(sampleId,' CUP Report',sep=''), gp=gpar(fontface="bold",fontsize=16))

if(separateFeaturePlot)
{
  grid.arrange(plot_grid(title,summaryPlot,genderPlot,sigPlot,svTraitsPlot,
            ncol=1,nrow=7,rel_heights=plotHeights,align='v',axis='l'))
  
  featurePlot = featurePlot +
    scale_x_discrete(position = "top") +
    theme(axis.text.x.top=element_text(angle=90,hjust=0,size=10,face='bold',family=font))
    
  grid.arrange(plot_grid(featurePlot,ncol=1,nrow=1),newpage = T)
} else
{
  plot_grid(title,summaryPlot,genderPlot,sigPlot,svTraitsPlot,featurePlot,
            ncol=1,nrow=8,rel_heights=plotHeights,align='v',axis='l')
}

dev.off()

