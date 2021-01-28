library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(stringi)
library(gtable)

# Parse and check inputs
args <- commandArgs(trailingOnly = TRUE)
print(args)

if(length(args) < 3)
{
  print("Requires arguments 1=SampleId, 2=Input directory, 3=Output directory 4=Nearest neighbour")
  stop()
}

sampleId = args[1]
inputDir = args[2]
outputDir = args[3]
runNN = ifelse(length(args) == 4 && args[4]=='true',T,F)
    
cupDataFile = paste0(inputDir,sampleId,'.cup.data.csv')
cupSimFile = paste0(inputDir,sampleId,'.cup.similarities.csv')

if(!file.exists(cupDataFile))
{
  print(sprintf('Missing CUP sample data file: %s', cupDataFile))
  stop()
}

cupSampleResults = read.csv(cupDataFile)

cupSampleSims = data.frame()

if(runNN & file.exists(cupSimFile))
{
  cupSampleSims = read.csv(cupSimFile)
}

print(sprintf('sample(%s) loaded %d results, %d similarities',sampleId,nrow(cupSampleResults),nrow(cupSampleSims)))

# Data preparation
cupPlotData = cupSampleResults %>% select(Category,ResultType,DataType,Value,RefCancerType,RefValue)

cupPlotData = cupPlotData %>% mutate(RefValueLabel=sprintf('%.0f%%',RefValue*100),
                                     DataType=stri_replace_all_fixed(DataType,'_',' '))

cupClassData = cupPlotData %>% filter(Category=='CLASSIFIER') %>% mutate(DataLabel=DataType)

cupGender = cupPlotData  %>% filter(DataType=='GENDER') %>% 
  mutate(DataLabel=sprintf('SEX (%s)',Value),
         PrevColour=ifelse(RefValue==0,'high',ifelse(RefValue<=0.02,'low','norm')))

# ensure features are shown in alphabetical order
cupFeatures = cupPlotData %>% filter(Category=='FEATURE'&ResultType!='LIKELIHOOD')
featureOrder = cupFeatures %>% group_by(DataType,Value) %>% count %>% arrange(DataType,Value) %>% ungroup()
rowIndex = data.frame(as.numeric(as.character(rownames(featureOrder))))
colnames(rowIndex) = c("FeatureIndex")
featureOrder = cbind(featureOrder,rowIndex)
cupFeatures = merge(cupFeatures,featureOrder %>% select(Value,FeatureIndex),by='Value',all.x=T)
cupFeatures = cupFeatures %>% mutate(DataLabel=Value)

cupOtherData = cupPlotData %>% filter(Category!='CLASSIFIER'&ResultType!='LIKELIHOOD'&DataType!='GENDER'&Category!='FEATURE') %>%
  mutate(DataLabel=ifelse(Category=='SNV_SIG'|Category=='SV',sprintf('%s (%.0f)',DataType,as.numeric(as.character(Value))),
                          ifelse(DataType %in% c('PURITY','PLOIDY','MS INDELS TMB','CHORD HRD'),sprintf('%s (%.2f)',DataType,as.numeric(as.character(Value))),
                                 sprintf('%s (%s)',DataType,Value))),
         PercColour=ifelse(RefValue<(-2)|RefValue>2,'high',ifelse(RefValue<0|RefValue>1,'medium',ifelse(RefValue<=0.02|RefValue>=0.98,'low','norm'))),
         DataTypeOrder=ifelse(Category=='SAMPLE_TRAIT',0,ifelse(Category=='SV',1,2)))


# Common themes for plots
font = 'sans'
defaultFontSize = 10

theme_set(theme_bw() + theme(axis.text=element_text(size=defaultFontSize),
                             axis.text.y=element_text(size=10,face='bold',family=font),
                             panel.grid=element_blank(),
                             panel.border=element_blank(),
                             axis.text.x.bottom=element_blank(),
                             axis.ticks.x=element_blank(),
                             axis.ticks.y=element_blank(),
                             legend.position='none'))

gradColourMin='white'
gradColourMax='seagreen3'
prevColours = c('high'='indianred3','medium'='salmon2','low'='peachpuff','norm'='white')

# Plotting functions
generate_css_plot<-function(simData,titleStr)
{
  rowIndex = data.frame(as.numeric(as.character(rownames(simData))))
  colnames(rowIndex) = c("RowIndex")
  simData = cbind(rowIndex,simData)
  simData = simData %>% mutate(CSS=sprintf('%d: %s',RowIndex,CSS))
  
  simPlotData = simData %>% gather('Field','Value',5:ncol(simData))
  
  simPlotData = simPlotData %>% mutate(TypeLabel=ifelse(Field=='CuppaCategory','Cuppa Category',ifelse(Field=='PrimaryLocation','Primary Location','Primary Sub-location')),
                                       TypeIndex=ifelse(Field=='CuppaCategory',1,ifelse(Field=='PrimaryLocation',2,3)))
  
  simPlot = ggplot(simPlotData,aes(x=reorder(TypeLabel,TypeIndex),y=reorder(CSS,-RowIndex))) +
    geom_tile(aes(fill='white'),colour='grey',stat="identity",position="identity") +
    scale_fill_manual(values=c('white')) +
    geom_text(aes(label=Value),size=3) +
    theme(axis.text.x=element_text(size=defaultFontSize,face='bold',family=font),
          axis.text.y=element_text(size=defaultFontSize,face='bold',family=font),
          panel.border=element_blank(),
          panel.grid=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          plot.margin=margin(t=0,r=10,b=0,l=0,unit='cm'),
          legend.position='none') +
    scale_x_discrete(position = "top") +
    labs(x='',y='',title=titleStr)
  
  return(simPlot)
}


# Generate Plots
summaryPlot = ggplot(cupClassData,aes(x=RefCancerType,y=DataLabel)) +
  geom_tile(aes(fill=RefValue),colour="grey",stat="identity",position="identity") +
  geom_text(aes(label=RefValueLabel),size=3) +
  scale_fill_gradient(low=gradColourMin,high=gradColourMax) +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top=element_text(angle=90,hjust=0,size=10,face='bold',family=font)) +
  labs(x='',y='',title='')

genderPlot = ggplot(cupGender,aes(x=RefCancerType,y=DataLabel)) +
  geom_tile(aes(fill=PrevColour),colour="grey",stat="identity",position="identity") +
  geom_text(aes(label=RefValueLabel),size=3) +
  scale_colour_manual(values=prevColours) +
  scale_fill_manual(values=prevColours,limits=names(prevColours)) +
  labs(x='',y='',title='')

featurePlot = ggplot(cupFeatures,aes(x=RefCancerType,y=reorder(DataLabel,-FeatureIndex))) +
  geom_tile(aes(fill=RefValue),colour="grey",stat="identity",position="identity") + 
  geom_text(aes(label=RefValueLabel),size=3) +
  scale_fill_gradient(low=gradColourMin,high=gradColourMax) +
  labs(x='',y='',title='FEATURES')

svTraitsPlot = ggplot(cupOtherData %>% filter(Category=='SV'|Category=='SAMPLE_TRAIT'),aes(x=RefCancerType,y=reorder(DataLabel,-DataTypeOrder))) +
  geom_tile(aes(fill=PercColour),colour="grey",stat="identity",position="identity") +
  geom_text(aes(label=RefValueLabel),size=3) +
  scale_colour_manual(values=prevColours) +
  scale_fill_manual(values=prevColours,limits=names(prevColours)) +
  labs(x='',y='',title='PERCENTILES')

sigPlot = ggplot(cupOtherData %>% filter(Category=='SNV'),aes(x=RefCancerType,y=DataLabel)) +
  geom_tile(aes(fill=PercColour),colour="grey",stat="identity",position="identity") +
  geom_text(aes(label=RefValueLabel),size=3) +
  scale_colour_manual(values=prevColours) +
  scale_fill_manual(values=prevColours,limits=names(prevColours)) +
  labs(x='',y='',title='SNV SIGNATURES')

outputFile = paste(outputDir,sampleId,'_cup_report.pdf', sep='')
print(paste("writing output to file: ", outputFile, sep=''))

featureLimit = 15
featureCount = nrow(cupFeatures %>% group_by(Value) %>% count)
titleHeight=10
summaryHeight=205
genderHeight=45
sigHeight=110
percHeight=90
disclaimer1Height=12
disclaimer2Height=14

if(featureCount > featureLimit)
{
  separateFeaturePlot=T
  print(sprintf('features(%d) print separately',featureCount))
  plotHeights = c(titleHeight,disclaimer1Height,disclaimer2Height,summaryHeight,genderHeight,sigHeight,percHeight)
} else
{
  separateFeaturePlot=F
  featureHeight = 45 + (featureCount - 1) * 10
  print(sprintf('features(%d) featureHeight(%d)',featureCount,featureHeight))
  plotHeights = c(titleHeight,disclaimer1Height,disclaimer2Height,summaryHeight,genderHeight,sigHeight,percHeight,featureHeight)
}

pdf(file=outputFile,height=14,width=20)

par(mar=c(1,1,1,1))

title = textGrob(paste(sampleId,' CUP Report',sep=''), gp=gpar(fontface="bold",fontsize=16))
disclaimer1 = textGrob(paste('All results and data described in this report are for research-use-only and have not been generated using a clinically validated and controlled procedure.'), gp=gpar(fontface="bold",fontsize=13))
disclaimer2 = textGrob(paste('These results should not be used for clinical decision making.'), gp=gpar(fontface="bold",fontsize=13))

if(separateFeaturePlot)
{
  grid.arrange(plot_grid(title,disclaimer1,disclaimer2,summaryPlot,genderPlot,sigPlot,svTraitsPlot,
                         ncol=1,nrow=7,rel_heights=plotHeights,align='v',axis='l'))
  
  featurePlot = featurePlot +
    scale_x_discrete(position = "top") +
    theme(axis.text.x.top=element_text(angle=90,hjust=0,size=10,face='bold',family=font))
  
  grid.arrange(plot_grid(featurePlot,ncol=1,nrow=1),newpage=T)
} else
{
  plot_grid(title,disclaimer1,disclaimer2,summaryPlot,genderPlot,sigPlot,svTraitsPlot,featurePlot,
            ncol=1,nrow=8,rel_heights=plotHeights,align='v',axis='l')
}

#if(nrow(cupSampleSims) > 0)
#{
#  cupSampleSims = cupSampleSims %>% 
#    mutate(CSS=sprintf('  %.1f%%', Score * 100),
#           CuppaCategory=MatchCancerType,
#           PrimaryType=ifelse(!is.na(MatchPrimarySubtype)&MatchPrimarySubtype!='',paste(MatchPrimaryType,MatchPrimarySubtype,sep=' - '),as.character(MatchPrimaryType)),
#           PrimaryLocation=ifelse(!is.na(MatchSubLocation)&MatchSubLocation!='',paste(MatchLocation,MatchSubLocation,sep=' - '),as.character(MatchLocation))) %>% 
#    select(MatchSampleId,MatchType,CSS,Score,CuppaCategory,PrimaryType,PrimaryLocation)
  
#  snv96TableData = cupSampleSims %>% filter(MatchType=='SNV_96_PAIRWISE_SIMILARITY') %>% select(-MatchType)
#  snvPlot = generate_css_plot(snv96TableData,'SNV 96 PAIRWISE SIMILARITY')
  
  # genPosTableData = cupSampleSims %>% filter(MatchType=='GENOMIC_POSITION_SIMILARITY') %>% select(-MatchType)
  # genPosPlot = generate_css_plot(genPosTableData,'GENOMIC POSITION SIMILARITY')
  
#  title = textGrob('Nearest Neighbours', gp=gpar(fontface="bold",fontsize=16))
  
#  plot_grid(title,snvPlot,ncol=1,nrow=2,rel_heights=c(1,5),align='v',axis='l',scale=0.9)
  # plot_grid(title,snvPlot,genPosPlot,ncol=1,nrow=3,rel_heights=c(1,5,5),align='v',axis='l',scale=0.9)
#}

dev.off()
