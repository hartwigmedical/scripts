library(tidyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(stringr)

localPath = '~/hmf/analyses/SVAnalysis/'
sharedPath = '~/Dropbox/HMF Australia team folder/RData/'

##### LOAD DATA ##### 
load(paste0(sharedPath,"highestPurityCohort.RData"))
svDrivers = read.csv(file = paste0(localPath,"LNX_DRIVERS.csv"))
svGermline = read.csv(paste(localPath,'SVGermline.csv',sep=''), header = T, stringsAsFactors = F)
svDrivers=rbind(svDrivers %>% select(Gene,SampleId,DriverType,DriverLikelihood=Likelihood) %>% filter(DriverLikelihood>0.8), svGermline %>% mutate(DriverLikelihood=1,DriverType=ifelse(biallelic==1,'Germline Biallleic','Germline Monoallelic')) %>% select(Gene=gene,SampleId=sampleId,DriverType,DriverLikelihood))

svData = read.csv(file = paste0(localPath,"LNX_SVS.csv"))
svCluster = read.csv(file = paste0(localPath,"LNX_CLUSTERS.csv"))
svEnriched = left_join(svData,svCluster %>% separate('Annotations',c('SynLength1','SynLength2','SynGapLength'),sep=';') %>% 
  select(SampleId, ClusterId, SuperType,SynLength1,SynLength2,SynGapLength), by = c("SampleId", "ClusterId")) %>%
  filter(SuperType != 'ARTIFACT') %>% 
  mutate(IsFragile = FSStart == 'true'|FSEnd == 'true',
  IsLineElement = LEStart != 'false'| LEEnd !='false',
  IsFoldback = Type == 'INV' & FoldbackLnkStart > 0, 
  Length = PosEnd - PosStart + 1,
  HomLength = pmin(10,nchar(as.character(HomologyStart)))) %>%
  left_join(highestPurityCohort %>% select(SampleId = sampleId, CancerType = cancerType), by = "SampleId")

# Samples with key drivers
sample_filter_by_driver<-function(driverGene){
  return((svDrivers %>% filter(DriverLikelihood>0.8,Gene == driverGene) %>% .$SampleId))
}
CDK12samples=sample_filter_by_driver('CDK12')
CCNE1samples=sample_filter_by_driver('CCNE1')
BRCA2samples=sample_filter_by_driver('BRCA2')
BRCA1samples=sample_filter_by_driver('BRCA1')
TP53samples=sample_filter_by_driver('T')
######################### 

##### FUNCTIONS ######
plot_eight_violins <- function(svEnriched,complexType = 'COMPLEX',violinScale = 'area') {
  
  svSimple = svEnriched %>% filter(!IsFoldback, Type %in% c("DUP","DEL", "INV"), ClusterCount == 1|ResolvedType=='SIMPLE_GRP') %>% mutate(Feature = paste0("Simple ",Type)) %>% select(Feature, Length)
  svComplex = svEnriched %>% filter(!IsFoldback, Type %in% c("DUP","DEL", "INV"), ResolvedType == complexType) %>% mutate(Feature = paste0(complexType," ",Type)) %>% select(Feature, Length)
  svRecipInv = svEnriched %>% filter(!IsFoldback, ResolvedType == "RECIP_INV") %>% mutate(Feature = "Recip INV", Length = pmax(PosEnd, PosStart) - pmin(PosEnd, PosStart) + 1) %>% select(Feature, Length)
  
  fbStart = svEnriched %>% filter(FoldbackLenStart>=0) %>% select(SampleId,Id,ClusterId,Chr=ChrStart,Arm=ArmStart,FoldbackLength=FoldbackLenStart,OtherId=FoldbackLnkStart)
  fbEnd = svEnriched %>% filter(FoldbackLenEnd>=0) %>% select(SampleId,Id,ClusterId,Chr=ChrEnd,Arm=ArmEnd,FoldbackLength=FoldbackLenEnd,OtherId=FoldbackLnkEnd)
  svFoldbacks = rbind(fbStart,fbEnd) %>% mutate(Feature='Foldback',Length = FoldbackLength+1) %>% select(Feature,Length)  #+1 allows dsiplay of 0 length
  
  svComplete = bind_rows(svSimple,svComplex) %>% bind_rows(svRecipInv) %>% bind_rows(svFoldbacks) 
  svFeatureLevels = unique(svComplete$Feature)
  svComplete = svComplete %>% mutate(Feature = factor(Feature, svFeatureLevels))
  
  ggplot(svComplete, aes(Feature, Length)) + 
    geom_violin(scale = violinScale,fill='light blue',color=NA) + 
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 90,size=8),axis.title.x=element_blank(),
          axis.text.y = element_text(angle = 90,size=8),axis.title.y = element_text(size=8),
          panel.grid.major.y=element_line(linetype = 8,size=0.1))
}
create_violin_plot_cancer_type <- function(x,feature='Length',scaleLogY=T,violinScale = 'area') {
  
  cancerTypeCounts = highestPurityCohort %>% group_by(cancerType) %>% count() %>% arrange(-n) %>%  ungroup()  %>% 
    mutate(CancerType = cancerType, weight = 1.0/n,Label = paste0(CancerType," (n=", n,")")) %>% select(CancerType, weight,Label)
  plotDF = x %>%
    filter(!is.na(CancerType)) %>%
    left_join(cancerTypeCounts, by = "CancerType") %>%
    mutate(Label = factor(Label, cancerTypeCounts$Label, ordered = T))
  
  p1 = ggplot(plotDF, aes_string('Label', feature)) + 
    geom_violin(scale = violinScale, aes(weight = weight),fill='light blue',color=NA)+
    theme(axis.text.x = element_text(angle = 90,size=8),axis.title.x=element_blank(),
          axis.text.y = element_text(angle = 90,size=8),axis.title.y = element_text(size=8),
          panel.grid.major.y=element_line(linetype = 8,size=0.5))
  
  if (scaleLogY==T) {
    p1 = p1+ scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x))
  }
  return (p1)
}
create_top_n_violin_plot <- function(x, topN=40,minVariants = 0) {
  xSummary = x %>% group_by(SampleId) %>% 
    count() %>% 
    filter(n > minVariants) %>%
    arrange(-n) %>% 
    ungroup() %>% 
    top_n(topN, n) %>%
    left_join(highestPurityCohort %>% select(SampleId = sampleId, CancerType = cancerType), by = "SampleId") %>%
    arrange(CancerType, -n) %>% 
    mutate(Label = str_wrap(paste0(SampleId,' ', ifelse(is.na(CancerType),"Unknown",CancerType),"(n=", n,")"),16)) %>% select(SampleId, Label)
  
  plotDF = x %>% filter(SampleId %in% xSummary$SampleId) %>% left_join(xSummary, by = "SampleId") %>% mutate(Label = factor(Label, xSummary$Label, ordered = T))
  
  p1 = ggplot(plotDF, aes(Label, Length)) + 
    geom_violin(scale = 'area',fill='light blue',color=NA) + 
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 90,size=8),axis.title.x=element_blank(),
          axis.text.y = element_text(angle = 90,size=8),axis.title.y = element_text(size=8),
          panel.grid.major.y=element_line(linetype = 8,size=0.25))
  
  return (p1)
}
plot_synthetics <- function(svEnriched,violinScale = 'area') {
  svSimple = svEnriched %>% filter(Type %in% c("DUP","DEL"), ClusterCount == 1) %>% mutate(Feature = paste0("Simple ",Type)) %>% select(Feature, Length) %>% arrange(Feature)
  synStart = svEnriched %>% filter(ClusterCount>1,ResolvedType %in% c('DEL','DUP'),is.na(LnkSvStart)) %>% select(SampleId,Id,ClusterId,Chr=ChrStart,Arm=ArmStart,Pos=PosStart,ResolvedType)
  synEnd = svEnriched %>% filter(ClusterCount>1,ResolvedType %in% c('DEL','DUP'),is.na(LnkSvEnd)) %>% select(SampleId,Id,ClusterId,Chr=ChrEnd,Arm=ArmEnd,Pos=PosEnd,ResolvedType)
  svSyn = bind_rows(synStart,synEnd) %>% 
    group_by(SampleId,ClusterId,ResolvedType) %>%
    summarise(Length=max(Pos)-min(Pos)) %>% ungroup() %>% mutate(Feature = paste0("Synthetic ",ResolvedType)) %>% select(Feature, Length) %>% arrange(Feature)
  svSyn= bind_rows(svSyn,svSyn)
  svSyn= bind_rows(svSyn,svSyn)
  svSyn= bind_rows(svSyn,svSyn)
  svSyn= bind_rows(svSyn,svSyn)
  svSyn= bind_rows(svSyn,svSyn)
  svRecipDups = bind_rows(svEnriched %>% filter(ResolvedType=='RECIP_INV_DUPS'|ResolvedType=='RECIP_TRANS_DUPS', ClusterCount == 2) %>% mutate(Length=as.numeric(SynLength1),Feature = "RECIP_DUPS") %>% select(Feature, Length) %>% arrange(Feature),
                         svEnriched %>% filter(ResolvedType=='RECIP_INV_DUPS'|ResolvedType=='RECIP_TRANS_DUPS', ClusterCount == 2) %>% mutate(Length=as.numeric(SynLength2),Feature = "RECIP_DUPS") %>% select(Feature, Length) %>% arrange(Feature))
  svRecipDups= bind_rows(svRecipDups,svRecipDups)
  svRecipDups= bind_rows(svRecipDups,svRecipDups)
  
  fbStart = svEnriched %>% filter(FoldbackLenStart>=0) %>% select(SampleId,Id,ClusterId,Chr=ChrStart,Arm=ArmStart,FoldbackLength=FoldbackLenStart,OtherId=FoldbackLnkStart)
  fbEnd = svEnriched %>% filter(FoldbackLenEnd>=0) %>% select(SampleId,Id,ClusterId,Chr=ChrEnd,Arm=ArmEnd,FoldbackLength=FoldbackLenEnd,OtherId=FoldbackLnkEnd)
  
  svComplete = bind_rows(svSimple,svSyn) %>% bind_rows(svRecipDups) 
  #svComplete=svFoldbacks
  svFeatureLevels = unique(svComplete$Feature)
  svComplete = svComplete %>% mutate(Feature = factor(Feature, svFeatureLevels))
  
  ggplot(svComplete, aes(Feature, Length)) + 
    geom_violin(scale = violinScale,fill='light blue',color=NA) + 
    scale_y_log10()  +
    theme(axis.text.x = element_text(angle = 90,size=8),axis.title.x=element_blank(),
          axis.text.y = element_text(angle = 90,size=8),axis.title.y = element_text(size=8),
          panel.grid.major.y=element_line(linetype = 8,size=0.1),panel.grid.minor.y=element_line(linetype = 8,size=0.1))
}
create_violin_plot_resolved_type <- function(x,feature='Length',scaleLogY=T) {
  
  plotDF = x %>% mutate(Label = factor(ResolvedType, ordered = T)) #%>%left_join(cancerTypeCounts)
  
  p1 = ggplot(plotDF, aes_string('Label', feature)) + 
    geom_violin(bw=0.1,scale = "area",fill='light blue',color=NA)  +theme(axis.text.x = element_text(angle = 90)) 
  
  if (scaleLogY==T) {
    p1 = p1+ scale_y_log10()
  }
  return (p1)
}
######################
##### LENGTH ANALYSES #######
### 1. High level INV,DEL and DUP Lengths ###
p1 = plot_eight_violins(svEnriched,'COMPLEX','count') + ggtitle("ALL simple & complex variants + foldback & reciprocal INV")
p2 = plot_eight_violins(svEnriched %>% filter(IsFragile),'COMPLEX','count') + ggtitle("Fragile site only")
plot_grid(p1, p2, ncol = 1)

### 2. Simple DELS and DUPS  by Cancer Type ###
pCancerTypeDups = create_violin_plot_cancer_type(svEnriched %>% filter(Type %in% c("DUP"), ClusterCount == 1)) + ggtitle("Length Distribution: Simple DUP by CancerType (count per sample)")
pCancerTypeDels = create_violin_plot_cancer_type(svEnriched %>% filter(Type %in% c("DEL"), ClusterCount == 1)) + ggtitle("Simple DEL by CancerType (count per sample)")
plot_grid(pCancerTypeDups, pCancerTypeDels, ncol = 1)

### 3. Simple DELS and DUPS top 50 samples ####
#TO DO: can we colour violin by enriched driver genes?
p1 = create_top_n_violin_plot(svEnriched %>% filter(Type %in% c("DUP"), ClusterCount == 1),50) + ggtitle("Simple Top 50 Dups")
p2 = create_top_n_violin_plot(svEnriched %>% filter(Type %in% c("DEL"), ClusterCount == 1) ,50) + ggtitle("Simple Top 50 Dels")
plot_grid(p1,p2, ncol = 1)

### 4.Simple DUP top N by enriched driver Gene ###
p3 = create_top_n_violin_plot(svEnriched %>% filter(Type %in% c("DUP"), ClusterCount == 1,(SampleId %in% CDK12samples))) + ggtitle("Simple Top CDK12 Dups")
p4 = create_top_n_violin_plot(svEnriched %>% filter(Type %in% c("DUP"), ClusterCount == 1,(SampleId %in% CCNE1samples))) + ggtitle("Simple Top CCNE1 Dups")
p5 = create_top_n_violin_plot(svEnriched %>% filter(Type %in% c("DUP"), ClusterCount == 1,(SampleId %in% BRCA1samples))) + ggtitle("Simple Top BRCA1 Dups")
plot_grid(p3,p4,p5, ncol = 1)
print(p3)

### 5.Simple DEL top N by enriched driver Gene ###
p6 = create_top_n_violin_plot(svEnriched %>% filter(Type %in% c("DUP"), ClusterCount == 1,(SampleId %in% BRCA1samples))) + ggtitle("Simple Top BRCA1 DEL")
p7 = create_top_n_violin_plot(svEnriched %>% filter(Type %in% c("DUP"), ClusterCount == 1,(SampleId %in% BRCA2samples))) + ggtitle("Simple Top BRCA2 DEL")
plot_grid(p6,p7, ncol = 1)

### 4.Simple DUP top N for enriched cancer types ###
p10 = create_top_n_violin_plot(svEnriched %>% filter(Type=="DUP", ClusterCount == 1,CancerType=='Ovary',20)) + ggtitle("Simple top Ovary DUP")
p11 = create_top_n_violin_plot(svEnriched %>% filter(Type=="DUP", ClusterCount == 1,CancerType %in% c('Esophagus','Stomach')),20) + ggtitle("Simple Top Esophagus & Stomach DUP")
p12 = create_top_n_violin_plot(svEnriched %>% filter(Type=="DUP", ClusterCount == 1,CancerType=='Prostate'),20) + ggtitle("Simple Top Prostate DUP")
p13 = create_top_n_violin_plot(svEnriched %>% filter(Type=="DUP", ClusterCount == 1,CancerType=='Colon/Rectum'),20) + ggtitle("Simple Top CRC DUP")
p14 = create_top_n_violin_plot(svEnriched %>% filter(Type=="DUP", ClusterCount == 1,CancerType=='Breast'),20) + ggtitle("Simple Top Breast DUP")
plot_grid(p10,p11,p12,p13,p14, ncol = 1)

### 5. Simple vs synthetic lengths overall ###
plot_synthetics(svEnriched,'count') + ggtitle("Synthetics")

### 6. Simple vs synthetic lengths by driver gene ###
p0= plot_synthetics(svEnriched,'count') + ggtitle("All Samples")
p1=plot_synthetics(svEnriched %>% filter(SampleId %in% CDK12samples),'count') + ggtitle("Samples with  CDK12 drivers")
p2=plot_synthetics(svEnriched %>% filter(SampleId %in% CCNE1samples),'count') + ggtitle("Samples with CCNE1 drivers")
p3=plot_synthetics(svEnriched %>% filter(SampleId %in% BRCA1samples),'count') + ggtitle("Samples with BRCA1 drivers")
p4=plot_synthetics(svEnriched %>% filter(SampleId %in% BRCA2samples),'count') + ggtitle("Samples with BRCA2 drivers")
plot_grid(p0,p1,p2,p3,p4,ncol=1)

### 7. Reciprocal DUPS ###
pairClusters = svCluster %>% filter(ClusterCount==2,ResolvedType=='RECIP_TRANS_DUPS'|ResolvedType=='RECIP_INV_DUPS') %>% separate('Annotations',c('SynLength1','SynLength2','SynGapLength'),sep=';') %>% 
  mutate(SynGapLength=as.numeric(as.character(SynGapLength)),SynLength1=as.numeric(as.character(SynLength1)),SynLength2=as.numeric(as.character(SynLength2)))
p0=ggplot(data=pairClusters,aes(SynLength1,SynLength2)) + geom_hex(bins=30) + facet_wrap(~ResolvedType) + scale_x_log10() + scale_y_log10()+
  scale_fill_gradient(low = 'light blue', high = 'dark blue') + ggtitle("All")
p1=ggplot(data=pairClusters %>% filter(SampleId %in% BRCA1samples),aes(SynLength1,SynLength2)) + geom_hex(bins=30) + facet_wrap(~ResolvedType) + scale_x_log10() + scale_y_log10()+
  scale_fill_gradient(low = 'light blue', high = 'dark blue') + ggtitle("BRCA1")
p2=ggplot(data=pairClusters %>% filter(SampleId %in% CCNE1samples),aes(SynLength1,SynLength2)) + geom_hex(bins=30) + facet_wrap(~ResolvedType) + scale_x_log10() + scale_y_log10()+
  scale_fill_gradient(low = 'light blue', high = 'dark blue') + ggtitle("CCNE1")
p3=ggplot(data=pairClusters %>% filter(SampleId %in% CDK12samples),aes(SynLength1,SynLength2)) + geom_hex(bins=30) + facet_wrap(~ResolvedType) + scale_x_log10() + scale_y_log10()+
  scale_fill_gradient(low = 'light blue', high = 'dark blue') + ggtitle("CDK12")
plot_grid(p0,p1,p2,p3,ncol=1)

################################

##### Variant Counts  #######
#1. By cancer Type and Resolved Type
topResolvedTypes = svEnriched  %>% filter(SuperType!='INCOMPLETE') %>% group_by(ResolvedType) %>% count %>% filter(n>2500) %>% .$ResolvedType
allCounts = svEnriched  %>% group_by(SampleId,CancerType,ResolvedType) %>% summarise(count=n()) %>% tidyr::complete(SampleId,CancerType,ResolvedType,fill=list(count=0.5))
create_violin_plot_cancer_type(allCounts %>% filter(ResolvedType %in% topResolvedTypes) ,'count',T,'width') +  ggtitle("Counts of Variants by Resolved Type") + facet_wrap(~ResolvedType)
create_violin_plot_cancer_type(allCounts %>% filter(ResolvedType=='RECIP_TRANS') ,'count',T,'area') +  ggtitle("Counts of Variants by Resolved Type") + facet_wrap(~ResolvedType)
View(svEnriched %>% group_by(SuperType,ResolvedType) %>% count)
################################

##### DB LENGTH ANALYSES #######
### 1. Short DB lengths by Resolved Type ###
#TODO: switch to breakend based
print(create_violin_plot_resolved_type(svEnriched %>% filter(ClusterCount>1,DBLenStart<100,DBLenStart>-100,SuperType!='INCOMPLETE',
                                                             !ResolvedType %in% c('DOUBLE_MINUTE')),'DBLenStart',F) + ggtitle(" DB by ResolvedType (abs counts)"))
################################

##### HOM LENGTH ANALYSES ######
### 1. HOM engths by Resolved Type ###
# TODO: convert to bar chart.  
# LINE has least homology, followed by COMPLEX and RECIP events.  PAIR_OTHER has a very long upper tail.  DELS and DUPS have nuances by length
View(svEnriched %>% filter(SuperType!='INCOMPLETE',Type!='SGL',Type!='INF') %>% group_by(ResolvedType,HL=pmin(20,nchar(as.character(Homology)))) %>% count %>% spread(HL,n))

################################

########## SCRATCH##############
# 4 main genes are almost mutually exclusive
View(svDrivers %>% filter(Gene %in% c('CDK12','BRCA1','BRCA2','CCNE1')) %>% group_by(Gene,SampleId) %>% count %>% spread(Gene,n))
plot_eight_violins(svEnriched %>% filter(SampleId=='CPCT02060120T'),'COMPLEX','count') + ggtitle("CPCT02060120T")

################################

p1 = ggplot(plotDF, aes_string('Label','SampleId')) + 
  geom_violin(scale = "area", aes(weight = weight),fill='light blue',color=NA)  +theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text.x = element_text(angle = 90,size=8),axis.title.x=element_blank(),
        axis.text.y = element_text(angle = 90,size=8),axis.title.y = element_text(size=8),
        panel.grid.major.y=element_line(linetype = 8,size=0.1))
print(p1)
View(svEnriched %>%  filter(SampleId %in% (svDrivers %>% filter(Gene != 'BRCA1') %>% .$SampleId)) %>% filter(Type %in% c("INV"), ClusterCount == 2,ResolvedType=='RECIP_INV') %>% filter(CancerType!='OOvary') %>% group_by(SampleId,CancerType) %>% count)

View(svDrivers %>% filter(SampleId %in% c('CPCT02080092T','CPCT02020659T','CPCT02070410T','CPCT02100117T','DRUP01070063T','DRUP01180001T','CPCT02070145T','CPCT02020630T','CPCT02060056T')) %>% group_by(Gene) %>% count)
print(plot_eight_violins(svEnriched %>% filter(SampleId=='CPCT02020683T'),'COMPLEX','count') + ggtitle("ALL simple & complex variants + foldback & reciprocal INV"))# control
print(plot_eight_violins(svEnriched %>% filter(SampleId=='CPCT02100161T'),'COMPLEX','count') + ggtitle("ALL simple & complex variants + foldback & reciprocal INV"))# control
print(plot_eight_violins(svEnriched %>% filter(SampleId=='CPCT02250003T'),'COMPLEX','count') + ggtitle("ALL simple & complex variants + foldback & reciprocal INV"))# control
print(plot_eight_violins(svEnriched %>% filter(SampleId=='CPCT02060056T'),'COMPLEX','count') + ggtitle("ALL simple & complex variants + foldback & reciprocal INV"))
print(plot_eight_violins(svEnriched %>% filter(SampleId=='CPCT02020422T'),'COMPLEX','count') + ggtitle("ALL simple & complex variants + foldback & reciprocal INV"))
print(plot_eight_violins(svEnriched %>% filter(SampleId=='CPCT02020867T'),'COMPLEX','count') + ggtitle("ALL simple & complex variants + foldback & reciprocal INV"))
print(plot_eight_violins(svEnriched %>% filter(SampleId=='CPCT02020630T'),'COMPLEX','count') + ggtitle("ALL simple & complex variants + foldback & reciprocal INV"))
print(plot_eight_violins(svEnriched %>% filter(SampleId=='CPCT02070145T'),'COMPLEX','count') + ggtitle("ALL simple & complex variants + foldback & reciprocal INV"))
print(plot_eight_violins(svEnriched %>% filter(SampleId=='CPCT02080092T'),'COMPLEX','count') + ggtitle("ALL simple & complex variants + foldback & reciprocal INV"))
print(plot_eight_violins(svEnriched %>% filter(SampleId=='CPCT02020659T'),'COMPLEX','count') + ggtitle("ALL simple & complex variants + foldback & reciprocal INV"))
print(plot_eight_violins(svEnriched %>% filter(SampleId=='CPCT02070410T'),'COMPLEX','count') + ggtitle("ALL simple & complex variants + foldback & reciprocal INV"))
print(plot_eight_violins(svEnriched %>% filter(SampleId=='CPCT02100117T'),'COMPLEX','count') + ggtitle("ALL simple & complex variants + foldback & reciprocal INV"))
print(plot_eight_violins(svEnriched %>% filter(SampleId=='DRUP01070063T'),'COMPLEX','count') + ggtitle("ALL simple & complex variants + foldback & reciprocal INV"))
print(plot_eight_violins(svEnriched %>% filter(SampleId=='DRUP01180001T'),'COMPLEX','count') + ggtitle("ALL simple & complex variants + foldback & reciprocal INV"))
print(plot_eight_violins(svEnriched %>% filter(SampleId=='CPCT02230050T'),'COMPLEX','count') + ggtitle("ALL simple & complex variants + foldback & reciprocal INV"))
print(plot_eight_violins(svEnriched %>% filter(SampleId=='CPCT02450007T'),'COMPLEX','count') + ggtitle("ALL simple & complex variants + foldback & reciprocal INV"))
View(svEnriched %>% filter(SampleId=='CPCT02030381T') %>% select(Length,ResolvedType,Type,IsFoldback))
View(svEnriched %>% filter(SampleId=='CPCT02450007T') %>% group_by(ResolvedType,Type,IsFoldback) %>% count)
create_violin_plot_cancer_type(svEnriched %>%  filter(SampleId %in% (svDrivers %>% filter(Gene != 'CCNE1') %>% .$SampleId)) %>% filter(Type %in% c("INV"), ClusterCount == 2,ResolvedType=='RECIP_INV')) + ggtitle("Length Distribution: Simple DUP by CancerType (count per sample)")
create_top_n_violin_plot(svEnriched  %>% filter(Type %in% c("INV"), ClusterCount == 2,ResolvedType=='RECIP_INV') %>% filter(CancerType!='OOvary'),40) + ggtitle("Simple top Ovary DUP")

View(svEnriched  %>% filter(Type %in% c("INV"), ClusterCount == 2,ResolvedType=='RECIP_INV') %>% group_by(SampleId,CancerType,long=Length<2e5) %>% count %>% spread(long,n))
View(svEnriched %>% filter(SampleId=='CPCT02070410T') %>% group_by(ResolvedType,Type,IsFoldback) %>% count)

p10 = create_top_n_violin_plot(svEnriched %>% filter(Type=="DEL", ClusterCount == 1,CancerType=='Ovary',20)) + ggtitle("Simple top Ovary DUP")
p11 = create_top_n_violin_plot(svEnriched %>% filter(Type=="DEL", ClusterCount == 1,CancerType %in% c('Esophagus','Stomach')),20) + ggtitle("Simple Top Esophagus & Stomach DUP")
p12 = create_top_n_violin_plot(svEnriched %>% filter(Type=="DEL", ClusterCount == 1,CancerType=='Prostate'),20) + ggtitle("Simple Top Prostate DUP")
p13 = create_top_n_violin_plot(svEnriched %>% filter(Type=="DEL", ClusterCount == 1,CancerType=='Colon/Rectum'),20) + ggtitle("Simple Top CRC DUP")
p14 = create_top_n_violin_plot(svEnriched %>% filter(Type=="DEL", ClusterCount == 1,CancerType=='Breast'),20) + ggtitle("Simple Top Breast DUP")
plot_grid(p10,p11,p12,p13,p14, ncol = 1)

create_violin_plot_cancer_type(svEnriched %>% filter(SampleId %in% CDK12samples,ClusterCount == 2,ResolvedType=='RECIP_TRANS_DUPS') %>%mutate(SynLength1=as.numeric(SynLength1)),'SynLength1',violinScale = 'count') + ggtitle("RECIP_INV")

View(svEnriched %>% filter(SampleId %in% CDK12samples,ClusterCount == 2,ResolvedType=='RECIP_TRANS_DUPS') %>% mutate(SynLength1=as.numeric(SynLength1)))
p10 = create_top_n_violin_plot(svEnriched %>% filter(Type=="DUP", ClusterCount == 1,CancerType=='Ovary',20)) + ggtitle("Simple top Ovary DUP")
