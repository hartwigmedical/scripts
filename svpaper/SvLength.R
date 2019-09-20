library(tidyr)
library(dplyr)
library(ggplot2)

localPath = '~/hmf/analyses/SVAnalysis/'
sharedPath = '~/Dropbox/HMF Australia team folder/RData/'

load(paste0(sharedPath,"highestPurityCohort.RData"))

##### LOAD ##########
svDrivers = read.csv(file = paste0(localPath,"SVA_DRIVERS.csv"))
svData = read.csv(file = paste0(localPath,"SVA_SVS.csv"))
svCluster = read.csv(file = paste0(localPath,"SVA_CLUSTERS.csv"))
svCombined = left_join(svData,svCluster %>% select(SampleId, ClusterId, SuperType), by = c("SampleId", "ClusterId")) %>%
  filter(SuperType != 'ARTIFACT')

### SIMPLE TOP 50
create_top_50_violin_plot <- function(x, filter = 100) {
  xSummary = x %>% group_by(SampleId) %>% 
    count() %>% 
    arrange(-n) %>% 
    ungroup() %>% 
    top_n(50, n) %>%
    left_join(highestPurityCohort %>% select(SampleId = sampleId, CancerType = cancerType), by = "SampleId") %>%
    arrange(CancerType, -n)
  
  if (filter > 0) {
    xSummary = xSummary %>% filter(n > filter)
  }
  
  xSummary[is.na(xSummary)] <- "Unknown"
  xSummary = xSummary %>% mutate(Label = paste(SampleId, CancerType, n, sep = "-")) %>% select(SampleId, Label)
  
  plotDF = x %>% filter(SampleId %in% xSummary$SampleId) %>% left_join(xSummary, by = "SampleId") %>% mutate(Label = factor(Label, xSummary$Label, ordered = T))
  
  p1 = ggplot(plotDF, aes(Label, length)) + 
    geom_violin(scale = "area") + 
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 90))
  
  return (p1)
}


simpleTop50DupsDF = svCombined %>% filter(Type %in% c("DUP"), ClusterCount == 1) %>% mutate(length = PosEnd - PosStart + 1)
p1 = create_top_50_violin_plot(simpleTop50DupsDF) + ggtitle("Simple Top 50 Dups")

simpleTop50DelsDF = svCombined %>% filter(Type %in% c("DEL"), ClusterCount == 1) %>% mutate(length = PosEnd - PosStart + 1)
p2 = create_top_50_violin_plot(simpleTop50DelsDF) + ggtitle("Simple Top 50 Dels")

cdk12Drivers = svDrivers %>% filter(Gene == 'CDK12') %>% select(SampleId)
simpleTop50DupsCDK12 = svCombined %>% filter(Type %in% c("DUP"), ClusterCount == 1) %>% mutate(length = PosEnd - PosStart + 1) %>% filter(SampleId %in% cdk12Drivers$SampleId)
p3 = create_top_50_violin_plot(simpleTop50DupsCDK12) + ggtitle("Simple Top CDK12 Dups")

simpleTop50DelsCDK12 = svCombined %>% filter(Type %in% c("DEL"), ClusterCount == 1) %>% mutate(length = PosEnd - PosStart + 1) %>% filter(SampleId %in% cdk12Drivers$SampleId)
p4 = create_top_50_violin_plot(simpleTop50DelsCDK12) + ggtitle("Simple Top CDK12 Dels")

ccne1Drivers = svDrivers %>% filter(Gene == 'CCNE1') %>% select(SampleId)
simpleTop50DupsCCNE1 = svCombined %>% filter(Type %in% c("DUP"), ClusterCount == 1) %>% mutate(length = PosEnd - PosStart + 1) %>% filter(SampleId %in% ccne1Drivers$SampleId)
p5 = create_top_50_violin_plot(simpleTop50DupsCCNE1) + ggtitle("Simple Top CCNE1 Dups")

simpleTop50DelsCCNE1 = svCombined %>% filter(Type %in% c("DEL"), ClusterCount == 1) %>% mutate(length = PosEnd - PosStart + 1) %>% filter(SampleId %in% ccne1Drivers$SampleId)
p6 = create_top_50_violin_plot(simpleTop50DelsCCNE1) + ggtitle("Simple Top CCNE1 Dels")

plot_grid(p1,p2, ncol = 1)
plot_grid(p3,p4, ncol = 1)
plot_grid(p5,p6, ncol = 1)


### Cancer type TOP 50

create_top_50_violin_plot_cancer_type <- function(x) {
  xSummary = x %>% 
    left_join(highestPurityCohort %>% select(SampleId = sampleId, CancerType = cancerType), by = "SampleId") %>%
    group_by(CancerType) %>% 
    count() %>% 
    arrange(-n)
  
  cancerTypeCounts = highestPurityCohort %>% group_by(cancerType) %>% count() %>% ungroup() %>% mutate(CancerType = cancerType, weight = 1.0/n) %>% select(CancerType, weight)
  
  xSummary[is.na(xSummary)] <- "Unknown"
  xSummary = xSummary %>% mutate(Label = paste(CancerType, n, sep = "-")) %>% select(CancerType, Label)
  
  plotDF = x %>% 
    left_join(highestPurityCohort %>% select(SampleId = sampleId, CancerType = cancerType), by = "SampleId") %>%
    mutate(CancerType = ifelse(is.na(CancerType), "Unknown", CancerType)) %>%
    left_join(xSummary, by = "CancerType") %>%
    mutate(Label = factor(Label, xSummary$Label, ordered = T)) %>%
    left_join(cancerTypeCounts)
  
  p1 = ggplot(plotDF, aes(Label, length)) + 
    geom_violin(scale = "count", aes(weight = weight)) + 
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 90))
  
  return (p1)
}

cancerTypeDups = svCombined %>% filter(Type %in% c("DUP"), ClusterCount == 1) %>% mutate(length = PosEnd - PosStart + 1)
pCancerTypeDups = create_top_50_violin_plot_cancer_type(cancerTypeDups) + ggtitle("DUPS by CancerType")

cancerTypeDels = svCombined %>% filter(Type %in% c("DEL"), ClusterCount == 1) %>% mutate(length = PosEnd - PosStart + 1)
pCancerTypeDels = create_top_50_violin_plot_cancer_type(cancerTypeDels) + ggtitle("DELS by CancerType")

plot_grid(pCancerTypeDups, pCancerTypeDels, ncol = 1)


### DRIVER SIMPLE DELS / DUPS
plot_nine_violins <- function(svEnriched) {

  svSimple = svEnriched %>% filter(!IsFoldback, Type %in% c("DUP","DEL", "INV"), ClusterCount == 1) %>% mutate(Feature = paste0("Simple ",Type)) %>% select(Feature, Length)
  svComplex = svEnriched %>% filter(!IsFoldback, Type %in% c("DUP","DEL", "INV"), ResolvedType == 'COMPLEX') %>% mutate(Feature = paste0("Complex ",Type)) %>% select(Feature, Length)
  svRecipInv = svEnriched %>% filter(!IsFoldback, ResolvedType == "RECIP_INV") %>% mutate(Feature = "Recip Inv", Length = pmax(PosEnd, PosStart) - pmin(PosEnd, PosStart) + 1) %>% select(Feature, Length)
  
  fbStart = svEnriched %>% filter(FoldbackLenStart>=0) %>% select(SampleId,Id,ClusterId,Chr=ChrStart,Arm=ArmStart,FoldbackLength=FoldbackLenStart,OtherId=FoldbackLnkStart)
  fbEnd = svEnriched %>% filter(FoldbackLenEnd>=0) %>% select(SampleId,Id,ClusterId,Chr=ChrEnd,Arm=ArmEnd,FoldbackLength=FoldbackLenEnd,OtherId=FoldbackLnkEnd)
  svFoldbacks = rbind(fbStart,fbEnd)
  svFoldbacks = svFoldbacks %>% mutate(
      IsChained=(OtherId!=Id),                               
      SingleBreakend=(OtherId==Id&FoldbackLength==0),
      FoldbackId=ifelse(Id<OtherId,Id,OtherId),
      Feature = ifelse(IsChained, "Chained Foldback", "Simple Foldback")) %>%
    group_by(SampleId,ClusterId,FoldbackId,IsChained, Feature) %>% 
    summarise(Chr=first(Chr),Arm=first(Arm),FoldbackLength=first(FoldbackLength)) %>%
    ungroup() %>%
    select(Feature, Length = FoldbackLength)
  
  svComplete = bind_rows(svSimple, svRecipInv) %>% bind_rows(svComplex) %>% bind_rows(svFoldbacks)
  svFeatureLevels = unique(svComplete$Feature)
  svComplete = svComplete %>% mutate(Feature = factor(Feature, svFeatureLevels))
  
  ggplot(svComplete, aes(Feature, Length)) + 
    geom_violin(scale = "area") + 
    scale_y_log10() 
}

svEnriched = svCombined %>%
  mutate(
    IsFragile = FSStart == 'true'|FSEnd == 'true',
    IsLineElement = LEStart != 'false'| LEEnd !='false',
    IsFoldback = Type == 'INV' & FoldbackLnkStart > 0, 
    Length = PosEnd - PosStart + 1) 

plot_nine_violins(svEnriched) + ggtitle("No Filter")
plot_nine_violins(svEnriched %>% filter(IsFragile)) + ggtitle("So Very Fragile")
plot_nine_violins(svEnriched %>% filter(IsLineElement)) + ggtitle("Line Element")
