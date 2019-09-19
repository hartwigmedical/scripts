load("~/hmf/analysis/cohort/highestPurityCohort.RData")


svDrivers = read.csv(file = "/Users/jon/hmf/analysis/fusions/SV")
svCluster = read.csv(file = "/Users/jon/hmf/analysis/fusions/SVA_CLUSTERS.csv")
svData = read.csv(file = "/Users/jon/hmf/analysis/fusions/SVA_SVS.csv")
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

### DRIVER SIMPLE DELS / DUPS

x = svCombined %>% filter(Type %in% c("DEL"), ClusterCount == 1) %>% mutate(length = PosEnd - PosStart + 1)








### COMPLEX V SIMPLE
complexVSimpleDF = svCombined %>% filter(Type %in% c("DUP","DEL", "INV")) %>%
  mutate(IsFoldback = Type == 'INV' & FoldbackLnkStart > 0) %>%
  filter(ClusterCount == 1 | ResolvedType == 'COMPLEX', !IsFoldback) %>%
  mutate(
    length = PosEnd - PosStart + 1, 
    ResolvedType = ifelse(ResolvedType == 'COMPLEX', 'COMPLEX', 'SIMPLE'))

ggplot(complexVSimpleDF, aes(Type, length)) + 
  geom_violin(scale = "area") + 
  scale_y_log10() + 
  facet_wrap(~ResolvedType)


#### FOLDBACK
fbStart = svCombined %>% filter(FoldbackLenStart>=0) %>% select(SampleId,Id,ClusterId,Chr=ChrStart,Arm=ArmStart,FoldbackLength=FoldbackLenStart,OtherId=FoldbackLnkStart)
fbEnd = svCombined %>% filter(FoldbackLenEnd>=0) %>% select(SampleId,Id,ClusterId,Chr=ChrEnd,Arm=ArmEnd,FoldbackLength=FoldbackLenEnd,OtherId=FoldbackLnkEnd)
foldbacks = rbind(fbStart,fbEnd)
foldbacks = foldbacks %>% mutate(
  IsChained=(OtherId!=Id),                               
  SingleBreakend=(OtherId==Id&FoldbackLength==0),
  FoldbackId=ifelse(Id<OtherId,Id,OtherId))

# unique foldback data
foldbackSummary = foldbacks %>% group_by(SampleId,ClusterId,FoldbackId,IsChained) %>% summarise(Chr=first(Chr),Arm=first(Arm),FoldbackLength=first(FoldbackLength))
ggplot(foldbackSummary, aes(NA, FoldbackLength)) + 
  geom_violin(scale = "area") + 
  scale_y_log10() + 
  facet_wrap(~IsChained)


