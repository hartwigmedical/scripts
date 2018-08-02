library(tidyr)
library(dplyr)
library(GenomicRanges)
library(MutationalPatterns)
detach("package:purple", unload=TRUE)
library(purple)

combined_signature <- function(GT1, GT2, MT1, MT2) {
  sampleFactors = c("GT1","GT2","MT1","MT2", 
                    "GT1,GT2","GT1,MT1","GT1,MT2","GT2,MT1","GT2,MT2","MT1,MT2",
                    "GT1,GT2,MT1","GT1,GT2,MT2","GT1,MT1,MT2","GT2,MT1,MT2",
                    "GT1,GT2,MT1,MT2")
  
  GT1Contrib = GT1 %>% 
    unite(scope, GT1, GT2, MT1, MT2, sep = ",") %>%
    mutate(scope = gsub(',NA',"", scope)) %>%
    select(sampleId = scope, startPosition, endPosition, type)
  
  GT2Contrib = GT2 %>% 
    filter(is.na(GT1)) %>%
    unite(scope, GT2, MT1, MT2, sep = ",") %>%
    mutate(scope = gsub(',NA',"", scope)) %>%
    select(sampleId = scope, startPosition, endPosition, type)
  
  MT1Contrib = MT1 %>% 
    filter(is.na(GT1), is.na(GT2)) %>%
    unite(scope, MT1, MT2, sep = ",") %>%
    mutate(scope = gsub(',NA',"", scope)) %>%
    select(sampleId = scope, startPosition, endPosition, type)
  
  MT2Contrib = MT2 %>% 
    filter(is.na(GT1), is.na(GT2), is.na(MT1)) %>%
    mutate(scope = "MT2") %>%
    select(sampleId = scope, startPosition, endPosition, type)
  
  total = bind_rows(GT1Contrib, GT2Contrib) %>% bind_rows(MT1Contrib) %>% bind_rows(MT2Contrib) %>%  mutate(sampleId = factor(sampleId, rev(sampleFactors)))
  return (sv_signature(total))
}


combined_scope <- function(GT1, GT2, MT1, MT2) {
  result = list()

  GT1 = GT1 %>% mutate(GT1 = "GT1", GT2 = NA, MT1 = NA, MT2 = NA, GT1ClusterCount = ClusterCount, GT2ClusterCount = NA, MT1ClusterCount = NA, MT2ClusterCount = NA)
  GT2 = GT2 %>% mutate(GT1 = NA, GT2 = "GT2", MT1 = NA, MT2 = NA, GT1ClusterCount = NA, GT2ClusterCount = ClusterCount, MT1ClusterCount = NA, MT2ClusterCount = NA)
  MT1 = MT1 %>% mutate(GT1 = NA, GT2 = NA, MT1 = "MT1", MT2 = NA, GT1ClusterCount = NA, GT2ClusterCount = NA, MT1ClusterCount = ClusterCount, MT2ClusterCount = NA)
  MT2 = MT2 %>% mutate(GT1 = NA, GT2 = NA, MT1 = NA, MT2 = "MT2", GT1ClusterCount = NA, GT2ClusterCount = NA, MT1ClusterCount = NA, MT2ClusterCount = ClusterCount)
  
  scope = append_scope("GT1", "GT2", GT1, GT2); GT1 = scope[["S1"]]; GT2 = scope[["S2"]]
  scope = append_scope("GT1", "MT1", GT1, MT1); GT1 = scope[["S1"]]; MT1 = scope[["S2"]]
  scope = append_scope("GT1", "MT2", GT1, MT2); GT1 = scope[["S1"]]; MT2 = scope[["S2"]]
  scope = append_scope("GT2", "MT1", GT2, MT1); GT2 = scope[["S1"]]; MT1 = scope[["S2"]]
  scope = append_scope("GT2", "MT2", GT2, MT2); GT2 = scope[["S1"]]; MT2 = scope[["S2"]]
  scope = append_scope("MT1", "MT2", MT1, MT2); MT1 = scope[["S1"]]; MT2 = scope[["S2"]]

  result[["GT1"]] <- GT1
  result[["GT2"]] <- GT2
  result[["MT1"]] <- MT1
  result[["MT2"]] <- MT2
  
  return (result)
}

append_scope <- function(sample1Name, sample2Name, S1, S2) {
  S1_S2 = sv_overlaps(S1, S2, maxgap = 100)
  
  S1[S1_S2$queryHits, sample2Name] <- sample2Name
  S1[S1_S2$queryHits, paste0(sample2Name, "ClusterCount")] <- S2[S1_S2$subjectHits, "ClusterCount"]
  
  S2[S1_S2$subjectHits, sample1Name] <- sample1Name
  S2[S1_S2$subjectHits, paste0(sample1Name, "ClusterCount")] <- S1[S1_S2$queryHits, "ClusterCount"]
  
  result = list()
  result[["S1"]] <- S1
  result[["S2"]] <- S2
  
  return (result)
}


old_column_names <- function(variants) {
  variants %>% 
    mutate(
      sampleId = SampleId,
      startChromosome = ChrStart,  
      endChromosome = ChrEnd, 
      startPosition = PosStart,
      endPosition = PosEnd,
      startOrientation = OrientStart,
      endOrientation = OrientEnd,
      type = Type) %>%
    select(-SampleId, -ChrStart, - ChrEnd, -PosStart, -PosEnd, -OrientStart, -OrientEnd) %>%
    select(sampleId, startChromosome, endChromosome, startPosition,endPosition,startOrientation,endOrientation,type, everything())
}

gridssCohortVariants = old_column_names(read.csv('/Users/jon/hmf/gridss/clustering/CLUSTER_GRIDSS.csv', header = T, stringsAsFactors = F))
gridssCohort = gridssCohortVariants %>% 
  select(sampleId) %>% 
  distinct() %>% 
  mutate(patientId = substring(sampleId, 1, 12)) %>% 
  group_by(patientId) %>% 
  filter(n() > 1)

prodCohortVariants = old_column_names(read.csv('/Users/jon/hmf/gridss/clustering/CLUSTER_V23.csv', header = T, stringsAsFactors = F)) %>% 
  filter(sampleId %in% gridssCohort$sampleId)
#save(prodCohortVariants, file = "/Users/jon/hmf/gridss/RData/mantaVariants.RData")
#load(file = "/Users/jon/hmf/gridss/RData/mantaVariants.RData")
#prodCohortVariants = old_column_names(read.csv('/Users/jon/hmf/gridss/clustering/CLUSTER_V24_COMPACT.csv', header = T, stringsAsFactors = F)) %>% 
#  filter(sampleId %in% gridssCohort$sampleId)
prodCohort = prodCohortVariants %>% 
  select(sampleId) %>% 
  distinct() %>% 
  mutate(patientId = substring(sampleId, 1, 12)) %>% 
  group_by(patientId) %>% 
  filter(n() > 1)

cohort = inner_join(gridssCohort, prodCohort, by = c("patientId", "sampleId")) %>%
  group_by(patientId) %>% 
  filter(n() == 2) %>%
  arrange(sampleId) %>%
  mutate(scope = paste0("Sample", row_number())) %>%
  spread(scope, sampleId)

rm(gridssCohort, prodCohort)

######## START PROCESSING SINGLE PATIENT
annotatedGT1 = data.frame(stringsAsFactors = F)
annotatedGT2 = data.frame(stringsAsFactors = F)
annotatedMT1 = data.frame(stringsAsFactors = F)
annotatedMT2 = data.frame(stringsAsFactors = F)

for (i in 1:nrow(cohort)) {
  patient = cohort[i, ]
  GT1 = gridssCohortVariants %>% filter(sampleId == patient$Sample1)
  GT2 = gridssCohortVariants %>% filter(sampleId == patient$Sample2)
  MT1 = prodCohortVariants %>% filter(sampleId == patient$Sample1)
  MT2 = prodCohortVariants %>% filter(sampleId == patient$Sample2)
  
  combinedScope = combined_scope(GT1, GT2, MT1, MT2)
  GT1 = combinedScope[["GT1"]]
  GT2 = combinedScope[["GT2"]]
  MT1 = combinedScope[["MT1"]]
  MT2 = combinedScope[["MT2"]]
  
  #signature = combined_signature(GT1, GT2, MT1, MT2)
  #sigPlot = plot_sv_signature(signature) + ggtitle(patient$patientId)
  #save_plot(paste0("/Users/jon/hmf/gridss/vsManta/", patient$patientId, ".png"), sigPlot, base_width = 18, base_height = 6)
  
  annotatedGT1 = bind_rows(GT1, annotatedGT1)
  annotatedGT2 = bind_rows(GT2, annotatedGT2)
  annotatedMT1 = bind_rows(MT1, annotatedMT1)
  annotatedMT2 = bind_rows(MT2, annotatedMT2)
  
}


save(annotatedGT1, annotatedGT2, annotatedMT1, annotatedMT2, file = "/Users/jon/hmf/gridss/RData/annotatedVariants.RData")

