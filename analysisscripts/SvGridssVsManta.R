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
  
  GT1 = GT1 %>% mutate(GT1 = "GT1", GT2 = NA, MT1 = NA, MT2 = NA)
  GT2 = GT2 %>% mutate(GT1 = NA, GT2 = "GT2", MT1 = NA, MT2 = NA)
  MT1 = MT1 %>% mutate(GT1 = NA, GT2 = NA, MT1 = "MT1", MT2 = NA)
  MT2 = MT2 %>% mutate(GT1 = NA, GT2 = NA, MT1 = NA, MT2 = "MT2")
  
  GT1_GT2 = sv_overlaps(GT1, GT2, maxgap = 100)
  GT1[GT1_GT2$queryHits, "GT2"] <- "GT2"
  GT2[GT1_GT2$subjectHits, "GT1"] <- "GT1"
  
  GT1_MT1 = sv_overlaps(GT1, MT1, maxgap = 100)
  GT1[GT1_MT1$queryHits, "MT1"] <- "MT1"
  MT1[GT1_MT1$subjectHits, "GT1"] <- "GT1"
  
  GT1_MT2 = sv_overlaps(GT1, MT2, maxgap = 100)
  GT1[GT1_MT2$queryHits, "MT2"] <- "MT2"
  MT2[GT1_MT2$subjectHits, "GT1"] <- "GT1"
  
  GT2_MT1 = sv_overlaps(GT2, MT1, maxgap = 100)
  GT2[GT2_MT1$queryHits, "MT1"] <- "MT1"
  MT1[GT2_MT1$subjectHits, "GT2"] <- "GT2"
  
  GT2_MT2 = sv_overlaps(GT2, MT2, maxgap = 100)
  GT2[GT2_MT2$queryHits, "MT2"] <- "MT2"
  MT2[GT2_MT2$subjectHits, "GT2"] <- "GT2"
  
  MT1_MT2 = sv_overlaps(MT1, MT2, maxgap = 100)
  MT1[MT1_MT2$queryHits, "MT2"] <- "MT2"
  MT2[MT1_MT2$subjectHits, "MT1"] <- "MT1"
  
  result[["GT1"]] <- GT1
  result[["GT2"]] <- GT2
  result[["MT1"]] <- MT1
  result[["MT2"]] <- MT2
  
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

### RUN THIS TO SAVE LOTS OF LOADING TIME
#prodVariants = read.csv('/Users/jon/hmf/gridss/clustering/CLUSTER_V24.csv', header = T, stringsAsFactors = F)
#prodVariants = prodVariants %>% filter(sampleId %in% gridssCohort$sampleId)
#write.csv(prodVariants, file = '/Users/jon/hmf/gridss/clustering/CLUSTER_V24_COMPACT.csv')

prodCohortVariants = old_column_names(read.csv('/Users/jon/hmf/gridss/clustering/CLUSTER_V24_COMPACT.csv', header = T, stringsAsFactors = F)) %>% 
  filter(sampleId %in% gridssCohort$sampleId)
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

for (i in 1:1) {
  patient = cohort[i, ]
  GT1 = gridssCohortVariants %>% filter(sampleId == patient$Sample1)
  GT2 = gridssCohortVariants %>% filter(sampleId == patient$Sample2)
  MT1 = prodCohortVariants %>% filter(sampleId == patient$Sample1)
  MT2 = prodCohortVariants %>% filter(sampleId == patient$Sample2)
  
  combinedScope = combined_scope(GT1, GT2, MT1, MT2)
  signature = combined_signature(combinedScope[["GT1"]], combinedScope[["GT2"]], combinedScope[["MT1"]], combinedScope[["MT1"]])
  plot_sv_signature(signature) + ggtitle(patient$patientId)
}


