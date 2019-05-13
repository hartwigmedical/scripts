library(dplyr)
library(tidyr)
library(ggplot2)
library(purple)
library(ggrepel)
library(scales)
library(cowplot)
theme_set(theme_bw())

#### DRIVER TYPE CLASSIFICATION
load("~/hmf/analysis/genepanel/genePanel.RData")
genePanel = genePanel %>% filter(martincorena | hmf | cosmicCurated)

load("~/hmf/analysis/genepanel/HmfRefCDSCv.RData")
genePanelCv = HmfRefCDSCv %>% filter(cancerType == "All", gene_name %in% genePanel$gene_name) %>% 
  mutate(
    d_ind =  ifelse(n_ind>0, n_ind * pmax(0,(wind_cv-1)/wind_cv), 0),
    d_mis =  ifelse(n_mis>0, n_mis * pmax(0,(wmis_cv-1)/wmis_cv), 0)) %>%
  select(gene_name, wmis_cv, wnon_cv, d_ind, d_mis)
genePanelCv = left_join(genePanelCv, genePanel[, c("gene_name", "cosmicTsg", "cosmicOncogene", "hmf", "martincorena")], by = "gene_name")

trainTsg = genePanelCv %>% filter(cosmicTsg, !cosmicOncogene) %>% mutate(tsg = 1) %>% select(wmis_cv, wnon_cv,  tsg)
trainOnco = genePanelCv %>% filter(cosmicOncogene, !cosmicTsg) %>% mutate(tsg = 0) %>% select(wmis_cv, wnon_cv,  tsg)
trainData = rbind(trainTsg, trainOnco)
model <- glm(tsg ~.,family=binomial(link='logit'),data=trainData)

summary(model)
anova(model, test="Chisq")

genePanelCv$continuousClassification <- predict(model,newdata = genePanelCv %>% select(wmis_cv, wnon_cv),type='response')
genePanelCv[is.na(genePanelCv)] <- F
genePanelCv$unambigiousCosmic = xor(genePanelCv$cosmicTsg, genePanelCv$cosmicOncogene)
genePanelCv$classification = ifelse(genePanelCv$continuousClassification > 0.5,"tsg","onco")
genePanelCv$classification = ifelse(!genePanelCv$hmf & !genePanelCv$martincorena & !genePanelCv$cosmicTsg & genePanelCv$cosmicOncogene, "onco", genePanelCv$classification)
genePanelCv$classification = ifelse(!genePanelCv$hmf & !genePanelCv$martincorena & genePanelCv$cosmicTsg & !genePanelCv$cosmicOncogene, "tsg", genePanelCv$classification)
genePanelCv$classification = ifelse(genePanelCv$d_ind > 2 * genePanelCv$d_mis & (genePanelCv$hmf | !genePanelCv$unambigiousCosmic), "tsg", genePanelCv$classification)
genePanelCv[genePanelCv$gene_name == "TERT", "classification"] <- "onco"

genePanel = left_join(genePanel, genePanelCv[, c("gene_name","classification")], by = "gene_name")
genePanel[genePanel$gene_name == "TERT", "classification"] <- "onco"

rm(trainTsg, trainOnco, trainData, model)

tsGenes = genePanel %>% filter(classification == "tsg")
oncoGenes = genePanel %>% filter(classification == "onco")
save(tsGenes, oncoGenes, file = "~/hmf/analysis/genepanel/driverGenes.RData")


###### DIFFREENCES WITH OLD
load(file = '~/hmf/analysis/genepanel/priorDriverGenes.RData')
old = bind_rows(priorTsGenes, priorOncoGenes) %>% select(gene_name, classification) 
new = bind_rows(tsGenes, oncoGenes) %>% select(gene_name, classification)
compareClassification = full_join(old, new, by = "gene_name", suffix = c(".old", ".new"))
View(compareClassification %>% filter(classification.old != classification.new | is.na(classification.old) | is.na(classification.new) ))
save(compareClassification, file = "~/hmf/analysis/genepanel/compareClassification.RData")

load(file = "~/hmf/analysis/genepanel/compareClassification.RData")
old = HmfRefCDSCv %>% filter(gene_name == "RACGAP1")
new = HmfRefCDSCv

load(file = "~/hmf/RData/Processed/HmfRefCDSCv.RData")



##### PRIOR CLASSICAITION
priorOncoGenes = oncoGenes
priorTsGenes = tsGenes
save(priorOncoGenes, priorTsGenes, file = '~/hmf/analysis/genepanel/priorDriverGenes.RData')


