library(dplyr)
library(tidyr)
library(ggplot2)
library(purple)
library(ggrepel)
library(scales)
library(cowplot)
theme_set(theme_bw())

#### DRIVER TYPE CLASSIFICATION
load("~/hmf/analysis/cohort/processed/completeGenePanel.RData")
dndsGenePanel = completeGenePanel %>% filter(martincorenaDnds | hmfDnds | cosmicCurated | actionableVariant | actionableDrup)

load("~/hmf/analysis/cohort/processed/HmfRefCDSCv.RData")
genePanelCv = HmfRefCDSCv %>% filter(cancerType == "All", gene_name %in% dndsGenePanel$gene) %>% 
  mutate(
    d_ind =  ifelse(n_ind>0, n_ind * pmax(0,(wind_cv-1)/wind_cv), 0),
    d_mis =  ifelse(n_mis>0, n_mis * pmax(0,(wmis_cv-1)/wmis_cv), 0),
    d_spl =  ifelse(n_spl>0, n_spl * pmax(0,(wspl_cv-1)/wspl_cv), 0),
    d_non =  ifelse(n_non>0, n_non * pmax(0,(wnon_cv-1)/wnon_cv), 0)
    ) %>%
  select(gene = gene_name, wmis_cv, wnon_cv, wspl_cv, wnon_cv, d_ind, d_mis, d_spl, d_non, n_ind, n_mis, n_spl, n_non)
genePanelCv = left_join(genePanelCv, dndsGenePanel[, c("gene", "cosmicTsg", "cosmicOncogene", "hmfDnds", "martincorenaDnds", "actionableVariant", "actionableDrup", "actionableDrupCategory", "actionableAmplification", "actionableDeletion", "hmfAmplification", "hmfDeletion")], by = "gene")

trainTsg = genePanelCv %>% filter(cosmicTsg, !cosmicOncogene) %>% mutate(tsg = 1) %>% select(wmis_cv, wnon_cv,  tsg)
trainOnco = genePanelCv %>% filter(cosmicOncogene, !cosmicTsg) %>% mutate(tsg = 0) %>% select(wmis_cv, wnon_cv,  tsg)
trainData = rbind(trainTsg, trainOnco)
model <- glm(tsg ~.,family=binomial(link='logit'),data=trainData)

summary(model)
anova(model, test="Chisq")

genePanelCv$continuousClassification <- predict(model,newdata = genePanelCv %>% select(wmis_cv, wnon_cv),type='response')
genePanelCv[is.na(genePanelCv)] <- F
genePanelCv$inDnds <- genePanelCv$hmfDnds | genePanelCv$martincorenaDnds
genePanelCv$unambigiousCosmic = xor(genePanelCv$cosmicTsg, genePanelCv$cosmicOncogene)

genePanelCv$classification = ifelse(genePanelCv$continuousClassification > 0.5,"tsg","onco")
genePanelCv$classification = ifelse(!genePanelCv$inDnds & !genePanelCv$cosmicTsg & genePanelCv$cosmicOncogene, "onco", genePanelCv$classification)
genePanelCv$classification = ifelse(!genePanelCv$inDnds & genePanelCv$cosmicTsg & !genePanelCv$cosmicOncogene, "tsg", genePanelCv$classification)
genePanelCv$classification = ifelse(!genePanelCv$inDnds & !genePanelCv$unambigiousCosmic & (genePanelCv$hmfAmplification | genePanelCv$actionableAmplification), "onco", genePanelCv$classification)
genePanelCv$classification = ifelse(!genePanelCv$inDnds & !genePanelCv$unambigiousCosmic & (genePanelCv$hmfDeletion | genePanelCv$actionableDeletion), "tsg", genePanelCv$classification)
genePanelCv$classification = ifelse(genePanelCv$d_ind > 2 * genePanelCv$d_mis & (genePanelCv$hmfDnds | !genePanelCv$unambigiousCosmic), "tsg", genePanelCv$classification)

save(genePanelCv, file = "~/hmf/analysis/cohort/processed/genePanelCv.RData")

#genePanelCv[genePanelCv$gene == "TERT", "classification"] <- "onco"


load("~/hmf/analysis/cohort/processed/genePanelCv.RData")
load("~/hmf/analysis/cohort/processed/completeGenePanel.RData")
completeGenePanel = left_join(completeGenePanel, genePanelCv %>% select(gene, reportablePointMutation = classification), by = "gene") 

completeGenePanel = completeGenePanel %>%
  mutate(reportableAmp = reportablePointMutation == "onco" | actionableAmplification | hmfAmplification,
         reportableDel = reportablePointMutation == "tsg" | actionableDeletion | hmfDeletion)

save(completeGenePanel, file = "~/hmf/analysis/cohort/processed/completeGenePanel.RData")


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


