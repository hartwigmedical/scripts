#library(IRanges)
#library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purple)

#### DRIVER TYPE CLASSIFICATION
load("~/hmf/RData/processed/genePanel.RData")
genePanel = genePanel %>% filter(martincorena | hmf | cosmicCurated)

load("~/hmf/RData/processed/HmfRefCDSCv.RData")
HmfRefCDSCv = purple::dnds_excess(HmfRefCDSCv)

genePanelCv = HmfRefCDSCv %>% filter(cancerType == "All", gene_name %in% genePanel$gene_name) %>% select(gene_name, wmis_cv, wnon_cv, prob_mis, prob_non, excess_mis, excess_non)
genePanelCv = left_join(genePanelCv, genePanel[, c("gene_name", "cosmicTsg", "cosmicOncogene", "hmf", "martincorena")], by = "gene_name")

trainTsg = genePanelCv %>% filter(cosmicTsg, !cosmicOncogene) %>% mutate(tsg = 1) %>% select(wmis_cv, wnon_cv,  tsg)
trainOnco = genePanelCv %>% filter(cosmicOncogene, !cosmicTsg) %>% mutate(tsg = 0) %>% select(wmis_cv, wnon_cv,  tsg)
trainData = rbind(trainTsg, trainOnco)
model <- glm(tsg ~.,family=binomial(link='logit'),data=trainData)

summary(model)
anova(model, test="Chisq")

continuousClassification <- predict(model,newdata = genePanelCv %>% select(wmis_cv, wnon_cv),type='response')
genePanelCv[is.na(genePanelCv)] <- F
genePanelCv$classification = ifelse(continuousClassification > 0.5,"tsg","onco")
genePanelCv$classification = ifelse(!hmf & !martincorena & !genePanelCv$cosmicTsg & genePanelCv$cosmicOncogene, "onco", genePanelCv$classification)
genePanelCv$classification = ifelse(!hmf & !martincorena & genePanelCv$cosmicTsg & !genePanelCv$cosmicOncogene, "tsg", genePanelCv$classification)

genePanel = left_join(genePanel, genePanelCv[, c("gene_name","classification")], by = "gene_name")
genePanel[genePanel$gene_name == "TERT", "classification"] <- "onco"

rm(trainTsg, trainOnco, trainData, model)

tsGenes = genePanel %>% filter(classification == "tsg")
oncoGenes = genePanel %>% filter(classification == "onco")
save(tsGenes, oncoGenes, file = "~/hmf/RData/processed/driverGenes.RData")

ggplot(data=genePanelCv,aes(prob_mis,prob_non,label=gene_name))+
  geom_point(aes(colour = factor(classification)))+
  geom_text(size=2,hjust = 0, nudge_x = 0.01)

ggplot(data=genePanelCv,aes(excess_mis+0.1,excess_non+0.1,label=gene_name))+
  geom_point(aes(colour = factor(classification)))+
  geom_text(size=2,hjust = 0, nudge_x = 0.01)+
  scale_x_log10()+
  scale_y_log10() 

