library(dplyr)
library(tidyr)
library(ggplot2)
library(purple)
load("~/hmf/RData/HmfRefCDSCv.RData")
load("~/hmf/RData/HmfRefCDSCvNullPcawg.RData")
load("~/hmf/RData/PcawgRefCDSCv.RData")
#load("~/hmf/RData/HmfRefCDSCvNullGlobal.RData")

absSignificance = 0.02
relSignificance = 0.05

HmfRefCDSCv$gene_name <- as.character(HmfRefCDSCv$gene_name)
PcawgRefCDSCv$gene_name <- as.character(PcawgRefCDSCv$gene_name)
HmfRefCDSCvNullPcawg$gene_name <- as.character(HmfRefCDSCvNullPcawg$gene_name)

HmfRefCDSCv$prob_mis = ifelse(HmfRefCDSCv$n_mis>0,pmax(0,(HmfRefCDSCv$wmis_cv-1)/HmfRefCDSCv$wmis_cv),0)
HmfRefCDSCv$prob_non = ifelse(HmfRefCDSCv$n_non,pmax(0,(HmfRefCDSCv$wnon_cv-1)/HmfRefCDSCv$wnon_cv),0)
HmfRefCDSCv$prob_spl = ifelse(HmfRefCDSCv$n_spl>0,pmax(0,(HmfRefCDSCv$wspl_cv-1)/HmfRefCDSCv$wspl_cv),0)
HmfRefCDSCv$prob_ind = ifelse(HmfRefCDSCv$n_ind>0,pmax(0,(HmfRefCDSCv$wind_cv-1)/HmfRefCDSCv$wind_cv),0)
HmfRefCDSCv$excess_mis = HmfRefCDSCv$prob_mis*HmfRefCDSCv$n_mis
HmfRefCDSCv$excess_non = HmfRefCDSCv$prob_non*HmfRefCDSCv$n_non
HmfRefCDSCv$excess_spl = HmfRefCDSCv$prob_spl*HmfRefCDSCv$n_spl
HmfRefCDSCv$excess_ind = HmfRefCDSCv$prob_ind*HmfRefCDSCv$n_ind

HmfRefCDSCv$hmf = HmfRefCDSCv$qglobal < absSignificance
PcawgRefCDSCv$martincorena = PcawgRefCDSCv$qglobal < absSignificance
HmfRefCDSCvNullPcawg$martincorena_null = 
  (HmfRefCDSCvNullPcawg$qmis_cv< relSignificance & HmfRefCDSCvNullPcawg$null_wmis_cv>0)|
  (HmfRefCDSCvNullPcawg$qtrunc_cv< relSignificance & HmfRefCDSCvNullPcawg$null_wnon_cv+HmfRefCDSCvNullPcawg$null_wspl_cv>0)|
  (HmfRefCDSCvNullPcawg$qind_cv< relSignificance & HmfRefCDSCvNullPcawg$null_wind_cv>0)

combinedCv = left_join(HmfRefCDSCv, PcawgRefCDSCv[, c("gene_name", "cancerType", "martincorena")], by = c("gene_name", "cancerType"))
combinedCv = left_join(combinedCv, HmfRefCDSCvNullPcawg[, c("gene_name", "cancerType", "martincorena_null")], by = c("gene_name", "cancerType"))

significantCv = combinedCv %>% filter(hmf | martincorena)
significantCv$status = "Unchanged"
significantCv$status = ifelse(!significantCv$martincorena & significantCv$hmf, "Added", significantCv$status)
significantCv$status = ifelse(significantCv$martincorena & ! significantCv$hmf, "Removed", significantCv$status)
rm(combinedCv)

ggplot(significantCv, aes(x = cancerType, y = gene_name, alpha = martincorena_null))+
  geom_tile( aes(fill=factor(status))) + theme(axis.text.y = element_text(size = 4), axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_discrete(range=c(0.1,1)) + 
  scale_fill_manual(values=c("green", "red", "blue"), na.value = c("grey")) + xlab("") + ylab("") +
  guides(fill=guide_legend(title=NULL))

ggplot(significantCv[!is.na(significantCv$martincorena_null) & significantCv$martincorena_null, ], aes(x = cancerType, y = gene_name))+
  geom_tile( aes(fill=factor(status))) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_discrete(range=c(0.4,1)) + 
  scale_fill_manual(values=c("green", "red", "blue"), na.value = c("grey"))+ xlab("") + ylab("") + 
  guides(fill=guide_legend(title=NULL))

##################### GENES

load("~/hmf/RData/cosmicGenes.RData")

hmfGenes =  HmfRefCDSCv %>% filter(qglobal_cv < absSignificance) %>% distinct(gene_name) %>% mutate(hmf = T)
martincorenaGenes =  PcawgRefCDSCv %>% filter(qglobal < absSignificance) %>% distinct(gene_name)  %>% mutate(martincorena = T)
genePanel = merge(hmfGenes, martincorenaGenes, by = "gene_name", all = T)
genePanel = merge(genePanel, cosmicGenes, by = "gene_name", all = T)
genePanel[is.na(genePanel)] <- FALSE
genePanel = genePanel %>% filter(hmf | martincorena | cosmicCurated)
rm(hmfGenes, martincorenaGenes, cosmicGenes)



##################### dNdS Classification
genePanelCv = HmfRefCDSCv %>% filter(cancerType == 'All', gene_name %in% genePanel$gene_name) %>% select(gene_name, cancerType, wmis_cv, wnon_cv, prob_mis, prob_non, excess_mis, excess_non)
genePanelCv = left_join(genePanelCv, genePanel[, c("gene_name", "cosmicTsg", "cosmicOncogene")], by = "gene_name")

trainTsg = genePanelCv %>% filter(cosmicTsg, !cosmicOncogene) %>% mutate(tsg = 1) %>% select(wmis_cv, wnon_cv,  tsg)
trainOnco = genePanelCv %>% filter(!cosmicTsg, cosmicOncogene) %>% mutate(tsg = 0) %>% select(wmis_cv, wnon_cv,  tsg)

trainData = rbind(trainTsg, trainOnco)
model <- glm(tsg ~.,family=binomial(link='logit'),data=trainData)

summary(model)
anova(model, test="Chisq")

continuousClassification <- predict(model,newdata= genePanelCv %>% select(wmis_cv, wnon_cv),type='response')
genePanelCv$classification = ifelse(continuousClassification > 0.5,"tsg","onco")
genePanel = left_join(genePanel, genePanelCv[, c("gene_name","classification")], by = "gene_name")

rm(trainTsg, trainOnco, trainData, model)

ggplot(data=genePanelCv,aes(prob_mis,prob_non,label=gene_name))+
  geom_point(aes(colour = factor(classification)))+
  geom_text(size=2,hjust = 0, nudge_x = 0.01)

ggplot(data=genePanelCv,aes(excess_mis+0.1,excess_non+0.1,label=gene_name))+
  geom_point(aes(colour = factor(classification)))+
  geom_text(size=2,hjust = 0, nudge_x = 0.01)+
  scale_x_log10()+
  scale_y_log10() 


load("~/hmf/RData/cancerTypeColours.RData")

oncoMutationsPerGene = HmfRefCDSCv %>% filter(cancerType != 'All', gene_name %in% oncoGenes$gene_name) %>% select(gene = gene_name, N = n_mis, cancerType)
oncoMutationsPerGeneLevels = oncoMutationsPerGene %>% group_by(gene) %>% summarise(N = sum(N)) %>% arrange(-N)
oncoMutationsPerGene$gene <- factor(oncoMutationsPerGene$gene, levels = oncoMutationsPerGeneLevels$gene)
ggplot(data=oncoMutationsPerGene, aes(gene, N)) +
  geom_bar(aes(fill = cancerType), stat = "identity") + 
  ggtitle("Oncogene mutations per cancer type") + xlab("Genes") + ylab("Number of variants") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual( values= cancerTypeColours)



oncoMutationsPerGene = oncoSomatics %>% filter(cancerType != 'All', gene_name %in% oncoGenes$gene_name) %>% select(gene = gene_name, N = n_mis, cancerType)
oncoMutationsPerGeneLevels = oncoMutationsPerGene %>% group_by(gene) %>% summarise(N = sum(N)) %>% arrange(-N)
oncoMutationsPerGene$gene <- factor(oncoMutationsPerGene$gene, levels = oncoMutationsPerGeneLevels$gene)
ggplot(data=oncoMutationsPerGene, aes(gene, N)) +
  geom_bar(aes(fill = cancerType), stat = "identity") + 
  ggtitle("Oncogene mutations per cancer type") + xlab("Genes") + ylab("Number of variants") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual( values= cancerTypeColours)




plot_mutations_by_cancer_type( HmfRefCDSCv %>% filter(cancerType != 'All', gene_name %in% tsGenes$gene_name))


