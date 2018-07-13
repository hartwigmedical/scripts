#library(IRanges)
#library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purple)
library(ggrepel)

#### DRIVER TYPE CLASSIFICATION
load("~/hmf/RData/processed/genePanel.RData")
genePanel = genePanel %>% filter(martincorena | hmf | cosmicCurated)

load("~/hmf/RData/processed/HmfRefCDSCv.RData")
genePanelCv = HmfRefCDSCv %>% filter(cancerType == "All", gene_name %in% genePanel$gene_name) %>% select(gene_name, wmis_cv, wnon_cv)
genePanelCv = left_join(genePanelCv, genePanel[, c("gene_name", "cosmicTsg", "cosmicOncogene", "hmf", "martincorena")], by = "gene_name")

trainTsg = genePanelCv %>% filter(cosmicTsg, !cosmicOncogene) %>% mutate(tsg = 1) %>% select(wmis_cv, wnon_cv,  tsg)
trainOnco = genePanelCv %>% filter(cosmicOncogene, !cosmicTsg) %>% mutate(tsg = 0) %>% select(wmis_cv, wnon_cv,  tsg)
trainData = rbind(trainTsg, trainOnco)
model <- glm(tsg ~.,family=binomial(link='logit'),data=trainData)

summary(model)
anova(model, test="Chisq")

genePanelCv$continuousClassification <- predict(model,newdata = genePanelCv %>% select(wmis_cv, wnon_cv),type='response')
genePanelCv[is.na(genePanelCv)] <- F
genePanelCv$classification = ifelse(genePanelCv$continuousClassification > 0.5,"tsg","onco")
genePanelCv$classification = ifelse(!genePanelCv$hmf & !genePanelCv$martincorena & !genePanelCv$cosmicTsg & genePanelCv$cosmicOncogene, "onco", genePanelCv$classification)
genePanelCv$classification = ifelse(!genePanelCv$hmf & !genePanelCv$martincorena & genePanelCv$cosmicTsg & !genePanelCv$cosmicOncogene, "tsg", genePanelCv$classification)
genePanelCv[genePanel$gene_name == "TERT", "classification"] <- "onco"

genePanel = left_join(genePanel, genePanelCv[, c("gene_name","classification")], by = "gene_name")
genePanel[genePanel$gene_name == "TERT", "classification"] <- "onco"

rm(trainTsg, trainOnco, trainData, model)

tsGenes = genePanel %>% filter(classification == "tsg")
oncoGenes = genePanel %>% filter(classification == "onco")
save(tsGenes, oncoGenes, file = "~/hmf/RData/processed/driverGenes.RData")

ggplot(data=genePanelCv,aes(wmis_cv+0.1,wnon_cv+0.1,label=gene_name))+
  geom_point(aes(colour = factor(classification)))+
  geom_text_repel(size=2,hjust = 0, nudge_x = 0.01) +
  scale_x_log10()+
  scale_y_log10() 
  



singleBlue = "#6baed6"
singleOnco = "#d94701"
#hotspotColours = setNames(c("#d94701","#fd8d3c", "#fdbe85"), hotspotFactors) #E69F00

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 


genePanelCvNames = genePanelCv %>% filter(gene_name %in% c("KRAS","NRAS","PIK3CA","AR","KIT")| wnon_cv > 30 )

scatterPlot <- ggplot(genePanelCv, aes(wmis_cv+0.1, wnon_cv+0.1, color=factor(classification)))+ 
    geom_point() + 
    geom_text_repel(data = genePanelCvNames, aes(wmis_cv+0.1, wnon_cv+0.1, label = gene_name), color = "black", size=3, hjust = 1, nudge_x = 0.01,  direction = "both")+ 
    scale_x_log10()  +
    scale_y_log10()  +
    scale_color_manual(name = "Gene Type", values = c(singleOnco,singleBlue), labels = c("Oncogene", "Tumor Suppressor Gene")) +
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank()) +
    xlab("w_missense") + ylab("w_nonsense") + ggtitle("")
    #theme(legend.position=c(0,1), legend.justification=c(0,1))
scatterPlot

xdensity <- ggplot(genePanelCv, aes(wmis_cv+0.1, fill=factor(classification))) + 
    geom_density(alpha=.5, adjust = 6.5) + 
    stat_smooth(aes(y=continuousClassification), fill = singleBlue, method="glm", method.args=list(family="binomial"), se=T) +
    #scale_y_continuous(trans = "identity", sec.axis = sec_axis(~./2, name = "Relative humidity [%]")) +
    scale_y_continuous(labels = percent, limits = c(0, 1)) +
    scale_x_log10()  + ylab("TSG Likelihood") + xlab("") + ggtitle("") + 
    scale_fill_manual(values = c(singleOnco,singleBlue)) + 
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank()) +
    theme(legend.position = "none")
xdensity

ydensity <- ggplot(genePanelCv, aes(wnon_cv+0.1, fill=factor(classification))) + 
  geom_density(alpha=.5) + 
  stat_smooth(aes(y=continuousClassification), fill = singleBlue, method="glm", method.args=list(family="binomial"), se=T) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  scale_x_log10()  + ylab("TSG Likelihood") + xlab("") + ggtitle("") + 
  scale_fill_manual(values = c(singleOnco,singleBlue)) +
  theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank()) +
  theme(legend.position = "none") + coord_flip()
ydensity

myLegend = g_legend(scatterPlot)
scatterPlot = scatterPlot + theme(legend.position="none")
  
pLogistic = plot_grid(xdensity, myLegend, scatterPlot, ydensity,  ncol=2, nrow=2, rel_widths=c(4, 2), rel_heights=c(2, 4), labels = c("A","","B","C"))
pLogistic  
save_plot("~/hmf/RPlot/Methods Figure 4 - LogisticRegression.png", pLogistic, base_width = 6, base_height = 6)
