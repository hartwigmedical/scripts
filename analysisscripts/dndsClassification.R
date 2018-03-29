#library(IRanges)
#library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purple)

#### DRIVER TYPE CLASSIFICATION
load("~/hmf/RData/genePanel.RData")
genePanel = genePanel %>% filter(martincorena | hmf | cosmicCurated)

load("~/hmf/RData/HmfRefCDSCv.RData")
HmfRefCDSCv$probmis = ifelse(HmfRefCDSCv$n_mis>0,pmax(0,(HmfRefCDSCv$wmis_cv-1)/HmfRefCDSCv$wmis_cv),0)
HmfRefCDSCv$probnon = ifelse(HmfRefCDSCv$n_non,pmax(0,(HmfRefCDSCv$wnon_cv-1)/HmfRefCDSCv$wnon_cv),0)
HmfRefCDSCv$probspl = ifelse(HmfRefCDSCv$n_spl>0,pmax(0,(HmfRefCDSCv$wspl_cv-1)/HmfRefCDSCv$wspl_cv),0)
HmfRefCDSCv$probind = ifelse(HmfRefCDSCv$n_ind>0,pmax(0,(HmfRefCDSCv$wind_cv-1)/HmfRefCDSCv$wind_cv),0)
HmfRefCDSCv$excess_mis = HmfRefCDSCv$probmis*HmfRefCDSCv$n_mis
HmfRefCDSCv$excess_non = HmfRefCDSCv$probnon*HmfRefCDSCv$n_non
HmfRefCDSCv$excess_spl = HmfRefCDSCv$probspl*HmfRefCDSCv$n_spl
HmfRefCDSCv$excess_ind = HmfRefCDSCv$probind*HmfRefCDSCv$n_ind

HmfRefCDSCv$excess_trunc = HmfRefCDSCv$probnon*(HmfRefCDSCv$n_non+HmfRefCDSCv$n_spl)
HmfRefCDSCv$type = ifelse(HmfRefCDSCv$excess_trunc > 4 | HmfRefCDSCv$excess_trunc > 0.3 * HmfRefCDSCv$excess_mis, "tsg", "onco")
HmfRefCDSCv = left_join(HmfRefCDSCv, genePanel[, c("gene_name", "cosmicTsg", "cosmicOncogene")])

sig = 0.02
hmfSignificant =  HmfRefCDSCv %>% filter(qglobal_cv < sig, cancerType == 'All')
hmfSignificant =  HmfRefCDSCv %>% filter(gene_name %in% genePanel$gene_name, cancerType == 'All')


#Logistic Regression 3
tsg = hmfSignificant %>% filter(cosmicTsg) %>% mutate(tsg = 1) %>% select(wmis_cv, wnon_cv,  tsg)
onco = hmfSignificant %>% filter(cosmicOncogene) %>%  mutate(tsg = 0) %>%  select(wmis_cv, wnon_cv,  tsg)
train = rbind(tsg, onco)
test = hmfSignificant %>% select(wmis_cv, wnon_cv)
model <- glm(tsg ~.,family=binomial(link='logit'),data=train)
summary(model)
anova(model, test="Chisq")
fitted.results <- predict(model,newdata=test,type='response')
fitted.results <- ifelse(fitted.results > 0.5,"tsg","onco")
hmfSignificant$logistic3 <- fitted.results


#Logistic Regression 4
tsg = hmfSignificant %>% filter(cosmicTsg) %>% mutate(tsg = 1) %>% select(excess_mis, excess_trunc,  tsg)
onco = hmfSignificant %>% filter(cosmicOncogene) %>%  mutate(tsg = 0) %>%  select(excess_mis, excess_trunc,  tsg)
train = rbind(tsg, onco)
test = hmfSignificant %>% select(excess_mis, excess_trunc)
model <- glm(tsg ~.,family=binomial(link='logit'),data=train)
summary(model)
anova(model, test="Chisq")
fitted.results <- predict(model,newdata=test,type='response')
fitted.results <- ifelse(fitted.results > 0.5,"tsg","onco")
hmfSignificant$logistic4 <- fitted.results

tsGenes = hmfSignificant %>% filter(logistic3 == "tsg") %>% select(gene_name)
oncoGenes = hmfSignificant %>% filter(logistic3 == "onco") %>% select(gene_name)
save(tsGenes, oncoGenes, file = "~/hmf/RData/geneClassification.RData")


#differences
diff = hmfSignificant %>% filter(type != logistic1 | type != logistic2 | type != logistic3 | logistic1 != logistic2 | logistic1 != logistic3 | logistic2 != logistic3)
diff = hmfSignificant %>% filter(type != logistic1 | type != logistic2 | type != logistic3 | type != logistic4 | logistic1 != logistic2 | logistic1 != logistic3 | logistic1 != logistic4 | logistic2 != logistic3 | logistic2 != logistic4 | logistic3 != logistic4)
diff = hmfSignificant %>% filter(type != logistic3 | type != logistic4 | logistic3 != logistic4)
diff$unknown <- !diff$cosmicTsg & !diff$cosmicOncogene

l0a = ggplot(data=hmfSignificant %>% filter(cancerType=='All'),aes(probmis,probnon,label=gene_name))+geom_point(aes(colour = factor(type)))+geom_text(size=2,hjust = 0, nudge_x = 0.01)
l0b = (ggplot(data=hmfSignificant %>% filter(cancerType=='All'),aes(excess_mis+0.1,excess_trunc+0.1,label=gene_name))+geom_point(aes(colour = factor(type)))+geom_text(size=2,hjust = 0, nudge_x = 0.01)+scale_x_log10()  +scale_y_log10() )

l3a = ggplot(data=hmfSignificant %>% filter(cancerType=='All'),aes(probmis,probnon,label=gene_name))+geom_point(aes(colour = factor(logistic3)))+geom_text(size=2,hjust = 0, nudge_x = 0.01)
l3b = (ggplot(data=hmfSignificant %>% filter(cancerType=='All'),aes(excess_mis+0.1,excess_trunc+0.1,label=gene_name))+geom_point(aes(colour = factor(logistic3)))+geom_text(size=2,hjust = 0, nudge_x = 0.01)+scale_x_log10()  +scale_y_log10() )

l4a = ggplot(data=hmfSignificant %>% filter(cancerType=='All'),aes(probmis,probnon,label=gene_name))+geom_point(aes(colour = factor(logistic4)))+geom_text(size=2,hjust = 0, nudge_x = 0.01)
l4b = (ggplot(data=hmfSignificant %>% filter(cancerType=='All'),aes(excess_mis+0.1,excess_trunc+0.1,label=gene_name))+geom_point(aes(colour = factor(logistic4)))+geom_text(size=2,hjust = 0, nudge_x = 0.01)+scale_x_log10()  +scale_y_log10() )

multiplot(l3a, l3b, cols = 1)


