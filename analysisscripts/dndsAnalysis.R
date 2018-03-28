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
