library(RMySQL)
library(dndscv)
library(IRanges)
detach("package:purple", unload=TRUE); 
library(purple);
library(Biostrings)
library(GenomicRanges)
library(MASS)
library(seqinr)
library(dplyr)
library(ggplot2)

load("~/hmf/RData/HmfRefCDSCv.RData")
load("~/hmf/RData/HmfRefCDSCvNullGlobal.RData")
load("~/hmf/RData/HmfRefCDSCvNullPcawg.RData")
load("~/hmf/RData/PcawgRefCDSCv.RData")

HmfRefCDSCvNullPcawg = HmfRefCDSCvNullPcawg[HmfRefCDSCvNullPcawg$cancerType != "Neuroendocrine", ]
PcawgRefCDSCv = PcawgRefCDSCv[PcawgRefCDSCv$cancerType != "Neuroendocrine", ]

#save(HmfRefCDSCvNullPcawg, file = "~/hmf/RData/HmfRefCDSCvNullPcawg.RData")
#save(PcawgRefCDSCv, file = "~/hmf/RData/PcawgRefCDSCv.RData")

# Attach significance
HmfRefCDSCv$hmf_significant = HmfRefCDSCv$qglobal_cv < 0.05
PcawgRefCDSCv$pcawg_significant = PcawgRefCDSCv$qglobal < 0.05
HmfRefCDSCvNullPcawg$pcawg_null_significant = 
  (HmfRefCDSCvNullPcawg$qmis_cv<0.05&HmfRefCDSCvNullPcawg$null_wmis_cv>0)|
  (HmfRefCDSCvNullPcawg$qtrunc_cv<0.05&HmfRefCDSCvNullPcawg$null_wnon_cv+HmfRefCDSCvNullPcawg$null_wspl_cv>0)|
  (HmfRefCDSCvNullPcawg$qind_cv<0.05&HmfRefCDSCvNullPcawg$null_wind_cv>0)

# Merge
CombinedCv = left_join(HmfRefCDSCv, PcawgRefCDSCv[, c("gene_name", "cancerType", "pcawg_significant")], by = c("gene_name", "cancerType"))
CombinedCv = left_join(CombinedCv, HmfRefCDSCvNullPcawg[, c("gene_name", "cancerType", "pcawg_null_significant")], by = c("gene_name", "cancerType"))

SignificantCv = CombinedCv %>% filter(hmf_significant | pcawg_significant)
SignificantCv$gene_name <- as.character(SignificantCv$gene_name)
SignificantCv$status = "Unchanged"
SignificantCv$status = ifelse(!SignificantCv$pcawg_significant & SignificantCv$hmf_significant, "Added", SignificantCv$status)
SignificantCv$status = ifelse(SignificantCv$pcawg_significant & ! SignificantCv$hmf_significant, "Removed", SignificantCv$status)

ggplot(SignificantCv, aes(y = cancerType, x = gene_name, alpha = pcawg_null_significant))+
  geom_tile( aes(fill=factor(status))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_discrete(range=c(0.1,1)) + 
  scale_fill_manual(values=c("green", "red", "blue"), na.value = c("grey"))
#+ facet_wrap(~cancerType)


ggplot(SignificantCv[SignificantCv$pcawg_null_significant, ], aes(y = cancerType, x = gene_name))+
  geom_tile( aes(fill=status)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values=c("green", "red", "blue"))

