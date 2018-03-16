library(RMySQL)
library(dplyr)
library(tidyr)


pilotDB = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
genePanelDB = (dbGetQuery(pilotDB, "SELECT * FROM genePanel"))
dbDisconnect(pilotDB)
rm(pilotDB)
genePanelDB = genePanelDB %>% mutate(true = TRUE) %>% spread(panel, true) %>% select(-matches("HMF"))
colnames(genePanelDB) <- c("gene_name", "cosmic", "lawrence","martincorena")
save(genePanelDB, file = "~/hmf/RData/genePanelDB.RData")


load("~/hmf/RData/cosmicCensusGenes.RData")
load("~/hmf/RData/genePanelDB.RData")
cosmic = cosmicCensusGenes %>% select(-cosmic_type) %>% filter(!gene_name %in% c(38596, 38961, 40057))
colnames(cosmic) <- c("gene_name", "cosmicOncogene", "cosmicTsg")
GenePanel = merge(genePanelDB, cosmic, by = "gene_name", all = T)
save(GenePanel, file = "~/hmf/RData/GenePanel.RData")


load("~/hmf/RData/GenePanel.RData")
load("~/hmf/RData/HmfRefCDSCv.RData")
load("~/hmf/RData/PcawgRefCDSCv.RData")

PanPcawgRefCDSCv = PcawgRefCDSCv[PcawgRefCDSCv$cancerType == 'All', c("gene_name", "qglobal")]
colnames(PanPcawgRefCDSCv) <- c("gene_name", "pcawg_qglobal_sv")

PanHmfRefCDSCvWithGenePanel = HmfRefCDSCv[HmfRefCDSCv$cancerType == 'All', ]
PanHmfRefCDSCvWithGenePanel$gene_name <- as.character(PanHmfRefCDSCvWithGenePanel$gene_name)
PanHmfRefCDSCvWithGenePanel = left_join(PanHmfRefCDSCvWithGenePanel, PanPcawgRefCDSCv, by="gene_name")
PanHmfRefCDSCvWithGenePanel = left_join(PanHmfRefCDSCvWithGenePanel, GenePanel, by="gene_name")

rm(GenePanel)
rm(PanPcawgRefCDSCv)
rm(PcawgRefCDSCv)
rm(HmfRefCDSCv)

sigLevel = 0.05
View(PanHmfRefCDSCvWithGenePanel %>% filter(qglobal_cv < sigLevel | pcawg_qglobal_sv < sigLevel | cosmic | lawrence | martincorena))
