library(RMySQL)
library(dplyr)
library(tidyr)

#Cosmic
cosmicCurated = read.csv("~/Documents/CosmicCurated.csv", stringsAsFactors = F)
cosmicCurated$cosmicCurated <- TRUE
cosmicCurated = cosmicCurated[, c("Genes","cosmicCurated")]
colnames(cosmicCurated) <- c("gene_name", "cosmicCurated")

cosmicCensus = read.csv("~/Documents/CosmicCensus.csv", stringsAsFactors = F)
colnames(cosmicCensus) <- c("gene_name", "cosmicOncogene", "cosmicTsg")
cosmicCensus = cosmicCensus %>% filter(cosmicOncogene | cosmicTsg)

#dNdS results
load(file="~/hmf/RData/PcawgRefCDSCv.RData")
load(file="~/hmf/RData/HmfRefCDSCv.RData")
sig = 0.05
hmfSignificant =  HmfRefCDSCv %>% filter(qglobal_cv < sig) %>% distinct(gene_name)
hmfSignificant$hmf <- TRUE
martincorenaSignificant =  PcawgRefCDSCv %>% filter(qglobal < sig) %>% distinct(gene_name)
martincorenaSignificant$martincorena <- TRUE

genePanel = merge(martincorenaSignificant, hmfSignificant, by = "gene_name", all = T)
genePanel = merge(genePanel, cosmicCurated, by = "gene_name", all = T)
genePanel = merge(genePanel, cosmicCensus, by = "gene_name", all = T)
genePanel[is.na(genePanel)] <- FALSE

save(genePanel, file="~/hmf/RData/genePanel.RData")



#### FRAGILE SITES
library(GenomicRanges)
library(dplyr)
fragileSites = read.csv(file = "~/Downloads/FragileSite.csv")
fragileSites %>% select(chrom, start, end) %>% mutate(range = GRanges(chrom, IRanges(start, end)))
fragileSites$range = GRanges(fragileSites$chrom, IRanges(fragileSites$start, fragileSites$end))

prodDB = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
canonicalTranscripts = purple::query_canonical_transcript(prodDB)
save(canonicalTranscripts, file = "~/hmf/RData/canonicalTranscripts.RData")
dbDisconnect(prodDB)
rm(prodDB)

canonicalTranscripts$range = GRanges(canonicalTranscripts$chromosome, IRanges(canonicalTranscripts$geneStart, canonicalTranscripts$geneEnd))
ol = as.matrix(findOverlaps(fragileSites$range, canonicalTranscripts$range))
fragileGenes = canonicalTranscripts[ol[,2], c("gene", "chromosome")]
fragileGenes$fragile <- TRUE
fragileGenes$chromosome <- NULL
colnames(fragileGenes) <- c("gene_name", "fragile")
save(fragileGenes, file = "~/hmf/RData/fragileGenes.RData")







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




#### TRAVIS AMPS
library(GenomicRanges)
library(dplyr)
load(file="~/hmf/RData/canonicalTranscripts.RData")
travisAmps = read.csv(file = "~/Documents/TravisAmp.csv")
travisAmps$range = GRanges(travisAmps$chrom, IRanges(travisAmps$start, travisAmps$end))

jon = canonicalTranscripts
jon$range = GRanges(jon$chromosome, IRanges(jon$geneStart, jon$geneEnd))
ol = as.matrix(findOverlaps(travisAmps$range, jon$range))

travisAmps = jon[ol[,2], c("gene", "chromosome")]
travisAmps$travisAmp <- TRUE
travisAmps$chromosome <- NULL
colnames(travisAmps) <- c("gene_name", "travisAmp")
save(travisAmps, file = "~/hmf/RData/travisAmps.RData")


#### TRAVIS DELTS
library(GenomicRanges)
library(dplyr)
travisDels = read.csv(file = "~/Documents/TravisDel.csv")
travisDels$range = GRanges(travisDels$chrom, IRanges(travisDels$start, travisDels$end))

jon = canonicalTranscripts
jon$range = GRanges(jon$chromosome, IRanges(jon$geneStart, jon$geneEnd))
ol = as.matrix(findOverlaps(travisDels$range, jon$range))

travisDels = jon[ol[,2], c("gene", "chromosome")]
travisDels$travisDel <- TRUE
travisDels$chromosome <- NULL
colnames(travisDels) <- c("gene_name", "travisDel")
save(travisDels, file = "~/hmf/RData/travisDels.RData")


#### PCAWG AMPS
pcawgAmps = read.csv("~/Documents/PCAWGDriverCodingAmplifications.csv", stringsAsFactors = F) %>% filter(gene != "near_CCND1") %>% select(gene, PCAWGCount)
pcawgAmpsTelomereIndex = gregexpr("telomere", pcawgAmps$gene) >= 0
pcawgAmpsTelomere = pcawgAmps[pcawgAmpsTelomereIndex, ]
pcawgAmpsArmIndex = grepl("p", pcawgAmps$gene) | grepl("q", pcawgAmps$gene)
pcawgAmpsArm = pcawgAmps[pcawgAmpsArmIndex & !pcawgAmpsTelomereIndex, ]
pcawgAmpsArm$chromosomeBand = substr(pcawgAmpsArm$gene, 1, ifelse(regexpr("[abcdef]", pcawgAmpsArm$gene) == -1,  nchar(pcawgAmpsArm$gene), regexpr("[abcdef]", pcawgAmpsArm$gene) - 1))
pcawgAmps = pcawgAmps[!pcawgAmpsTelomereIndex & !pcawgAmpsArmIndex, ]
save(pcawgAmps, pcawgAmpsTelomere, pcawgAmpsArm, file = "~/hmf/RData/pcawgAmps.RData")


pcawgDels = read.csv("~/Documents/PCAWGDriverCodingDeletions.csv", stringsAsFactors = F) %>% select(gene, PCAWGCount)
pcawgDelsTelomereIndex = gregexpr("telomere", pcawgDels$gene) >= 0
pcawgDelsTelomere = pcawgDels[pcawgDelsTelomereIndex, ]
pcawgDelsArmIndex = grepl("p", pcawgDels$gene) | grepl("q", pcawgDels$gene)
pcawgDelsArm = pcawgDels[pcawgDelsArmIndex & !pcawgDelsTelomereIndex, ]
pcawgDelsArm$chromosomeBand = substr(pcawgDelsArm$gene, 1, ifelse(regexpr("[abcdef_]", pcawgDelsArm$gene) == -1,  nchar(pcawgDelsArm$gene), regexpr("[abcdef_]", pcawgDelsArm$gene) - 1))
pcawgDels = pcawgDels[!pcawgDelsTelomereIndex & !pcawgDelsArmIndex, ]
save(pcawgDels, pcawgDelsTelomere, pcawgDelsArm, file = "~/hmf/RData/pcawgDels.RData")


