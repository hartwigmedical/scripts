library(GenomicRanges)
library(tidyverse)
library(stringr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(scales)
library(RMySQL)

# Paper common
dropbox = "~/../Dropbox (HMF Australia)/HMF Australia team folder/Structural Variant Analysis"
figdir=paste0(dropbox, "/figures/")
if (!exists("cohort")) load(paste0(dropbox, "/paper/cohort.RData"))
#if (!exists("db")) db = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
if (!exists("sva_svs")) sva_svs = read_csv('../../sv/SVA_SVS.csv')
if (!exists("sva_links")) sva_links = read_csv('../../sv/SVA_LINKS.csv')
if (!exists("sva_clusters")) sva_clusters = read_csv('../../sv/SVA_CLUSTERS.csv')

cancerTypeOrder = sort(unique(cohort$cancerType))
cancerTypeOrder = c(cancerTypeOrder[!(cancerTypeOrder %in% c("Unknown", "Other"))], "Other", "Unknown")
cancerTypeOrdinal = seq_along(cancerTypeOrder)
names(cancerTypeOrdinal) = cancerTypeOrder

figsave = function(figureName, ...) {
  ggsave(paste0(figdir, figureName, ".pdf"), ...)
  ggsave(paste0(figdir, figureName, ".png"), ...)
}

theme_set(theme_bw())
