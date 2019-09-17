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
if (!file.exists(dropbox)) {
  dropbox = "~/Dropbox/HMF Australia team folder/Structural Variant Analysis"
}
figdir=paste0(dropbox, "/figures/")
if (!exists("cohort")) load(paste0(dropbox, "/paper/cohort.RData"))
if (!exists("db")) db = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
if (!exists("sva_svs")) sva_svs = read_csv('../../sv/SVA_SVS.csv')
if (!exists("sva_links")) sva_links = read_csv('../../sv/SVA_LINKS.csv')
if (!exists("sva_clusters")) sva_clusters = read_csv('../../sv/SVA_CLUSTERS.csv')
#if (!exists("sva_segments")) sva_segments = read_tsv('../../sv/SVA_VIS_SEGMENTS.tsv')
if (!exists("sample_rename_lookup")) {
  sample_lookup_df = DBI::dbGetQuery(db, "SELECT * FROM sampleMapping")
  sample_rename_lookup = sample_lookup_df$hmfId
  names(sample_rename_lookup) = sample_lookup_df$sampleId
}


cancerTypeOrder = sort(unique(cohort$cancerType))
cancerTypeOrder = c(cancerTypeOrder[!(cancerTypeOrder %in% c("Unknown", "Other"))], "Other", "Unknown")
cancerTypeOrdinal = seq_along(cancerTypeOrder)
names(cancerTypeOrdinal) = cancerTypeOrder

figsave = function(figureName, ...) {
  ggsave(paste0(figdir, figureName, ".pdf"), ...)
  ggsave(paste0(figdir, figureName, ".png"), ...)
}
get_legend = function(p) {
  grobs <- ggplotGrob(p)$grobs
  legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}

theme_set(theme_bw() + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()))
sva_sampleids = cohort %>% filter(hpc) %>% pull(sampleId)


gridss_fig_tp_colours = c("#6baed6", "#3182bd", "#08519c")
gridss_fig_fp_colours = c("#fb6a4a", "#de2d26", "#a50f15")
