detach("cowplot", unload=TRUE)
library(cowplot)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)
theme_set(theme_bw())

##### ACTUAL DATA
load(file = "~/hmf/RData/Processed/hpcDndsOncoDrivers.RData")
load(file = "~/hmf/RData/Processed/hpcDndsTsgMutations.RData")

load("~/hmf/RData/Reference/canonicalTranscripts.RData")
load(file = "~/hmf/RData/Reference/highestPurityCohort.RData")
load(file = "~/hmf/RData/Reference/cancerTypeColours.RData")

extract_position <- function(AAChange) {
  prot.spl = strsplit(x = as.character(AAChange), split = ".", fixed = TRUE)
  prot.conv = sapply(sapply(prot.spl, function(x) x[length(x)]), "[", 1)
  pos = gsub(pattern = "Ter.*", replacement = "", x = prot.conv)
  pos = gsub(pattern = "[[:alpha:]]", replacement = "", x = pos)
  pos = gsub(pattern = "\\*$", replacement = "", x = pos)
  pos = gsub(pattern = "^\\*", replacement = "", x = pos)
  pos = gsub(pattern = "\\*.*", replacement = "", x = pos)
  pos = as.numeric(sapply(X = strsplit(x = pos, split = "_", fixed = TRUE), FUN = function(x) x[1]))
  return (pos)
}

hotspot_category <- function(hotspot, nearHotspot) {
  result = ifelse(hotspot, "OnHotspot", "OffHotspot")
  result = ifelse(nearHotspot, "NearHotspot", result)
  return (result)
}

variantShape = c(23, 21, 22)
variantShape = setNames(variantShape, c("INDEL", "SNV", "MNV"))

hotspotColor = c("#de2d26","#fb6a4a","#3182bd")
hotspotColor = setNames(hotspotColor, c("OnHotspot","NearHotspot","OffHotspot"))

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

lollipopData = hpcDndsOncoDrivers %>%
  filter(pHGVS != "") %>%
  mutate(
  hotspot = hotspot_category(hotspot, nearHotspot),
  pos = extract_position(pHGVS)) %>%
  filter(!is.na(pos)) %>%
  group_by(gene, variant, pHGVS, pos, hotspot) %>%
  summarise(driverLikelihood = sum(driverLikelihoodAdjusted), n = n()) %>%
  group_by(gene) %>%
  mutate(pHGVS = ifelse( min_rank(-driverLikelihood) <= 8 & (driverLikelihood > 2 | (driverLikelihood > 1 & hotspot != 'OffHotspot')) , pHGVS, ""))

lollipopCancerTypeData = hpcDndsOncoDrivers %>%
  mutate(hotspot = ifelse(hotspot | nearHotspot, "Hotspot", "NonHotspot")) %>%
  left_join(highestPurityCohort %>% select(sampleId, cancerType), by = "sampleId") %>%
  group_by(gene, cancerType, hotspot) %>% count() %>%
  spread(hotspot, n, fill = 0) %>%
  gather(hotspot, n, 3,4)

lollipop <- function(selectedGene) {
  lollipopDataGene = lollipopData %>% filter(gene == selectedGene)
  totalCodons = canonicalTranscripts %>% filter(gene == selectedGene) %>% mutate(codons = codingBases / 3) %>% pull(codons)
  nudgeDistance = totalCodons / 100 + 5
  maxN = max(lollipopDataGene$n)

  maxY = totalCodons + 10

  p1 = ggplot(data = lollipopDataGene) +
    geom_segment(aes(x = pos, xend = pos, y = 0, yend = n, color = hotspot), linetype = "dotted", size = 0.5) +
    geom_segment(aes(x = pos, xend = pos, y = 0, yend = driverLikelihood, color = hotspot), size = 0.55) +
    geom_point(aes(x = pos, y = n, color = hotspot, fill = hotspot, shape = variant), size = 4, alpha = 1) +
    geom_text_repel(data = lollipopDataGene, aes(x = pos, y = n, label = pHGVS), hjust = 0, nudge_x = nudgeDistance, direction = "both", xlim = c(0, maxY)) +
    scale_shape_manual(name = "Type", values = variantShape) +
    scale_color_manual(name = "Location", values = hotspotColor)  + scale_fill_manual(name = "Location", values = hotspotColor) +
    ylab("") + xlab("Codon") + ggtitle(paste0(selectedGene, " Variants")) +
    scale_x_continuous(limits = c(0, maxY), expand = c(0,0)) +
    scale_y_continuous(breaks = c(0.5, 1,2,4,8,16,32,64,128), trans="log2", limits = c(0.5, max(4, max(maxN)))) +
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(legend.position = "bottom")

  p1_legend = g_legend(p1)
  p1 = p1+theme(legend.position="none")

  p2 = ggplot(data = lollipopCancerTypeData %>% filter(gene == selectedGene)) +
    geom_bar(aes(x = hotspot, y = n, fill = cancerType), stat = "identity", width = 0.8) +
    scale_fill_manual(values = cancerTypeColours) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(), legend.title = element_blank(), axis.ticks = element_blank())+
    ggtitle("") + xlab("") + ylab("") +
    scale_y_continuous( expand = c(0,0)) +
    scale_x_discrete( expand = c(0,0))

  p3 = plot_grid(p1, p2, p1_legend, nrow = 2, ncol = 2, rel_widths = c(4, 1), rel_heights = c(7,1))

  return (p3)
}




selectedGene = "BRAF"
pOncoLollipop = lollipop(selectedGene)
pOncoLollipop

for (selectedGene in unique(lollipopData$gene)) {
  pOncoLollipop = lollipop(selectedGene)
  save_plot(paste0("~/hmf/RPlot/oncoLollipop/", selectedGene, ".png"), pOncoLollipop, base_width = 20, base_height = 6)
}

genes = unique(lollipopData$gene)
plots = list()
for (i in c(1:length(genes))) {
  plots[[i]] = lollipop(genes[i])
}

length(plots)

for (i in seq(1, length(plots), 5)) {
  j = min(length(plots), i + 4)
  cat(i, ":", j, "\n")

  myPlotList = plots[i:j]
  myPlot = plot_grid(plotlist = myPlotList, ncol = 1)
  save_plot(paste0("~/hmf/RPlot/oncoLollipop/OncoLollipop", i, ".png"), myPlot, base_width = 20, base_height = 6 * (j-i + 1), limitsize = FALSE)
}




########################################################################### TSG

tsgData = hpcDndsTsgMutations %>%
  ungroup() %>%
  filter(pHGVS != '') %>%
  mutate(
  pos = extract_position(pHGVS),
  status = ifelse(biallelic, "Biallelic", NA),
  status = ifelse(is.na(status) & driverType == "MultiHit", "MultiHit", status),
  status = ifelse(is.na(status) & driverType == "SingleHit", "SingleHit", status),
  status = ifelse(is.na(status) & driverType == "Redundant", "Redundant", status),
  type = ifelse(type == 'SNP', 'SNV', type),
  type = ifelse(type == 'MNP', 'MNV', type)
  ) %>%
  group_by(gene, type, pHGVS, pos, status) %>%
  summarise(driverLikelihood = n())

selectedGene = "TP53"
tsgDataGene = tsgData %>% filter(gene == selectedGene)

ggplot(data = tsgDataGene) +
  geom_segment(aes(x = pos, xend = pos, y = 0, yend = driverLikelihood, color = status), size = 0.5) +
  geom_point(aes(x = pos, y = driverLikelihood, color = status, fill = status, shape = type), size = 4, alpha = 1) +
  scale_shape_manual(values = variantShape) + #scale_color_manual(values = hotspotColor)  + scale_fill_manual(values = hotspotColor) 
  ylab("Variants") + xlab("Codon") + ggtitle(selectedGene) +
  scale_x_continuous(limits = c(0, max(tsgDataGene$pos) + 10), expand = c(0,0)) +
  scale_y_continuous(breaks = c(1,2,4,8,16,32,64,128), trans="log2", limits = c(NA, max(1.1, max(tsgDataGene$driverLikelihood))))+
#scale_y_continuous(trans=log2_trans()) +
#scale_y_sqrt() +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), legend.position = "bottom") #, panel.border = element_blank())

tsg_lollipop <- function(selectedGene) {
  tsgDataGene = tsgData %>% filter(gene == selectedGene)

  ggplot(data = tsgDataGene) +
    geom_segment(aes(x = pos, xend = pos, y = 0, yend = driverLikelihood, color = status), size = 0.5) +
    geom_point(aes(x = pos, y = driverLikelihood, color = status, fill = status, shape = type), size = 4, alpha = 1) +
    scale_shape_manual(values = variantShape) + #scale_color_manual(values = hotspotColor)  + scale_fill_manual(values = hotspotColor) 
    ylab("Variants") + xlab("Codon") + ggtitle(selectedGene) +
    scale_x_continuous(limits = c(0, max(tsgDataGene$pos) + 10), expand = c(0,0)) +
    scale_y_continuous(breaks = c(1,2,4,8,16,32,64,128), trans="log2", limits = c(NA, max(1.1, max(tsgDataGene$driverLikelihood))))+
  #scale_y_continuous(trans=log2_trans()) +
  #scale_y_sqrt() +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), legend.position = "bottom") #, panel.border = element_blank())
}

tsg_lollipop("APC")

pdf(file='~/hmf/lollipopTSG.pdf', onefile=T, width = 20, height = 7)
for (selected in unique(tsgData$gene)) {
  print(tsg_lollipop(selected))
}
dev.off()
