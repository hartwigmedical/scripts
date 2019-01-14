detach("cowplot", unload=TRUE)
library(cowplot)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)
theme_set(theme_bw())

##### ACTUAL DATA
load(file = "~/hmf/RData/Processed/hpcDndsOncoDrivers.RData")
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
  result = ifelse(hotspot, "OnHotspot", "NonHotspot")
  result = ifelse(nearHotspot, "NearHotspot", result)
  return (result)
}

hotspotColor = c("#2c7bb6","#d7191c","#e66101")
hotspotColor = setNames(hotspotColor, c("NonHotspot", "OnHotspot","NearHotspot"))

variantShape = c(23, 21, 22)
variantShape = setNames(variantShape, c("INDEL", "SNV", "MNV"))

crossAlpha = c(1, 0)
crossAlpha = setNames(crossAlpha, c("Nonsense", "Missense"))

#variantColour = c("#2c7bb6","#5e3c99","#d7191c","#1a9641" ,"#e66101")
#variantColour = setNames(variantColour, c("Missense", "NonsenseOnHotspot", "MissenseOnHotspot","Nonsense", "MissenseNearHotspot"))

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

lollipopData = hpcDndsOncoDrivers %>%
  filter(pHGVS != "") %>%
  #filter(pHGVS != "", impact != 'Frameshift') %>%
  mutate(
    hotspot = hotspot_category(hotspot, nearHotspot),
    impact = ifelse(variant == 'INDEL' & impact == 'Inframe', "Missense", impact),
    impact = ifelse(variant == 'INDEL' & impact == 'Frameshift', "Nonsense", impact),
    pos = extract_position(pHGVS)) %>%
  filter(!is.na(pos)) %>%
  group_by(gene, variant, pHGVS, pos, impact) %>%
  summarise(hotspot = first(hotspot), driverLikelihood = sum(driverLikelihoodAdjusted), n = n()) %>%
  group_by(gene) %>%
  mutate(pHGVS = ifelse( min_rank(-driverLikelihood) <= 8 & (driverLikelihood > 2 | (driverLikelihood > 1 & impact != 'Missense')) , pHGVS, ""))

lollipopCancerTypeData = hpcDndsOncoDrivers %>%
  mutate(hotspot = ifelse(hotspot | nearHotspot, "Hotspot", "NonHotspot")) %>%
  left_join(highestPurityCohort %>% select(sampleId, cancerType), by = "sampleId") %>%
  mutate(driver = driverLikelihoodAdjusted, passenger = 1 - driver) %>%
  group_by(gene, cancerType) %>% summarise(Driver = sum(driver), Passenger = sum(passenger)) %>%
  gather(type, value, 3, 4)

selectedGene = "BRAF"

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
    geom_point(aes(x = pos, y = n, alpha = impact), color = "White", size = 3, shape = 4) + 
    #geom_text_repel(data = lollipopDataGene, aes(x = pos, y = n, label = pHGVS), hjust = 0, nudge_x = nudgeDistance, direction = "both", xlim = c(0, maxY)) +
    scale_shape_manual(name = "Type", values = variantShape) +
    ylab("") + xlab("Codon") + ggtitle(paste0(selectedGene, " Variants")) +
    scale_x_continuous(limits = c(0, maxY), expand = c(0,0)) +
    scale_y_continuous(labels = c(0.5, 1,2,4,8,16,32,64,128), breaks = c(0.5, 1,2,4,8,16,32,64,128), trans="log2", limits = c(0.5, max(4, max(maxN)))) +
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(legend.position = "none", legend.text=element_text(size=8), axis.text = element_text(size = 12)) +
    scale_color_manual(name = "Hotspot", values = hotspotColor) +
    scale_fill_manual(name = "Hotspot", values = hotspotColor) + 
    scale_alpha_manual(name = "Nonsense", values = crossAlpha)

  p2 = ggplot(data = lollipopCancerTypeData %>% filter(gene == selectedGene)) +
    geom_bar(aes(x = type, y = value, fill = cancerType), stat = "identity", width = 0.8) +
    scale_fill_manual(values = cancerTypeColours) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(), legend.title = element_blank(), axis.ticks = element_blank())+
    theme(legend.position = "right", legend.text=element_text(size=8), axis.text.y = element_text(size = 12)) + 
    ggtitle("") + xlab("") + ylab("") +
    scale_y_continuous( expand = c(0,0)) +
    scale_x_discrete( expand = c(0,0)) +
    guides(fill=guide_legend(ncol=1))

  #p2_legend = g_legend(p2)
  #p2 = p2+theme(legend.position="none")
  #p_legends = plot_grid(p1_legend, p2_legend, ncol = 1)

  p_graphs = plot_grid(p1, p2, ncol = 2, nrow = 1, rel_widths = c(4, 1))
  return (p_graphs)
}

lollipop("PREX2")
lollipop("BRAF")
lollipop("NRAS")

#for (selectedGene in unique(lollipopData$gene)) {
#  pOncoLollipop = lollipop(selectedGene)
#  save_plot(paste0("~/hmf/RPlot/oncoLollipop/", selectedGene, ".png"), pOncoLollipop, base_width = 20, base_height = 6)
#}

genes = unique(lollipopData$gene)
plots = list()
for (i in c(1:length(genes))) {
  plots[[i]] = lollipop(genes[i])
}

pdf(file="~/hmf/RPlot/Extended Figure 7a -  Onco Lollipops.pdf",width=15, height = 6)
for (i in 1:length(plots)) {
  print(plots[[i]])
}
dev.off()
