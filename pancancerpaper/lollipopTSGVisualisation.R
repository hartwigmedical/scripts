detach("cowplot", unload=TRUE)
library(cowplot)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)
theme_set(theme_bw())

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

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


variantShape = c(23, 21, 22)
variantShape = setNames(variantShape, c("INDEL", "SNV", "MNV"))

variantColour = c("#2c7bb6","#5e3c99","#d7191c","#1a9641" ,"#e66101")
variantColour = setNames(variantColour, c("Missense", "NonsenseOnHotspot", "MissenseOnHotspot","Nonsense", "MissenseNearHotspot"))

load(file = "~/hmf/RData/Reference/highestPurityCohort.RData")
load("~/hmf/RData/Reference/canonicalTranscripts.RData")
load(file = "~/hmf/RData/Processed/hpcDndsTsgDrivers.RData")
load(file = "~/hmf/RData/Processed/hpcDndsTsgMutations.RData")
load(file = "~/hmf/RData/Reference/cancerTypeColours.RData")

inactivationMechanism = c("SingleHit", "MultiHit", "Mixed", "Biallelic")
variantAlpha = c(0.1,0.3,0.7,1)
variantAlpha = setNames(variantAlpha, inactivationMechanism)

tsgDriverLikelihoods = hpcDndsTsgDrivers %>% select(sampleId, gene, driverLikelihoodAdjusted)
tsgDriverCombinedData = hpcDndsTsgMutations %>% 
  filter(pHGVS != "", redundant == F, impact != 'Splice') %>%
  select(sampleId, gene, pHGVS, biallelic, hotspot, type, inactivation = driverType, impact, canonicalCodingEffect) %>%
  mutate(
    impact = ifelse(type == "INDEL" & canonicalCodingEffect == "NONSENSE_OR_FRAMESHIFT", "Nonsense", impact),
    impact = ifelse(type == "INDEL" & canonicalCodingEffect == "MISSENSE", "Missense", impact)) %>%
  left_join(tsgDriverLikelihoods, by = c("sampleId", "gene"))

tsgDrivers = tsgDriverCombinedData %>%
  mutate(
    pos = extract_position(pHGVS),
    type = ifelse(type == 'SNP', 'SNV', type),
    type = ifelse(type == 'MNP', 'MNV', type),
    impact = paste0(impact, ifelse(hotspot, "OnHotspot", ""))
  ) %>%
  group_by(gene, type, pHGVS, pos, impact) %>%
  summarise(
    driverLikelihood = sum(driverLikelihoodAdjusted), 
    inactivation = ifelse(length(unique(inactivation)) == 1, unique(inactivation), "Mixed"),
    n = n()) %>%
  filter(!is.na(pos)) %>%
  group_by(gene) %>%
  mutate(
    pHGVS = ifelse( min_rank(-driverLikelihood) <= 8 & driverLikelihood > 2, pHGVS, ""),
    inactivation = factor(inactivation, levels = inactivationMechanism, ordered = T))


tsgDriversCancerData = tsgDriverCombinedData %>%
  left_join(highestPurityCohort %>% select(sampleId, cancerType), by = "sampleId") %>%
  mutate(driver = driverLikelihoodAdjusted, passenger = 1 - driver) %>%
  group_by(gene, cancerType) %>% summarise(Driver = sum(driver), Passenger = sum(passenger)) %>%
  gather(type, value, 3, 4)


tsg_lollipop <- function(selectedGene) {
  tsgDataGene = tsgDrivers %>% filter(gene == selectedGene)
  
  totalCodons = canonicalTranscripts %>% filter(gene == selectedGene) %>% mutate(codons = codingBases / 3) %>% pull(codons)
  maxY = totalCodons + 10
  maxN = max(tsgDataGene$n)
  nudgeDistance = totalCodons / 100 + 5
  
  p1 = ggplot(data = tsgDataGene) +
    geom_segment(aes(x = pos, xend = pos, y = 0, yend = n, color = impact), linetype = "dotted", size = 0.5) + 
    geom_segment(aes(x = pos, xend = pos, y = 0, yend = driverLikelihood, color = impact), size = 0.55) +
    geom_point(aes(x = pos, y = n, color = impact, shape = type), size = 4, alpha = 1, fill = "White") +
    geom_point(aes(x = pos, y = n, color = impact, fill = impact, shape = type, alpha = inactivation), size = 4) +
    geom_text_repel(aes(x = pos, y = n, label = pHGVS), hjust = 0, nudge_x = nudgeDistance, direction = "both", xlim = c(0, maxY)) +
    ylab("") + xlab("Codon") + ggtitle(paste0(selectedGene, " Variants")) +
    scale_x_continuous(limits = c(0, maxY), expand = c(0,0)) +
    scale_y_continuous(labels = c(0.5, 1,2,4,8,16,32,64,128), breaks = c(0.5, 1,2,4,8,16,32,64,128), trans="log2", limits = c(0.5, max(4, max(maxN)))) +
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(legend.position = "bottom") + 
    scale_shape_manual(name = "Type", values = variantShape) + 
    scale_color_manual(name = "Impact", values = variantColour) +
    scale_fill_manual(name = "Impact", values = variantColour) +
    scale_alpha_manual(name = "Inactivation", values= variantAlpha)
  
  p1_legend = g_legend(p1)
  p1 = p1+theme(legend.position="none")

  p2 = ggplot(data = tsgDriversCancerData %>% filter(gene == selectedGene)) +
    geom_bar(aes(x = type, y = value, fill = cancerType), stat = "identity", width = 0.8) +
    scale_fill_manual(values = cancerTypeColours) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(), legend.title = element_blank(), axis.ticks = element_blank())+
    theme(legend.position = "right", legend.text=element_text(size=8), axis.text.y = element_text(size = 12)) + 
    ggtitle("") + xlab("") + ylab("") +
    scale_y_continuous( expand = c(0,0)) +
    scale_x_discrete( expand = c(0,0)) +
    guides(fill=guide_legend(ncol=1))
  
  p_graphs = plot_grid(p1, p2, p1_legend, ncol = 2, nrow = 2, rel_widths = c(4, 1), rel_heights = c(7,1))
  return (p_graphs)
  
}

selectedGene = "TP53"
tsg_lollipop("TP53")


genes = unique(tsgDrivers$gene)
plots = list()
for (i in c(1:length(genes))) {
  plots[[i]] = tsg_lollipop(genes[i])
}

for (i in seq(1, length(plots), 4)) {
  j = min(length(plots), i + 3)
  cat(i, ":", j, "\n")
  
  myPlotList = plots[i:j]
  myPlot = plot_grid(plotlist = myPlotList, ncol = 1)
  save_plot(paste0("~/hmf/RPlot/tsgLollipop/TsgLollipop", i, ".png"), myPlot, base_width = 15, base_height = 6 * (j-i + 1), limitsize = FALSE)
}

