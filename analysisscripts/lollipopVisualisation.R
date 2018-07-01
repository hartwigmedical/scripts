
require(maftools)
detach("cowplot", unload=TRUE)
library(cowplot)


cancerTypeColours = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
                      "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#ff4791","#01837a",
                      "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
                      "#dea185","#a0729d","#8a392f")



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
  result = ifelse(hotspot, "OnHotspot", "OffHotspot")
  result = ifelse(nearHotspot, "NearHotspot", result)
  return (result)
}

lollipopData = hpcDndsOncoDrivers %>%
  mutate(
    hotspot = hotspot_category(hotspot, nearHotspot),
    pos = extract_position(pHGVS)) %>%
  filter(!is.na(pos)) %>%
  group_by(gene, variant, pHGVS, pos, hotspot) %>%
  summarise(driverLikelihood = n())


selectedGene = "ERBB2"
variantShape = c(23, 21, 22)
variantShape = setNames(variantShape, c("INDEL", "SNV", "MNV"))

hotspotColor = c("#de2d26","#fb6a4a","#3182bd")
hotspotColor = setNames(hotspotColor, c("OnHotspot","NearHotspot","OffHotspot"))

a1 = lollipop("PIK3CA")
a2 = lollipop("KRAS")
a3 = lollipop("BRAF")
a4 = lollipop("NRAS")
a5 = lollipop("ESR1")
a6 = lollipop("KIT")
a7 = lollipop("FGFR1")
a8 = lollipop("ERBB2")
a9 = lollipop("MUC6")

lollipop("NRAS")
genes = c("PIK3CA","KRAS","BRAF","NRAS","ESR1","KIT","FGFR1","MUC6","ERBB2")

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 

lollipop <- function(selectedGene) {
  lollipopDataGene = lollipopData %>% filter(gene == selectedGene)
  p1 = ggplot(data = lollipopDataGene) +
    geom_segment(aes(x = pos, xend = pos, y = 0, yend = driverLikelihood, color = hotspot), size = 0.5) + 
    geom_point(aes(x = pos, y = driverLikelihood, color = hotspot, fill = hotspot, shape = variant), size = 4, alpha = 1) +
    scale_shape_manual(values = variantShape) + scale_color_manual(values = hotspotColor)  + scale_fill_manual(values = hotspotColor) +
    ylab("Variants") + xlab("Codon") + ggtitle(selectedGene) + 
    scale_x_continuous(limits = c(0, max(lollipopDataGene$pos) + 10), expand = c(0,0)) +
    scale_y_continuous(breaks = c(1,2,4,8,16,32,64,128), trans="log2", limits = c(NA, max(1.1, max(lollipopDataGene$driverLikelihood))))+
    #scale_y_continuous(trans=log2_trans()) +
    #scale_y_sqrt() +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), legend.position = "bottom") #, panel.border = element_blank())
  
  p1_legend = g_legend(p1)
  p1 = p1+theme(legend.position="none")
  
  lollipopCancerTypeData = hpcDndsOncoDrivers %>%
    mutate(hotspot = ifelse(hotspot | nearHotspot, "Hotspot", "NonHotspot")) %>%
    left_join(highestPurityCohort %>% select(sampleId, cancerType), by = "sampleId") %>%
    group_by(gene, cancerType, hotspot) %>% count()
  
  p2 = ggplot(data = lollipopCancerTypeData %>% filter(gene == selectedGene)) +
    geom_bar(aes(x = hotspot, y = n, fill = cancerType), stat = "identity", width = 0.8) + 
    scale_fill_manual(values = cancerTypeColours) + 
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(), legend.title = element_blank())+
    ggtitle("") + xlab("") + ylab("Variants") +
    scale_y_continuous( expand = c(0,0)) +
    scale_x_discrete( expand = c(0,0))

  p3 = plot_grid(p1, p2, p1_legend, nrow = 2, ncol = 2, rel_widths = c(4, 1), rel_heights = c(7,1))
  
  return (p3)
}

lollipop("ACVR1")
lollipop("ERBB2")
lollipop("KRAS")
lollipop("BRAF")

genes = unique(lollipopData$gene)
jon = list()
for(selectedGene in genes) {
  jon[[selectedGene]] <- lollipop(selectedGene)
}


pdf(file='~/hmf/lollipop2.pdf', onefile=T, width = 20, height = 7) 
for(selectedGene in genes) {
  print(jon[[selectedGene]])
}
dev.off()
