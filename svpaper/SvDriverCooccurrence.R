library(cowplot)
library(ggrepel)

library("scales")

theme_set(theme_bw() + theme(
axis.text = element_text(size=5), axis.title = element_text(size=7), legend.title = element_text(size=5), legend.text = element_text(size=5), legend.key.size = unit(0.2, "cm")))

reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv,
    log_breaks(base = base),
    domain = c(1e-100, Inf))
}

load(file = "~/hmf/RData/Reference/cancerTypeColours.RData")
driverGeneCooccurrence = read.table(file = "~/hmf/analysis/svPaper/sv_driver_gene_cooccurrence.csv", sep = ",", header = T)


qValueThreshold = 0.05
cooccurenceData = driverGeneCooccurrence %>%
  mutate(
    PositivelyCorrelated = Count_GT_Exp,
    QValue = FDR) %>%
  filter(QValue<qValueThreshold, CancerType != 'All') %>%
  mutate(
    PositivelyCorrelated = Count_GT_Exp,
    correlation = ifelse(PositivelyCorrelated, 0.3, -0.3),
    label = paste(Gene, Category, sep = "|"),
    nudge = ifelse(PositivelyCorrelated, 0.1, -0.1),
    hjust = ifelse(correlation > 0, 0, 1),
    facet = ifelse(QValue < 3e-06, T, F)) 


left_plot <- function(df, breaks, limits) {
  p2 = ggplot(data = df, aes(x = QValue, y = correlation)) +
    geom_segment(aes(xend = QValue, y = 0, yend = correlation, color = CancerType), size = 0.5) +
    geom_point(aes(color = CancerType), size = 4, alpha = 1 ) +
    #geom_text(aes(label = label), size = 5 * 24.5/72, hjust=df$hjust, nudge_y = df$nudge) +
    geom_text_repel(aes(label = label), size = 5 * 24.5/72, hjust=df$hjust, direction = "y", nudge_y = -0.1) +
    scale_color_manual(values = cancerTypeColours) +
    scale_x_continuous(trans=reverselog_trans(10), breaks = breaks, limits = limits, position = "top") +
    scale_y_continuous(limits = c(-2,0), breaks = c(0), expand = c(0,0)) +
    theme(legend.position = "bottom", legend.title = element_blank()) + ggtitle(" ") + ylab("") +
    theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 6)) +
    theme(axis.text = element_blank(), axis.ticks.x = element_blank()) +
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(),  panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "black", size = 1)) +
    coord_flip()
  
   return(p2)
}

right_plot <- function(df, breaks, limits) {
  p1 = ggplot(data = df, aes(x = QValue, y = correlation)) +
    geom_segment(aes(xend = QValue, y = 0, yend = correlation, color = CancerType), size = 0.5) +
    geom_point(aes(color = CancerType), size = 4, alpha = 1 ) +
    #geom_text(aes(label = label), size = 5 * 24.5/72, hjust=df$hjust, nudge_y = df$nudge) +
    geom_text_repel(aes(label = label), size = 5 * 24.5/72, hjust=df$hjust, direction = "y", nudge_y = 0.1) +
    scale_color_manual(values = cancerTypeColours) +
    scale_x_continuous(trans=reverselog_trans(10), breaks = breaks, limits = limits, position = "bottom") +
    scale_y_continuous(limits = c(0,2), breaks = c(0), expand = c(0,0)) +
    theme(legend.position = "none", legend.title = element_blank()) + ggtitle(" ") + ylab("") +
    theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 6)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(margin = margin(r = 10))) +
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(),  panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "black", size = 1)) +
    coord_flip()
  
  return (p1)
}

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

plot_side_by_side <- function(df, title) {
  
  p1 = left_plot(df, breaks = c(1e-4,1e-6, 1e-8, 1e-10,1e-12, 1e-14, 1e-16), limits = c(1e-3, 1e-17)) +
    annotate("text", y = -0.3, x = 1e-17, label = "Negative Correlation" , size =  6 * 24.5/72, hjust = 1) 
  legend = g_legend(p1)
  p1 =  p1 + theme(legend.position = "none")
  
  p2 = right_plot(df, breaks = c(1e-4,1e-6, 1e-8, 1e-10,1e-12, 1e-14, 1e-16), limits = c(1e-3, 1e-17)) +
    annotate("text", y = 0.3, x = 1e-17, label = "Positive Correlation" , size =  6 * 24.5/72, hjust = 0) 
  
  p3 = left_plot(df, breaks = c(0.05, 0.0015, 0.0025, 0.003), limits = c(0.05, 0.003)) +
    theme(legend.position = "none") +
    annotate("text", y = -0.3, x = 0.003, label = "Negative Correlation" , size =  6 * 24.5/72, hjust = 1) 
  
  p4 = right_plot(df, breaks = c(0.05, 0.0015, 0.0025, 0.003), limits = c(0.05, 0.003)) +
    annotate("text", y = 0.3, x = 0.003, label = "Positive Correlation" , size =  6 * 24.5/72, hjust = 0) 
  
  p1234 = plot_grid(p1, p2, p3, p4, nrow = 1)
  pFinal = plot_grid(p1234, legend, ncol = 1, rel_heights = c(7,1), labels = c(title, ""), label_size = 8)
  
  return(pFinal)
}
pDel =  plot_side_by_side(cooccurenceData %>% filter(grepl("DEL", Category)), "DEL")
pDup = plot_side_by_side(cooccurenceData %>% filter(grepl("DUP", Category)), "DUP")
pLine =  plot_side_by_side(cooccurenceData %>% filter(grepl("LINE", Category)), "LINE")



ggplot2::ggsave("~/hmf/analysis/svPaper/plot/DupCooccurrence.pdf", pDup, width = 189, height = 140, units = "mm", dpi = 300)
ggplot2::ggsave("~/hmf/analysis/svPaper/plot/DelCooccurrence.pdf", pDel, width = 189, height = 140, units = "mm", dpi = 300)
ggplot2::ggsave("~/hmf/analysis/svPaper/plot/LineCooccurrence.pdf", pLine, width = 189, height = 140, units = "mm", dpi = 300)

ggplot2::ggsave("~/hmf/analysis/svPaper/plot/DupCooccurrence.png", pDup, width = 189, height = 140, units = "mm", dpi = 300)
ggplot2::ggsave("~/hmf/analysis/svPaper/plot/DelCooccurrence.png", pDel, width = 189, height = 140, units = "mm", dpi = 300)
ggplot2::ggsave("~/hmf/analysis/svPaper/plot/LineCooccurrence.png", pLine, width = 189, height = 140, units = "mm", dpi = 300)

