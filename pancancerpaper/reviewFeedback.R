library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(scales)
theme_set(theme_grey())
theme_set(theme_bw())
singleBlue = "#6baed6"


####### ASCAT V PURPLE COMPARE

load(file = "hmf/RData/Processed/highestPurityCohortSummary.RData")
ascat=read.table("~/hmf/resources/ascat_results", header = T, sep = "\t")
hmf=highestPurityCohortSummary %>% select(sampleId, purity, ploidy)

combined = inner_join(hmf, ascat, by = "sampleId", suffix = c(".hmf", ".ascat"))

purityCompare = ggplot(data=combined, aes(x = purity.hmf, y = purity.ascat)) + 
  geom_point(color = singleBlue) + 
  xlab("PURPLE") + ylab("ASCAT") + ggtitle("Purity Comparison - ASCAT Vs PURPLE") + 
  theme(panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks = element_blank()) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) + scale_x_continuous(labels = percent, limits = c(0, 1)) 
purityCompare

ploidyCompare = ggplot(data=combined, aes(x = ploidy.hmf, y = ploidy.ascat)) + 
  geom_point(color = singleBlue) + 
  xlab("PURPLE") + ylab("ASCAT") + ggtitle("Ploidy Comparison - ASCAT Vs PURPLE") + 
  theme(panel.border = element_blank(), axis.ticks = element_blank()) +
  scale_y_continuous(limits = c(1, 5)) + scale_x_continuous(limits = c(1, 5)) 
ploidyCompare

p = plot_grid(purityCompare, ploidyCompare, ncol = 2, labels = "AUTO")
p

ascatComparison = combined
save(ascatComparison, file = "~/hmf/RData/ascatComparison.RData")
save_plot("~/hmf/RPlot/Supplementary Figure - Ascat.png", p, base_width = 12, base_height = 6)





