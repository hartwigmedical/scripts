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


############# Y DELETES
load(file = "hmf/RData/Processed/highestPurityCohortSummary.RData")
load(file = "~/hmf/RData/reference/hpcCopyNumbers.RData")
clinical = highestPurityCohortSummary %>% select(sampleId, cancerType)
maleCopyNumberY = hpcCopyNumbers %>% 
  filter(chromosome == 'Y') %>% 
  left_join(clinical, by = "sampleId") %>% 
  mutate(status = ifelse(copyNumber < 0.5, "deleted", "normal")) %>%
  group_by(cancerType, status) %>%
  summarise(length = sum(end - start + 1)) %>%
  group_by(cancerType) %>%
  mutate(percentage = length / sum(length))

maleCopyNumberYTotal = maleCopyNumberY %>% group_by(status) %>% summarise(length = sum(length)) %>% spread(status, length) %>% mutate(percentage = deleted / (normal + deleted))
maleCopyNumberY = maleCopyNumberY %>%
  mutate(totalPercentage = maleCopyNumberYTotal$percentage) %>%
  filter(cancerType != "Other")

maleCopyNumberYLevels = maleCopyNumberY %>% filter(status == 'deleted') %>% arrange(-percentage)
maleCopyNumberY = maleCopyNumberY %>% ungroup() %>% 
  mutate(
    cancerType = factor(cancerType, maleCopyNumberYLevels$cancerType),
    status = factor(status, c("normal", "deleted"), ordered = T))
maleCopyNumberY[maleCopyNumberY$cancerType == "CNS", "totalPercentage"] <- NA
maleCopyNumberY[maleCopyNumberY$cancerType == "NET", "totalPercentage"] <- NA

ggplot(data = maleCopyNumberY, aes(x = cancerType, y = percentage)) +
  geom_bar(aes(fill=status), stat = "identity") + 
  geom_line(aes(x = as.numeric(cancerType), y = totalPercentage), linetype = 2) +
  ggtitle("") + xlab("") + ylab("% of males with somatic Y chromosome deletion") +
  coord_flip() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.ticks = element_blank(), legend.position="none") + 
  scale_y_continuous(labels = percent, expand=c(0.01, 0.01), limits = c(0, 1)) +
  scale_fill_manual(values = c("#deebf7", "#2171b5")) +
  annotate("text", x = 18, y = maleCopyNumberYTotal$percentage, label = "Pan Cancer", size = 3) +
  annotate("text", x = 17, y = maleCopyNumberYTotal$percentage, label = sprintf(fmt='(%.1f%%)', 100*maleCopyNumberYTotal$percentage), size = 3) 


###### SVs V Indels V Snvs etc 
load(file = "~/hmf/RData/Reference/cancerTypeColours.RData")
load(file = "hmf/RData/Processed/highestPurityCohortSummary.RData")
highestPurityCohortSummary = highestPurityCohortSummary %>% mutate(TOTAL_SV = DEL + DUP + INS + INV +TRL )

p0 = ggplot(data=highestPurityCohortSummary, aes(x = TOTAL_SNV, y = TOTAL_INDEL)) +
  geom_segment(aes(x = 1e2, xend=1e6, y = 12400, yend=12380), linetype = "dashed") + annotate("text", x = 1e2, y = 15000, label = "MSI Threshold", size = 3, hjust = 0) +
  geom_segment(aes(y = 1e2, yend=1e6, x = 30950, xend=30950), linetype = "dashed") + annotate("text", x = 32000, y = 1.1e2, label = "TMB High Threshold", size = 3, hjust = 0) +
  geom_point(aes(color = cancerType)) + 
  scale_color_manual(values = cancerTypeColours) + 
  scale_x_continuous(trans="log10", limits = c(1e2, 1e6)) + 
  scale_y_continuous(trans="log10", limits = c(1e2, 1e6)) + 
  theme(legend.position = "right", legend.title = element_blank()) + guides(colour = guide_legend(ncol = 1)) +
  xlab("SNVs") + ylab("Indels") + ggtitle("")

p1 = ggplot(data=highestPurityCohortSummary, aes(x = TOTAL_SV, y = TOTAL_INDEL)) +
  geom_segment(aes(x = 1e1, xend=1e4, y = 12400, yend=12380), linetype = "dashed") + annotate("text", x = 1e3, y = 15000, label = "MSI Threshold", size = 3, hjust = 0) +
  geom_point(aes(color = cancerType)) + 
  scale_color_manual(values = cancerTypeColours) + 
  scale_x_continuous(trans="log10", limits = c(1e1, 1e4)) + 
  scale_y_continuous(trans="log10", limits = c(1e2, 1e6)) + 
  theme(legend.position = "right", legend.title = element_blank()) + guides(colour = guide_legend(ncol = 1)) +
  xlab("SVs") + ylab("Indels") + ggtitle("")

p2 = ggplot(data=highestPurityCohortSummary, aes(x = TOTAL_SNV, y = TOTAL_SV)) +
  geom_segment(aes(y = 1e1, yend=1e4, x = 30950, xend=30950), linetype = "dashed") + annotate("text", x = 33000, y = 5000, label = "TMB High Threshold", size = 3, hjust = 0) +
  geom_point(aes(color = cancerType)) + 
  scale_color_manual(values = cancerTypeColours) + 
  scale_x_continuous(trans="log10", limits = c(1e2, 1e6)) + 
  scale_y_continuous(trans="log10", limits = c(1e1, 1e4)) + 
  theme(legend.position = "right", legend.title = element_blank()) + guides(colour = guide_legend(ncol = 1)) +
  xlab("SNVs") + ylab("SVs") + ggtitle("")

plot_grid(p0, p1, p2, ncol = 2, labels = "AUTO")


####  Actionable Pie
library(ggforce) 

load(file = "~/hmf/RData/reference/simplifiedDrivers.RData")
actionableDriverColours = simplifiedDriverColours[c("Amp", "Fusion", "Indel", "Missense", "Nonsense", "Del")]
names(actionableDriverColours) <- c("Amp", "Fusion", "Indel","SNV","MSI","MNV")

load(file = "~/hmf/RData/Processed/responsiveVariants.RData")
responseChartData = responsiveVariants %>% 
  group_by(eventType) %>% 
  dplyr::count() %>% 
  ungroup() %>% 
  mutate(share = n / sum(n)) %>%
  mutate(
    eventType = ifelse(eventType == "Amplification", "Amp", eventType),
    eventType = ifelse(eventType == "INDEL", "Indel", eventType)
    )
  
responseChartData = responseChartData %>%
  mutate(end = 2 * pi * cumsum(share)/sum(share),
       start = lag(end, default = 0),
       middle = 0.5 * (start + end),
       hjust = ifelse(middle > pi, 1, 0),
       vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))

ggplot(responseChartData) + 
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1, start = start, end = end, fill = eventType)) +
  geom_text(aes(x = 1.05 * sin(middle), y = 1.05 * cos(middle), label = paste0(round(responseChartData$share*100, 1), "%"), hjust = hjust, vjust = vjust)) +
  coord_fixed() +
  scale_x_continuous(limits = c(-1.5, 1.4), name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1, 1), name = "", breaks = NULL, labels = NULL) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  scale_fill_manual(values = actionableDriverColours, name = "Driver") + 
  ggtitle("Actionable Drivers") 

