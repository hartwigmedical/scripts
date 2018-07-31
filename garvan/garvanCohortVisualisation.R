detach("package:purple", unload=TRUE)
library(purple)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(scales)
theme_set(theme_grey())
theme_set(theme_bw())

outputDir = "~/garvan/RData/"
outputDir = "~/Documents/LKCGP_projects/RData/"
referenceDir = paste0(outputDir, "reference/")
processedDir = paste0(outputDir, "processed/")
plotDir = paste0(outputDir, "plots/")

load(paste0(processedDir, "highestPurityCohortSummary.RData"))
load(paste0(processedDir, "hpcCancerTypeCounts.RData"))
load(paste0(processedDir, "cancerTypeColours.RData"))
load(paste0(processedDir, "cohortSNPSummary.RData"))
load(paste0(processedDir, "cohortIndelSummary.RData"))
cancerTypeFactors =  factor(hpcCancerTypeCounts$cancerType, levels = hpcCancerTypeCounts$cancerType)


#################### SETUP ####################
somaticColours = c("#a6611a","#dfc27d","#80cdc1","#018571")
somaticColours = setNames(somaticColours, c("HMF SNV","PCAWG SNV", "PCAWG MNV", "HMF MNV"))
somaticLinetypes = c("solid","dashed","dashed","solid")
somaticLinetypes = setNames(somaticLinetypes, c("HMF SNV","PCAWG SNV", "PCAWG MNV", "HMF MNV"))

indelSVColours = c("#d01c8b","#f1b6da","#b8e186","#4dac26")
indelSVColours = setNames(indelSVColours, c("HMF INDEL","PCAWG INDEL", "PCAWG SV", "HMF SV"))
indelSVLinetypes = c("solid","dashed","dashed","solid")
indelSVLinetypes = setNames(indelSVLinetypes, c("HMF INDEL","PCAWG INDEL", "PCAWG SV", "HMF SV"))

singleSubstitutionColours = c("#14B0EF","#060809","#E00714","#BFBEBF","#90CA4B","#E9BBB8")
singleSubstitutionColours = setNames(singleSubstitutionColours, c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))

doubleSubstitutions = c("CC>TT","CC>AA","CC>NN","TC>NN","TT>NN","AC>NN","GC>NN","TG>NN","CT>NN","TA>NN","CG>NN","AT>NN","3+Substitutions")
doubleSubstitutionColours = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928", "#060809")
doubleSubstitutionColours = setNames(doubleSubstitutionColours, doubleSubstitutions)

indelColours = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")
indelColours = setNames(indelColours, c("INS-repeat","INS-other","DEL-repeat", "DEL-other", "DEL-MH"))

svTypes = c("DUP","DEL","BND","INS","INV")
svColours = c("#33a02c","#e31a1c","#1f78b4","#ffff33","#060809","#984ea3")
svColours = setNames(svColours, svTypes)

#################### PREPARE DATA ####################
agePlotData = highestPurityCohortSummary %>%
  select(sampleId, ageAtSample, cancerType) %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  arrange(cancerType, -ageAtSample)

cancerTypeData = highestPurityCohortSummary %>%
  group_by(cancerType) %>%
  count() %>%
  arrange(-n) %>%
  ungroup() %>%
  mutate(percentage = round(n / sum(n), 3))  %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors))

hmfMutationalLoad = highestPurityCohortSummary %>%
  select(sampleId, cancerType, ends_with("INDEL"), ends_with("SNV"), ends_with("MNV"), BND, DEL, INS, INV, DUP) %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  mutate(
    INDEL = TOTAL_INDEL,
    SNV = TOTAL_SNV,
    SV = BND + DEL + INS + INV + DUP)

combinedMutationalLoad =  hmfMutationalLoad %>% select(sampleId, cancerType, INDEL, SNV, SV) %>%
  mutate(source = "HMF") %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors))

combinedMutationalLoad = combinedMutationalLoad %>%
  group_by(cancerType) %>%
  mutate(
    medianSNV = median(SNV, na.rm = T),
    medianINDEL = median(INDEL, na.rm = T),
    medianSV = median(SV, na.rm = T)
  ) %>% ungroup()


hpcSNP = cohortSNPSummary %>%
  filter(sampleId %in% highestPurityCohortSummary$sampleId) %>%
  filter(nchar(ref) == 1, nchar(alt) == 1) %>%
  mutate(non_standard_type = paste(ref, alt, sep = '>')) %>%
  mutate(type = standard_mutation(non_standard_type)) %>%
  ungroup() %>%
  group_by(sampleId, type) %>%
  summarise(n = sum(n)) %>%
  left_join(highestPurityCohortSummary %>% select(sampleId, cancerType), by = "sampleId") %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  group_by(sampleId) %>%
  mutate(sampleMutationalLoad = sum(n), sampleRelativeN = n / sampleMutationalLoad) %>%
  ungroup() %>%
  arrange(sampleMutationalLoad)
hpcSNP$sampleId = factor(hpcSNP$sampleId, levels = unique(hpcSNP$sampleId))

hpcINDEL = cohortIndelSummary %>%
  filter(sampleId %in% highestPurityCohortSummary$sampleId) %>%
  group_by(sampleId, type = category) %>%
  summarise(n = sum(n)) %>%
  left_join(highestPurityCohortSummary %>% select(sampleId, cancerType), by = "sampleId") %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  group_by(sampleId) %>%
  mutate(sampleMutationalLoad = sum(n), sampleRelativeN = n / sampleMutationalLoad) %>%
  ungroup() %>%
  arrange(sampleMutationalLoad)
hpcINDEL$sampleId = factor(hpcINDEL$sampleId, levels = unique(hpcINDEL$sampleId))

hpcSV = highestPurityCohortSummary %>% select(sampleId, cancerType, BND, DEL, DUP, INS, INV) %>%
  gather(type, n, BND, DEL, DUP, INS, INV) %>%
  group_by(sampleId, cancerType, type) %>%
  summarise(n = sum(n)) %>%
  group_by(sampleId) %>%
  mutate(sampleMutationalLoad = sum(n), sampleRelativeN = n / sampleMutationalLoad) %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  ungroup() %>%
  arrange(sampleMutationalLoad)
hpcSV$sampleId = factor(hpcSV$sampleId, levels = unique(hpcSV$sampleId))

#################### Cancer Type Summary FACETED ####################
p1 = ggplot(data=cancerTypeData, aes(x = NA, y = n)) +
  geom_bar(aes(fill = cancerType), stat = "identity") +
  scale_fill_manual(values=cancerTypeColours, guide=FALSE) +
  geom_text(aes(label=paste0("(", percentage*100, "%)")), vjust=-0.5, size = 2) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text.x = element_text(size = 9, face = "bold")) +
  ylab("Samples") +
  facet_grid(~cancerType)

p2 = ggplot(agePlotData, aes(NA, ageAtSample)) +
  geom_violin(aes(fill=cancerType), draw_quantiles = c(0.25, 0.5, 0.75), scale = "area") +
  scale_fill_manual(values=cancerTypeColours, guide=FALSE) +
  scale_colour_manual(values=cancerTypeColours, guide=FALSE) +
  ylab("Age") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank(), strip.text = element_blank()) +
  facet_grid(~cancerType)

p3 = ggplot(data=combinedMutationalLoad) +
  stat_ecdf(aes(SNV,color='HMF SNV',linetype='HMF SNV'), geom = "step", pad = FALSE) + geom_segment(aes(x = medianSNV, xend = medianSNV, y = 0.25, yend = 0.75, color='HMF SNV'), show.legend = F) +
  scale_x_log10() + facet_grid(~cancerType) +
  scale_colour_manual(name = "Combined", values=somaticColours) +
  scale_linetype_manual(name = "Combined", values = somaticLinetypes) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(), strip.text = element_blank(), legend.position="top", legend.title = element_blank()) +
  xlab("Somatic Variants") +
  coord_flip()

p4 = ggplot(data=combinedMutationalLoad) +
  stat_ecdf(aes(INDEL, color='HMF INDEL', linetype = 'HMF INDEL'),geom = "step", pad = FALSE) + geom_segment(aes(x = medianINDEL, xend = medianINDEL, y = 0.25, yend = 0.75, color='HMF INDEL'), show.legend = F) +
  stat_ecdf(aes(SV,color='HMF SV',linetype='HMF SV'),geom = "step", pad = FALSE) + geom_segment(aes(x = medianSV, xend = medianSV, y = 0.25, yend = 0.75, color='HMF SV'), show.legend = F) +
  scale_x_log10() + facet_grid(~cancerType) +
  scale_colour_manual(name = "Combined", values=indelSVColours) +
  scale_linetype_manual(name = "Combined", values = indelSVLinetypes) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(), strip.text = element_blank(), legend.position="bottom", legend.title = element_blank()) +
  xlab("Somatic Variants") +
  xlab("INDELs & SVs") +
  coord_flip()

p5 = ggplot(data=hpcSNP, aes(x = sampleId, y = sampleRelativeN)) +
  geom_bar(aes(fill = type), stat = "identity", width=1) + ylab("") +
  scale_fill_manual(values=singleSubstitutionColours) +
  theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #strip.background = element_blank(),
    strip.text.x = element_text(size = 9),
    legend.position="bottom", legend.title = element_blank()) +
  facet_grid(~cancerType, scales = "free_x") +
  guides(fill = guide_legend(nrow = 1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1))


p7 = ggplot(data=hpcINDEL, aes(x = sampleId, y = sampleRelativeN)) +
  geom_bar(aes(fill = type), stat = "identity", width=1) + ylab("") +
  scale_fill_manual(values=indelColours) +
  theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), strip.text.x =element_blank(),
    legend.position="bottom", legend.title = element_blank()) +
  facet_grid(~cancerType, scales = "free_x") +
  guides(fill = guide_legend(nrow = 1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1))

p8 = ggplot(data=hpcSV, aes(x = sampleId, y = sampleRelativeN)) +
  geom_bar(aes(fill = type), stat = "identity", width=1) + ylab("") +
  scale_fill_manual(values=svColours) +
  theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), strip.text.x =element_blank(),
    legend.position="bottom", legend.title = element_blank()) +
  facet_grid(~cancerType, scales = "free_x") +
  guides(fill = guide_legend(nrow = 1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1))

pFigure1 = plot_grid(p1, p2, p3, p4, p5, p7, p8, ncol=1, align="v", rel_heights = c(1, 1, 3, 3, 2, 2, 2), labels = c("A", "B", "C", "D", "E", "F", "G"))
save_plot( paste0(plotDir, "Figure 1 - Overview.png"), pFigure1, base_width = 14, base_height = 20)


####################################
### Purity PLOT @@@@@@@@
purityData = highestPurityCohortSummary %>%
  select(sampleId, purity,cancerType ) %>%
  arrange(cancerType, -purity)

pPurity = ggplot(data=purityData)+
  stat_ecdf(aes(purity,color='Purity'),geom = "step", pad = FALSE) +
  coord_flip() +
  xlab("Purity") + ylab("Samples")+ 
  theme(axis.title.x =  element_blank()) +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(labels = percent) + theme(legend.position = "none")

save_plot(paste0(plotDir, "Cohort Purity.png"), pPurity, base_width = 5, base_height = 5)
