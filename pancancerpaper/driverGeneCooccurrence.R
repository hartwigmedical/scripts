library(data.table)
library(IRanges)
library(dplyr)
library(tidyr)
library(stringi)
library(devtools)


# Correlation between genes using a probability of co-occurence

getSampleGeneList<-function(dgData)
{
  # create a dataframe with a row for each sampleId and a string list of the genes involved
  sgList = data.frame(matrix(ncol = 3, nrow = 0))
  sgList = setNames(sgList, c("SampleId", "GeneList", "GeneCount"))

  curGeneList = ""
  curSample = ""
  curGeneCount = 0

  for(i in 1:nrow(dgData))
  {
    sgRow = dgData[i,]

    if(curSample != sgRow$sampleId)
    {
      if(nrow(sgList) > 0)
      {
        rowIndex = nrow(sgList)
        sgList[rowIndex,2] = paste(curGeneList, ",", sep='') # append a comma for the grep search below
        sgList[rowIndex,3] = curGeneCount
      }

      curSample = sgRow$sampleId
      curGeneList = paste(",", sgRow$gene, sep='')
      rowIndex = nrow(sgList)+1
      sgList[rowIndex,1] = curSample
      sgList[rowIndex,2] = curGeneList
      sgList[rowIndex,3] = 1
      curGeneCount = 1
    }
    else
    {
      curGeneList = paste(curGeneList, ",", sgRow$gene, sep='')
      curGeneCount = curGeneCount + 1
    }
  }

  return (sgList)
}

calc_gene_pairing_probs<-function(gsCounts, sgList, sampleCount, logCalcs = T)
{
  ggPPResults = data.frame(matrix(ncol = 9, nrow = 0))
  ggPPResults = setNames(ggPPResults, c("Gene1", "Gene2", "Gene1SC", "Gene2SC", "BothGenesSC", "BothGenesExpected", "Fisher", "GeneChr1", "GeneChr2"))

  geneCount = n_distinct(gsCounts$gene)

  for(i in 1:nrow(gsCounts))
  {
    geneRow = gsCounts[i,]
    gene1 = geneRow$gene
    scWithGene1 = geneRow$SampleCount
    gene1SamplesPerc = round(scWithGene1/sampleCount,4)

    # if(i > 4) # temp..
    #   break

    # now search each gene and filter on each signature about a significant threshold
    for(j in i+1:nrow(gsCounts))
    {
      if(j > nrow(gsCounts))
        break

      gene2Row = gsCounts[j,]
      gene2 = gene2Row$gene
      scWithGene2 = gene2Row$SampleCount
      gene2SamplesPerc = round(scWithGene2/sampleCount,4)

      # how many samples have both genes
      grep1 = paste(',', gene1, ',', sep='')
      grep2 = paste(',', gene2, ',', sep='')
      scWithGene1And2 = nrow(sgList %>% filter(grepl(grep1, GeneList)&grepl(grep2, GeneList)))

      # expected count of gene 2 in samples with gene 1
      ggExpectedCount = round(sampleCount *gene1SamplesPerc * gene2SamplesPerc,4)

      scWithGene1NoGene2 = scWithGene1 - scWithGene1And2
      scNoGene1WithGene2 = scWithGene2 - scWithGene1And2
      scNoGene2 = sampleCount - scWithGene2
      scNoGene1NoGene2 = scNoGene2 - scWithGene1NoGene2

      if(scWithGene1And2 < 0 | scNoGene1WithGene2 < 0 | scWithGene1NoGene2 < 0 | scNoGene1NoGene2 < 0)
      {
        stop(paste("Invalid inputs: gene1=", gene1, ", gene2=", gene2, ", sampleCount=", sampleCount, ", withGene1=", scWithGene1, ", withGene2=", scWithGene2,
                   ", withGene1AndGene2=", scWithGene1And2, ", noGene1WithGene2=", scNoGene1WithGene2, ", withGene1NoGene2=", scWithGene1NoGene2, ", noGene1NoGene2=", scNoGene1NoGene2, sep=''))
      }

      fishMatrix = rbind(c(scWithGene1And2,scNoGene1WithGene2), c(scWithGene1NoGene2,scNoGene1NoGene2))

      if(scWithGene1And2 < ggExpectedCount)
        ggFisherProb = fisher.test(fishMatrix, alternative="less")$p.value
      else
        ggFisherProb = fisher.test(fishMatrix, alternative="greater")$p.value

      # inclusionThreshold = 0.001
      inclusionThreshold = 1.01

      if(ggFisherProb < inclusionThreshold)
      {
        if(logCalcs)
        {
          print(paste("gene1=", gene1, ", gene2=", gene2, ", sampleCount=", sampleCount, ", withGene1=", scWithGene1, ", withGene2=", scWithGene2,
                      ", withGene1AndGene2=", scWithGene1And2, ", noGene1WithGene2=", scNoGene1WithGene2, ", withGene1NoGene2=", scWithGene1NoGene2, ", noGene1NoGene2=", scNoGene1NoGene2,
                      ", fetProb=", round(ggFisherProb,6), sep=''))
        }

        rowIndex = nrow(ggPPResults)+1
        ggPPResults[rowIndex,1] = gene1
        ggPPResults[rowIndex,2] = gene2
        ggPPResults[rowIndex,3] = scWithGene1
        ggPPResults[rowIndex,4] = scWithGene2
        ggPPResults[rowIndex,5] = scWithGene1And2
        ggPPResults[rowIndex,6] = ggExpectedCount
        ggPPResults[rowIndex,7] = ggFisherProb
        ggPPResults[rowIndex,8] = geneRow$Chromosome
        ggPPResults[rowIndex,9] = gene2Row$Chromosome
      }
    }
  }

  if(nrow(ggPPResults) > 0)
  {
    ggPPResults$PositivelyCorrelated = ggPPResults$BothGenesSC > ggPPResults$BothGenesExpected
    ggPPResults = ggPPResults %>% arrange(Fisher)

    # set ranking values
    rowIndex = data.frame(as.numeric(as.character(rownames(ggPPResults))))
    colnames(rowIndex) <- c("Rank")
    ggPPResults = cbind(ggPPResults, rowIndex)

    ggPPResults$TestCount = geneCount*geneCount/2
    ggPPResults$QValue = ggPPResults$Fisher*ggPPResults$TestCount/ggPPResults$Rank
  }

  return (ggPPResults)
}

calc_gene_cooccurence<-function(cancerTypesList, driverGeneList, logCalcs = T)
{
  ggAllProbs = data.frame(matrix(ncol = 16, nrow = 0))

  i = 1
  for(cancerTypeStr in cancerTypesList)
  {
    if(cancerTypeStr == "All")
      dgData = driverGeneList
    else
      dgData = driverGeneList %>% filter(cancerType==cancerTypeStr)

    sampleGeneList = getSampleGeneList(dgData)

    sgList = dgData %>% group_by(gene) %>% summarise(SampleCount=n(), Chromosome=first(chromosome))
    sampleCount = n_distinct(dgData$sampleId)

    print(paste(i, ": cancer=", cancerTypeStr, ", sampleCount=", sampleCount, ", geneRecords=", nrow(dgData), sep=''))

    ggProbs = calc_gene_pairing_probs(sgList, sampleGeneList, sampleCount, logCalcs)

    if(nrow(ggProbs) > 0)
    {
      ggProbs$CancerType = cancerTypeStr
      ggProbs$GeneCount = n_distinct(sgList)
      ggProbs$SampleCount = sampleCount

      ggAllProbs = rbind(ggAllProbs, ggProbs)
    }

    i = i + 1
  }

  return (ggAllProbs)
}

View(ggProbs)

# load("~/data/driversByGene.RData")
# View(driversByGene)

rm(driversByGene)
load("~/data/hpcDriversByGene_20180628.RData")
load("~/hmf/RData/processed/hpcDriversByGene.RData")

View(hpcDriversByGene)
driversByGene = hpcDriversByGene

View(driversByGene %>% group_by(type) %>% summarise(Count=n()))
View(driversByGene %>% group_by(gene) %>% summarise(Count=n()))

cancerTypes = driversByGene %>% group_by(cancerType) %>% summarise(Count=n())
cancerTypesList = cancerTypes$cancerType
View(cancerTypes)


# only consider TSGs and Oncogones of high confidence
driversNonFusions = driversByGene %>% filter(type!='FUSION'&driverLikelihood>=0.5) %>% arrange(sampleId)
# nrow(driversNonFusions)

# run the co-occurrence logic for each cancer type
allGenePairProbs = calc_gene_cooccurence(cancerTypes$cancerType, driversNonFusions, F)
save(allGenePairProbs, file = "~/hmf/RData/Processed/allGenePairProbs.RData")

nrow(allGenePairProbs)
View(allGenePairProbs)
write.csv(allGenePairProbs, "~/logs/r_output/genePairCo-occurence_20180702_all.csv", row.names=F, quote=F)

View(allGenePairProbs)

View(allGenePairProbs %>% filter(GeneChr1!=GeneChr2) %>% select(CancerType,everything()))

# adjust column names as required
colnames(allGenePairProbs) = c()

# set column names
# CancerType   (eg. Colon/Rectum)
# Driver1
# Driver2
# SamplesInCohort
# SamplesWithDriver1Only
# SamplesWithDriver2Only
# SamplesWithBothDrivers
# PositivelyCorrelated  (T/F)
# FisherExactTestPValue
# SignificanceRankInCohort
# QValue


# analysis of the results

# filter for significant results
pValueThreshold = 1e-5
qValueThreshold = 0.05
nrow(ggAllProbs %>% filter(PositivelyCorrelated==T&Fisher<pValueThreshold&GeneChr1!=GeneChr2))
nrow(ggAllProbs %>% filter(PositivelyCorrelated==F&Fisher<pValueThreshold))


# are there any positively correctly oncogenes and TSGs?
nrow(ggAllProbs2 %>% filter(PositivelyCorrelated==T&Fisher<0.01&Gene1Type!=Gene2Type))
View(ggAllProbs2 %>% filter(PositivelyCorrelated==T&Fisher<0.001&Gene1Type!=Gene2Type))


################################## PLOT RESULTS ##################################
library(cowplot)
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
load(file = "~/hmf/RData/Processed/allGenePairProbs.RData")
#old = allGenePairProbs  %>% filter(QValue<qValueThreshold, GeneChr1!=GeneChr2)
#new = allGenePairProbs %>% filter(QValue<qValueThreshold, GeneChr1!=GeneChr2)
#save(new, old, file = "~/hmf/RData/Processed/driverCoorcurrenceDiff.RData")

qValueThreshold = 0.05
cooccurenceData = allGenePairProbs %>% 
  filter(QValue<qValueThreshold, GeneChr1!=GeneChr2) %>%
  mutate(
    correlation = ifelse(PositivelyCorrelated, 0.3, -0.3),
    label = paste(Gene1, Gene2, sep = "|"),
    nudge = ifelse(PositivelyCorrelated, 0.1, -0.1),
    hjust = ifelse(correlation > 0, 0, 1),
    facet = ifelse(QValue < 3e-06, T, F),
    jon = log10(1 / QValue)
    )

p1 = ggplot(data = cooccurenceData, aes(x = QValue, y = correlation)) +
  geom_segment(aes(xend = QValue, y = 0, yend = correlation, color = CancerType), size = 0.5) +
  geom_point(aes(color = CancerType), size = 4, alpha = 1 ) +
  geom_text(aes(label = label), size = 5 * 24.5/72, hjust=cooccurenceData$hjust, nudge_y = cooccurenceData$nudge) +
  scale_color_manual(values = cancerTypeColours) +
  scale_x_continuous(trans=reverselog_trans(10), breaks = c(1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-11, 1e-12), limits = c(1e-1,1e-6), position = "bottom") +
  scale_y_continuous(limits = c(0,2), breaks = c(0), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank()) + ggtitle(" ") + ylab("") +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 6)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(margin = margin(r = 10))) +
  theme(panel.border = element_blank(), panel.grid.minor = element_blank(),  panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "black", size = 1)) +
  #theme(panel.background = element_rect(fill = "#f7fcf5"))+ panel.grid.major.y = element_blank(),
  coord_flip()

p2 = ggplot(data = cooccurenceData, aes(x = QValue, y = correlation)) +
  geom_segment(aes(xend = QValue, y = 0, yend = correlation, color = CancerType), size = 0.5) +
  geom_point(aes(color = CancerType), size = 4, alpha = 1 ) +
  geom_text(aes(label = label), size = 5 * 24.5/72, hjust=cooccurenceData$hjust, nudge_y = cooccurenceData$nudge) +
  scale_color_manual(values = cancerTypeColours) +
  scale_x_continuous(trans=reverselog_trans(10), breaks = c(1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-11, 1e-12), limits = c(1e-1,1e-6), position = "top") +
  scale_y_continuous(limits = c(-2,0), breaks = c(0), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank()) + ggtitle(" ") + ylab("") +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 6)) +
  theme(axis.text = element_blank(), axis.ticks.x = element_blank()) +
  theme(panel.border = element_blank(), panel.grid.minor = element_blank(),  panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "black", size = 1)) +
  #theme(panel.background = element_rect(fill = "#fff5eb"))+
  coord_flip()


g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend <- g_legend(p2)
p2 = p2 + theme(legend.position = "none")

p3 = ggplot(data = cooccurenceData, aes(x = QValue, y = correlation)) +
  geom_segment(aes(xend = QValue, y = 0, yend = correlation, color = CancerType), size = 0.5) +
  geom_point(aes(color = CancerType), size = 4, alpha = 1 ) +
  geom_text(aes(label = label), size = 5 * 24.5/72, hjust=cooccurenceData$hjust, nudge_y = cooccurenceData$nudge) +
  annotate("text", y = -0.3, x = 1.5e-12, label = "Negative Correlation" , size = 6 * 24.5/72, hjust = 1) +
  scale_color_manual(values = cancerTypeColours) +
  scale_x_continuous(trans=reverselog_trans(10), breaks = c(1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-10,1e-11,1e-12, 1e-13), limits = c(1e-11,1e-12), position = "top") +
  scale_y_continuous(limits = c(-1.8,0), breaks = c(-3,0,3), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank()) + ggtitle(" ") + ylab("") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text = element_blank(), axis.ticks.x = element_blank()) +
  theme(panel.border = element_blank(), panel.grid.minor = element_blank(),  panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "black", size = 1)) +
  coord_flip()


p4 = ggplot(data = cooccurenceData, aes(x = QValue, y = correlation)) +
  geom_segment(aes(xend = QValue, y = 0, yend = correlation, color = CancerType), size = 0.5) +
  geom_point(aes(color = CancerType), size = 4, alpha = 1 ) +
  geom_text(aes(label = label), size = 5 * 24.5/725, hjust=cooccurenceData$hjust, nudge_y = cooccurenceData$nudge) +
  annotate("text", y = 0.3, x = 1.5e-12, label = "Postive Correlation" , size =  6 * 24.5/72, hjust = 0) +
  scale_color_manual(values = cancerTypeColours) +
  scale_x_continuous(trans=reverselog_trans(10), breaks = c(1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-10,1e-11,1e-12, 1e-13), limits = c(1e-11,1e-12)) +
  scale_y_continuous(limits = c(0,2), breaks = c(-3,0,3), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank()) + ggtitle(" ") + ylab("") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(margin = margin(r = 10))) +
  theme(panel.border = element_blank(), panel.grid.minor = element_blank(),  panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "black", size = 1)) +
  coord_flip()

pMain = plot_grid(p3, p4, p2, p1, nrow = 2, ncol = 2, rel_heights = c(1, 6))
pCooccurence = plot_grid(pMain, legend, ncol = 1, rel_heights = c(7, 1), labels = c("a"), label_size = 8)
pCooccurence
ggplot2::ggsave("~/hmf/RPlot/Extended Figure 9.pdf", pCooccurence, width = 89, height = 180, units = "mm", dpi = 300)
ggplot2::ggsave("~/hmf/RPlot/Extended Figure 9.png", pCooccurence, width = 89, height = 180, units = "mm", dpi = 300)
ggplot2::ggsave("~/hmf/RPlot/Extended Figure 9.eps", pCooccurence, width = 89, height = 180, units = "mm", dpi = 300)
