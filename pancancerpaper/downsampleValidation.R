library(RMySQL)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)
theme_set(theme_bw())

##########  DATA COLLECTION
dbDownsample = dbConnect(MySQL(), dbname='downsample_validation', groups="RAnalysis")
downsampleCohort = dbGetQuery(dbDownsample, "SELECT * FROM purity")
downsampleSomatics = query_somatic_variants(dbDownsample, downsampleCohort)
downsampleSVs = query_structural_variant_summary(dbDownsample, downsampleCohort)
downsampleSomaticsSummary = cohort_somatic_summary(downsampleSomatics)
dbDisconnect(dbDownsample)
rm(dbDownsample)


dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
prodCohort = dbGetQuery(dbProd, "SELECT * FROM purity") %>% filter(sampleId %in% downsampleCohort$sampleId)
prodSomatics = query_somatic_variants(dbProd, downsampleCohort)
prodSVs = query_structural_variant_summary(dbProd, downsampleCohort)
prodSomaticsSummary = cohort_somatic_summary(prodSomatics)
dbDisconnect(dbProd)
rm(dbProd)

save(downsampleCohort, downsampleSVs, downsampleSomaticsSummary, prodCohort, prodSVs, prodSomaticsSummary, file = "~/hmf/analysis/downsample/data.RData")


##########  ANALYSIS
load(file = "~/hmf/analysis/downsample/data.RData")

downsampleSVs = downsampleSVs %>% mutate(SV_TOTAL = BND + INS + INV + DEL + DUP)
prodSVs = prodSVs %>% mutate(SV_TOTAL = BND + INS + INV + DEL + DUP)

downsampleSummary = left_join(downsampleSVs, downsampleSomaticsSummary, by = "sampleId") %>%
  left_join(downsampleCohort %>% select(sampleId, purity, ploidy), by = "sampleId")
prodSummary = left_join(prodSVs, prodSomaticsSummary, by = "sampleId") %>%
  left_join(prodCohort %>% select(sampleId, purity, ploidy), by = "sampleId")

downsamplingComparison = left_join(downsampleSummary, prodSummary, by = "sampleId", suffix = c(".down",".prod"))
downsamplingComparison = downsamplingComparison %>% 
  mutate(
    SV_DIFF = SV_TOTAL.prod - SV_TOTAL.down,  SV_DIFF_PERCENT = (SV_TOTAL.down - SV_TOTAL.prod) / SV_TOTAL.prod,
    SNV_DIFF = TOTAL_SNV.prod - TOTAL_SNV.down,
    MNV_DIFF = TOTAL_MNV.prod - TOTAL_MNV.down,
    INDEL_DIFF = TOTAL_INDEL.prod - TOTAL_INDEL.down
  ) 

#save(downsamplingComparison, file = "~/hmf/analysis/downsample/downsamplingComparison.RData")


plot_changes <- function(data, multiplier, percentageOffset, title) {
  colnames(data) <- c("sampleId", "prod", "down")
  totalAverage = data %>% ungroup() %>% summarise(prod = sum(prod), down = sum(down)) %>% mutate(change = (prod - down) / prod) %>% pull(change)
  tidyData = data %>% mutate(change = (prod - down) / prod) %>% gather(variable, value, 2, 3) %>% ungroup() %>% mutate(variable = ifelse(variable == "prod", "Normal Depth","Downsampled"))
  tidyData$average = totalAverage
  
  sampleFactor = tidyData %>% filter(variable == "Normal Depth") %>% arrange(-value) %>% pull(sampleId)
  variableFactor = c("Normal Depth","Downsampled")
  variableFactorColours = setNames(c("#6baed6", "#d94701"), variableFactor)
  tidyData = tidyData %>% mutate(sampleId = factor(sampleId, sampleFactor), variable = factor(variable, variableFactor))
  
  ggplot(tidyData, aes(x = sampleId)) + 
    geom_bar(alpha = 0.8, aes(x = sampleId, y = value, fill = variable), stat = "identity", position = "dodge") +
    geom_point(aes(y = multiplier * (percentageOffset + change) )) +
    #geom_line(aes(y = multiplier * (percentageOffset + change), group = 1)) +
    geom_line(aes(y = multiplier * (percentageOffset + average)), group = 1, linetype = "dashed") + 
    xlab("Samples") + ylab("Count") + ggtitle(title) + 
    scale_y_continuous(expand = c(0,0), sec.axis = sec_axis(~./multiplier - percentageOffset, name = " % difference", labels = percent)) +
    scale_fill_manual(name = "", values = variableFactorColours) +
    theme(panel.border = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), legend.position = c(0.7,0.9), panel.grid.major.x = element_blank())
}


#data = downsamplingComparison %>% select(sampleId, prod = TOTAL_INDEL.prod, downsample = TOTAL_INDEL.down)
#multiplier = 500000
#percentageOffset = 0.15
#title = "Structural Variants"

pSV = plot_changes(downsamplingComparison %>% select(sampleId, prod = SV_TOTAL.prod, downsample = SV_TOTAL.down), 1000, 0, "Structural Variants")
pSNV = plot_changes(downsamplingComparison %>% select(sampleId, prod = TOTAL_SNV.prod, downsample = TOTAL_SNV.down), 500000, 0, "SNV")
pMNV = plot_changes(downsamplingComparison %>% select(sampleId, prod = TOTAL_MNV.prod, downsample = TOTAL_MNV.down), 1000, 0, "MNV")
pIndel = plot_changes(downsamplingComparison %>% select(sampleId, prod = TOTAL_INDEL.prod, downsample = TOTAL_INDEL.down), 500000, 0.15, "Indel")

pPurity = ggplot(downsamplingComparison, aes(x = purity.prod, y = purity.down)) + geom_point() + 
  xlab("Normal Depth") + ylab("Downsampled") + ggtitle("Purity") + 
  scale_y_continuous(labels = percent) + scale_x_continuous(labels = percent) + 
  theme(panel.border = element_blank(), axis.ticks = element_blank())


pPloidy = ggplot(downsamplingComparison, aes(x = ploidy.prod, y = ploidy.down)) + geom_point() + 
  xlab("Normal Depth") + ylab("Downsampled") + ggtitle("Ploidy") + 
  #scale_y_continuous(labels = percent) + scale_x_continuous(labels = percent) + 
  theme(panel.border = element_blank(), axis.ticks = element_blank())

p = plot_grid(pPurity, pSNV, pSV, pPloidy, pMNV, pIndel, ncol = 3, labels = "AUTO")
save_plot("~/hmf/RPlot/Extended Figure 3 - Downsampling.png", p, base_width = 10, base_height = 6)
