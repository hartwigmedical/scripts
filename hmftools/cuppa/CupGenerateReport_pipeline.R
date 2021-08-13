library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(stringi)
library(gtable)

# Parse and check inputs
args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) < 3)
  {
  print("Requires arguments 1=SampleId, 2=Input directory, 3=Output directory")
  stop()
}

sampleId = args[1]
inputDir = args[2]
outputDir = args[3]

cupDataFile = paste0(inputDir, sampleId, '.cup.data.csv')

if (!file.exists(cupDataFile))
  {
  print(sprintf('Missing CUP sample data file: %s', cupDataFile))
  stop()
}

cupSampleResults = read.csv(cupDataFile)

print(sprintf('sample(%s) loaded %d results', sampleId, nrow(cupSampleResults)))

# Data preparation
cupPlotData <- cupSampleResults %>% select(Category, ResultType, DataType, Value, RefCancerType, RefValue)

cupPlotData <- cupPlotData %>% mutate(RefValueLabel = sprintf('%.0f%%', RefValue * 100),
                                      DataType = stri_replace_all_fixed(DataType, '_', ' '))

## Do not show GENDER CLASSIFIER on CUP report
cupClassData <- cupPlotData %>%
  filter(Category == 'CLASSIFIER' | Category == 'COMBINED') %>%
  filter(DataType != 'GENDER') %>%
  mutate(DataLabel = DataType)

## Make sure data in SEX is not overwritten with GENDER CLASSIFIER data
cupGender <- cupPlotData %>%
  filter(DataType == 'GENDER' & Category == 'SAMPLE_TRAIT') %>%
  mutate(DataLabel = sprintf('SEX (%s)', Value),
         PrevColour = ifelse(RefValue == 0, 'high', ifelse(RefValue <= 0.02, 'low', 'norm')))

# Ensure features are shown in alphabetical order
cupFeatures <- cupPlotData %>% filter(Category == 'FEATURE' & ResultType != 'LIKELIHOOD')
featureOrder <- cupFeatures %>%
  group_by(DataType, Value) %>%
  count %>%
  arrange(DataType, Value) %>%
  ungroup()
rowIndex <- data.frame(as.numeric(as.character(rownames(featureOrder))))
colnames(rowIndex) <- "FeatureIndex"
featureOrder <- cbind(featureOrder, rowIndex)
cupFeatures <- merge(cupFeatures, featureOrder %>% select(Value, FeatureIndex), by = 'Value', all.x = T)
cupFeatures <- cupFeatures %>% mutate(DataLabel = Value)

cupOtherData <- cupPlotData %>%
  filter(Category != 'CLASSIFIER' &
           ResultType != 'LIKELIHOOD' &
           DataType != 'GENDER' &
           Category != 'FEATURE') %>%
  mutate(DataLabel = ifelse(Category == 'SNV_SIG' | Category == 'SV', sprintf('%s (%.0f)', DataType, as.numeric(as.character(Value))),
                            ifelse(DataType %in% c('PURITY', 'PLOIDY', 'MS INDELS TMB', 'CHORD HRD'), sprintf('%s (%.2f)', DataType, as.numeric(as.character(Value))),
                                   sprintf('%s (%s)', DataType, Value))),
         PercColour = ifelse(RefValue < (-2) | RefValue > 2, 'high', ifelse(RefValue < 0 | RefValue > 1, 'medium', ifelse(RefValue <= 0.02 | RefValue >= 0.98, 'low', 'norm'))),
         DataTypeOrder = ifelse(Category == 'SAMPLE_TRAIT', 0, ifelse(Category == 'SV', 1, 2)))


# Common themes for plots
font <- 'sans'
defaultFontSize <- 10

theme_set(theme_bw() + theme(axis.text = element_text(size = defaultFontSize),
                             axis.text.y = element_text(size = 10, face = 'bold', family = font),
                             panel.grid = element_blank(),
                             panel.border = element_blank(),
                             axis.text.x.bottom = element_blank(),
                             axis.ticks.x = element_blank(),
                             axis.ticks.y = element_blank(),
                             legend.position = 'none'))

gradColourMin <- 'white'
gradColourMax <- 'seagreen3'
prevColours <- c('high' = 'indianred3', 'medium' = 'salmon2', 'low' = 'peachpuff', 'norm' = 'white')

# Generate Plots
summaryPlot <- ggplot(cupClassData, aes(x = RefCancerType, y = DataLabel)) +
  geom_tile(aes(fill = RefValue), colour = "grey", stat = "identity", position = "identity") +
  geom_text(aes(label = RefValueLabel), size = 3) +
  scale_fill_gradient(low = gradColourMin, high = gradColourMax) +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(angle = 90, hjust = 0, size = 10, face = 'bold', family = font)) +
  labs(x = '', y = '', title = '')

genderPlot <- ggplot(cupGender, aes(x = RefCancerType, y = DataLabel)) +
  geom_tile(aes(fill = PrevColour), colour = "grey", stat = "identity", position = "identity") +
  geom_text(aes(label = RefValueLabel), size = 3) +
  scale_colour_manual(values = prevColours) +
  scale_fill_manual(values = prevColours, limits = names(prevColours)) +
  labs(x = '', y = '', title = '')

featurePlot <- ggplot(cupFeatures, aes(x = RefCancerType, y = reorder(DataLabel, -FeatureIndex))) +
  geom_tile(aes(fill = RefValue), colour = "grey", stat = "identity", position = "identity") +
  geom_text(aes(label = RefValueLabel), size = 3) +
  scale_fill_gradient(low = gradColourMin, high = gradColourMax) +
  labs(x = '', y = '', title = 'FEATURES')

svTraitsPlot <- ggplot(cupOtherData %>% filter(Category == 'SV' | Category == 'SAMPLE_TRAIT'), aes(x = RefCancerType, y = reorder(DataLabel, -DataTypeOrder))) +
  geom_tile(aes(fill = PercColour), colour = "grey", stat = "identity", position = "identity") +
  geom_text(aes(label = RefValueLabel), size = 3) +
  scale_colour_manual(values = prevColours) +
  scale_fill_manual(values = prevColours, limits = names(prevColours)) +
  labs(x = '', y = '', title = 'PERCENTILES')

sigPlot <- ggplot(cupOtherData %>% filter(Category == 'SNV'), aes(x = RefCancerType, y = DataLabel)) +
  geom_tile(aes(fill = PercColour), colour = "grey", stat = "identity", position = "identity") +
  geom_text(aes(label = RefValueLabel), size = 3) +
  scale_colour_manual(values = prevColours) +
  scale_fill_manual(values = prevColours, limits = names(prevColours)) +
  labs(x = '', y = '', title = 'SNV SIGNATURES')

outputFile <- paste0(outputDir, sampleId, '_cup_report.pdf')
print(paste0("writing output to file: ", outputFile))

featureLimit <- 15
featureCount <- nrow(cupFeatures %>% group_by(Value) %>% count)
titleHeight <- 10
summaryHeight <- 205
genderHeight <- 45
sigHeight <- 110
percHeight <- 90

if (featureCount > featureLimit)
  {
  separateFeaturePlot <- T
  print(sprintf('features(%d) print separately', featureCount))
  plotHeights <- c(titleHeight, summaryHeight, genderHeight, sigHeight, percHeight)
} else
  {
  separateFeaturePlot <- F
  featureHeight <- 45 + (featureCount - 1) * 10
  print(sprintf('features(%d) featureHeight(%d)', featureCount, featureHeight))
  plotHeights <- c(titleHeight, summaryHeight, genderHeight, sigHeight, percHeight, featureHeight)
}

pdf(file = outputFile, height = 14, width = 20)

par(mar = c(1, 1, 1, 1))

title <- textGrob(paste0(sampleId, ' CUP Report'), gp = gpar(fontface = "bold", fontsize = 16))

if (separateFeaturePlot)
  {
  grid.arrange(plot_grid(title, summaryPlot, genderPlot, sigPlot, svTraitsPlot,
                         ncol = 1, nrow = 5, rel_heights = plotHeights, align = 'v', axis = 'l'))

  featurePlot = featurePlot +
    scale_x_discrete(position = "top") +
    theme(axis.text.x.top = element_text(angle = 90, hjust = 0, size = 10, face = 'bold', family = font))

  grid.arrange(plot_grid(featurePlot, ncol = 1, nrow = 1), newpage = T)
} else
  {
  plot_grid(title, summaryPlot, genderPlot, sigPlot, svTraitsPlot, featurePlot,
            ncol = 1, nrow = 6, rel_heights = plotHeights, align = 'v', axis = 'l')
}

dev.off()
