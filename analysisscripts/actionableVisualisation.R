library(dplyr)
library(tidyr)
library(scales)

load(file = '~/hmf/RData/Reference/hpcCancerTypeCounts.RData')
load(file = '~/hmf/RData/Processed/highestPurityCohortSummary.RData')
actionableVariants = read.table('~/hmf/resources/actionableVariantsPerSample.tsv',header=T,sep = '\t')
treatmentTypeFactors = rev(c("A_OnLabel","A_OffLabel","B_OnLabel","B_OffLabel"))
treatmentTypeColors = setNames(c("#fff7bc","#fee391","#fe9929","#993404"), treatmentTypeFactors)
treatmentTypeColors = setNames(c("#f0f9e8","#bae4bc","#43a2ca","#0868ac"), treatmentTypeFactors)

unique(actionableVariants$source)
#actionableVariants = actionableVariants %>% filter(source == 'oncoKb')

msiResponsive = highestPurityCohortSummary %>% 
  filter(msiStatus == 'MSI', cancerType != 'Other') %>% 
  select(sampleId, patientCancerType = cancerType) %>% 
  mutate(drug = "MSI", treatmentType = factor("A_OnLabel", treatmentTypeFactors))

responsiveDrugs = actionableVariants %>% 
  filter(
   !gene %in% c('PTEN','KRAS'),
    hmfLevel %in% c('A','B'),
    !level %in% c('Early trials'),
    hmfResponse %in% c('Responsive')) %>%
  mutate(
    treatmentType = ifelse(treatmentType == "On-label", "OnLabel", "OffLabel")) %>%
  group_by(sampleId,patientCancerType, drug, treatmentType, hmfLevel) %>%
  count() %>% 
  unite(treatmentType, hmfLevel, treatmentType) %>%
  mutate(treatmentType = factor(treatmentType, treatmentTypeFactors)) %>%
  group_by(sampleId,patientCancerType, drug) %>%
  top_n(1, treatmentType) %>% 
  select(-n) %>%
  bind_rows(msiResponsive)

responsiveDrugGenes = actionableVariants %>% 
  filter(
    !gene %in% c('PTEN','KRAS'),
    hmfLevel %in% c('A','B'),
    !level %in% c('Early trials'),
    hmfResponse %in% c('Responsive')) %>%
  mutate(
    treatmentType = ifelse(treatmentType == "On-label", "OnLabel", "OffLabel")) %>%
  group_by(sampleId,patientCancerType, drug, gene, treatmentType, hmfLevel) %>%
  count() %>% 
  unite(treatmentType, hmfLevel, treatmentType) %>%
  mutate(treatmentType = factor(treatmentType, treatmentTypeFactors)) %>%
  group_by(sampleId,patientCancerType, drug, gene) %>%
  top_n(1, treatmentType) %>% 
  select(-n) %>%
  bind_rows(msiResponsive)

resistantDrugs = actionableVariants %>% 
  filter(
    #!gene %in% c('TP53','RB1'),
    hmfLevel %in% c('A','B'),
    !level %in% c('Early trials'),
    hmfResponse %in% c('Resistant')) %>%
  mutate(
    treatmentType = ifelse(treatmentType == "On-label","OnLabel", "OffLabel")) %>%
  group_by(sampleId,patientCancerType, drug, treatmentType, hmfLevel) %>%
  count() %>% 
  unite(treatmentType, hmfLevel, treatmentType) %>%
  mutate(treatmentType = factor(treatmentType, treatmentTypeFactors)) %>%
  group_by(sampleId,patientCancerType, drug) %>%
  top_n(1, treatmentType) %>% 
  select(-n)

actionableDrugs = merge(responsiveDrugs,resistantDrugs,by=c('sampleId','patientCancerType','drug'),all=T,suffixes=c('_Response','_Resistance'), fill=0) %>%
  mutate(
    responsive = !is.na(treatmentType_Response) & is.na(treatmentType_Resistance),
    resistant = is.na(treatmentType_Response) & !is.na(treatmentType_Resistance),  
    inconsistent = !is.na(treatmentType_Response) & !is.na(treatmentType_Resistance),
    status = ifelse(responsive, "Responsive", "Resistance"),
    status = ifelse(inconsistent, "Inconsistent", status)
    ) 

actionableDrugsSummary = actionableDrugs %>%
  group_by(sampleId, patientCancerType, drug, status) %>% count() %>%
  spread(status, n) 

actionableSampleSummary = actionableDrugs %>%
  group_by(sampleId, patientCancerType, status) %>% count() %>%
  spread(status, n) %>% arrange(Inconsistent, Responsive)
  

drugResponseSummary = actionableDrugs %>% filter(responsive) %>% distinct(sampleId, drug) %>% group_by(drug) %>% count()
geneResponseSummary = inner_join(responsiveDrugGenes, actionableDrugs %>% filter(responsive) %>% select(sampleId, drug, responsive), by = c("sampleId","drug")) %>%
  ungroup() %>%
  distinct(sampleId, gene) %>%
  group_by(gene) %>%
  count()

actionablePlotData = actionableDrugs %>%
  filter(responsive) %>%
  group_by(patientCancerType, sampleId) %>%
  top_n(1, treatmentType_Response) %>%
  summarise(treatmentType_Response = dplyr::first(treatmentType_Response)) %>%
  group_by(patientCancerType, treatmentType_Response) %>% 
  count() %>%
  left_join(hpcCancerTypeCounts %>% select(patientCancerType = cancerType, N), by = "patientCancerType" ) %>%
  mutate(percentage = n/N) %>%
  filter(!is.na(percentage)) %>% 
  arrange(percentage)

actionablePlotDataFactors= actionablePlotData %>% 
  select(patientCancerType, treatmentType_Response, n, N) %>%
  group_by(patientCancerType) %>% 
  spread(treatmentType_Response, n, fill = 0) %>%
  mutate(
    APercent = (A_OffLabel + A_OnLabel) / N,
    BPercent = B_OnLabel / N) %>%
  arrange(APercent, BPercent)

actionablePlotData = actionablePlotData %>% 
  ungroup() %>%
  mutate(patientCancerType = factor(patientCancerType, actionablePlotDataFactors$patientCancerType))

ggplot(data = actionablePlotData, aes(x = patientCancerType, y = percentage)) +
  geom_bar(stat = "identity", aes(fill = treatmentType_Response)) + 
  scale_fill_manual(values = treatmentTypeColors) +
  ggtitle("Actionability") + xlab("Cancer Type") + ylab("% samples with treatment options") +
  scale_y_continuous(labels = percent, limits = c(0, 1), expand = c(0.02,0)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(axis.ticks = element_blank()) +
  coord_flip()


plot_grid(pAll, pCivic, pCgi, pOnco)


 

