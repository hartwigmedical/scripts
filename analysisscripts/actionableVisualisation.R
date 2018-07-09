library(dplyr)
library(tidyr)
library(scales)

########################################### Prepare Data
alphabetical_drug <- function(drugs) {
  result = c()
  for (i in 1:length(drugs)) {
    drug = drugs[i]
    sortedDrugList = sort(unlist(strsplit(drug, " + ", fixed = T)))
    result[i] <- paste(sortedDrugList, collapse = " + ")
  }
  return (result)
}

levelTreatmentFactors = rev(c("A_OnLabel","A_OffLabel","B_OnLabel","B_OffLabel"))
levelTreatmentColors = setNames(rev(c("#2171b5","#6baed6","#bdd7e7","#eff3ff")), levelTreatmentFactors)

load(file = '~/hmf/RData/Processed/highestPurityCohortSummary.RData')
actionableVariants = read.table('~/hmf/resources/actionableVariantsPerSample.tsv',header=T,sep = '\t', stringsAsFactors = F) %>%
  filter(!startsWith(sampleEvent, "Synonymous"), 
         !gene %in% c('PTEN','KRAS'),
         !level %in% c('Early trials'),
         hmfLevel %in% c('A','B')) %>%
  mutate(
    drug = ifelse(drug == "Fluvestrant", "Fulvestrant", drug),
    drug = ifelse(drug == "Ado-trastuzumab Emtansine", "Ado-Trastuzumab Emtansine", drug),
    drug = ifelse(drug == "AZD4547", "AZD-4547", drug),
    drug = ifelse(drug == "AZD5363", "AZD-5363", drug),
    drug = ifelse(drug == "BGJ398", "BGJ-398", drug),
    drug = alphabetical_drug(drug),
    treatmentType = ifelse(treatmentType == "On-label", "OnLabel", "OffLabel"),
    levelTreatment = factor(paste(hmfLevel, treatmentType, sep = "_"), levelTreatmentFactors))

########################################### Supplementary Data
drugResponse <- function(actionableVariants, response) {
  actionableVariants %>% 
    filter(hmfResponse == response) %>%
    group_by(sampleId,patientCancerType, drug, levelTreatment) %>%
    count() %>% 
    group_by(sampleId,patientCancerType, drug) %>%
    top_n(1, levelTreatment) %>% 
    select(-n) 
}

msiResponsive = highestPurityCohortSummary %>% 
  filter(msiStatus == 'MSI', cancerType != 'Other') %>% 
  select(sampleId, patientCancerType = cancerType) %>% 
  mutate(drug = "MSI", levelTreatment = factor("A_OnLabel", levelTreatmentFactors))

responsiveDrugs = drugResponse(actionableVariants, 'Responsive') %>% bind_rows(msiResponsive)
resistantDrugs = drugResponse(actionableVariants, 'Resistant')

actionableDrugs = merge(responsiveDrugs,resistantDrugs,by=c('sampleId','patientCancerType','drug'),all=T,suffixes=c('_Response','_Resistance'), fill=0) %>%
  mutate(
    responsive = !is.na(levelTreatment_Response) & is.na(levelTreatment_Resistance),
    resistant = is.na(levelTreatment_Response) & !is.na(levelTreatment_Resistance),  
    inconsistent = !is.na(levelTreatment_Response) & !is.na(levelTreatment_Resistance),
    status = ifelse(responsive, "Responsive", "Resistance"),
    status = ifelse(inconsistent, "Inconsistent", status)
  ) %>%
  filter(responsive)

responsiveVariants = actionableDrugs %>% select(sampleId, drug) %>%
  left_join(actionableVariants, by = c("sampleId","drug")) %>%
  select(sampleId, gene, drug, eventType, sampleEvent, levelTreatment, hmfLevel, treatmentType) %>%
  group_by(sampleId, gene, drug, eventType, sampleEvent, levelTreatment, hmfLevel, treatmentType) %>%
  distinct() %>%
  group_by(sampleId, gene, drug, eventType, sampleEvent) %>%
  top_n(1, levelTreatment) %>%
  select(-levelTreatment) %>%
  group_by(sampleId, gene, eventType, sampleEvent) %>%
  mutate(variantTreatmentOptions = n()) %>%
  ungroup()
save(responsiveVariants, file = "~/hmf/RData/Processed/responsiveVariants.RData")


geneResponseSummary = responsiveVariants %>% ungroup() %>% distinct(sampleId, gene) %>% group_by(gene) %>% count()
drugResponseSummary = responsiveVariants %>% ungroup() %>% distinct(sampleId, drug) %>% group_by(drug) %>% count()
#PLATINUM v Platinum Agent ??
#Binimetinib v Binimetinib (MEK162) v Binimetinib + Ribociclib
#Alpelisib v Alpelisib + Fulvestrant
#Buparlisib v Buparlisib + Fulvestrant
#Dabrafenib v Dabrafenib + Trametinib

########################################### Visualise
load(file = '~/hmf/RData/Reference/hpcCancerTypeCounts.RData')

actionablePlotData = actionableDrugs %>%
  group_by(patientCancerType, sampleId) %>%
  top_n(1, levelTreatment_Response) %>%
  summarise(levelTreatment_Response = dplyr::first(levelTreatment_Response)) %>%
  group_by(patientCancerType, levelTreatment_Response) %>% 
  count() %>%
  left_join(hpcCancerTypeCounts %>% select(patientCancerType = cancerType, N), by = "patientCancerType" ) %>%
  mutate(percentage = n/N) %>%
  filter(!is.na(percentage)) %>% 
  arrange(percentage)

actionablePlotDataFactors = actionablePlotData %>% 
  select(patientCancerType, levelTreatment_Response, n, N) %>%
  group_by(patientCancerType) %>% 
  spread(levelTreatment_Response, n, fill = 0) %>%
  mutate(
    APercent = (A_OffLabel + A_OnLabel) / N,
    BPercent = B_OnLabel / N) %>%
  arrange(APercent, BPercent)

actionablePlotData = actionablePlotData %>% 
  ungroup() %>%
  mutate(patientCancerType = factor(patientCancerType, actionablePlotDataFactors$patientCancerType))

p1 = ggplot(data = actionablePlotData, aes(x = patientCancerType, y = percentage)) +
  geom_bar(stat = "identity", aes(fill = levelTreatment_Response)) + 
  scale_fill_manual(values = levelTreatmentColors) +
  ggtitle("") + xlab("") + ylab("% Samples with treatment options") +
  scale_y_continuous(labels = percent, limits = c(0, 1), expand = c(0.02,0)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(axis.ticks = element_blank()) +
  coord_flip()

pActionability = plot_grid(p1, labels = "AUTO")
pActionability

save_plot("~/hmf/RPlot/Figure 8 - Actionable.png", pActionability, base_width = 6, base_height = 6)