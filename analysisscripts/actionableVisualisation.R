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

actionableVariantsPerSample = read.csv('~/hmf/resources/actionableVariantsPerSample.tsv',header=T,sep = '\t', stringsAsFactors = F) %>%
  filter(!gene %in% c('PTEN','KRAS'),
         hmfLevel %in% c('A','B')) %>%
  mutate(
    treatmentType = ifelse(treatmentType == "Unknown", "OffLabel", treatmentType),
    drug = ifelse(drug == "Fluvestrant", "Fulvestrant", drug),
    drug = ifelse(drug == "Ado-trastuzumab Emtansine", "Ado-Trastuzumab Emtansine", drug),
    drug = ifelse(drug == "AZD4547", "AZD-4547", drug),
    drug = ifelse(drug == "AZD5363", "AZD-5363", drug),
    drug = ifelse(drug == "BGJ398", "BGJ-398", drug),
    drug = alphabetical_drug(drug),
    levelTreatment = factor(paste(hmfLevel, treatmentType, sep = "_"), levelTreatmentFactors))

load(file = '~/hmf/RData/Processed/highestPurityCohortSummary.RData')
pembrolizumabVariants = highestPurityCohortSummary %>% 
  filter(msiStatus == 'MSI', cancerType != 'Other') %>% 
  select(sampleId, patientCancerType = cancerType) %>% 
  mutate(
    drug = "Pembrolizumab", 
    levelTreatment = factor("A_OnLabel", levelTreatmentFactors),
    gene = "",eventType = "MSI", pHgvs = "", hmfResponse= "Responsive")

nivolumabVariants = highestPurityCohortSummary %>% 
  filter(msiStatus == 'MSI', cancerType != 'Other') %>% 
  select(sampleId, patientCancerType = cancerType) %>% 
  mutate(
    drug = "Nivolumab", 
    levelTreatment = factor(ifelse(patientCancerType == "Colon/Rectum", "A_OnLabel", "A_OffLabel"), levelTreatmentFactors),
    gene = "",eventType = "MSI", pHgvs = "", hmfResponse= "Responsive")

actionableVariants = actionableVariantsPerSample %>% select(sampleId, patientCancerType, drug, levelTreatment, gene, eventType, pHgvs, hmfResponse) %>%
  bind_rows(pembrolizumabVariants) %>%
  bind_rows(nivolumabVariants)
  
#actionableGenes = actionableVariants %>% select(gene) %>% distinct()
#save(actionableGenes, file = "~/hmf/resources/actionableGenes.RData")

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

responsiveDrugs = drugResponse(actionableVariants, 'Responsive')
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

responsiveVariants = actionableDrugs %>% select(sampleId, patientCancerType, drug) %>%
  left_join(actionableVariants, by = c("sampleId", "patientCancerType","drug")) %>%
  select(sampleId, cancerType = patientCancerType, gene, drug, eventType, pHgvs, levelTreatment) %>%
  group_by(sampleId, cancerType, gene, drug, eventType) %>% top_n(1, levelTreatment) %>%
  group_by(sampleId, cancerType, gene, drug, eventType, pHgvs, levelTreatment) %>%
  distinct() %>%
  ungroup() %>%
  mutate(levelTreatment = factor(levelTreatment, rev(levelTreatmentFactors))) %>%
  group_by(sampleId, cancerType, gene, eventType, pHgvs, levelTreatment) %>%
  summarise(drug = paste(drug, collapse = ";")) %>%
  spread(levelTreatment, drug, fill = "") %>%
  ungroup()
save(responsiveVariants, file = "~/hmf/RData/Processed/responsiveVariants.RData")

sampleIdMap = read.csv(file = "/Users/jon/hmf/secure/SampleIdMap.csv", stringsAsFactors = F)
actionability = responsiveVariants %>% left_join(sampleIdMap, by = "sampleId") %>%
  select(-sampleId) %>%
  select(sampleId = hmfSampleId, everything())
write.csv(actionability, file = "~/hmf/RData/Actionability.csv", row.names = F) 

geneResponseSummary = responsiveVariants %>% ungroup() %>% distinct(sampleId, gene) %>% group_by(gene) %>% count()
drugResponseSummary = responsiveVariants %>% ungroup() %>% distinct(sampleId, drug) %>% group_by(drug) %>% count()
#PLATINUM v Platinum Agent ??
#Binimetinib v Binimetinib (MEK162) v Binimetinib + Ribociclib
#Alpelisib v Alpelisib + Fulvestrant
#Buparlisib v Buparlisib + Fulvestrant
#Dabrafenib v Dabrafenib + Trametinib

########################################### Visualise
load(file = '~/hmf/RData/Reference/hpcCancerTypeCounts.RData')

actionablePlotData = responsiveVariants %>% 
  mutate(
    response = ifelse(B_OnLabel != "", "B_OnLabel", "B_OffLabel"),
    response = ifelse(A_OffLabel != "", "A_OffLabel", response),
    response = ifelse(A_OnLabel != "", "A_OnLabel", response),
    response = factor(response, levelTreatmentFactors)) %>%
  group_by(sampleId, cancerType) %>% arrange(response) %>% summarise(response = last(response)) %>%
  group_by(cancerType, response) %>% count() %>% arrange(cancerType, response) %>%
  left_join(hpcCancerTypeCounts %>% select(cancerType, N), by = "cancerType" ) %>% 
  mutate(percentage = n/N) %>%
  arrange(cancerType, response)

actionablePlotDataFactors = actionablePlotData %>% 
  select(cancerType, response, n, N) %>%
  group_by(cancerType) %>% 
  spread(response, n, fill = 0) %>%
  mutate(
    APercent = (A_OffLabel + A_OnLabel) / N,
    BPercent = B_OnLabel / N) %>%
  arrange(APercent, BPercent)

actionablePlotData = actionablePlotData %>% 
  ungroup() %>%
  mutate(cancerType = factor(cancerType, actionablePlotDataFactors$cancerType))

p1 = ggplot(data = actionablePlotData, aes(x = cancerType, y = percentage)) +
  geom_bar(stat = "identity", aes(fill = response)) + 
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
