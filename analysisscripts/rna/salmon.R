library(ggplot2)
library(RMySQL)
library(random)

#### AHR Data
ahrData = read.table(file = "~/hmf/analysis/salmon/ahr/ahr.tsv", sep = "\t", header = T, stringsAsFactors = F) %>%
  mutate(
    AHR = ifelse(Amplification == "High Amp", "AMP", "NONE"), 
    AHR = ifelse(!is.na(Hotspot) & Hotspot, "HOTSPOT", AHR),
    AHR = ifelse(!is.na(HasExonDeletion) & HasExonDeletion, "DEL", AHR)) %>%
  filter(AHR != "NONE") %>%
  select(sampleId, AHR)
names(ahrData) <- c("Sample", "AHR")

#### GENE Panel
load(file = "~/hmf/analysis/cohort/reference/canonicalTranscripts.RData")
load(file = "~/hmf/analysis/cohort/processed/genePanel.RData")
genePanel = genePanel %>% filter(!is.na(reportablePointMutation) | !is.na(reportableAmp) | !is.na(reportableDel)) %>% 
  select(gene) %>%
  bind_rows(data.frame(gene = c("AHR"), stringsAsFactors = F)) %>%
  left_join(canonicalTranscripts %>% select(gene, geneId), by = "gene") 
names(genePanel) <- c("GeneName", "GeneId")

### Raw ISOFOX Data
fileNames <- Sys.glob("/Users/jon/hmf/analysis/salmon/isofox/*.isf.transcript_data.csv")
rawIsofoxData = data.frame(stringsAsFactors = F)
for (fileName in fileNames) {
  sample <- gsub("/Users/jon/hmf/analysis/salmon/isofox/", "", fileName)
  sample <- gsub(".isf.transcript_data.csv", "", sample)
  cat ("Processing", sample, "\n")
  data = read.table(fileName, sep = ",", header = T, stringsAsFactors = F) %>% 
    mutate(NumReads = FitAllocation) %>%
    select(GeneId, GeneName, TransId = TransName, EffectiveLength, NumReads) %>%
    mutate(normalisationFactor =  sum(ifelse(EffectiveLength > 0, NumReads / EffectiveLength, 0))) %>%
    mutate(TPM = NumReads / EffectiveLength / normalisationFactor * 1000000) %>%
    filter(GeneName %in% genePanel$GeneName)
  data$Sample <- sample
  rawIsofoxData = bind_rows(rawIsofoxData, data)
}

#### Raw Salmon Data
ensemblTransExonData = read.csv('~/hmf/analysis/RNA/ensembl_trans_exon_data.csv', stringsAsFactors = F)
ensemblTransExonData = ensemblTransExonData %>% select(GeneId, TransId = Trans) %>% distinct()

fileNames <- Sys.glob("/Users/jon/hmf/analysis/salmon/data/*.sf")
rawSalmonData = data.frame(stringsAsFactors = F)
for (fileName in fileNames) {
  sample <- gsub("/Users/jon/hmf/analysis/salmon/data/", "", fileName)
  sample <- gsub("-quant.sf", "", sample)
  cat ("Processing", sample, "\n")
  data = read.table(fileName, sep = "\t", header = T, stringsAsFactors = F) %>%
    mutate(TransId = substr(Name, 1, 15)) %>% select(-Name) %>%
    left_join(ensemblTransExonData, by = "TransId") %>%
    filter(GeneId %in% genePanel$GeneId) %>%
    left_join(genePanel, by = "GeneId")
  data$Sample <- sample
  rawSalmonData = bind_rows(rawSalmonData, data)
}

#fileName = "/Users/jon/hmf/analysis/salmon/data/CPCT02010944T-quant.sf"
#ENST00000284881


#### ADD PATIENT DATA
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis", host = "127.0.0.1")
purityTable = dbGetQuery(dbProd, "SELECT * from purity")
clinicalTable = dbGetQuery(dbProd, "SELECT * from clinical")
dbDisconnect(dbProd)
rm(dbProd)

patientData = purityTable %>% filter(sampleId %in% rawIsofoxData$Sample | sampleId %in% rawSalmonData$Sample) %>%
  left_join(clinicalTable, by = "sampleId") %>%
  select(sampleId, purity, primaryTumorLocation, status, qcStatus)
names(patientData) <- c("Sample", "Purity", "CancerType", "Status", "QCStatus")
patientData = patientData %>% left_join(ahrData, by = "Sample")
patientData[is.na(patientData)] <- "NONE"


#### GET DRIVER DATA
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis", host = "127.0.0.1")
sampleIdString = paste("'", patientData$Sample, "'", collapse = ",", sep = "")
driverQuery = paste0("select * from driverCatalog where sampleId in (",sampleIdString, ")")
driverData = dbGetQuery(dbProd, driverQuery)
dbDisconnect(dbProd)
rm(dbProd)
driverData = driverData %>% group_by(sampleId, gene) %>% arrange(gene, driver) %>% mutate(n = row_number()) %>% filter(n == 1) %>% ungroup() %>%
  select(sampleId, gene, driver, driverLikelihood)
names(driverData) <- c("Sample", "GeneName", "Driver", "DriverLikelihood")


#### SAVE DATA
save(genePanel, driverData, rawIsofoxData, rawSalmonData, patientData, file = "/Users/jon/hmf/analysis/salmon/expression.RData")


salmonData = result %>% group_by(Name) %>% mutate(TotalNumReads = sum(NumReads)) %>% filter(TotalNumReads > 0) %>% ungroup() %>%
  filter(Name %in% genePanel$transcriptId) %>%
  select(-TotalNumReads)

################ PROCESS
load(file = "/Users/jon/hmf/analysis/salmon/expression.RData")

sampleIntersection = rawSalmonData %>% filter(Sample %in% rawIsofoxData$Sample, Sample != "CPCT02020618T") %>% distinct(Sample) %>% pull(Sample)
sampleIntersection

genePanel = genePanel %>% arrange(GeneName) %>% mutate( Group = pmin(17,row_number() %/% 25))

driverData = driverData %>%
  mutate(
    Impact = ifelse(DriverLikelihood > 0.2, "MED", "LOW"), 
    Impact = ifelse(DriverLikelihood > 0.8, "HIGH", Impact),
    Driver = ifelse(Driver == "MUTATION", paste(Driver, Impact, sep = "_"), Driver)) %>%
  select(-Impact, -DriverLikelihood)

combinedData = bind_rows(
  rawIsofoxData %>% select(Sample, GeneId, GeneName, TransId, EffectiveLength, NumReads, TPM) %>% mutate(Source = "ISOFOX"),
  rawSalmonData %>% select(Sample, GeneId, GeneName, TransId, EffectiveLength, NumReads, TPM) %>% mutate(Source = "SALMON")) %>%
  left_join(patientData %>% select(-ends_with("Status")), by = "Sample") %>%
  left_join(driverData, by = c("GeneName", "Sample")) %>%
  left_join(genePanel %>% select(GeneName, Group), by = "GeneName") %>%
  group_by(Sample, Source) %>%
  mutate(normalisationFactor = sum(ifelse(TPM > 5000,0, NumReads) / EffectiveLength)) %>%
  ungroup() %>%
  mutate(TPM = NumReads / EffectiveLength / normalisationFactor * 1000000) %>%
  select(-normalisationFactor) %>%
  mutate(
    DriverCatalog = ifelse(is.na(Driver), "NONE", Driver),
    DriverCategory = ifelse(is.na(Driver), "NONE", "DRIVER")
  ) %>%
  select(-Driver)

head(combinedData)
compareData = combinedData %>% 
#  filter(grepl("CPCT02010944T", Sample)) %>%
  select(-NumReads) %>% select(-EffectiveLength) %>% spread(Source, TPM) 



  
geneData = combinedData %>% 
  group_by(Source, Sample, GeneName, Purity, CancerType, Group, DriverCatalog, DriverCategory) %>% 
  summarise(TPM = sum(TPM)) %>% ungroup()

geneData$random <- runif(nrow(geneData), 0.02, 0.1)

cohortData = geneData %>%
  #filter( grepl("CPCT02010944", Sample)) %>%
  filter(CancerType == "Urinary tract", Sample != "CPCT02010944T", Sample %in% sampleIntersection, grepl("CPCT02010944", Sample)) %>%
  group_by(GeneName, Source) %>%
  mutate(
    logTPM =log(pmax(TPM,random),2),
    logMedian = median(logTPM),
    logMax = logMedian + 2, 
    logMin = logMedian - 2) 

cohortDataSummary = cohortData %>% group_by(GeneName, Group, Source) %>%
  summarise(
    logMedian = median(logTPM),
    logMax = logMedian + 2, 
    logMin = logMedian - 2) 

sampleData = geneData %>% filter(Sample == "CPCT02010944T") %>%
  mutate(logTPM =log(pmax(TPM,0.1),2))


ggplot(geneData %>% filter(GeneName == "ACVR1B"), aes(Purity, TPM)) + geom_point(aes(color = DriverCatalog)) + ggtitle("ACVR1B") + facet_grid(~Source)
ggplot(compareData %>% filter(GeneName == "FHIT"), aes(ISOFOX, SALMON)) + geom_point()  + scale_y_log10() + scale_x_log10()


driverCatalogColours = c("red","blue","orange")
driverCatalogColours = setNames(driverCatalogColours, c("DRIVER", "NONE", "SPLICE"))
sourceColours = c("red", "blue")
sourceColours = setNames(sourceColours, c("ISOFOX", "SALMON"))

targetGroup = 16
pdf(file="/Users/jon/hmf/analysis/salmon/expression.pdf",width=20, height = 6)
for (targetGroup in unique(genePanel$Group)) {
  cat("Processing", targetGroup, "\n")
  
  plot = ggplot() +
    geom_tile(aes(x=GeneName, width = 0.8, y = (logMin + logMax) /2, height = logMax - logMin), data = cohortDataSummary %>% filter(Group == targetGroup),  stat="identity", fill = "green") + 
    geom_violin(aes(x=GeneName,y=logTPM, fill = Source), data = cohortData %>% filter(Group == targetGroup), scale='count',draw_quantiles = c(0.05,0.5,0.95),bw=0.7, size = 0.1) + 
    geom_point(aes(x = GeneName, y = logTPM, shape = Source), data = sampleData %>% filter(Group == targetGroup)) +
    scale_fill_manual(values=sourceColours) 
  
  print(plot)
}
dev.off()




pdf(file="/Users/jon/hmf/analysis/salmon/singleSample.pdf",width=20, height = 6)
for (targetGroup in unique(genePanel$Group)) {
  cat("Processing", targetGroup, "\n")
  plot = ggplot(cohortData %>% filter(Group == targetGroup)) + geom_point(aes(x = GeneName, y = logTPM, shape = Source, color = Sample))
  print(plot)
}
dev.off()



