wd <- paste0(Sys.getenv("HOME"), "/hmf/R/2023/Q4/")
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(DBI)
library(RMySQL)
dbConnect <- dbConnect(MySQL(), dbname='actin_pilot', groups="RAnalysis")

queryGeneCount <- "select 
gene, 
count(distinct m.patientId) as gene_count,
count(distinct case when isHotspot then m.patientId end) as hotspot_count,
count(distinct case when driverLikelihood = 'High' then m.patientId end) as high_driver_likelihood_count
from tumor t 
inner join molecular m on m.patientId=t.patientId 
left join variant v on v.sampleId=m.sampleId
where containsTumorCells and isReportable and primaryTumorSubLocation is not null
group by gene"
queryGeneCountResult <- dbGetQuery(dbConnect, queryGeneCount)

queryGeneCountCUP <- "select 
gene, 
count(distinct m.patientId) as gene_count,
count(distinct case when isHotspot then m.patientId end) as hotspot_count,
count(distinct case when driverLikelihood = 'High' then m.patientId end) as high_driver_likelihood_count
from tumor t 
inner join molecular m on m.patientId=t.patientId 
left join variant v on v.sampleId=m.sampleId
where containsTumorCells and isReportable and primaryTumorSubLocation='CUP'
group by gene"
queryGeneCountCUPResult <- dbGetQuery(dbConnect, queryGeneCountCUP)

queryGeneCountLR <- "select
gene, 
count(distinct m.patientId) as gene_count,
count(distinct case when isHotspot then m.patientId end) as hotspot_count,
count(distinct case when driverLikelihood = 'High' then m.patientId end) as high_driver_likelihood_count
from tumor t 
inner join molecular m on m.patientId=t.patientId 
left join variant v on v.sampleId=m.sampleId
where containsTumorCells and isReportable and primaryTumorSubLocation != 'CUP' and primaryTumorSubLocation is not null
group by gene"
queryGeneCountLRResult <- dbGetQuery(dbConnect, queryGeneCountLR)

queryTMLCount <- "select hasHighTumorMutationalLoad, count(distinct sampleId) as TML_count
from molecular
group by hasHighTumorMutationalLoad"
queryTMLCountResult <- dbGetQuery(dbConnect, queryTMLCount)

queryTMBCount <- "select hasHighTumorMutationalBurden, count(distinct sampleId) as TMB_count
from molecular
group by hasHighTumorMutationalBurden"
queryTMBCountResult <- dbGetQuery(dbConnect, queryTMBCount)

queryHRDCount <- "select isHomologousRepairDeficient, count(distinct sampleId) as HRD_count
from molecular
group by isHomologousRepairDeficient"
queryHRDCountResult <- dbGetQuery(dbConnect, queryHRDCount)

queryPurity <- "select purity from molecular"
queryPurityResult <- dbGetQuery(dbConnect, queryPurity)

dbDisconnect(dbConnect)
View(queryGeneCountResult)
View(queryGeneCountCUPResult)
View(queryGeneCountLRResult)
View(queryTMLCountResult)
View(queryTMBCountResult)
View(queryHRDCountResult)
View(queryPurityResult)

# Predominantly Altered Genes in Last Resort Cases (Reportable) ----------------
LR_sorted <- queryGeneCountLRResult %>% arrange(desc(gene_count))
LR_gene_10 <- head(LR_sorted, 10)

pdf(file= paste0(wd,"Genes_LR_Reportable.pdf"), width = 10, height = 7)
LR_table <- tableGrob(LR_gene_10, rows = NULL)
grid.draw(LR_table)
invisible(dev.off())

# Predominantly Altered Genes in CUP Cases (Reportable) -----------------------
CUP_sorted <- queryGeneCountCUPResult %>% arrange(desc(gene_count))
CUP_gene_10 <- head(CUP_sorted, 10)

pdf(file= paste0(wd,"Genes_CUP_Reportable.pdf"), width = 10, height = 7)
CUP_table <- tableGrob(CUP_gene_10, rows = NULL)
grid.draw(CUP_table)
invisible(dev.off())

# Tumor Mutational Load - entire cohort ----------------------------------------------
TMLhigh <- queryTMLCountResult$TML_count[queryTMLCountResult$hasHighTumorMutationalLoad==1]
TMLlow <- queryTMLCountResult$TML_count[queryTMLCountResult$hasHighTumorMutationalLoad==0]
slices <- c(TMLhigh, TMLlow)
lbls <- c("Yes:", "No:")
pct <- round(slices/sum(slices)*100, 1)
lbls <- c(paste0("Yes: ", pct[1], "% (N=", TMLhigh, ")"),
          paste0("No: ", pct[2], "% (N=", TMLlow, ")"))

pdf(file= paste0(wd,"Tumor_Mutational_Load.pdf"), width = 10, height = 7)
pie(slices, labels = lbls, col=c("red","blue"),main="High Tumor Mutational Load")
invisible(dev.off())

noquote(paste0(pct[1], "% of patients had a high tumor mutational load"))

# Tumor Mutational Burden - entire cohort ----------------------------------------------
TMBhigh <- queryTMBCountResult$TMB_count[queryTMBCountResult$hasHighTumorMutationalBurden==1]
TMBlow <- queryTMBCountResult$TMB_count[queryTMBCountResult$hasHighTumorMutationalBurden==0]
slices <- c(TMBhigh, TMBlow)
lbls <- c("Yes:", "No:")
pct <- round(slices/sum(slices)*100, 1)
lbls <- c(paste0("Yes: ", pct[1], "% (N=", TMBhigh, ")"),
          paste0("No: ", pct[2], "% (N=", TMBlow, ")"))

pdf(file= paste0(wd,"Tumor_Mutational_Burden.pdf"), width = 10, height = 7)
pie(slices, labels = lbls, col=c("red","blue"),main="High Tumor Mutational Burden")
invisible(dev.off())

noquote(paste0(pct[1], "% of patients had a high tumor mutational burden"))

# HRD - entire cohort ----------------------------------------------
HRD <- na.omit(queryHRDCountResult$HRD_count[queryHRDCountResult$isHomologousRepairDeficient==1])
HRP <- na.omit(queryHRDCountResult$HRD_count[queryHRDCountResult$isHomologousRepairDeficient==0])
Undetermined <- queryHRDCountResult$HRD_count[is.na(queryHRDCountResult$isHomologousRepairDeficient)]

slices <- c(HRD, HRP)
lbls <- c("Yes:", "No:", "Undetermined:")
pct <- round(slices/sum(slices)*100, 1)
lbls <- c(paste0("Yes: ", pct[1], "% (N=", HRD, ")"),
          paste0("No: ", pct[2], "% (N=", HRP, ")"))

pdf(file= paste0(wd,"Homologous_Repair_Deficiency.pdf"), width = 10, height = 7)
pie(slices, labels = lbls, col=c("red","blue"),main="Homologous Repair Deficiency")
invisible(dev.off())

noquote(paste0(pct[1], "% of patients were homologous repair deficient (NA for n=", Undetermined,")"))

# Tumor purity distribution - entrire cohort ----------------------------------------------
pdf(file= paste0(wd,"Tumor_Purity.pdf"), width = 10, height = 7)

data <- queryPurityResult
ggplot(data, aes(x = purity)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "darkblue", alpha = 0.7) +
  geom_vline(xintercept = median(data$purity), color = "red", linetype = "dashed", linewidth = 1.5) +
  geom_segment(aes(x = median(purity) - IQR(purity) / 2, xend = median(purity) + IQR(purity) / 2,
                   y = 0, yend = 0), color = "blue", linetype = "dashed", linewidth = 1.5) +
  labs(title = "Histogram of Tumor Purity with Median and IQR",
       x = "Tumor Purity",
       y = "Frequency") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme_minimal()
invisible(dev.off())

noquote(paste0("Median Tumor purity ", median(data$purity), "% (IQR: ", IQR(data$purity),")"))
