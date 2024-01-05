wd <- paste0(Sys.getenv("HOME"), "/hmf/tmp/")
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(DBI)
library(RMySQL)
dbConnect <- dbConnect(MySQL(), dbname='actin_pilot', groups="RAnalysis")

## Count of CUP patients
queryCUPCount <- "
select distinct(m.patientId)
from tumor t
inner join molecular m on m.patientId=t.patientId
left join variant v on v.sampleId=m.sampleId
where primaryTumorSubLocation = 'CUP' and containsTumorCells"
queryCUPCountResult <- dbGetQuery(dbConnect, queryCUPCount)
CUPpatients=as.integer(count(queryCUPCountResult))

## Count of LR patients
queryLRCount <- "
select distinct(m.patientId)
from tumor t
inner join molecular m on m.patientId=t.patientId
left join variant v on v.sampleId=m.sampleId
where primaryTumorSubLocation != 'CUP' and primaryTumorSubLocation is not null and containsTumorCells"
queryLRCountResult <- dbGetQuery(dbConnect, queryLRCount)
LRpatients=as.integer(count(queryLRCountResult))

## Count of Altered Genes in entire cohort (Reportable)
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

## Count of Altered Genes in CUP Cases (Reportable)
queryGeneCountCUP <- "select 
gene, 
count(distinct case when driverLikelihood = 'High' then m.patientId end) as high_driver_likelihood_count,
count(distinct m.patientId) as gene_count,
count(distinct case when isHotspot then m.patientId end) as hotspot_count
from tumor t 
inner join molecular m on m.patientId=t.patientId 
left join variant v on v.sampleId=m.sampleId
where containsTumorCells and isReportable and primaryTumorSubLocation='CUP'
group by gene"
queryGeneCountCUPResult <- dbGetQuery(dbConnect, queryGeneCountCUP)

## Count of Altered Genes in LR Cases (Reportable)
queryGeneCountLR <- "select
gene, 
count(distinct case when driverLikelihood = 'High' then m.patientId end) as high_driver_likelihood_count,
count(distinct m.patientId) as gene_count,
count(distinct case when isHotspot then m.patientId end) as hotspot_count
from tumor t 
inner join molecular m on m.patientId=t.patientId 
left join variant v on v.sampleId=m.sampleId
where containsTumorCells and isReportable and primaryTumorSubLocation != 'CUP' and primaryTumorSubLocation is not null
group by gene"
queryGeneCountLRResult <- dbGetQuery(dbConnect, queryGeneCountLR)

## Overview of Copy Number alterations - CUP
queryCopyNumberCUP <- "select 
event, count(distinct m.patientId) as Count
from tumor t 
inner join molecular m on m.patientId=t.patientId 
left join copyNumber c on c.sampleId=m.sampleId
where containsTumorCells and isReportable and primaryTumorSubLocation = 'CUP'
group by event;"
queryCopyNumberCUPResult <- dbGetQuery(dbConnect, queryCopyNumberCUP)

## Overview of Copy Number alterations - LR
queryCopyNumberLR <- "select 
event, count(distinct m.patientId) as Count
from tumor t 
inner join molecular m on m.patientId=t.patientId 
left join copyNumber c on c.sampleId=m.sampleId
where containsTumorCells and isReportable and primaryTumorSubLocation != 'CUP' and primaryTumorSubLocation is not null
group by event;"
queryCopyNumberLRResult <- dbGetQuery(dbConnect, queryCopyNumberLR)

## Overview of variant treatment evidence in CUP Cases (Reportable and High driver)
queryTreatmentCountCUP <- "select e.type, Count(distinct m.patientId) as treatment_type_count
from variant v
inner join variantEvidence e on e.variantId=v.id 
inner join molecular m on m.sampleId=v.sampleId
inner join tumor t on t.patientId = m.patientId
where containsTumorCells and isReportable and driverLikelihood = 'High' and primaryTumorSubLocation = 'CUP'
group by e.type"
queryTreatmentCountCUPResult <- dbGetQuery(dbConnect, queryTreatmentCountCUP)

## Overview of variant treatment evidence in LR Cases (Reportable and High driver)
queryTreatmentCountLR <- "select e.type, Count(distinct m.patientId) as treatment_type_count
from variant v
inner join variantEvidence e on e.variantId=v.id 
inner join molecular m on m.sampleId=v.sampleId
inner join tumor t on t.patientId = m.patientId
where containsTumorCells and isReportable and driverLikelihood = 'High' and primaryTumorSubLocation != 'CUP' and primaryTumorSubLocation is not null
group by e.type"
queryTreatmentCountLRResult <- dbGetQuery(dbConnect, queryTreatmentCountLR)

## Average number of high driver likelihood variants per LR sampleId
queryHighDriverPerLRSample <- "select avg(highDriverCount)
from
(select m.sampleId, count(case when driverLikelihood='High' then m.sampleId end) as highDriverCount
from molecular m
inner join tumor t on m.patientId=t.patientId
left join variant v on m.sampleId=v.sampleId
where containsTumorCells and primaryTumorSubLocation != 'CUP' and primaryTumorSubLocation is not null
group by m.sampleId) as subquery"
queryHighDriverPerLRSampleResult <- dbGetQuery(dbConnect, queryHighDriverPerLRSample)

## Average number of high driver likelihood variants per CUP sampleId
queryHighDriverPerCUPSample <- "select avg(highDriverCount)
from
(select m.sampleId, count(case when driverLikelihood='High' then m.sampleId end) as highDriverCount
from molecular m
inner join tumor t on m.patientId=t.patientId
left join variant v on m.sampleId=v.sampleId
where containsTumorCells and primaryTumorSubLocation = 'CUP'
group by m.sampleId) as subquery"
queryHighDriverPerCUPSampleResult <- dbGetQuery(dbConnect, queryHighDriverPerCUPSample)


## Average number of variant treatment evidence (for high driver likelihood variants) per sampleId
queryVariantEvidencePerSample <- "
select round(avg(variantEvidenceCount), 1) as average_evidence_per_sampleId
    from ( 
		select m.sampleId, count(distinct treatment) as variantEvidenceCount
		from molecular m
		left join variant v on v.sampleId = m.sampleId
		left join variantEvidence e on e.variantId = v.id
		where containsTumorCells
		group by m.sampleId
	) as subquery;
"
queryVariantEvidencePerSampleResult <- dbGetQuery(dbConnect, queryVariantEvidencePerSample)

## Tumor mutational load - entire cohort
queryTMLCount <- "select hasHighTumorMutationalLoad, count(distinct patientId) as TML_count
from molecular
where containsTumorCells
group by hasHighTumorMutationalLoad"
queryTMLCountResult <- dbGetQuery(dbConnect, queryTMLCount)

## Tumor mutational burden - entire cohort
queryTMBCount <- "select hasHighTumorMutationalBurden, count(distinct patientId) as TMB_count
from molecular
where containsTumorCells
group by hasHighTumorMutationalBurden"
queryTMBCountResult <- dbGetQuery(dbConnect, queryTMBCount)

## Tumor homologous repair deficiency - entire cohort
queryHRDCount <- "select isHomologousRepairDeficient, count(distinct patientId) as HRD_count
from molecular
where containsTumorCells
group by isHomologousRepairDeficient"
queryHRDCountResult <- dbGetQuery(dbConnect, queryHRDCount)

## Tumor purity overview
queryPurity <- "select purity
from molecular
where containsTumorCells"
queryPurityResult <- dbGetQuery(dbConnect, queryPurity)

dbDisconnect(dbConnect)
#View(queryGeneCountResult)
#View(queryGeneCountCUPResult)
#View(queryGeneCountLRResult)
#View(queryCopyNumberCUPResult)
#View(queryCopyNumberLRResult)
#View(queryTreatmentCountCUPResult)
#View(queryTreatmentCountLRResult)
#View(queryHighDriverPerLRSampleResult)
#View(queryHighDriverPerCUPSampleResult)
#View(queryVariantEvidencePerSampleResult)
#View(queryTMLCountResult)
#View(queryTMBCountResult)
#View(queryHRDCountResult)
#View(queryPurityResult)

# Predominantly Altered Genes in CUP Cases (Reportable) -----------------------
CUP_gene_sorted <- queryGeneCountCUPResult %>% arrange(desc(high_driver_likelihood_count)) %>% rename("Gene" = gene, "Gene Count" = gene_count, "Hotspot Count" = hotspot_count, "High Driver Likelihood Count" = high_driver_likelihood_count)
CUP_gene_10 <- head(CUP_gene_sorted, 10)

CUP_gene_10_with_percentage <- CUP_gene_10 %>%
  mutate(across(where(is.numeric), 
                ~ paste0(round(.), " (", round((. / CUPpatients) * 100, 0), "%)")))

pdf(file= paste0(wd,"Genes_CUP_Reportable.pdf"), width = 10, height = 7)
CUP_gene_table <- tableGrob(CUP_gene_10_with_percentage, rows = NULL)
grid.draw(CUP_gene_table)
invisible(dev.off())

# Predominantly Altered Genes in Last Resort Cases (Reportable) ----------------
LR_gene_sorted <- queryGeneCountLRResult %>% arrange(desc(high_driver_likelihood_count)) %>% rename("Gene" = gene, "Gene Count" = gene_count, "Hotspot Count" = hotspot_count, "High Driver Likelihood Count" = high_driver_likelihood_count)
LR_gene_10 <- head(LR_gene_sorted, 10)

LR_gene_10_with_percentage <- LR_gene_10 %>%
  mutate(across(where(is.numeric), 
                ~ paste0(round(.), " (", round((. / LRpatients) * 100, 0), "%)")))

pdf(file= paste0(wd,"Genes_LR_Reportable.pdf"), width = 10, height = 7)
LR_gene_table <- tableGrob(LR_gene_10_with_percentage, rows = NULL)
grid.draw(LR_gene_table)
invisible(dev.off())

# Predominant Copy Number Alterations in CUP Cases (Reportable) ----------------
CUP_CN_sorted <- queryCopyNumberCUPResult %>% arrange(desc(Count)) %>% rename("Event" = event)
CUP_CN_10 <- head(CUP_CN_sorted, 10)

CUP_CN_10_with_percentage <- CUP_CN_10 %>%
  mutate(across(where(is.numeric), 
                ~ paste0(round(.), " (", round((. / CUPpatients) * 100, 0), "%)")))

pdf(file= paste0(wd,"Copy_Number_CUP.pdf"), width = 10, height = 7)
CUP_CN_table <- tableGrob(CUP_CN_10_with_percentage, rows = NULL)
grid.draw(CUP_CN_table)
invisible(dev.off())

# Predominantly Copy Number Alterations in Last Resort Cases (Reportable) ----------------
LR_CN_sorted <- queryCopyNumberLRResult %>% arrange(desc(Count)) %>% rename("Event" = event)
LR_CN_10 <- head(LR_CN_sorted, 10)

LR_CN_10_with_percentage <- LR_CN_10 %>%
  mutate(across(where(is.numeric), 
                ~ paste0(round(.), " (", round((. / LRpatients) * 100, 0), "%)")))

pdf(file= paste0(wd,"Copy_Number_LR.pdf"), width = 10, height = 7)
LR_CN_table <- tableGrob(LR_CN_10_with_percentage, rows = NULL)
grid.draw(LR_CN_table)
invisible(dev.off())

# Overview of treatment evidence type in CUP Cases (Reportable and High driver) ----------------
CUP_treatment_sorted <- queryTreatmentCountCUPResult %>% arrange(desc(treatment_type_count)) %>% rename("Treatment Category" = type, "Count (%)" = treatment_type_count)
CUP_treatment_10 <- head(CUP_treatment_sorted, 10)

CUP_treatment_10_with_percentage <- CUP_treatment_10 %>%
  mutate(across(where(is.numeric), 
                ~ paste0(round(.), " (", round((. / CUPpatients) * 100, 0), "%)")))

pdf(file= paste0(wd,"Variant_treatment_Evidence_CUP.pdf"), width = 10, height = 7)
CUP_treatment_table <- tableGrob(CUP_treatment_10_with_percentage, rows = NULL)
grid.draw(CUP_treatment_table)
invisible(dev.off())

# Overview of treatment evidence type in LR Cases (Reportable and High driver) ----------------
LR_treatment_sorted <- queryTreatmentCountLRResult %>% arrange(desc(treatment_type_count)) %>% rename("Treatment Category" = type, "Count (%)" = treatment_type_count)
LR_treatment_10 <- head(LR_treatment_sorted, 10)

LR_treatment_10_with_percentage <- LR_treatment_10 %>%
  mutate(across(where(is.numeric), 
                ~ paste0(round(.), " (", round((. / LRpatients) * 100, 0), "%)")))

pdf(file= paste0(wd,"Variant_treatment_Evidence_LR.pdf"), width = 10, height = 7)
LR_treatment_table <- tableGrob(LR_treatment_10_with_percentage, rows = NULL)
grid.draw(LR_treatment_table)
invisible(dev.off())

## Average amount of high driver likelihood variants per sampleId 
cat("Average found variants per LR sample:", as.double(queryHighDriverPerLRSampleResult), "\n")

## Average amount of high driver likelihood variants per sampleId
cat("Average found variants per CUP sample:", as.double(queryHighDriverPerCUPSampleResult), "\n")

## Average amount of variant treatment evidence (for high driver likelihood variants) per sampleId
cat("Average found evidence-based treatments per sample:", as.double(queryVariantEvidencePerSampleResult), "\n")

# Tumor Mutational Load - entire cohort ----------------------------------------------
TMLhigh <- queryTMLCountResult$TML_count[queryTMLCountResult$hasHighTumorMutationalLoad==1]
TMLlow <- queryTMLCountResult$TML_count[queryTMLCountResult$hasHighTumorMutationalLoad==0]
slices <- c(TMLhigh, TMLlow)
lbls <- c("Yes:", "No:")
pct <- round(slices/sum(slices)*100, 1)
lbls <- c(paste0("Yes: ", pct[1], "% (N=", TMLhigh, ")"),
          paste0("No: ", pct[2], "% (N=", TMLlow, ")"))

pdf(file= paste0(wd,"Tumor_Mutational_Load.pdf"), width = 10, height = 7)
pie(slices, labels = lbls, col=c("green4","red"),main="High Tumor Mutational Load")
invisible(dev.off())

noquote(paste0(pct[1], "% of patients had a high tumor mutational load sample"))

# Tumor Mutational Burden - entire cohort ----------------------------------------------
TMBhigh <- queryTMBCountResult$TMB_count[queryTMBCountResult$hasHighTumorMutationalBurden==1]
TMBlow <- queryTMBCountResult$TMB_count[queryTMBCountResult$hasHighTumorMutationalBurden==0]
slices <- c(TMBhigh, TMBlow)
lbls <- c("Yes:", "No:")
pct <- round(slices/sum(slices)*100, 1)
lbls <- c(paste0("Yes: ", pct[1], "% (N=", TMBhigh, ")"),
          paste0("No: ", pct[2], "% (N=", TMBlow, ")"))

pdf(file= paste0(wd,"Tumor_Mutational_Burden.pdf"), width = 10, height = 7)
pie(slices, labels = lbls, col=c("green4","red"),main="High Tumor Mutational Burden")
invisible(dev.off())

noquote(paste0(pct[1], "% of patients had a high tumor mutational burden sample"))

# HRD - entire cohort ----------------------------------------------
HRD <- na.omit(queryHRDCountResult$HRD_count[queryHRDCountResult$isHomologousRepairDeficient==1])
HRP <- na.omit(queryHRDCountResult$HRD_count[queryHRDCountResult$isHomologousRepairDeficient==0])
Undetermined <- queryHRDCountResult$HRD_count[is.na(queryHRDCountResult$isHomologousRepairDeficient)]

slices <- c(HRD, HRP, Undetermined)
lbls <- c("Yes", "No", "NA")
pct <- round(slices/sum(slices)*100, 1)
lbls <- paste0(lbls, ": ", pct, "% (N=", slices, ")")

pdf(file= paste0(wd,"Homologous_Repair_Deficiency.pdf"), width = 10, height = 7)
pie(slices, labels = lbls, col=c("green4","red", "gray"),main="Homologous Repair Deficiency")
invisible(dev.off())

noquote(paste0(pct[1], "% of patients had a homologous repair deficient sample"))

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

noquote(paste0("Median tumor purity for samples with detected tumor ", median(data$purity)*100, "% (IQR: ", IQR(data$purity),")"))

