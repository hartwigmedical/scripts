library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(DBI)
library(RMySQL)

wd <- paste0(Sys.getenv("HOME"), "/hmf/tmp/")
dbConnect <- dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")

## Misallocation
v2MisallocationQuery <- "
select percent from hmfpatients.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.clinical c on c.sampleId = s.sampleId
where s.signature = 'MISALLOC';
"

v3MisallocationQuery <- "
select percent from hmfpatients_pilot.signature s
inner join hmfpatients.hpc h on h.sampleId = s.sampleId
inner join hmfpatients.clinical c on c.sampleId = s.sampleId
where s.signature = 'MISALLOC';
"

v2MisallocationResult <- dbGetQuery(dbConnect, v2MisallocationQuery)
v3MisallocationResult <- dbGetQuery(dbConnect, v3MisallocationQuery)

v2_mean_misalloc_percent <- mean(v2MisallocationResult$percent, na.rm = TRUE)
v2_median_misalloc_percent <- median(v2MisallocationResult$percent, na.rm = TRUE)
v3_mean_misalloc_percent <- mean(v3MisallocationResult$percent, na.rm = TRUE)
v3_median_misalloc_percent <- median(v3MisallocationResult$percent, na.rm = TRUE)

print(paste("Mean Percent v2:", v2_mean_misalloc_percent))
print(paste("Median Percent v2:", v2_median_misalloc_percent))
print(paste("Mean Percent v3:", v3_mean_misalloc_percent))
print(paste("Median Percent v3:", v3_median_misalloc_percent))

## v2 signatures means for entire cohort
v2_mean_table <- data.frame(Signature = character(), MeanPercent = numeric(), stringsAsFactors = FALSE)

for (i in 1:30) {
  v2_query <- sprintf("
    select percent from hmfpatients.signature s
    inner join hmfpatients.hpc h on h.sampleId = s.sampleId
    where s.signature = 'Sig%d';
  ", i)

  result <- dbGetQuery(dbConnect, v2_query)
  mean_percent <- mean(result$percent, na.rm = TRUE)
  v2_mean_table <- rbind(v2_mean_table, data.frame(Signature = paste0(i), MeanPercent = mean_percent))
}
print(v2_mean_table)

dbDisconnect(dbConnect)
dbConnect <- dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")

## v3 signatures means for entire cohort
v3_mean_table <- data.frame(Signature = character(), MeanPercent = numeric(), stringsAsFactors = FALSE)

signatures <- c(
  "BI_COMPOSITE_SNV_SBS1_P", "BI_COMPOSITE_SNV_SBS2_P", "BI_COMPOSITE_SNV_SBS3_P", "BI_COMPOSITE_SNV_SBS4_P",
  "BI_COMPOSITE_SNV_SBS5_P", "BI_COMPOSITE_SNV_SBS6_S", "BI_COMPOSITE_SNV_SBS7a_S", "BI_COMPOSITE_SNV_SBS7b_S",
  "BI_COMPOSITE_SNV_SBS7c_S", "BI_COMPOSITE_SNV_SBS8_P", "BI_COMPOSITE_SNV_SBS9_P", "BI_COMPOSITE_SNV_SBS10a_S",
  "BI_COMPOSITE_SNV_SBS11_S", "BI_COMPOSITE_SNV_SBS12_P", "BI_COMPOSITE_SNV_SBS13_P", "BI_COMPOSITE_SNV_SBS14_S",
  "BI_COMPOSITE_SNV_SBS15_S", "BI_COMPOSITE_SNV_SBS16_P", "BI_COMPOSITE_SNV_SBS17a_P", "BI_COMPOSITE_SNV_SBS17b_P",
  "BI_COMPOSITE_SNV_SBS18_P", "BI_COMPOSITE_SNV_SBS19_P", "BI_COMPOSITE_SNV_SBS21_S", "BI_COMPOSITE_SNV_SBS22_P",
  "BI_COMPOSITE_SNV_SBS26_S", "BI_COMPOSITE_SNV_SBS28_P", "BI_COMPOSITE_SNV_SBS30_P"
)

for (sig in signatures) {
  v3_query <- sprintf("
    select percent from hmfpatients_pilot.signature s
    inner join hmfpatients.hpc h on h.sampleId = s.sampleId
    where s.signature like '%%%s%%';
  ", sig)

  result <- dbGetQuery(dbConnect, v3_query)
  mean_percent <- mean(result$percent, na.rm = TRUE)
  signature_number = gsub("\\D", "", sig)

  if (signature_number %in% v3_mean_table$Signature) {
    v3_mean_table$MeanPercent[v3_mean_table$Signature == signature_number] <-
      v3_mean_table$MeanPercent[v3_mean_table$Signature == signature_number] + mean_percent
  } else {
    v3_mean_table <- rbind(v3_mean_table, data.frame(Signature = signature_number, MeanPercent = mean_percent))
  }
}

# Merge v2_mean_table and v3_mean_table
final_table <- merge(v2_mean_table, v3_mean_table, by = "Signature", all = TRUE)
colnames(final_table) <- c("SBS_Signature", "Mean_v2", "Mean_v3")
final_table <- final_table[order(as.numeric(final_table$SBS_Signature)), ]

print(v2_mean_table)
print(v3_mean_table)
print(final_table)

dbDisconnect(dbConnect)
