library(tidyr)
library(dplyr)
library(RMySQL)
library(GenomicRanges)
library(ggplot2)
####################################### RAW DATA COLLECTION

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
dbDisconnect(dbProd)
rm(dbProd)

query_indels <- function(dbConnect, cohort) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste0("select * from somaticVariant where type = 'INDEL' and sampleId in(", sampleIdString, ")")
  return (dbGetQuery(dbConnect, query))
}

query_svs <- function(dbConnect, cohort) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste0("select * from structuralVariant where sampleId in (", sampleIdString, ")")
  return (dbGetQuery(dbConnect, query))
}

load(file = "~/hmf/paper2/RData/Reference/highestPurityCohort.RData")

hpcIndels = query_indels(dbProd, highestPurityCohort)
save(hpcIndels, file = "~/hmf/paper2/RData/Reference/hpcIndels.RData")

hpcStructuralVariants = query_svs(dbProd, highestPurityCohort)
save(hpcStructuralVariants, file= "~/hmf/paper2/RData/Reference/hpcStructuralVariants.RData")

####################################### Quality Checks
indelSampleIds = unique(hpcIndels$sampleId)
svSampleIds = unique(hpcStructuralVariants$sampleId)

missingIndels = highestPurityCohort %>% ungroup() %>% select(sampleId) %>% filter(!sampleId %in% indelSampleIds)
missingSvs = highestPurityCohort %>% ungroup() %>% select(sampleId) %>% filter(!sampleId %in% svSampleIds)
rm(missingSvs, missingIndels, indelSampleIds, svSampleIds)

####################################### ANALYSIS
load(file = "~/hmf/paper2/RData/Reference/highestPurityCohort.RData")
load(file = "~/hmf/paper2/RData/Reference/hpcStructuralVariants.RData")
load(file = "~/hmf/paper2/RData/Reference/hpcIndels.RData")

scope_summary <- function(svs, indels) {
  indelRanges = GRanges(paste0(indels$chromosome, indels$sampleId), IRanges(indels$position, indels$position + indels$length))
  svRanges = GRanges(paste0(svs$startChromosome, svs$sampleId), IRanges(svs$startPosition, svs$endPosition))
  ol = as.matrix(findOverlaps(indelRanges, svRanges, type="any", select="all", maxgap = 20))
  
  indels$scope <- "Private"
  indels[ol[, 1], "scope"]<- "Shared"
  indels[ol[, 1], "svId"] <- svs[ol[,2], "id"]
  
  svs$scope <- "Private"
  svs[ol[, 2], "scope"]<- "Shared"
  svs[ol[, 2], "indelId"] <- indels[ol[, 1], "id"]
  
  svSummary = svs %>% group_by(length, scope) %>% count() %>% spread(scope, n) %>% select(length, sv = Private, svShared = Shared)
  indelSummary = indels %>% group_by(length, scope) %>% count() %>% spread(scope, n) %>% select(length, indel = Private, indelShared = Shared)
  
  summary = full_join(svSummary, indelSummary, by = "length") 
  summary[is.na(summary)] <- 0
  summary = summary %>% mutate(shared = round((svShared + indelShared) / 2))
  
  return (summary)
}

scope_plot <- function(summary) {
  df = summary %>% select(-svShared) %>% gather(scope, count, sv, indel, shared) %>%
    mutate(scope = factor(scope, levels = c("sv", "indel", "shared"), ordered = T))
  
  ggplot(data = df) +
    geom_bar(stat = "identity", aes(x = length, y = count, fill = scope))
}


indels = hpcIndels %>% 
  filter(nchar(ref) > nchar(alt), filter == 'PASS') %>%
  mutate(length = abs(nchar(ref) - nchar(alt))) %>%
  filter(length >= 10, length <= 100) 

svs = hpcStructuralVariants %>%
  filter(type == 'DEL', filter == '') %>%
  mutate(length = endPosition - startPosition - 1, oddHomology = nchar(startHomologySequence) %% 2) %>%
  mutate(length = ifelse(oddHomology == 1, length + 1, length)) %>%
  filter(length >= 10, length <= 100) 

dels = scope_summary(svs, indels)
delPlot = scope_plot(dels) + xlim(24, 76) + ggtitle("DELS")


indels = hpcIndels %>% 
  filter(nchar(ref) < nchar(alt), filter == 'PASS') %>%
  mutate(length = abs(nchar(ref) - nchar(alt))) %>%
  filter(length >= 10, length <= 100) 

svs = hpcStructuralVariants %>%
  filter(type %in% c('DUP', 'INS'), filter == '') %>%
  mutate(length = endPosition - startPosition + nchar(insertSequence), oddHomology = nchar(startHomologySequence) %% 2) %>%
#  mutate(length = ifelse(oddHomology == 1, length + 1, length)) %>%
  filter(length >= 10, length <= 100) 

dups = scope_summary(svs, indels)
dupPlot = scope_plot(dups) + xlim(24, 76) + ggtitle("DUPS")

cowplot::plot_grid(delPlot, dupPlot, nrow = 1)


############################################# GRIDSS COMPARISON
svLengthLabels = c("low", "high")
svLengthBreaks = c(32, 64, 128)

svLengthsPerSample = hpcStructuralVariants %>%
  filter(type == 'DEL', filter == '') %>%
  mutate(length = endPosition - startPosition - 1, oddHomology = nchar(startHomologySequence) %% 2) %>%
  mutate(length = ifelse(oddHomology == 1, length + 1, length)) %>%
  filter(length >=32, length <= 127) %>%
  mutate(svLength = cut(length, breaks = svLengthBreaks, labels = svLengthLabels, include.lowest = T, right = F)) %>%
  group_by(sampleId, svLength) %>% count() %>% spread(svLength, n, fill = 0)

ggplot(svLengthsPerSample) + 
  geom_point(aes(x = low, y = high)) + ggtitle("Gridss Comparison") + xlab("Lengths 32-63") + ylab("Length 64-127")

############################################# Scatter DELS
svLengthsPerSample = hpcStructuralVariants %>%
  filter(type == 'DEL', filter == '') %>%
  mutate(length = endPosition - startPosition - 1, oddHomology = nchar(startHomologySequence) %% 2) %>%
  mutate(length = ifelse(oddHomology == 1, length + 1, length)) %>%
  filter(length >= 33, length <= 100) %>%
  group_by(sampleId) %>% summarise(svCount = n()) 

indelLengthLabels = c("1-3", "4-7", "8-15", "16-31", "32+")
indelLengthBreaks = c(1, 4, 8, 16, 32, 100000000)

indelLengthsPerSample = hpcIndels %>% 
  filter(nchar(ref) > nchar(alt), filter == 'PASS') %>%
  mutate(length = abs(nchar(ref) - nchar(alt)), repeating = repeatCount >= 4) %>%
  filter(length >= 1, length <= 31) %>%
  mutate(indelLength = cut(length, breaks = indelLengthBreaks, labels = indelLengthLabels, include.lowest = T, right = F)) %>%
  group_by(sampleId, repeating, indelLength) %>% summarise(indelCount = n()) 


combined = full_join(svLengthsPerSample, indelLengthsPerSample, by = "sampleId") %>% mutate(repeating = ifelse(repeating, "RepeatCount >= 4", "RepeatCount < 4"))
combined[is.na(combined)] <- 0

ggplot(combined ) +
#ggplot(combined %>% filter(indelLength == '1-3')) + 
  ylim(0, 1000) +
  #scale_y_log10() +
  geom_point(aes(x = svCount, y = indelCount)) + facet_grid(repeating~indelLength, scales = "free_y") + xlab("Count SV DEL < 100 bases") 



############################################# Scatter DUPS
svLengthsPerSample = hpcStructuralVariants %>%
  filter(type %in% c('DUP', 'INS'), filter == '') %>%
  mutate(length = endPosition - startPosition + nchar(insertSequence), oddHomology = nchar(startHomologySequence) %% 2) %>%
  filter(length >= 33, length <= 100) %>%
  group_by(sampleId) %>% summarise(svCount = n()) 

indelLengthLabels = c("1-3", "4-7", "8-15", "16-31", "32+")
indelLengthBreaks = c(1, 4, 8, 16, 32, 100000000)

indelLengthsPerSample = hpcIndels %>% 
  filter(nchar(ref) < nchar(alt), filter == 'PASS') %>%
  mutate(length = abs(nchar(ref) - nchar(alt)), repeating = repeatCount >= 4) %>%
  filter(length >= 1, length <= 31) %>%
  mutate(indelLength = cut(length, breaks = indelLengthBreaks, labels = indelLengthLabels, include.lowest = T, right = F)) %>%
  group_by(sampleId, repeating, indelLength) %>% summarise(indelCount = n()) 


combined = full_join(svLengthsPerSample, indelLengthsPerSample, by = "sampleId") %>% mutate(repeating = ifelse(repeating, "RepeatCount >= 4", "RepeatCount < 4"))
combined[is.na(combined)] <- 0

ggplot(combined ) +
  #ggplot(combined %>% filter(indelLength == '1-3')) + 
  ylim(0, 200) +
  #scale_y_log10() +
  geom_point(aes(x = svCount, y = indelCount)) + facet_grid(repeating~indelLength, scales = "free_y") + xlab("Count SV (DUP, INS) < 100 bases") 




#################################### VALIDATION
validationCohort = c("CPCT02030516T","CPCT02450014T","CPCT02130091T","CPCT02050327T","CPCT02070386T","CPCT02330102T","CPCT02370037T","CPCT02030461T","CPCT02160052T","CPCT02120143T","CPCT02150016T","DRUP01330008T","DRUP01010096T")
indels = hpcIndels %>% filter(sampleId %in% validationCohort) %>%
  filter(nchar(ref) > nchar(alt), filter == 'PASS') %>%
  mutate(length = abs(nchar(ref) - nchar(alt))) %>%
  filter(length >= 15, length <= 75)
  
svs = hpcStructuralVariants %>% 
  filter(sampleId %in% validationCohort) %>%
  filter(type == 'DEL', filter == '') %>%
  mutate(
    oddHomology = nchar(startHomologySequence) %% 2, 
     endPosition = ifelse(oddHomology == 1,  endPosition + 1, endPosition),
    length = endPosition - startPosition - 1) %>%
  filter(length >= 15, length <= 75)
  

write.csv(indels, file = "~/hmf/paper2/indelValidationDel.csv")
write.csv(svs, file = "~/hmf/paper2/svValidationDel.csv")



indels = hpcIndels %>% filter(sampleId %in% validationCohort) %>%
  filter(nchar(ref) < nchar(alt), filter == 'PASS') %>%
  mutate(length = abs(nchar(ref) - nchar(alt))) %>%
  filter(length >= 15, length <= 75)

svs = hpcStructuralVariants %>% 
  filter(sampleId %in% validationCohort) %>%
  filter(type %in% c('INS', 'DUP'), filter == '') %>%
  mutate(
    oddHomology = nchar(startHomologySequence) %% 2, 
    endPosition = ifelse(oddHomology == 1,  endPosition + 1, endPosition),
    length = endPosition - startPosition - 1 + nchar(insertSequence)) %>%
  filter(length >= 15, length <= 75)

write.csv(indels, file = "~/hmf/paper2/indelValidationDup.csv")
write.csv(svs, file = "~/hmf/paper2/svValidationDup.csv")



