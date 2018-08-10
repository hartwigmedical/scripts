detach("package:purple", unload=TRUE)
library(purple)
library(RMySQL)
library(dplyr)
library(tidyr)
library(cowplot)
theme_set(theme_bw())

#### SAVE COHORT
#save(cohortToRun, cohort, file = "~/hmf/analysis/fit/cohort.RData")
load(file = "~/hmf/analysis/fit/cohort.RData")

load(file = "~/hmf/RData/Reference/cancerTypeColours.RData")
chromosomeColours = c(cancerTypeColours, "black")
chromosomeColours = setNames(chromosomeColours, c(1:22))

query_regions <- function(dbConnect, sampleId) {
  query = paste0("select  * from copyNumber where sampleId ='", sampleId, "'");
  return (dbGetQuery(dbConnect, query))
}

query_somatics <- function(dbConnect, sampleId) {
  query = paste0("select  * from somaticVariant where sampleId ='", sampleId, "' and filter = 'PASS'");
  return (dbGetQuery(dbConnect, query))
}

create_title <- function(env, purity, ploidy, qcStatus) {
  paste0(sample, " ", env, " (", round(100*purity,2), "%, ", round(ploidy,2), ", ", qcStatus, ")")
}

fitted_region_graph <- function(title, copyNumberRegions) {
  totalBafCount = sum(copyNumberRegions$bafCount)
  simpleCopyNumberRegions = copyNumberRegions %>%
    filter(bafCount > 10) %>%
    transmute(
      sampleId,
      chromosome,
      start,
      end,
      bafCount,
      observedBaf = round(observedBaf,2),
      copyNumber,
      majorAllele = round(copyNumber * actualBaf, 2),
      minorAllele = round(copyNumber - copyNumber * actualBaf, 2)
    ) 
  
  p = ggplot(simpleCopyNumberRegions, aes(x=majorAllele,y=minorAllele)) +
    geom_point(aes(size = bafCount), alpha = 0.4) +
    #geom_point(aes(size = bafCount)) +
    xlab("Major Allele") + ylab("Minor Allele") + ggtitle(title) +
    scale_x_continuous(breaks = c(-200:200), limits = c(-0.1, 4)) +
    scale_y_continuous(breaks = c(-200:200), limits = c(-0.1, 3)) +
    scale_color_gradientn(colours=rev(rainbow(1000, start=0, end=0.75))) +
    theme(panel.grid.minor = element_blank()) + 
    scale_size(range = c(1,10)) 
  
  return(p)
}

somatic_enrichment <- function(purity, somatics) {
  somatics$samplePurity <- purity
  somatics = somatics %>%
    filter(!chromosome %in% c('X','Y'), filter == 'PASS') %>%
    mutate(
      majorAllelePloidy = adjustedCopyNumber - minorAllelePloidy,
      somaticPloidy = adjustedVaf * adjustedCopyNumber,
      subclonalProbability = samplePurity / (adjustedCopyNumber * samplePurity + 2 * (1 - samplePurity)),
      majorAlleleProbability = pmax(majorAllelePloidy, 0) *  subclonalProbability,
      minConceivableAlleleCount = qbinom(0.01, 100000, totalReadCount * subclonalProbability / 100000, T),
      maxConceivableAlleleCount = qbinom(0.99, 100000, totalReadCount * majorAlleleProbability / 100000, T),
      inconsistent = alleleReadCount > maxConceivableAlleleCount,
      subclonal = alleleReadCount < minConceivableAlleleCount,
      clonality = ifelse(subclonal, "SUBCLONAL", "CLONAL"),
      clonality = ifelse(inconsistent, "INCONSISTENT", clonality),
      CN = round(adjustedCopyNumber),
      CN = ifelse(CN > 5, "6+", CN),
      CN = factor(CN, levels = c("6+", rev(-1:5)), ordered = T))
  
  return (somatics)
}

somatic_summary <- function(somatics) {
  
  emptySomaticSummary = data.frame(sampleId = "", SUBCLONAL = 0, CLONAL = 0, INCONSISTENT = 0, stringsAsFactors = F)
  
  somatics = somatics %>% 
    group_by(sampleId, clonality) %>% 
    count() %>% 
    spread(clonality, n) 
  
  result = bind_rows(somatics, emptySomaticSummary)
  result[is.na(result)] <- 0
  result = result %>%
    ungroup() %>%
    summarise(sampleId = first(sampleId), SUBCLONAL = sum(SUBCLONAL), CLONAL = sum(CLONAL), INCONSISTENT = sum(INCONSISTENT)) %>%
    mutate(inconsistentSomaticPercentage = round(INCONSISTENT / (CLONAL + SUBCLONAL + INCONSISTENT), 2))
  
  return (result)
}

find_peak_rubbish <- function(d) {
  peak_y = max(d$y)
  peak_x = first(d$x[d$y == peak_y])
  return (c(x = peak_x, y = peak_y))
}

model_peak_rubbish <- function(histo) {
  peak = find_peak_from_histogram(histo)
  peak_y = peak["y"]
  peak_x = peak["x"]
  result = dpois(round(domain * depth * expectedVaf,0), peak_x * depth * expectedVaf)
  result = result * peak_y / max(result)
  return (result)
}

find_peak_from_histogram <- function(histo) {
  histo = histo %>% filter(n > 0)
  somatics =  rep(histo$bucket, histo$n)
  d = density(somatics, bw = 0.1, from = 0, to = 10, n = 101)
  peak_y = max(d$y)
  peak_x = first(d$x[d$y == peak_y])
  
  peak_y =  peak_y * 0.1 * length(somatics)
  return (c(x = peak_x, y = peak_y))
}

create_somatic_model <- function(purity, ploidy, depth, somatics, binwidth = 0.1) {
  domain = seq(0, 10, binwidth)
  expectedVaf = purity / (ploidy * purity + 2 * (1-purity))
  myHist = create_initial_histogram(somatics)
  initialSomatics = sum(myHist$n)
  result = data.frame(x = domain, stringsAsFactors = F)
  
  i = 0
  exitLoop = F
  while (!exitLoop) { 
    modelPeak = find_peak_from_histogram(myHist)
    modelSomatics = dpois(round(domain * depth * expectedVaf,0), modelPeak["x"] * depth * expectedVaf)
    modelSomatics = round(modelSomatics * modelPeak["y"] / max(modelSomatics), 2)
    alreadyExists = modelPeak["x"] %in% colnames(result)
    if (!alreadyExists) {
      result[, paste0(modelPeak["x"])] <- modelSomatics
      myHist$n = pmax(0, myHist$n - modelSomatics)
      unexplainedSomatics = sum(myHist)
      cat(i, " ", modelPeak["x"], "\n")
    }
    i  = i + 1
    exitLoop = unexplainedSomatics / initialSomatics < 0.05 || i > 10 || alreadyExists
  }
  
  result$total = rowSums(result[, -1])
  return (result)
}

create_initial_histogram <- function(somatics) {
  myInitialSomatics = somatics %>% filter(somaticPloidy > 0, somaticPloidy < 10)
  myEmptyHist = data.frame(bucket = seq(0, 10, 0.1), n = 0)
  
  myInitialHist = myInitialSomatics %>%
    mutate(
      bucketFactor = cut(somaticPloidy, breaks = seq(0, 10, 0.1)),
      bucket = (as.numeric(bucketFactor) - 1) * 0.1
    ) %>%
    group_by(bucket) %>% 
    count() %>%
    bind_rows(myEmptyHist) %>%
    group_by(bucket) %>%
    summarise(n = sum(n)) %>% 
    arrange(bucket)
  
  return (myInitialHist)
}


create_somatic_model_rubbish <- function(purity, ploidy, depth, somatics, binwidth = 0.1) {
  
  somaticPloidyBuckets = somatics %>% filter(clonality == "CLONAL") %>% mutate(ploidyBucket = round(somaticPloidy, 0)) %>% group_by(ploidyBucket) %>% count()
  somaticPloidyBuckets

  expectedVaf = purity / (ploidy * purity + 2 * (1-purity))
  domain = seq(0, 10, binwidth)
  y1 = dpois(round(domain * depth * expectedVaf,0), depth * expectedVaf) * expectedVaf * binwidth * depth * somaticPloidyBuckets %>% filter(ploidyBucket == 1) %>% pull(n)
  if (length(y1) == 0) {
    y1 = rep(0, length(domain))
  }
  y2 = dpois(round(domain * depth * expectedVaf,0), 2 * depth * expectedVaf) * expectedVaf * binwidth * depth * somaticPloidyBuckets %>% filter(ploidyBucket == 2) %>% pull(n)
  if (length(y2) == 0) {
    y2 = rep(0, length(domain))
  }
  y3 = dpois(round(domain * depth * expectedVaf,0), 3 * depth * expectedVaf) * expectedVaf * binwidth * depth * somaticPloidyBuckets %>% filter(ploidyBucket == 3) %>% pull(n)
  if (length(y3) == 0) {
    y3 = rep(0, length(domain))
  }
  y4 = dpois(round(domain * depth * expectedVaf,0), 4 * depth * expectedVaf) * expectedVaf * binwidth * depth * somaticPloidyBuckets %>% filter(ploidyBucket == 4) %>% pull(n)
  if (length(y4) == 0) {
    y4 = rep(0, length(domain))
  }
  
  result = data.frame(x = domain, y1 = round(y1, 2), y2 = round(y2,2), y3 = round(y3,2), y4 = round(y4,2), stringsAsFactors = F)
  result$yTotal =  result$y1 +  result$y2 +  result$y3 + result$y4

  return (result)
}

somatic_graph <- function(purity, ploidy, depth, somatics, binwidth = 0.1) {

  clonalityColoursColours = setNames(c("#66c2a5","#8da0cb", "#fc8d62"), c("CLONAL", "SUBCLONAL","INCONSISTENT"))
  clonalityAlphas = setNames(c(0.8, 0.8, 0.8), c("CLONAL", "SUBCLONAL","INCONSISTENT"))
  
  somaticModel = create_somatic_model(purity, ploidy, depth, somatics, binwidth)
  tidySomaticModel = somaticModel %>% gather(peak, value, -1) %>% mutate(peak = factor(peak, colnames(somaticModel), ordered = T))
  
  p1 = ggplot(somatics) +
    geom_histogram(aes(x = somaticPloidy), alpha = 0.4, binwidth = binwidth, position = "stack") +
    xlab("Ploidy") + ggtitle("Somatic CNV") + ylab("") +
    theme(legend.position = c(0.9,0.77)) + 
    geom_line(data = tidySomaticModel, aes(x, value, color = peak)) + 
    theme(panel.border = element_blank(), panel.grid.minor.y = element_blank(), axis.ticks = element_blank()) +
    scale_x_continuous(limits = c(0, 10), breaks = c(1:10))
  
  #+ theme(legend.position = "none")
  
  p2 = ggplot(somatics) +
    geom_histogram(aes(x = somaticPloidy, fill = clonality, alpha = clonality), binwidth = 0.1, position = "identity") +
    #geom_density(bw = 0.03, aes(x = somaticPloidy,  group = clonality, y=0.1 * ..count..), position = "identity") + 
    scale_fill_manual(name = "Clonality", values = clonalityColoursColours) +
    scale_color_manual(name = "Clonality", values = clonalityColoursColours) +
    scale_alpha_manual(values = clonalityAlphas) +
    xlim(0, 4) + xlab("Ploidy") + ggtitle("Somatic Clonality") + ylab("") +
    theme(legend.position = "none") +
    theme(panel.border = element_blank(), panel.grid.minor.y = element_blank(), axis.ticks = element_blank())
  
  p3 = plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 1))
  return (p3)
}

#dev.off()
#somatic_graph(sampleDetails$purity.prod, sampleDetails$ploidy.prod,  sampleDetails$tumorCoverage, enrichedProdSomatics)
#somatic_graph(sampleDetails$purity.prod, sampleDetails$ploidy.prod, sampleDetails$tumorCoverage, enrichedProdSomatics)
#somatic_graph(sampleDetails$purity.pilot, sampleDetails$ploidy.pilot, sampleDetails$tumorCoverage, enrichedPilotSomatics)
#sum(model$y1)


copynumber_graph <- function(copyNumberRegions) {
  
  distance = copyNumberRegions %>% summarise(distance = sum(as.numeric(end) - as.numeric(start))) %>% pull(distance)
  totalBafCount = sum(copyNumberRegions$bafCount)
  
  copyNumberRegions = copyNumberRegions %>%
    filter(!chromosome %in% c('X','Y'), copyNumber < 10) %>%
    mutate(
      chromosome = factor(chromosome, levels= c(1:22), ordered = T),
      # weight = (end - start)/distance ) 
      weight = bafCount/totalBafCount ) 
  
  ggplot(copyNumberRegions) +
    geom_histogram(aes(x = copyNumber, fill = chromosome, weight = weight), binwidth = 0.1, position = "stack") +
    scale_fill_manual(values = chromosomeColours) + 
    scale_x_continuous(breaks = c(0:10)) + 
    theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
    xlab("Copy Number") + ylab("Baf Weighed Count")
}

env_graph <- function(env, purity, ploidy, qcstatus, depth, copyNumberRegions, somatics) {
  
  title = create_title(env, purity, ploidy, qcstatus)
  
  fitted_graph = fitted_region_graph(title, copyNumberRegions)
  somatic_graph = somatic_graph(purity, ploidy, depth, somatics)
  copynumber_graph = copynumber_graph(copyNumberRegions)
  return (plot_grid(fitted_graph, somatic_graph, copynumber_graph, rel_widths = c(2,2, 1.5),  nrow = 1))
}

### PRODUCTION COHORT
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
prodCohort = dbGetQuery(dbProd, "SELECT sampleId, purity, ploidy, status, qcStatus FROM purity")
prodMetrics = dbGetQuery(dbProd, "SELECT sampleId, round(tumorMeanCoverage,0) as tumorCoverage FROM metric")
dbDisconnect(dbProd); rm(dbProd)

### PILOT COHORT
dbPilot = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
pilotCohort = dbGetQuery(dbPilot, "SELECT sampleId, purity, ploidy, status, qcStatus FROM purity WHERE version = 2.15 and modified > '2018-08-04'")
dbDisconnect(dbPilot); rm(dbPilot)

#### COMPARE COHORTS
cohort = left_join(pilotCohort, prodCohort, by = "sampleId", suffix = c(".pilot", ".prod")) %>%
  mutate(
    changeInPurity = abs(purity.prod - purity.pilot),
    changeInPloidy = abs(ploidy.prod - ploidy.pilot)) %>%
  arrange(-changeInPloidy)
cohort = left_join(cohort, prodMetrics, by = "sampleId")  

cohortToRun = cohort %>% top_n(100, changeInPloidy) %>% mutate(patientId = substr(sampleId, 1, 12)) #%>% pull(patientId)

sample = "CPCT02050047T"
sample = "CPCT02060129T"
sample = "CPCT02030352T"
sample = "CPCT02040286T"
sample = "CPCT02290041T"
sample = "CPCT02080227T"
sample = "CPCT02010722TII"
sample = "CPCT02020483T"

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
dbPilot = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")

# Initialistaion
enrichedPilotSomaticsSummary = data.frame()
enrichedProdSomaticsSummary = data.frame()

for (sample in cohortToRun$sampleId) {
  sampleDetails = cohort %>% filter(sampleId == sample)
  
  #Collect Data
  prodSomatics = query_somatics(dbProd, sample)
  prodRegions = query_regions(dbProd, sample)
  pilotSomatics = query_somatics(dbPilot, sample)
  pilotRegions = query_regions(dbPilot, sample)
  
  # Enrich
  enrichedProdSomatics = somatic_enrichment(sampleDetails$purity.prod, prodSomatics)
  enrichedPilotSomatics = somatic_enrichment(sampleDetails$purity.pilot, pilotSomatics)
  
  # Somatic Summary
  enrichedProdSomaticsSummary = bind_rows(somatic_summary(enrichedProdSomatics), enrichedProdSomaticsSummary)
  enrichedPilotSomaticsSummary = bind_rows(somatic_summary(enrichedPilotSomatics), enrichedPilotSomaticsSummary)

  # GRAPHS
  prod_graphs = env_graph("Production ", sampleDetails$purity.prod, sampleDetails$ploidy.prod, sampleDetails$qcStatus.prod, sampleDetails$tumorCoverage, prodRegions, enrichedProdSomatics)
  pilot_graphs = env_graph("Pilot", sampleDetails$purity.pilot, sampleDetails$ploidy.pilot, sampleDetails$qcStatus.pilot, sampleDetails$tumorCoverage, pilotRegions, enrichedPilotSomatics)
  complete = plot_grid(prod_graphs, pilot_graphs, ncol = 1)
  save_plot(paste0("/Users/jon/hmf/analysis/fit/", sample, ".png"), complete, base_height = 10, base_width = 20)
}

dbDisconnect(dbProd); rm(dbProd)
dbDisconnect(dbPilot); rm(dbPilot)

#View(enrichedProdSomaticsSummary)
#View(enrichedPilotSomaticsSummary)

#inconsistentSomaticSummary = left_join(enrichedProdSomaticsSummary, enrichedPilotSomaticsSummary, by = "sampleId", suffix = c(".prod",".pilot"))
#save(inconsistentSomaticSummary, file = "~/hmf/analysis/fit/inconsistentSomaticSummary.RData")
#View(enrichedPilotSomatics)
somatic_graph(enrichedPilotSomatics)
somatic_graph(enrichedProdSomatics)
#enrichedPilotSomaticSummary = somatic_summary(enrichedPilotSomatics)


dev.off()
ggplot(enrichedPilotSomatics) + 
  stat_density(aes(x = somaticPloidy, fill = clonality, alpha = clonality), position = 'identity', bw = 0.05, n = 512) +
  scale_fill_manual(values = clonalityColoursColours) + 
  scale_alpha_manual(values = clonalityAlphas) 

#load(file = "~/hmf/analysis/fit/inconsistentSomaticSummary.RData")



rbinom(100, 1000, 0.5)


dbinom(15, 100, 0.1428)


clonalityColoursColours = setNames(c("#66c2a5","#8da0cb", "#fc8d62"), c("CLONAL", "SUBCLONAL","INCONSISTENT"))
clonalityAlphas = setNames(c(0.8, 0.8, 0.8), c("CLONAL", "SUBCLONAL","INCONSISTENT"))


complete

enrichedPilotSomatics %>% 
  mutate(ploidyBucket = round(somaticPloidy,0)) %>%
  group_by(ploidyBucket) %>% count()

sampleDetails

depth = 115

# Pilot
purity = 0.4
ploidy = 2.1327375
somatics = enrichedPilotSomatics
somatics %>% filter(somaticPloidy > 0.3, somaticPloidy < 2.8) %>%count()

#Prod
purity = 0.67
ploidy = 4.13
somatics = enrichedProdSomatics
somatics %>% mutate(roundPloidy = round(somaticPloidy, 0)) %>% group_by(roundPloidy) %>% count()

expectedVaf = purity / (ploidy * purity + 2 * (1-purity))
domain = seq(0, 4, 0.05)
range1 = dpois(round(domain * depth * expectedVaf,0), depth * expectedVaf) * expectedVaf * 0.05 *  depth * 1008
range2 = dpois(round(domain * depth * expectedVaf,0), 2 * depth * expectedVaf) * expectedVaf * 0.05 * depth * 1900
sum(range1) + sum(range2)

jon = data.frame(x = domain, y1 = range1, y2 = range2, y3 = range1 + range2) 

ggplot(somatics) +
  geom_histogram(aes(x = somaticPloidy), binwidth = 0.05) + 
  geom_point(data = jon, aes(x, y1), color  = "red") + 
  geom_point(data = jon, aes(x, y2), color  = "green") + 
  geom_point(data = jon, aes(x, y3), color  = "blue") + ggtitle("Prod") +
  xlim(0, 4)

sampleDetails
purity = 0.83
ploidy = 3.4
depth = 105
somatics = enrichedPilotSomatics

model = create_model(purity, ploidy, enrichedPilotSomatics)
sum(model$yTotal)



ggplot(somatics) +
  geom_histogram(aes(x = somaticPloidy), binwidth = 0.05) + 
  geom_point(data = model, aes(x, y1), color  = "red") + 
  geom_point(data = model, aes(x, y2), color  = "green") + 
  geom_point(data = model, aes(x, y3), color  = "orange") + 
  geom_point(data = model, aes(x, yTotal), color  = "blue") + ggtitle("Prod") +
  xlim(0, 4)





##### PILOT
purity = sampleDetails$purity.pilot
ploidy = sampleDetails$ploidy.pilot
depth = sampleDetails$tumorCoverage
somatics = enrichedPilotSomatics
binwidth = 0.1
somatics %>% count()



##### PROD
sampleDetails
purity = sampleDetails$purity.prod
ploidy = sampleDetails$ploidy.prod
depth = sampleDetails$tumorCoverage
somatics = enrichedProdSomatics
binwidth = 0.1
somatics %>% count()



