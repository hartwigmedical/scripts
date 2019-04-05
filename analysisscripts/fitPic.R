detach("package:purple", unload=TRUE)
library(purple)
library(RMySQL)
library(dplyr)
library(tidyr)
library(cowplot)
theme_set(theme_bw())


## Purity Range
range = read.table(file = "/Users/jon/hmf/analysis/COLO829T/purple/COLO829T.purple.purity.range", sep = "\t", header = T, comment.char = "!") %>%
  select(Purity = X.Purity, everything()) %>%
  #filter(Score < 3) %>%
  mutate(Score = pmin(4, Score))

bestPurity = range[1, "Purity"]
bestPloidy = range[1, "Ploidy"]
bestNormFactor = range[1, "NormFactor"]


range = range %>% 
#  mutate(Ploidy = round(Ploidy, 2)) %>% 
#  group_by(Purity) %>% 
  arrange(Purity, Ploidy) %>% 
  mutate(
    prevPloidy = lag(Ploidy), 
    width = abs((prevPloidy - Ploidy)) + 0.01,
    PloidyUpper = lag(Ploidy + width)
    ) %>%
  filter(Ploidy < 5)
range[is.na(range)] <- 0  
  
yMin = min(range$Purity)
yMax = max(range$Purity)
xMin = min(range$Ploidy)
xMax = max(range$Ploidy) 

ggplot(range, aes( x = Ploidy, y = Purity)) +
  geom_tile(aes(fill=Score, width = width), height = 0.01) + 
  scale_fill_gradientn(colours=rev(rainbow(1000, start=0, end=0.2))) + 
  #coord_cartesian(xlim = c(1.2, 5), ylim = c(0.08, 1)) + 
  #scale_y_continuous(expand=c(0.2, 0.2)) +
  #ylim(yMin, yMax) + 
  #geom_segment(aes(x = xMin, xend = xMax, y=bestPloidy, yend = bestPloidy), linetype = "dashed") +
  #geom_segment(aes(x = bestPurity, xend = bestPurity, y=yMin, yend = yMax), linetype = "dashed") + 
  #geom_label(data = data.frame(), aes(x = bestPurity, y = 1.2, label = bestPurity), size = 3) +
  #geom_label(data = data.frame(), aes(x = xMin, y = bestPloidy, label = round(bestPloidy,2)), size = 3, hjust = 0) + 
  theme(panel.grid.minor = element_blank()) + 
  ggtitle("Purity Range") 

xMin = min(range$Purity)
xMax = max(range$Purity) + 0.05
yMin = min(range$NormFactor) - 0.05
yMax = max(range$NormFactor) 

ggplot(range, aes(x = Purity, y = NormFactor))+
  geom_raster(aes(fill=Score)) +
  geom_segment(aes(x = xMin, xend = xMax, y=bestNormFactor, yend = bestNormFactor), linetype = "dashed") +
  geom_segment(aes(x = bestPurity, xend = bestPurity, y=yMin, yend = yMax), linetype = "dashed") +
#geom_text(aes(x = xMin, y = bestNormFactor, label = paste0("Ploidy: ", round(bestPloidy, 2)), hjust = 0, vjust = 0), nudge_y = 0.05) +
  #annotate("sdf", label = "sdf", aes(x = bestPurity, y = yMin)) +
  #scale_fill_gradientn(
  #  colours=c("#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"),
  #  values = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1),
  #  breaks = c(0, 0.5, 1, 2, 3), labels = c(0, 0.5, 1, 2, "3+")) +
  #scale_fill_gradientn(colours=rev(heat.colors(1000))) +
  scale_fill_gradientn(colours=rev(rainbow(1000, start=0, end=0.2))) +
  ggtitle(paste0("Purity Range")) + ylab("Norm Factor") +
  geom_label(data = data.frame(), aes(x = bestPurity, y = yMin, label = bestPurity), size = 3) +
  geom_label(data = data.frame(), aes(x = xMax, y = bestNormFactor, label = round(bestPloidy,2)), size = 3) + 
  theme(panel.grid.minor = element_blank())
  


### FITTED REGIONS

#copyNumberRegions = read.table(file = "/Users/jon/hmf/analysis/COLO829T/purple/COLO829T.purple.cnv", sep = "\t", header = T, comment.char = "!")

fittedRegions = read.table(file = "/Users/jon/hmf/analysis/COLO829T/purple/COLO829T.purple.fitted", sep = "\t", header = T, comment.char = "!") %>% 
  filter(germlineStatus == "DIPLOID", bafCount > 0) %>%
  mutate(Score = deviationPenalty * eventPenalty, Weight = bafCount)

ggplot(fittedRegions, aes(x=majorAllelePloidy,y=minorAllelePloidy)) + 
  geom_point(aes(size = Weight, color = Score), alpha = 0.8) +
  xlab("Major Allele") + ylab("Minor Allele") + ggtitle("Fitted Segments") +
  scale_x_continuous(breaks = c(-200:200), limits = c(-0.1, 4)) +
  scale_y_continuous(breaks = c(-200:200), limits = c(-0.1, 3)) +
  #scale_color_gradientn(colours=rev(heat.colors(10))) +
  scale_color_gradientn(colours=rev(rainbow(1000, start=0, end=0.2))) +
  theme(panel.grid.minor = element_blank())  + 
  scale_size(range = c(1,10), guide = "none")

  p = ggplot(simpleCopyNumberRegions, aes(x=majorAllele,y=minorAllele)) +
  geom_point(aes(size = bafCount), alpha = 0.4) +
  #geom_point(aes(size = bafCount)) +
  xlab("Major Allele") + ylab("Minor Allele") + ggtitle(title) +
  scale_x_continuous(breaks = c(-200:200), limits = c(-0.1, 4)) +
  scale_y_continuous(breaks = c(-200:200), limits = c(-0.1, 3)) +
  scale_color_gradientn(colours=rev(rainbow(1000, start=0, end=0.75))) +
  theme(panel.grid.minor = element_blank()) +
  scale_size(range = c(1,10))


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
    majorAlleleVaf = pmax(majorAllelePloidy, 0) *  subclonalProbability,
    minConceivableAlleleCount = qbinom(0.01, 100000, totalReadCount * subclonalProbability / 100000, T),
    maxConceivableAlleleCount = qbinom(0.99, 100000, totalReadCount * majorAlleleVaf / 100000, T),
    maxConceivableVaf = pmin(1, maxConceivableAlleleCount / totalReadCount),
    maxConceivablePloidy = (maxConceivableVaf * 2 * (1 - purity)) / (purity * (1 - maxConceivableVaf)),
    somaticPenalty = pmax(0, somaticPloidy - maxConceivablePloidy),
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

histo = myHist
find_peak_from_histogram <- function(histo, binwidth) {
  histo = histo %>% filter(n > 0)

  peak_y = max(histo$n)
  peak_x = first(histo$bucket[histo$n == peak_y])

  histoSomatics =  rep(histo$bucket, histo$n)
  d = density(histoSomatics, bw = binwidth, from = (binwidth/2), to = (10 + binwidth/2), n = (10 / binwidth + 1))
  peak_y =  max(d$y)
  peak_x = first(d$x[d$y == peak_y])
  peak_y =  max(d$y) * binwidth * sum(histo$n)
  return (c(x = peak_x, y = peak_y))
}

create_initial_histogram <- function(somatics, binwidth) {
  domain = seq(0 + binwidth / 2, 10 + binwidth / 2, binwidth)

  myInitialSomatics = somatics %>% filter(somaticPloidy > (binwidth / 2), somaticPloidy < (10 + binwidth / 2))
  myEmptyHist = data.frame(bucket = domain, n = 0)

  myInitialHist = myInitialSomatics %>%
    mutate(
    bucketFactor = cut(somaticPloidy, breaks = domain),
    bucket = (as.numeric(bucketFactor) - 1) * binwidth + binwidth / 2
    ) %>%
    group_by(bucket) %>%
    count() %>%
    bind_rows(myEmptyHist) %>%
    group_by(bucket) %>%
    summarise(n = sum(n)) %>%
    arrange(bucket)

  return (myInitialHist)
}

create_somatic_model <- function(purity, ploidy, depth, somatics, binwidth) {
  domain = seq(0 + binwidth / 2, 10 + binwidth / 2, binwidth)

  expectedVaf = purity / (ploidy * purity + 2 * (1-purity))
  myInitialHist = create_initial_histogram(somatics, binwidth)
  myHist = myInitialHist

  initialSomatics = sum(myHist$n)
  result = data.frame(x = (myInitialHist$bucket + binwidth/2), stringsAsFactors = F)

  i = 1
  exitLoop = F
  while (!exitLoop) {
    modelPeak = find_peak_from_histogram(myHist,  binwidth)

    peak_x = modelPeak["x"] + binwidth / 2
    modelSomatics = dpois(round(domain * depth * expectedVaf,0), peak_x * depth * expectedVaf)
    modelSomatics = round(modelSomatics * modelPeak["y"] / max(modelSomatics), 2)
    peakLabel = paste0(peak_x)
    alreadyExists = peakLabel %in% colnames(result)

    if (!alreadyExists) {
      initialHeightOfPeak = myInitialHist[myInitialHist$bucket == modelPeak["x"], ]$n
      cat(modelPeak["y"] / initialHeightOfPeak)
      if (initialHeightOfPeak > 0 & modelPeak["y"] / initialHeightOfPeak > 0.2) {
        result[, peakLabel] <- modelSomatics
      }

      myHist$n = pmax(0, myHist$n - modelSomatics)
      unexplainedSomatics = sum(myHist$n)
      cat(i, " ", modelPeak["x"], " ", unexplainedSomatics / initialSomatics, " ",  modelPeak["y"] / initialHeightOfPeak, "\n")
    }
    i  = i + 1
    exitLoop = unexplainedSomatics / initialSomatics < 0.05 || i > 10 || alreadyExists || !any(myHist$n > 1)
  }

  if (ncol(result) > 3) {
    result$total = rowSums(result[, -1])
  }

  return (result)
}

somatic_graph <- function(purity, ploidy, depth, somatics, binwidth = 0.05) {

  clonalityColoursColours = setNames(c("#66c2a5","#8da0cb", "#fc8d62"), c("CLONAL", "SUBCLONAL","INCONSISTENT"))
  clonalityAlphas = setNames(c(0.8, 0.8, 0.8), c("CLONAL", "SUBCLONAL","INCONSISTENT"))

  somaticModel = create_somatic_model(purity, ploidy, depth, somatics, binwidth)
  tidySomaticModel = somaticModel %>% gather(peak, value, -1) %>% mutate(peak = factor(peak, colnames(somaticModel), ordered = T))

  p1 = ggplot(somatics) +
    geom_histogram(aes(x = somaticPloidy, fill = CN), alpha = 0.4, binwidth = binwidth, position = "stack") +
    xlab("Ploidy") + ggtitle("Somatic CNV") + ylab("") +
    theme(legend.position = "right") +
  #theme(legend.position = c(0.9,0.77)) +
    geom_line(data = tidySomaticModel, aes(x, value, color = peak)) +
    theme(panel.border = element_blank(), panel.grid.minor.y = element_blank(), axis.ticks = element_blank()) +
    scale_x_continuous(limits = c(0, 5), breaks = c(1:10))

  #+ theme(legend.position = "none")

  p2 = ggplot(somatics) +
    geom_histogram(aes(x = somaticPloidy, fill = clonality, color = clonality, alpha = clonality), binwidth = 0.1, position = "identity") +
  #geom_density(bw = 0.03, aes(x = somaticPloidy,  group = clonality, y=0.1 * ..count..), position = "identity") +
    scale_fill_manual(name = "Clonality", values = clonalityColoursColours) +
    scale_color_manual(name = "Clonality", values = clonalityColoursColours) +
    scale_alpha_manual(values = clonalityAlphas) +
    xlim(0, 4) + xlab("Ploidy") + ggtitle("Somatic Clonality") + ylab("") +
    theme(legend.position = "none") +
    theme(panel.border = element_blank(), panel.grid.minor.y = element_blank(), axis.ticks = element_blank())

  p3 = plot_grid(p1, p2, nrow = 1, rel_widths = c(2, 1))
  return (p3)
}

#dev.off()
#somatic_graph(sampleDetails$purity.prod, sampleDetails$ploidy.prod,  sampleDetails$tumorCoverage, enrichedProdSomatics)
#somatic_graph(sampleDetails$purity.prod, sampleDetails$ploidy.prod, sampleDetails$tumorCoverage, enrichedProdSomatics)
#somatic_graph(sampleDetails$purity.pilot, sampleDetails$ploidy.pilot, sampleDetails$tumorCoverage, enrichedPilotSomatics)
#sum(model$y1)


copynumber_graph <- function(copyNumberRegions) {
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

##### ENVIRONMENT SETTINGS
#outputDir = "/Users/jon/hmf/analysis/lowcoverage/"
#pilotDbName = "low_coverages"
outputDir = "/Users/jon/hmf/analysis/fit/"
pilotDbName = "hmfpatients_pilot"

### PRODUCTION COHORT
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
prodCohort = dbGetQuery(dbProd, "SELECT sampleId, purity, ploidy, status, qcStatus FROM purity")
prodMetrics = dbGetQuery(dbProd, "SELECT sampleId, round(tumorMeanCoverage,0) as tumorCoverage FROM metric")
dbDisconnect(dbProd); rm(dbProd)

### PILOT COHORT
dbPilot = dbConnect(MySQL(), dbname=pilotDbName, groups="RAnalysisWrite")
pilotCohort = dbGetQuery(dbPilot, "SELECT sampleId, purity, ploidy, status, qcStatus FROM purity WHERE version = 2.15 and modified > '2018-08-23 09:19:04'")
pilotCohort = dbGetQuery(dbPilot, "SELECT sampleId, purity, ploidy, status, qcStatus FROM purity")
pilotCohort = dbGetQuery(dbPilot, "SELECT pilot.sampleId, pilot.purity, pilot.ploidy, pilot.status, pilot.qcStatus from hmfpatients_pilot.purity pilot, hmfpatients.purity prod
                         WHERE pilot.sampleId = prod.sampleId 
                         AND (pilot.qcStatus = 'PASS' OR prod.qcStatus = 'PASS')
                         and (
                         (pilot.qcstatus = 'FAIL_DELETED_GENES' OR prod.qcstatus = 'FAIL_DELETED_GENES') 
                         OR (pilot.maxDiploidProportion >= 0.95 AND pilot.status <> prod.status)
                         OR (pilot.purity > prod.purity + 0.1 AND pilot.ploidy > prod.ploidy + 0.2) 
                         OR (pilot.purity > 0.97 and pilot.status <> 'NO_TUMOR') 	
                         )
                         order by abs(prod.ploidy - pilot.ploidy) desc")

dbDisconnect(dbPilot); rm(dbPilot)

#### COMPARE COHORTS
cohort = left_join(pilotCohort, prodCohort, by = "sampleId", suffix = c(".pilot", ".prod")) %>%
  mutate(
  changeInPurity = abs(purity.prod - purity.pilot),
  changeInPloidy = abs(ploidy.prod - ploidy.pilot)) %>%
  arrange(-changeInPloidy)
cohort = left_join(cohort, prodMetrics, by = "sampleId")
cohortToRun = cohort
#cohortToRun = cohort %>% filter(changeInPloidy > 0.5)
#cohortToRun = cohort %>% top_n(100, changeInPloidy) %>% mutate(patientId = substr(sampleId, 1, 12)) #%>% pull(patientId)
save(cohort, cohortToRun, file = paste0(outputDir, "cohort.RData"))

###### EXECUTE FROM HERE
load(file = paste0(outputDir, "cohort.RData"))
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
dbPilot = dbConnect(MySQL(), dbname=pilotDbName, groups="RAnalysisWrite")

# Initialistaion
enrichedPilotSomaticsSummary = data.frame()
enrichedProdSomaticsSummary = data.frame()

#sample = "CPCT02380013T"

for (sample in cohortToRun$sampleId) {
  cat("Processing", sample, "\n")
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
  prodStatus = paste0(sampleDetails$qcStatus.prod,"|",sampleDetails$status.prod)
  pilotStatus = paste0(sampleDetails$qcStatus.pilot,"|",sampleDetails$status.pilot)

  prod_graphs = env_graph("Production ", sampleDetails$purity.prod, sampleDetails$ploidy.prod, prodStatus, sampleDetails$tumorCoverage, prodRegions, enrichedProdSomatics)
  pilot_graphs = env_graph("Pilot", sampleDetails$purity.pilot, sampleDetails$ploidy.pilot, pilotStatus, sampleDetails$tumorCoverage, pilotRegions, enrichedPilotSomatics)
  complete = plot_grid(prod_graphs, pilot_graphs, ncol = 1)
  save_plot(paste0(outputDir, sample, ".png"), complete, base_height = 10, base_width = 20)
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


#### SAVE COHORT
#save(cohortToRun, cohort, file = "~/hmf/analysis/fit/cohort.RData")
#load(file = "~/hmf/analysis/fit/cohort.RData")


##### PILOT
purity = sampleDetails$purity.pilot
ploidy = sampleDetails$ploidy.pilot
depth = sampleDetails$tumorCoverage
somatics = enrichedPilotSomatics
binwidth = 0.05
somatics %>% count()



##### PROD
sampleDetails
purity = sampleDetails$purity.prod
ploidy = sampleDetails$ploidy.prod
depth = sampleDetails$tumorCoverage
somatics = enrichedProdSomatics
binwidth = 0.05
somatics %>% count()

jon = enrichedPilotSomatics %>% filter(inconsistent) %>%
  select(sampleId, chromosome, position, alleleReadCount, totalReadCount, majorAllelePloidy, somaticPloidy, adjustedCopyNumber, maxConceivableAlleleCount, majorAlleleProbability)


jon = somatic_enrichment(sampleDetails$purity.pilot, pilotSomatics)

sampleDetails$purity.pilot
sampleDetails

jon[1, c("chromosome", "adjustedCopyNumber", "alleleReadCount","totalReadCount", "somaticPloidy", "majorAllelePloidy", "majorAlleleVaf", "maxConceivableAlleleCount", "maxConceivableVaf", "maxConceivablePloidy")]

head(jon)
sum(jon$somaticPenalty) / nrow(jon)


jon %>% filter(bafC)


