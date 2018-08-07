detach("package:purple", unload=TRUE)
library(purple)
library(RMySQL)
library(dplyr)
library(tidyr)
library(cowplot)
theme_set(theme_bw())


#library(GenomicRanges)

query_regions <- function(dbConnect, sampleId) {
  query = paste0("select  * from copyNumberRegion where sampleId ='", sampleId, "'");
  return (dbGetQuery(dbConnect, query))
}

query_somatics <- function(dbConnect, sampleId) {
  query = paste0("select  * from somaticVariant where sampleId ='", sampleId, "'");
  return (dbGetQuery(dbConnect, query))
}



### PRODUCTION COHORT
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
prodCohort = dbGetQuery(dbProd, "SELECT sampleId, purity, ploidy, status, qcStatus FROM purity")
dbDisconnect(dbProd)
rm(dbProd)

### PILOT COHORT
dbPilot = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
pilotCohort = dbGetQuery(dbPilot, "SELECT sampleId, purity, ploidy, status, qcStatus FROM purity WHERE version = 2.15 and modified > '2018-08-04'")
dbDisconnect(dbPilot)
rm(dbPilot)

#### COMPARE COHORTS
cohort = left_join(pilotCohort, prodCohort, by = "sampleId", suffix = c(".pilot", ".prod")) %>%
  mutate(
  changeInPurity = abs(purity.prod - purity.pilot),
  changeInPloidy = abs(ploidy.prod - ploidy.pilot)) %>%
  arrange(-changeInPloidy)


dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
dbPilot = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")

sample = "CPCT02050047T"
sampleDetails = cohort %>% filter(sampleId == sample)

prodSomatics = query_somatics(dbProd, sample)
prodRegions = query_regions(dbProd, sample)
pilotSomatics = query_somatics(dbPilot, sample)
pilotRegions = query_regions(dbPilot, sample)

dbDisconnect(dbProd);rm(dbProd)
dbDisconnect(dbPilot); rm(dbPilot)

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
    observedTumorRatio = round(observedTumorRatio,2),
    observedNormalRatio = round(observedNormalRatio,2),
    actualTumorCopyNumber,
    majorAllele = round(actualTumorCopyNumber * actualTumorBaf, 2),
    minorAllele = round(actualTumorCopyNumber - actualTumorCopyNumber * actualTumorBaf, 2),
    ploidyPenalty,
    deviation = (cnvDeviation + bafDeviation),
    score = (cnvDeviation + bafDeviation) * observedBaf * ploidyPenalty * bafCount / totalBafCount,
    totalDeviation
    ) %>%
    arrange(score)

  p = ggplot(simpleCopyNumberRegions, aes(x=majorAllele,y=minorAllele,z=deviation, color=score)) +
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

fitted_region_graph("Jon", pilotRegions)

somatic_graph <- function(somatics) {
  somatics = somatics %>%
    filter(!chromosome %in% c('X','Y')) %>%
    mutate(
    majorAllelePloidy = adjustedCopyNumber - minorAllelePloidy,
    somaticPloidy = adjustedVaf * adjustedCopyNumber,
    distanceFromMajorAllele = somaticPloidy - majorAllelePloidy,
    distanceFromMinorAllele = somaticPloidy - minorAllelePloidy,
    distanceFromOne = somaticPloidy - 1,
    distanceFromAllele = ifelse(abs(distanceFromMajorAllele) < abs(distanceFromMinorAllele), distanceFromMajorAllele, distanceFromMinorAllele),
    distanceFromAllele = ifelse(abs(distanceFromOne) < abs(distanceFromMajorAllele), distanceFromOne, distanceFromAllele),
    distanceFromWholePloidy = somaticPloidy - round(somaticPloidy),
    CN = round(adjustedCopyNumber),
    CN = ifelse(CN > 5, "6+", CN),
    CN = factor(CN, levels = c("6+", rev(-1:5)), ordered = T))

  p1 = ggplot(somatics) +
    geom_histogram(aes(x = somaticPloidy, fill = CN), binwidth = 0.1, position = "stack") +
    xlim(0, 4) + xlab("Somatic Ploidy") + theme(legend.position = "none")

  p2 = ggplot(somatics) +
    geom_histogram(aes(x = distanceFromAllele, fill = CN), binwidth = 0.1, position = "stack") +
    xlim(-5, 5) + xlab("Distance from Major or Minor Allele")


  p3 = plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 1.2))
  return (p3)
}


env_graph <- function(title, copyNumberRegions, somatics) {
  fitted_graph = fitted_region_graph(title, copyNumberRegions)
  somatic_graph = somatic_graph(somatics)
  return (plot_grid(fitted_graph, somatic_graph, nrow = 1))
}

prod_graphs = env_graph(paste0(sample, " Production (Purity:", round(100*sampleDetails$purity.prod,2), "% Ploidy:", round(sampleDetails$ploidy.prod,2), ")"), prodRegions, prodSomatics)
pilot_graphs = env_graph(paste0(sample, " Pilot (Purity:", round(100*sampleDetails$purity.pilot,2), "% Ploidy:", round(sampleDetails$ploidy.pilot,2), ")"), pilotRegions, pilotSomatics)
complete = plot_grid(prod_graphs, pilot_graphs, ncol = 1)
complete

somatics = prodSomatics
jon = somatics %>% filter(CN == 4)

save_plot(paste0("/Users/jon/hmf/analysis/fit/", sample, ".png"), complete, base_height = 10, base_width = 20)


####### TESTING
copyNumberRegions = prodRegions
jon = somatics 


library(plotly)
p <- plot_ly(z = ~simpleCopyNumberRegions) %>% add_surface()
p

jon = simpleCopyNumberRegions %>% filter(majorAllele < 5, minorAllele < 3)

kd <- with(jon, MASS::kde2d (majorAllele, minorAllele, n = 100))
p <- plot_ly(x = kd$x, y = kd$y, z = kd$z) %>% add_surface()
p


jon =  simpleCopyNumberRegions %>% filter(minorAllele > 1.5,minorAllele < 1.9, majorAllele < 2.5, chromosome == 6 )
