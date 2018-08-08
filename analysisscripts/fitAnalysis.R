detach("package:purple", unload=TRUE)
library(purple)
library(RMySQL)
library(dplyr)
library(tidyr)
library(cowplot)
theme_set(theme_bw())

load(file = "~/hmf/RData/Reference/cancerTypeColours.RData")
chromosomeColours = c(cancerTypeColours, "black")
chromosomeColours = setNames(chromosomeColours, c(1:22))


query_regions <- function(dbConnect, sampleId) {
  query = paste0("select  * from copyNumber where sampleId ='", sampleId, "'");
  return (dbGetQuery(dbConnect, query))
}

query_somatics <- function(dbConnect, sampleId) {
  query = paste0("select  * from somaticVariant where sampleId ='", sampleId, "'");
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

somatic_graph <- function(somatics) {
  
  inconsistentColours = setNames(c("#6baed6", "#cb181d"), c(F,T))
  inconsistentAlphas = setNames(c(0.4, 1), c(F,T))
  
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
      inconsistent = somaticPloidy > majorAllelePloidy,
      CN = round(adjustedCopyNumber),
      CN = ifelse(CN > 5, "6+", CN),
      CN = factor(CN, levels = c("6+", rev(-1:5)), ordered = T))
  
  p1 = ggplot(somatics) +
    geom_histogram(aes(x = somaticPloidy, fill = CN), binwidth = 0.1, position = "stack") +
    xlim(0, 4) + xlab("Ploidy") + ggtitle("Somatics by CN") + ylab("") +
    theme(legend.position = c(0.9,0.77)) + 
    theme(panel.border = element_blank(), panel.grid.minor.y = element_blank(), axis.ticks = element_blank())
  
  #+ theme(legend.position = "none")
  
  p2 = ggplot(somatics) +
    geom_histogram(aes(x = somaticPloidy, fill = inconsistent, alpha = inconsistent), binwidth = 0.1, position = "stack") +
    scale_fill_manual(values = inconsistentColours) +
    scale_alpha_manual(values = inconsistentAlphas) +
    xlim(0, 4) + xlab("Ploidy") + ggtitle("Inconsistent Somatics") + ylab("") +
    theme(legend.position = "none", axis.text.y = element_blank()) +
    theme(panel.border = element_blank(), panel.grid.minor.y = element_blank(), axis.ticks = element_blank())
   
  p3 = plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 1))
  return (p3)
}

copynumber_graph <- function(copyNumberRegions) {
  
  totalBafCount = sum(copyNumberRegions$bafCount)
  
  copyNumberRegions = copyNumberRegions %>%
    filter(!chromosome %in% c('X','Y'), copyNumber < 10) %>%
    mutate(
      chromosome = factor(chromosome, levels= c(1:22), ordered = T),
      weight = bafCount/totalBafCount ) 
  
  ggplot(copyNumberRegions) +
    geom_histogram(aes(x = copyNumber, fill = chromosome, weight = weight), binwidth = 0.1, position = "stack") +
    scale_fill_manual(values = chromosomeColours) + 
    scale_x_continuous(breaks = c(0:10)) + 
    theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
    xlab("Copy Number") + ylab("Baf Weighed Count")
}

env_graph <- function(title, copyNumberRegions, somatics) {
  fitted_graph = fitted_region_graph(title, copyNumberRegions)
  somatic_graph = somatic_graph(somatics)
  copynumber_graph = copynumber_graph(copyNumberRegions)
  return (plot_grid(fitted_graph, somatic_graph, copynumber_graph, rel_widths = c(2,2, 1.5),  nrow = 1))
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

cohortToRun = cohort %>% top_n(100, changeInPloidy) %>% mutate(patientId = substr(sampleId, 1, 12)) #%>% pull(patientId)
sample = "CPCT02050047T"

for (sample in cohortToRun$sampleId) {
  sampleDetails = cohort %>% filter(sampleId == sample)
  
  prodSomatics = query_somatics(dbProd, sample)
  prodRegions = query_regions(dbProd, sample)
  pilotSomatics = query_somatics(dbPilot, sample)
  pilotRegions = query_regions(dbPilot, sample)
  
  prod_graphs = env_graph(create_title("Production", sampleDetails$purity.prod, sampleDetails$ploidy.prod, sampleDetails$qcStatus.prod), prodRegions, prodSomatics)
  pilot_graphs = env_graph(create_title("Pilot", sampleDetails$purity.pilot, sampleDetails$ploidy.pilot, sampleDetails$qcStatus.pilot), pilotRegions, pilotSomatics)
  complete = plot_grid(prod_graphs, pilot_graphs, ncol = 1)
  complete
  
  save_plot(paste0("/Users/jon/hmf/analysis/fit/", sample, ".png"), complete, base_height = 10, base_width = 20)
}



dbDisconnect(dbProd);rm(dbProd)
dbDisconnect(dbPilot); rm(dbPilot)









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
