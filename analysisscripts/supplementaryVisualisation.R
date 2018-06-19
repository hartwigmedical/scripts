detach("package:purple", unload=TRUE)
library(purple)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(scales)
theme_set(theme_bw())



load("~/hmf/RData/processed/driverGenes.RData")
load(file = '~/hmf/RData/Reference/hpcCancerTypeCounts.RData')
load("~/hmf/RData/processed/hpcDriversByGene.RData")
load(file = "~/hmf/RData/Processed/highestPurityCohortSummary.RData")
highestPurityCohortSummary[is.na(highestPurityCohortSummary)] <- 0


########## Driver rates by type in WGD vs non WGD
wgdDriverSampleCount = highestPurityCohortSummary %>% group_by(WGD) %>% summarise(total=n())
wgdDriverRates = hpcDriversByGene %>%
  filter(driverLikelihood > 0) %>%
  mutate(
    driver = as.character(driver),
    driver = ifelse(substr(driver, 1, 6) == "Fusion", "Fusion", driver),
    driver = ifelse(driver %in% c("Frameshift", "Inframe"), "Indel", driver),
    driver = ifelse(driver %in% c("FragileDel"), "Del", driver),
    driver = ifelse(driver %in% c("Missense","Nonsense","Splice"), "SNV/MNV", driver)
    ) %>%
  group_by(sampleId, driver) %>%
  summarise(n = sum(driverLikelihood)) %>%
  #count() %>%
  left_join(highestPurityCohortSummary %>% select(sampleId, WGD), by = "sampleId") %>%
  group_by(driver, WGD) %>%
  summarise(n = sum(n)) %>%
  left_join(wgdDriverSampleCount, by = "WGD") %>%
  mutate(n = n / total) %>%
  ungroup()
wgdDriverLevels = wgdDriverRates %>% filter(WGD) %>% arrange(-n)
wgdDriverRates = mutate(wgdDriverRates, driver = factor(driver, wgdDriverLevels$driver))

p1 = ggplot(data = wgdDriverRates, aes(x = driver, y = n)) +
  geom_bar(aes(fill = WGD), stat = "identity", position="dodge") +
  scale_fill_manual(values = c("#f1a340", "#998ec3")) +
  ggtitle("Drivers Per Sample by WGD") + 
  xlab("Driver") + ylab("Drivers") +
  coord_flip() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.ticks = element_blank(), legend.position="bottom") + 
  scale_y_continuous(expand=c(0.01, 0.01))

plot_grid(p1, labels="AUTO")


