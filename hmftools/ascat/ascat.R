#library(devtools)
#install_github("Crick-CancerGenomics/ascat/ASCAT")
library(ASCAT)
library(dplyr)

args <- commandArgs(trailing = T)
runDir <- args[1]
outDir <- args[2]
tumor <- args[3]
normal <- args[4]

#tumor = "COLO829"
#normal = "COLO829R"
#runDir = "/Users/jon/hmf/analysis/COLO829/"
#outDir = "/Users/jon/hmf/analysis/ascat/"

cobaltPath = paste0(runDir, "/cobalt/", tumor, ".cobalt")
amberPath = paste0(runDir, "/amber/", tumor, ".amber.baf")
cobalt = read.csv(file = cobaltPath, sep = '\t')
amber = read.csv(file = amberPath, sep = '\t')
rm(cobaltPath, amberPath)

gender = ifelse(nrow(amber[amber$Chromosome == 'X' & amber$Position > 2699520 & amber$Position < 155260560, ]) < 1000, "XY", "XX")
cat("Processing", tumor, "with gender", gender, "\n")

cobaltCleaned = cobalt %>% filter(ReferenceGCDiploidRatio > 0, TumorGCRatio > 0);
amberCleaned = amber %>% select(Chromosome, Position, NormalBAF, TumorBAF) %>% 
  mutate(Position = floor(Position / 1000)*1000 + 1) %>%
  group_by(Chromosome, Position) %>%
  summarise(NormalBAF = sum(NormalBAF)/n(), TumorBAF = sum(TumorBAF)/n())
rm(cobalt, amber)

## Array (needs at least one homozygous location)
randomArray = anti_join(cobaltCleaned, amberCleaned,by = c("Chromosome", "Position")) %>% select(Chromosome, Position) %>% sample_n(1)
amberArray = amberCleaned %>% select(Chromosome, Position) %>% inner_join(cobaltCleaned %>% select(Chromosome, Position), by = c("Chromosome", "Position")) %>% ungroup()
array = full_join(randomArray, amberArray,  by = c("Chromosome", "Position")) %>% 
  mutate(Chromosome = factor(Chromosome, levels = c(1:22, 'X', 'Y'))) %>% 
  arrange(Chromosome, Position) %>%
  mutate(Chromosome = as.character(Chromosome))
rm(randomArray, amberArray)

## Generate Input
germlineBAF = array %>% left_join(amberCleaned, by = c("Chromosome", "Position")) %>% select(Chromosome, Position, NormalBAF)
germlineBAF[is.na(germlineBAF)] <- 0

tumorBAF = array %>% left_join(amberCleaned, by = c("Chromosome", "Position")) %>% select(Chromosome, Position, TumorBAF)
tumorBAF[is.na(tumorBAF)] <- 0

germlineRatio = cobaltCleaned %>% select(Chromosome, Position, ReferenceGCDiploidRatio) %>% 
  mutate(ReferenceGCDiploidRatio = log2(pmax(ReferenceGCDiploidRatio, 0.001))) %>%
  filter(!is.na(ReferenceGCDiploidRatio), !is.infinite(ReferenceGCDiploidRatio)) %>% 
  inner_join(array, by = c("Chromosome", "Position")) 
tumorRatio = cobaltCleaned %>% select(Chromosome, Position, TumorGCRatio) %>%   
  mutate(TumorGCRatio = log2(pmax(TumorGCRatio, 0.0001))) %>%
  filter(!is.na(TumorGCRatio), !is.infinite(TumorGCRatio)) %>% 
  inner_join(array, by = c("Chromosome", "Position")) 

names(germlineBAF) <- c("Chromosome", "Position", tumor)
names(germlineRatio) <- c("Chromosome", "Position", tumor)
names(tumorBAF) <- c("Chromosome", "Position", tumor)
names(tumorRatio) <- c("Chromosome", "Position", tumor)

## Persist Input
germlineRatioPath = paste0(outDir, "/input/", normal, "_LogR.txt")
tumorRatioPath = paste0(outDir, "/input/", tumor, "_LogR.txt")
germlineBAFPath = paste0(outDir, "/input/", normal, "_BAF.txt")
tumorBAFPath = paste0(outDir, "/input/", tumor, "_BAF.txt")

write.table(germlineRatio, file = germlineRatioPath, sep = '\t', quote = F, row.names = T, col.names = T)
write.table(tumorRatio, file = tumorRatioPath, sep = '\t', quote = F, row.names = T, col.names = T)
write.table(germlineBAF, file = germlineBAFPath, sep = '\t', quote = F, row.names = T, col.names = T)
write.table(tumorBAF, file = tumorBAFPath, sep = '\t', quote = F, row.names = T, col.names = T)

## Execute
ascat.bc = ascat.loadData(tumorRatioPath, tumorBAFPath, germlineRatioPath, germlineBAFPath, gender = gender)
ascat.bc = ascat.aspcf(ascat.bc,  out.dir = paste0(outDir, "/input/"))
ascat.output = ascat.runAscat(ascat.bc, gamma = 1, img.dir = paste0(outDir, "/output/"))
output = data.frame(tumor = tumor, purity = ascat.output$aberrantcellfraction, psi = ascat.output$psi, ploidy = ascat.output$ploidy, goodnessofFit = ascat.output$goodnessOfFit, stringsAsFactors = F)
write.table(output, file = paste0(outDir, "/output/", tumor, ".results.txt"), sep = '\t', quote = F, row.names = F, col.names = T )
write.table(ascat.output$segments, file = paste0(outDir, "/output/", tumor, ".segments.txt"), sep = '\t', quote = F, row.names = F, col.names = T )
write.table(ascat.output$segments_raw, file = paste0(outDir, "/output/", tumor, ".segments_raw.txt"), sep = '\t', quote = F, row.names = F, col.names = T )
