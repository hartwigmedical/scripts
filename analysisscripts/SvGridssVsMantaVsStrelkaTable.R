library(RMySQL)

#probes = read.csv(file = "/Users/jon/hmf/analysis/mantaVgridss/probes.csv")

#load(file = "~/hmf/analysis/probes/allOverlaps.RData")
#load(file = "~/hmf/analysis/probes/gridssSingles.RData")

load(file = "~/hmf/analysis/probes/probeResult.RData")
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis", host = "127.0.0.1")
patientMapping = dbGetQuery(dbProd, "SELECT sampleId, hmfId from sampleMapping")
dbDisconnect(dbProd)
rm(dbProd)

supplementaryTable =  left_join(probeResult, patientMapping, by = "sampleId") %>%
  mutate(
    gridss = callset == 'Gridss', 
    manta = callset == 'Manta' | scope == "SharedManta" | scope == "SharedBoth",
    strelka = callset == 'Strelka' | scope == "SharedStrelka" | scope == "SharedBoth") %>%
  ungroup() %>%
  select(sample = hmfId,gridss, manta, strelka, type,
         startChromosome, startPosition, startOrientation, 
         endChromosome, endPosition, endOrientation, insertSequence, 
         sampleRefDepth = sourceRefDepth, sampleAltDepth = sourceAltDepth, 
         sumControlRefDepth = otherRefDepth, sumControlAltDepth = otherAltDepth, maxControlAltDepth = otherAltDepthMax, 
         probeQuality, p, supported)


write.csv(supplementaryTable, file = "/Users/jon/hmf/analysis/mantaVgridss/probeSupplementaryTable.csv", row.names = F)




load(file = "/Users/jon/hmf/analysis/mantaVgridss/scopedData.RData")
library(BSgenome.Hsapiens.UCSC.hg19)
library(seqinr)
dups = gridss %>% filter(type == 'DUP', sampleId == 'XXXXX') %>%
  mutate(len = endPosition - startPosition + 1) %>% filter (len <= 200 ) %>%
  mutate(
    dup = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, paste0("chr", startChromosome), startPosition, endPosition)),
    repeatContext = substr(dup, 1, 2), 
    repeatCount = floor(nchar(dup) / 2),
    tmp = strrep(repeatContext, 12),
    repeatCount = ifelse(substr(dup, 1, nchar(tmp)) == tmp, repeatCount, 0),
    msi = repeatCount > 0) %>%
  select(-tmp, -dup) %>%
  filter(msi)

dupSummary = dups %>% mutate(isRepeated = repeatCount > 0) %>% group_by(sampleId, isRepeated) %>% count() %>% spread(isRepeated, n)
dupSummary

probeMsi = dups
save(probeMsi, file = "~/hmf/tmp/probeMsi.RData")
