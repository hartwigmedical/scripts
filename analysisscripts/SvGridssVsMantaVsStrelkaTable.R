library(RMySQL)

#probes = read.csv(file = "/Users/jon/hmf/analysis/mantaVgridss/probes.csv")

#load(file = "~/hmf/analysis/probes/allOverlaps.RData")
#load(file = "~/hmf/analysis/probes/gridssSingles.RData")

load(file = "~/hmf/analysis/probes/probeResult.RData")
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis", host = "127.0.0.1")
patientMapping = dbGetQuery(dbProd, "SELECT sampleId, hmfId from sampleMapping")
dbDisconnect(dbProd)
rm(dbProd)

unique(probeResult$callset)
unique(probeResult$scope)

probeResult %>% filter(callset == 'Gridss') %>% group_by(scope)  %>% summarise(n = n()) 
probeResult %>% filter(callset == 'Strelka') %>% group_by(scope)  %>% summarise(n = n()) 

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
