library(dplyr)

allDbFusions <- read.csv("~/data/r/allDbFusions.csv", stringsAsFactors=FALSE)
allKnownFusions <- read.csv("~/data/r/allKnownFusions.csv",stringsAsFactors=FALSE)
knownPromiscuousH <- allKnownFusions %>% filter(is.na(T_gene)) %>% select(H_gene)
knownPromiscuousT <- allKnownFusions %>% filter(is.na(H_gene)) %>% select(T_gene)

dbFusionPairs <- allDbFusions %>% mutate(H_gene = upGene, T_gene = downGene, H_exon = upExonRank, T_exon = downExonRank) %>% 
  mutate (H_transcript = upTranscript, T_transcript = downTranscript, H_canonical = upIsCanonical, T_canonical = downIsCanonical) %>% 
  group_by(sampleId, H_gene, T_gene) %>% arrange(-H_canonical, -T_canonical, -upExonMax, -downExonMax) %>% filter(row_number() == 1) %>%
  select (sampleId, H_gene, T_gene, H_exon, T_exon, H_transcript, T_transcript, H_canonical, T_canonical)

exactReportedFusions <- dbFusionPairs %>% merge(allKnownFusions, by = c("H_gene", "T_gene")) %>% 
  select (sampleId, H_gene, T_gene, H_exon, T_exon, H_transcript, T_transcript, H_canonical, T_canonical) %>% mutate(Source = "known")

intergenicPromT <- dbFusionPairs %>% filter(H_gene != T_gene) %>% merge(knownPromiscuousT, by = "T_gene") %>% mutate(Source = "promiscuous_T")

intergenicPromH <- dbFusionPairs %>% filter(H_gene != T_gene) %>% merge(knownPromiscuousH, by = "H_gene") %>% mutate(Source = "promiscuous_H")

intragenicProm <- dbFusionPairs %>% filter(H_gene == T_gene) %>% merge(knownPromiscuousT, by = "T_gene") %>% filter(T_exon - H_exon > 1) %>% mutate(Source = "intragenic_T")

reportedFusions <- rbind(exactReportedFusions, intergenicPromT, intergenicPromH, intragenicProm) %>% mutate(Value = TRUE) %>% spread(Source, Value)

write.csv(reportedFusions, "~/data/r/reportedFusions.csv", row.names = FALSE)

topReportedFusions <- reportedFusions %>% group_by(H_gene, T_gene) %>% count() %>% arrange(-n)
write.csv(topReportedFusions, "~/data/r/topReportedFusions.csv", row.names = FALSE)

samplesMultiH <- reportedFusions %>% group_by(sampleId, H_gene) %>% summarise(countH = n()) %>% filter(countH > 1)
samplesMultiT <- reportedFusions %>% group_by(sampleId, T_gene) %>% summarise(countT = n()) %>% filter(countT > 1)
samplesMulti <- samplesMultiH %>% merge(samplesMultiT, by = "sampleId", all = TRUE)
