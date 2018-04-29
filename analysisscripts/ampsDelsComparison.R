####### COMBINE WITH TRAVIS AND PCAWG - DELETES
travisDels = read.csv(file = "~/Documents/TravisDeletes.csv", stringsAsFactors = F)
colnames(travisDels) <- c("travisPeak","range","chromosome","start","end")
travisDels$range = GRanges(travisDels$chromosome, IRanges(travisDels$start, travisDels$end))

load("~/hmf/RData/geneCopyNumberDeleteTargets.RData")
dels = geneCopyNumberDeleteTargets %>% select(gene_name = target, chromosome, start = superCandidatesStart, end = superCandidatesEnd, chromosomeBand, telomere, hmfCount = N)
dels$chromosomeBand <- paste0(dels$chromosome, dels$chromosomeBand)
dels$range = GRanges(dels$chromosome, IRanges(dels$start, dels$end))

# Travis
ol = as.matrix(findOverlaps(dels$range, travisDels$range, select="all"))
dels[ol[,1], "travisPeak"] <- travisDels[ol[,2], "travisPeak"]
missingTravisDels = travisDels[!travisDels$travisPeak %in% dels$travisPeak, ]
rm(ol)

# PCAWG
load(file = "~/hmf/RData/canonicalTranscripts.RData")
load("~/hmf/RData/pcawgDels.RData")
colnames(pcawgDels) <- c("gene","pcawgGeneCount")
dels = merge(dels, pcawgDels, by.x = "gene_name", by.y = "gene", all.x = T)

pcawgDels = pcawgDels[!pcawgDels$gene %in% pcawgDelsArm$gene, ]
pcawgDels = pcawgDels[!pcawgDels$gene %in% pcawgDelsTelomere$gene, ]
found = dels[dels$pcawgGeneCount > 0, "gene_name"]
pcawgDels = pcawgDels[!pcawgDels$gene %in% found, ]

pcawgDelTranscripts = canonicalTranscripts[canonicalTranscripts$gene %in% pcawgDels$gene, c("gene", "chromosome","geneStart","geneEnd")]
pcawgDels = left_join(pcawgDels, pcawgDelTranscripts, by = "gene")
pcawgDels$range = GRanges(pcawgDels$chromosome, IRanges(pcawgDels$geneStart, pcawgDels$geneEnd))

ol = as.matrix(findOverlaps(dels$range, pcawgDels$range, select="all"))
dels[ol[,1], "pcawgGeneRegion"] <- pcawgDels[ol[,2], "gene"]
dels[ol[,1], "pcawgGeneRegionCount"] <- pcawgDels[ol[,2], "pcawgGeneCount"]
dels$range <- NULL
missingPcawgDels = pcawgDels[!pcawgDels$gene %in% dels$pcawgGeneRegion, ]

colnames(pcawgDelsArm) <- c("gene","pcawgBandCount","chromosomeBand")
pcawgDelsArm = pcawgDelsArm %>% group_by(chromosomeBand) %>% summarise(pcawgBandCount=sum(pcawgBandCount), n = n()) #### NOTE THE GROUP BY AND SUM HERE!
dels = left_join(dels, pcawgDelsArm[, c("chromosomeBand", "pcawgBandCount")], by = "chromosomeBand")
found = dels %>% filter(pcawgBandCount > 0) %>% select(chromosomeBand)
missingPcawgDelsArm = pcawgDelsArm[!pcawgDelsArm$chromosomeBand %in% found$chromosomeBand, ]

colnames(pcawgDelsTelomere) <- c("telomere","pcawgTelomereCount")
dels = left_join(dels, pcawgDelsTelomere, by = "telomere")
found = dels %>% filter(pcawgTelomereCount > 0) %>% select(telomere)
missingPcawgDelsTelomere = pcawgDelsTelomere[!pcawgDelsTelomere$telomere %in% found$telomere, ]

rm(found)
rm(ol)




####### AMPS
travisAmps = read.csv(file = "~/Documents/TravisAmplifications.csv", stringsAsFactors = F)
colnames(travisAmps) <- c("travisPeak","range","chromosome","start","end")
travisAmps$range = GRanges(travisAmps$chromosome, IRanges(travisAmps$start, travisAmps$end))

load("~/hmf/RData/geneCopyNumberAmplificationTargets.RData")
amps = geneCopyNumberAmplificationTargets %>% select(gene_name = target, chromosome, start = superCandidatesStart, end = superCandidatesEnd, chromosomeBand, telomere, hmfCount = N)
amps$range = GRanges(amps$chromosome, IRanges(amps$start, amps$end))
amps$chromosomeBand <- paste0(amps$chromosome, amps$chromosomeBand)

# Travis
ol = as.matrix(findOverlaps(amps$range, travisAmps$range, select="all"))
amps[ol[,1], "travisPeak"] <- travisAmps[ol[,2], "travisPeak"]
missingTravisAmps = travisAmps[!travisAmps$travisPeak %in% amps$travisPeak, ]
rm(ol)

# PCAWG
load(file = "~/hmf/RData/canonicalTranscripts.RData")
load("~/hmf/RData/pcawgAmps.RData")
colnames(pcawgAmps) <- c("gene","pcawgGeneCount")
amps = merge(amps, pcawgAmps, by.x = "gene_name", by.y = "gene", all.x = T)

pcawgAmps = pcawgAmps[!pcawgAmps$gene %in% pcawgAmpsArm$gene, ]
pcawgAmps = pcawgAmps[!pcawgAmps$gene %in% pcawgAmpsTelomere$gene, ]
found = amps[amps$pcawgGeneCount > 0, "gene_name"]
pcawgAmps = pcawgAmps[!pcawgAmps$gene %in% found, ]

pcawgAmpTranscripts = canonicalTranscripts[canonicalTranscripts$gene %in% pcawgAmps$gene, c("gene", "chromosome","geneStart","geneEnd")]
pcawgAmps = left_join(pcawgAmps, pcawgAmpTranscripts, by = "gene")
pcawgAmps$range = GRanges(pcawgAmps$chromosome, IRanges(pcawgAmps$geneStart, pcawgAmps$geneEnd))

ol = as.matrix(findOverlaps(amps$range, pcawgAmps$range, select="all"))
amps[ol[,1], "pcawgGeneRegion"] <- pcawgAmps[ol[,2], "gene"]
amps[ol[,1], "pcawgGeneRegionCount"] <- pcawgAmps[ol[,2], "pcawgGeneCount"]
amps$range <- NULL
missingPcawgAmps = pcawgAmps[!pcawgAmps$gene %in% amps$pcawgGeneRegion, ]

colnames(pcawgAmpsArm) <- c("gene","pcawgBandCount","chromosomeBand")
pcawgAmpsArm = pcawgAmpsArm %>% group_by(chromosomeBand) %>% summarise(pcawgBandCount=sum(pcawgBandCount), n = n()) #### NOTE THE GROUP BY AND SUM HERE!
amps = left_join(amps, pcawgAmpsArm[, c("chromosomeBand", "pcawgBandCount")], by = "chromosomeBand")
found = amps %>% filter(pcawgBandCount > 0) %>% select(chromosomeBand)
missingPcawgAmpsArm = pcawgAmpsArm[!pcawgAmpsArm$chromosomeBand %in% found$chromosomeBand, ]

colnames(pcawgAmpsTelomere) <- c("telomere","pcawgTelomereCount")
amps = left_join(amps, pcawgAmpsTelomere, by = "telomere")
found = amps %>% filter(pcawgTelomereCount > 0) %>% select(telomere)
missingPcawgAmpsTelomere = pcawgAmpsTelomere[!pcawgAmpsTelomere$telomere %in% found$telomere, ]

rm(ol)

save(amps, file = "~/hmf/RData/amps.RData")
save(dels, file = "~/hmf/RData/dels.RData")
