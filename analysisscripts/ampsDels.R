library(dplyr)
library(tidyr)
library(GenomicRanges)

collapseCandidates <- function(x) {
  x = x[!is.na(x)]
  if (length(x) == 0) {
    return (NA)
  }
  return (paste(x, collapse = ","))
}


categoriseCandidates <-function(gene, candidates, genePanel) {
 
  candidateVector = unlist(strsplit(candidates, split = ","))
  candidatePanel = genePanel[genePanel$gene_name %in% candidateVector, ]
  
  result = candidatePanel %>% mutate_all(funs(ifelse(. != F, gene_name, NA))) %>% select(-gene_name) %>% summarise_all(funs(collapseCandidates))
  
  result$remainders <- collapseCandidates(candidateVector[!candidateVector %in% genePanel$gene_name])
  result$gene = gene
  
  return(result)
}

superSizedCandidates <-function(candidates, superGenes) {

  candidateVector = unlist(strsplit(candidates, split = ","))
  supers = superGenes %>% filter(sub %in% candidateVector) %>% select(super)
  union = sort(unique(c(supers$super, candidateVector)))
  
  return(collapseCandidates(union))
}

candidatesRange <- function(gene, candidates, canonicalTranscripts) {
  candidateVector = unlist(strsplit(candidates, split = ","))
  candidateTranscripts = canonicalTranscripts[canonicalTranscripts$gene %in% candidateVector, ]
  return (data.frame(gene, superCandidatesStart = as.numeric(min(candidateTranscripts$geneStart)), superCandidatesEnd = max(candidateTranscripts$geneEnd)))
}

#### SUPER GENES
load(file = "~/hmf/RData/canonicalTranscripts.RData")
canonicalTranscripts$range = GRanges(canonicalTranscripts$chromosome, IRanges(canonicalTranscripts$geneStart, canonicalTranscripts$geneEnd))
canonicalTranscriptsOverlaps = data.frame(findOverlaps(canonicalTranscripts$range, canonicalTranscripts$range, type="within")) %>% filter(queryHits != subjectHits)
superGenes = data.frame(sub = canonicalTranscripts[canonicalTranscriptsOverlaps[, 1], c("gene")], super = canonicalTranscripts[canonicalTranscriptsOverlaps[, 2], c("gene")], stringsAsFactors = F) 

#### ADD FRAGILE SITES TO GENE PANEL
load("~/hmf/RData/fragileGenes.RData")
load("~/hmf/RData/genePanel.RData")
genePanel = merge(genePanel, fragileGenes, by = "gene_name", all=T)
rm(fragileGenes)

#### DELETIONS
load(file = "~/hmf/RData/geneCopyNumberDeletesDriverSummary.RData")
dels = data.frame(geneCopyNumberDeletesDriverSummary)
dels = dels[dels$unsupported < dels$N / 2, ]
#dels = dels[dels$N > 3, ]

dels$candidatesCount = sapply(dels[, c("candidates")], function(x) {length(unlist(strsplit(x, split = ",")))})
dels$superCandidates = sapply(dels[, c("candidates")], function(x) {superSizedCandidates(x, superGenes)})
dels$superCandidatesCount = sapply(dels[, c("superCandidates")], function(x) {length(unlist(strsplit(x, split = ",")))})

delCandidatesRange = apply(dels[, c("gene","superCandidates")], 1, function(x) {candidatesRange(x[1], x[2], canonicalTranscripts)})
dels = merge(dels, do.call(rbind, delCandidatesRange), by = "gene")

delCandidates = apply(dels[, c("gene","superCandidates")], 1, function(x) {categoriseCandidates(x[1], x[2], genePanel)})
dels = merge(dels, do.call(rbind, delCandidates), by = "gene")

dels$method <- "panel"
dels$target <- dels$hmf
dels$target <- ifelse(is.na(dels$target), dels$martincorena, dels$target)
dels$target <- ifelse(is.na(dels$target), dels$cosmicCurated, dels$target)
dels$method <- ifelse(is.na(dels$target), "cosmic", dels$method)
dels$target <- ifelse(is.na(dels$target), dels$cosmicTsg, dels$target)
#dels$method <- ifelse(is.na(dels$target), "telomere", dels$method)
#dels$target <- ifelse(is.na(dels$target) & dels$telomereSupported > dels$N / 2, paste0(dels$chromosome, substr(dels$chromosomeBand,1,1), "_telomere"), dels$target)
dels$method <- ifelse(is.na(dels$target), "highest", dels$method)
dels$target <- ifelse(is.na(dels$target), dels$gene, dels$target)
dels$telomere <- ifelse(dels$method == "highest" & dels$telomereSupported > dels$N / 2, paste0(dels$chromosome, substr(dels$chromosomeBand,1,1), "_telomere"), NA)

geneCopyNumberDeleteTargets = dels
save(geneCopyNumberDeleteTargets, file = "~/hmf/RData/geneCopyNumberDeleteTargets.RData")

rm(dels)
rm(delCandidates)
rm(geneCopyNumberDeletesDriverSummary)
rm(canonicalTranscriptsOverlaps)
rm(delCandidatesRange)

#### AMPLIFICAIONS
load(file = "~/hmf/RData/geneCopyNumberAmplificationSummary.RData")
amps = data.frame(geneCopyNumberAmplificationSummary)
amps$candidatesCount = sapply(amps[, c("candidates")], function(x) {length(unlist(strsplit(x, split = ",")))})
amps$superCandidates = sapply(amps[, c("candidates")], function(x) {superSizedCandidates(x, superGenes)})
amps$superCandidatesCount = sapply(amps[, c("superCandidates")], function(x) {length(unlist(strsplit(x, split = ",")))})

ampCandidatesRange = apply(amps[, c("gene","superCandidates")], 1, function(x) {candidatesRange(x[1], x[2], canonicalTranscripts)})
amps = merge(amps, do.call(rbind, ampCandidatesRange), by = "gene")

ampCandidates = apply(amps[, c("gene","superCandidates")], 1, function(x) {categoriseCandidates(x[1], x[2], genePanel)})
ampCandidates = do.call(rbind, ampCandidates)
amps = merge(amps, ampCandidates, by = "gene")

amps$method <- "panel"
amps$target <- amps$hmf
amps$target <- ifelse(is.na(amps$target), amps$martincorena, amps$target)
amps$target <- ifelse(is.na(amps$target), amps$cosmicCurated, amps$target)
amps$method <- ifelse(is.na(amps$target), "cosmic", amps$method)
amps$target <- ifelse(is.na(amps$target), amps$cosmicOncogene, amps$target)
#amps$method <- ifelse(is.na(amps$target), "telomere", amps$method)
#amps$target <- ifelse(is.na(amps$target) & amps$telomereSupported > amps$N / 2, paste0(amps$chromosome, substr(amps$chromosomeBand,1,1), "_telomere"), amps$target)
amps$method <- ifelse(is.na(amps$target), "highest", amps$method)
amps$target <- ifelse(is.na(amps$target), amps$gene, amps$target)
amps$telomere <- ifelse(amps$method == "highest" & amps$telomereSupported > amps$N / 2, paste0(amps$chromosome, substr(amps$chromosomeBand,1,1), "_telomere"), NA)

geneCopyNumberAmplificationTargets = amps
save(geneCopyNumberAmplificationTargets, file = "~/hmf/RData/geneCopyNumberAmplificationTargets.RData")

rm(amps)
rm(ampCandidates)
rm(ampCandidatesRange)
rm(geneCopyNumberAmplificationSummary)


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
