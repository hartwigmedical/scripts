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

firstCandidate <- function(candidates) {
  candidateVector = unlist(strsplit(candidates, split = ","))
  return (candidateVector[1])
}

longestCandidate <- function(candidates, canonicalTranscripts) {
  candidateVector = data.frame(gene = unlist(strsplit(candidates, split = ",")), stringsAsFactors = F)
  candidateVector = left_join(candidateVector, canonicalTranscripts %>% select(gene, codingBases), by = "gene") %>%
    filter(codingBases == max(codingBases), codingBases > 0) 
  return (collapseCandidates(candidateVector$gene))
}

highestScoringCodingCandidate <- function(candidates, canonicalTranscripts) {
  candidateVector = data.frame(gene = unlist(strsplit(candidates, split = ",")), stringsAsFactors = F) %>% mutate(rank = row_number())
  candidateVector = left_join(candidateVector, canonicalTranscripts %>% select(gene, codingBases), by = "gene") %>% 
    mutate(isCoding = codingBases > 0) %>%
    arrange(-isCoding, rank)%>% dplyr::slice(1)
  return (collapseCandidates(candidateVector$gene))
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
load(file = "~/hmf/RData/reference/canonicalTranscripts.RData")
canonicalTranscripts$range = GRanges(canonicalTranscripts$chromosome, IRanges(canonicalTranscripts$geneStart, canonicalTranscripts$geneEnd))
canonicalTranscriptsOverlaps = data.frame(findOverlaps(canonicalTranscripts$range, canonicalTranscripts$range, type="within")) %>% filter(queryHits != subjectHits)
superGenes = data.frame(sub = canonicalTranscripts[canonicalTranscriptsOverlaps[, 1], c("gene")], super = canonicalTranscripts[canonicalTranscriptsOverlaps[, 2], c("gene")], stringsAsFactors = F) 

#### Load Amps and Dels Gene Panel
load("~/hmf/RData/processed/genePanel.RData")
genePanel[is.na(genePanel)] <- F
delGenePanel  = genePanel %>% select(-cosmicOncogene, -amplification) %>%
  filter(martincorena| hmf | cosmicCurated | cosmicTsg | deletion)
ampGenePanel  = genePanel %>% select(-cosmicTsg, -deletion) %>%
  filter(martincorena| hmf | cosmicCurated | cosmicOncogene | amplification)
rm(genePanel)

#### DELETIONS
load(file = "~/hmf/RData/processed/geneCopyNumberDeletionsSummary.RData")
dels = geneCopyNumberDeletionsSummary %>% 
  group_by(gene) %>%
  filter(score > 5, unsupported < N / 2) %>%
  mutate(candidatesCount = length(unlist(strsplit(candidates, split = ",")))) %>%
  mutate(superCandidates = superSizedCandidates(candidates, superGenes)) %>%
  mutate(superCandidatesCount = length(unlist(strsplit(superCandidates, split = ",")))) 

delCandidatesRange = apply(dels[, c("gene","superCandidates")], 1, function(x) {candidatesRange(x[1], x[2], canonicalTranscripts)})
dels = merge(dels, do.call(rbind, delCandidatesRange), by = "gene")

delCandidates = apply(dels[, c("gene","superCandidates")], 1, function(x) {categoriseCandidates(x[1], x[2], delGenePanel)})
dels = merge(dels, do.call(rbind, delCandidates), by = "gene")

dels =  dels %>% 
  group_by(gene) %>%
  mutate(highestScoring = highestScoringCodingCandidate(remainders, canonicalTranscripts), longest = longestCandidate(remainders, canonicalTranscripts )) 

dels$method <- "panel"
dels$target <- dels$hmf
dels$target <- ifelse(is.na(dels$target), dels$martincorena, dels$target)
dels$target <- ifelse(is.na(dels$target), dels$cosmicCurated, dels$target)
dels$method <- ifelse(is.na(dels$target), "known", dels$method)
dels$target <- ifelse(is.na(dels$target), dels$deletion, dels$target)
dels$method <- ifelse(is.na(dels$target), "cosmic", dels$method)
dels$target <- ifelse(is.na(dels$target), dels$cosmicTsg, dels$target)
dels$method <- ifelse(is.na(dels$target), "highest", dels$method)
dels$target <- ifelse(is.na(dels$target), dels$highestScoring, dels$target)
dels$telomere <- ifelse(dels$method == "highest" & dels$telomereSupported > dels$N / 2, paste0(dels$chromosome, substr(dels$chromosomeBand,1,1), "_telomere"), NA)
dels$centromere <- ifelse(dels$method == "highest" & dels$centromereSupported > dels$N / 2, paste0(dels$chromosome, substr(dels$chromosomeBand,1,1), "_centromere"), NA)

delsWithMultipleTargets = dels %>% group_by(gene) %>% mutate(n = length(unlist(strsplit(target, split = ",")))) %>% filter(n > 1)
dels = dels %>% group_by(gene) %>% mutate(target = firstCandidate(target))

geneCopyNumberDeleteTargets = dels
save(geneCopyNumberDeleteTargets, file = "~/hmf/RData/processed/geneCopyNumberDeleteTargets.RData")

rm(dels)
rm(delCandidates)
rm(canonicalTranscriptsOverlaps)
rm(delCandidatesRange)
rm(geneCopyNumberDeletionsSummary)

#### AMPLIFICAIONS
load(file = "~/hmf/RData/processed/geneCopyNumberAmplificationSummary.RData")
amps = geneCopyNumberAmplificationSummary %>% 
  group_by(gene) %>%
  filter(score > 20) %>%
  mutate(candidatesCount = length(unlist(strsplit(candidates, split = ",")))) %>%
  mutate(superCandidates = superSizedCandidates(candidates, superGenes)) %>%
  mutate(superCandidatesCount = length(unlist(strsplit(superCandidates, split = ",")))) 

ampCandidatesRange = apply(amps[, c("gene","superCandidates")], 1, function(x) {candidatesRange(x[1], x[2], canonicalTranscripts)})
amps = merge(amps, do.call(rbind, ampCandidatesRange), by = "gene")

ampCandidates = apply(amps[, c("gene","superCandidates")], 1, function(x) {categoriseCandidates(x[1], x[2], ampGenePanel)})
amps = merge(amps, do.call(rbind, ampCandidates), by = "gene")

amps = amps %>% 
  group_by(gene) %>%
  mutate(highestScoring = highestScoringCodingCandidate(remainders, canonicalTranscripts), longest = longestCandidate(remainders, canonicalTranscripts ))

amps$method <- "panel"
amps$target <- amps$hmf
amps$target <- ifelse(is.na(amps$target), amps$martincorena, amps$target)
amps$target <- ifelse(is.na(amps$target), amps$cosmicCurated, amps$target)
amps$method <- ifelse(is.na(amps$target), "known", amps$method)
amps$target <- ifelse(is.na(amps$target), amps$amplification, amps$target)
amps$method <- ifelse(is.na(amps$target), "cosmic", amps$method)
amps$target <- ifelse(is.na(amps$target), amps$cosmicOncogene, amps$target)
amps$method <- ifelse(is.na(amps$target), "longest", amps$method)
amps$target <- ifelse(is.na(amps$target), amps$longest, amps$target)
amps$method <- ifelse(is.na(amps$target), "highest", amps$method)
amps$target <- ifelse(is.na(amps$target), amps$gene, amps$target)
amps$telomere <- ifelse(amps$method == "highest" & amps$telomereSupported > amps$N / 2, paste0(amps$chromosome, substr(amps$chromosomeBand,1,1), "_telomere"), NA)

ampsWithMultipleTargets = amps %>% group_by(gene) %>% mutate(n = length(unlist(strsplit(target, split = ",")))) %>% filter(n > 1)
amps = amps %>% group_by(gene) %>% mutate(target = firstCandidate(target))

geneCopyNumberAmplificationTargets = amps
save(geneCopyNumberAmplificationTargets, file = "~/hmf/RData/processed/geneCopyNumberAmplificationTargets.RData")

rm(amps)
rm(ampCandidates)
rm(ampCandidatesRange)
rm(geneCopyNumberAmplificationSummary)

