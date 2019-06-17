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

longestCandidate <- function(candidates, canonicalTranscripts) {
  candidateVector = data.frame(gene = unlist(strsplit(candidates, split = ",")), stringsAsFactors = F)
  candidateVector = left_join(candidateVector, canonicalTranscripts %>% select(gene, codingBases), by = "gene") %>%
    mutate(isCoding = codingBases > 0, rank = row_number()) %>%
    arrange(!isCoding, -codingBases, rank) %>% 
    dplyr::slice(1)
  return (collapseCandidates(candidateVector$gene))
}

#longestCandidate("MIR5708,TPD52,FABP5,ZBTB10", canonicalTranscripts)
#longestCandidate("WASH7P,MIR5708", canonicalTranscripts)

highestScoringCodingCandidate <- function(candidates, canonicalTranscripts) {
  candidateVector = data.frame(gene = unlist(strsplit(candidates, split = ",")), stringsAsFactors = F) %>% mutate(rank = row_number())
  candidateVector = left_join(candidateVector, canonicalTranscripts %>% select(gene, codingBases), by = "gene") %>% 
    mutate(isCoding = codingBases > 0) %>%
    arrange(!isCoding, rank) %>% dplyr::slice(1)
  return (collapseCandidates(candidateVector$gene))
}

categoriseCandidates <-function(gene, candidates, genePanel) {
 
  candidateVector = unlist(strsplit(candidates, split = ","))
  candidatePanel = genePanel[genePanel$gene %in% candidateVector, ]
  
  result = candidatePanel %>% mutate_all(funs(ifelse(. != F, gene, NA))) %>% select(-gene) %>% summarise_all(funs(collapseCandidates))
  
  result$remainders <- collapseCandidates(candidateVector[!candidateVector %in% genePanel$gene])
  result$gene = gene
  
  return(result)
}

superSizedCandidates <-function(candidates, superGenes) {
  candidateVector = data.frame(gene = unlist(strsplit(candidates, split = ",")), stringsAsFactors = F) %>% 
    mutate(rank = row_number())
  candidateSupers = inner_join(candidateVector, superGenes, by = c("gene" = "sub")) %>% select(gene = super, rank)
  result = bind_rows(candidateVector, candidateSupers) %>% arrange(rank)
  return(collapseCandidates(unique(result$gene)))
}

candidatesRange <- function(gene, candidates, canonicalTranscripts) {
  candidateVector = unlist(strsplit(candidates, split = ","))
  candidateTranscripts = canonicalTranscripts[canonicalTranscripts$gene %in% candidateVector, ]
  return (data.frame(gene, superCandidatesStart = as.numeric(min(candidateTranscripts$geneStart)), superCandidatesEnd = max(candidateTranscripts$geneEnd)))
}


singleTarget <- function(targets, superCandidates) {
  targetVector = unlist(strsplit(targets, split = ","))
  superCandidateDF = data.frame(gene = unlist(strsplit(superCandidates, split = ",")), stringsAsFactors = F)  %>% mutate(rank = row_number()) %>%
    filter(gene %in% targetVector) %>% arrange(rank)
  return (superCandidateDF[1, "gene"])
}

#copyNumberSummary= dels
#cosmic = "cosmicTsg"; known = "knownDeletion";  actionable = "actionableDeletion";
  
attach_target <- function(copyNumberSummary, cosmic, known, actionable) {
  copyNumberSummary$method <- "dnds"
  copyNumberSummary$target <- copyNumberSummary$hmfDnds
  copyNumberSummary$target <- ifelse(is.na(copyNumberSummary$target), copyNumberSummary$martincorenaDnds, copyNumberSummary$target)
  copyNumberSummary$target <- ifelse(is.na(copyNumberSummary$target), copyNumberSummary$cosmicCurated, copyNumberSummary$target)
  copyNumberSummary$method <- ifelse(is.na(copyNumberSummary$target), "actionable", copyNumberSummary$method)
  copyNumberSummary$target <- ifelse(is.na(copyNumberSummary$target), copyNumberSummary[[actionable]], copyNumberSummary$target)
  copyNumberSummary$method <- ifelse(is.na(copyNumberSummary$target), "known", copyNumberSummary$method)
  copyNumberSummary$target <- ifelse(is.na(copyNumberSummary$target), copyNumberSummary[[known]], copyNumberSummary$target)
  copyNumberSummary$method <- ifelse(is.na(copyNumberSummary$target), "cosmic", copyNumberSummary$method)
  copyNumberSummary$target <- ifelse(is.na(copyNumberSummary$target), copyNumberSummary[[cosmic]], copyNumberSummary$target)
  copyNumberSummary$method <- ifelse(is.na(copyNumberSummary$target), "drup", copyNumberSummary$method)
  copyNumberSummary$target <- ifelse(is.na(copyNumberSummary$target), copyNumberSummary$actionableDrup, copyNumberSummary$target)
  copyNumberSummary$method <- ifelse(is.na(copyNumberSummary$target), "longest", copyNumberSummary$method)
  copyNumberSummary$target <- ifelse(is.na(copyNumberSummary$target), copyNumberSummary$longest, copyNumberSummary$target)
  copyNumberSummary$telomere <- ifelse(copyNumberSummary$method == "longest" & copyNumberSummary$telomereSupported > copyNumberSummary$N / 2, paste0(copyNumberSummary$chromosome, substr(copyNumberSummary$chromosomeBand,1,1), "_telomere"), NA)
  copyNumberSummary$centromere <- ifelse(copyNumberSummary$method == "longest" & copyNumberSummary$centromereSupported > copyNumberSummary$N / 2, paste0(copyNumberSummary$chromosome, substr(copyNumberSummary$chromosomeBand,1,1), "_centromere"), NA)  
  
  copyNumberSummary = copyNumberSummary %>% group_by(gene) %>% mutate(target = singleTarget(target, superCandidates))
  return(copyNumberSummary)
}

#### SUPER GENES
load(file = "~/hmf/analysis/cohort/reference/canonicalTranscripts.RData")
canonicalTranscripts$range = GRanges(canonicalTranscripts$chromosome, IRanges(canonicalTranscripts$geneStart, canonicalTranscripts$geneEnd))
canonicalTranscriptsOverlaps = data.frame(findOverlaps(canonicalTranscripts$range, canonicalTranscripts$range, type="within")) %>% filter(queryHits != subjectHits)
superGenes = data.frame(sub = canonicalTranscripts[canonicalTranscriptsOverlaps[, 1], c("gene")], super = canonicalTranscripts[canonicalTranscriptsOverlaps[, 2], c("gene")], stringsAsFactors = F) 

#### Load Amps and Dels Gene Panel
load("~/hmf/analysis/cohort/processed/genePanel.RData")
genePanel[is.na(genePanel)] <- F
delGenePanel  = genePanel %>% select(gene, martincorenaDnds, hmfDnds, cosmicCurated, cosmicTsg, knownDeletion, actionableDeletion, actionableDrup) %>%
  filter(martincorenaDnds| hmfDnds | cosmicCurated | cosmicTsg | knownDeletion | actionableDeletion | actionableDrup)
ampGenePanel  = genePanel %>% select(gene, martincorenaDnds, hmfDnds, cosmicCurated, cosmicOncogene, knownAmplification, actionableAmplification, actionableDrup) %>%
  filter(martincorenaDnds| hmfDnds | cosmicCurated | cosmicOncogene | knownAmplification | actionableAmplification | actionableDrup)
rm(genePanel)

#### DELETIONS
load(file = "~/hmf/analysis/cohort/processed/geneCopyNumberDeletionsSummary.RData")
dels = geneCopyNumberDeletionsSummary %>% 
  group_by(gene) %>%
  filter(score > 7, unsupported < N / 2) %>%
  mutate(candidatesCount = length(unlist(strsplit(candidates, split = ",")))) %>%
  mutate(superCandidates = superSizedCandidates(candidates, superGenes)) %>%
  mutate(superCandidatesCount = length(unlist(strsplit(superCandidates, split = ",")))) 

delCandidatesRange = apply(dels[, c("gene","superCandidates")], 1, function(x) {candidatesRange(x[1], x[2], canonicalTranscripts)})
dels = merge(dels, do.call(rbind, delCandidatesRange), by = "gene")

delCandidates = apply(dels[, c("gene","superCandidates")], 1, function(x) {categoriseCandidates(x[1], x[2], delGenePanel)})
dels = merge(dels, do.call(rbind, delCandidates), by = "gene")

dels =  dels %>% 
  group_by(gene) %>%
  mutate(longest = longestCandidate(remainders, canonicalTranscripts)) 

dels = attach_target(dels, "cosmicTsg", "knownDeletion", "actionableDeletion")

geneCopyNumberDeleteTargets = dels
save(geneCopyNumberDeleteTargets, file = "~/hmf/analysis/cohort/processed/geneCopyNumberDeleteTargets.RData")

rm(dels)
rm(delCandidates)
rm(canonicalTranscriptsOverlaps)
rm(delCandidatesRange)
rm(geneCopyNumberDeletionsSummary)

#### AMPLIFICAIONS
load(file = "~/hmf/analysis/cohort/processed/geneCopyNumberAmplificationSummary.RData")
amps = geneCopyNumberAmplificationSummary %>% 
  group_by(gene) %>%
  filter(score > 35) %>%
  mutate(candidatesCount = length(unlist(strsplit(candidates, split = ",")))) %>%
  mutate(superCandidates = superSizedCandidates(candidates, superGenes)) %>%
  mutate(superCandidatesCount = length(unlist(strsplit(superCandidates, split = ",")))) 

ampCandidatesRange = apply(amps[, c("gene","superCandidates")], 1, function(x) {candidatesRange(x[1], x[2], canonicalTranscripts)})
amps = merge(amps, do.call(rbind, ampCandidatesRange), by = "gene")

ampCandidates = apply(amps[, c("gene","superCandidates")], 1, function(x) {categoriseCandidates(x[1], x[2], ampGenePanel)})
amps = merge(amps, do.call(rbind, ampCandidates), by = "gene")

amps = amps %>% 
  group_by(gene) %>%
  mutate(longest = longestCandidate(remainders, canonicalTranscripts))

amps = attach_target(amps, "cosmicOncogene", "knownAmplification", "actionableAmplification")

geneCopyNumberAmplificationTargets = amps
save(geneCopyNumberAmplificationTargets, file = "~/hmf/analysis/cohort/processed/geneCopyNumberAmplificationTargets.RData")
