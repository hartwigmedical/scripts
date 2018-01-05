
adjacent_to_gene<-function(topGene, geneCopyNumbers, genes) {

  adjacent_to_deleted<-function(previous_distance, copyNumbers) {
    previouslyAdjacentSamples = copyNumbers[geneDistance == previous_distance & adjacent, .N, by = .(sampleId)]
    result = !is.na(match(copyNumbers$sampleId, previouslyAdjacentSamples$sampleId))
    return (result)
  }

  geneCopyNumbers$geneDistance <- ifelse(geneCopyNumbers$gene == topGene, 0, NA)
  geneCopyNumbers$adjacent <- ifelse(geneCopyNumbers$gene == topGene, TRUE, NA)

  topGeneIndex = match(topGene, genes$gene)
  if (topGeneIndex != nrow(genes)) {
    for(neigbouringGeneIndex in (topGeneIndex+1):nrow(genes)) {

      neigbouringGene=genes[neigbouringGeneIndex, ]$gene
      neighbouringGeneDistance=gene_distance_index(topGeneIndex, neigbouringGeneIndex, genes)
      geneCopyNumbers$geneDistance <- ifelse(geneCopyNumbers$gene == neigbouringGene, neighbouringGeneDistance, geneCopyNumbers$geneDistance)

      geneCopyNumbers$adjacent <- ifelse(geneCopyNumbers$geneDistance == neighbouringGeneDistance, adjacent_to_deleted(neighbouringGeneDistance-1, geneCopyNumbers), geneCopyNumbers$adjacent)
      if (nrow(geneCopyNumbers[geneDistance == neighbouringGeneDistance & adjacent == TRUE]) == 0) {
        break
      }
    }
  }

  if (topGeneIndex != 1) {
    for(neigbouringGeneIndex in (topGeneIndex-1):1) {

      neigbouringGene=genes[neigbouringGeneIndex, ]$gene
      neighbouringGeneDistance=gene_distance_index(neigbouringGeneIndex, topGeneIndex, genes)
      geneCopyNumbers$geneDistance <- ifelse(geneCopyNumbers$gene == neigbouringGene, -neighbouringGeneDistance, geneCopyNumbers$geneDistance)

      geneCopyNumbers$adjacent <- ifelse(geneCopyNumbers$geneDistance == -neighbouringGeneDistance, adjacent_to_deleted(-neighbouringGeneDistance+1, geneCopyNumbers), geneCopyNumbers$adjacent)
      if (nrow(geneCopyNumbers[geneDistance == -neighbouringGeneDistance & adjacent == TRUE]) == 0) {
        break
      }
    }
  }

  return (geneCopyNumbers$adjacent %in% TRUE)
}


aggregate_gene_copy_numbers_by_cancer_type<-function(geneCopyNumbers) {
  cancerTypes = unique(geneCopyNumbers$cancerType)
  geneCopyNumbersByCancerType = dcast(geneCopyNumbers, ... ~ cancerType, fun.aggregate = length, value.var = c("cancerType"))
  geneCopyNumbersByCancerTypeAggregation = geneCopyNumbersByCancerType[, c(lapply(.SD, sum), list(sd=sd(minCopyNumber), score=sum(score), N=length(score)))  , by=list(gene, chromosome, start, end, chromosomeBand), .SDcols = as.character(cancerTypes)]
  return (geneCopyNumbersByCancerTypeAggregation[order(-score)])
}

aggregate_gene_copy_numbers<-function(geneCopyNumbers) {
  return (geneCopyNumbers[, .(score=sum(score), sd = sd(minCopyNumber)), by=list(gene, chromosome, start, end, chromosomeBand)][order(start)])
}

candidates<-function(topGene, adjacentAggregate) {
  candidates = c(adjacentAggregate[score >= topGene$score -2 | score >= 0.95*topGene$score]$gene)
  candidatesDF = data.frame(gene=topGene$gene, candidates=I(list(candidates)), stringsAsFactors = FALSE)
  return (candidatesDF)
}

copy_number_drivers<-function(allGenes, allGeneCopyNumbers, maxDriversPerChromosome = 20, chromosomes = c(1:22, "X")) {
  allGenes = data.table(allGenes)
  allGeneCopyNumbers = data.table(allGeneCopyNumbers)
  allGeneCopyNumbers$cancerType = ifelse(is.na(allGeneCopyNumbers$cancerType), "NA", allGeneCopyNumbers$cancerType)

  copyNumberDeletions = list()
  topGeneRemovals <- data.frame(sampleId=character(),
                   chromosome=character(),
                   start=integer(),
                   end=integer(),
                   gene=character(),
                   minCopyNumber=double(),
                   score=double(),
                   chromosomeBand=character(),
                   cancerType=character())

  allCanditates = data.frame(gene=character(), candidates=list(character()))

  for (currentChromosome in chromosomes) {
    #currentChromosome = 1
    cat("Processing chromosome:", currentChromosome, "\n")

    chromosomeGenes = allGenes[chromosome == currentChromosome][order(start)]
    chromosomeGeneCopyNumbers = allGeneCopyNumbers[chromosome == currentChromosome]
    chromosomeSummary = aggregate_gene_copy_numbers(chromosomeGeneCopyNumbers)

    localDeletes = chromosomeGeneCopyNumbers


    chromosomeResult = list()
    for (i in 1:maxDriversPerChromosome) {
      peakSummary = aggregate_gene_copy_numbers(localDeletes)
      if (nrow(peakSummary) == 0) {
        break
      }

      topGene = head(peakSummary[order(-score)], 1)
      genesWithSimilarScore = peakSummary[score >= topGene$score -1]
      cat("       Isolating gene:", topGene$gene, "\n")

      topGeneIndexLocalPeak = match(topGene$gene, peakSummary$gene)
      localPeak = peakSummary[max(1,topGeneIndexLocalPeak-7):min(topGeneIndexLocalPeak+7, nrow(peakSummary)),]
      #localPeakLeft = peakSummary[max(1, topGeneIndexLocalPeak-1)]
      #localPeakRight = peakSummary[min(nrow(peakSummary), topGeneIndexLocalPeak+1)]

      topGeneIndexTotalPeak = match(topGene$gene, chromosomeSummary$gene)
      totalPeak = chromosomeSummary[max(1,topGeneIndexTotalPeak-7):min(topGeneIndexTotalPeak+7, nrow(chromosomeSummary)),]

      isAdjacent = adjacent_to_gene(topGene$gene,  localDeletes, chromosomeGenes)

      adjacent = localDeletes[c(isAdjacent)]
      adjacentSummary = aggregate_gene_copy_numbers_by_cancer_type(adjacent)

      topGeneRemovals = rbind(topGeneRemovals, adjacent[gene==topGene$gene])


      geneCandidates = candidates(topGene, adjacentSummary)
      allCanditates = rbind(allCanditates, geneCandidates)

      localDeletes = localDeletes[!c(isAdjacent)]

      peak = list(topGene = topGene, candidates = geneCandidates$candidates, localPeak = localPeak, totalPeak = totalPeak, adjacent = adjacentSummary)
      chromosomeResult[[i]] <- peak
    }

    copyNumberDeletions[[currentChromosome]] <- chromosomeResult
  }

  summary <- aggregate_gene_copy_numbers_by_cancer_type(topGeneRemovals)[order(chromosome, -N)]
  summary$candidates <- allCanditates[match(summary$gene, allCanditates$gene), c("candidates")]
  copyNumberDeletions[["summary"]] <- summary
  return (copyNumberDeletions)
}


#load("~/hmf/geneCopyNumber.RData")
#allGenes = genes
#allGeneCopyNumbers = geneCopyNumberDeletes
#allGeneCopyNumbers = geneCopyNumberAmplifactions
#copyNumberDeletions = copy_number_drivers(allGenes, allGeneCopyNumbers, maxDriversPerChromosome = 3, chromosomes = c(1:2))
#summary = copyNumberDeletions$summary
