#topGene = "RB1"
#chromosomeCopyNumbers = localDeletes
#chromosomeGenes[chromosomeGenes$gene == topGene]



adjacent_to_gene<-function(topGene, chromosomeCopyNumbers, chromosomeGenes) {

  adjacent_to_deleted<-function(previous_distance, copyNumbers) {
    previouslyAdjacentSamples = copyNumbers[geneDistance == previous_distance & adjacent, .N, by = .(sampleId)]
    result = !is.na(match(copyNumbers$sampleId, previouslyAdjacentSamples$sampleId))
    return (result)
  }

  chromosomeCopyNumbers$geneDistance <- ifelse(chromosomeCopyNumbers$gene == topGene, 0, NA)
  chromosomeCopyNumbers$adjacent <- ifelse(chromosomeCopyNumbers$gene == topGene, TRUE, NA)

  topGeneIndex = match(topGene, chromosomeGenes$gene)
  if (topGeneIndex != nrow(chromosomeGenes)) {
    for(neigbouringGeneIndex in (topGeneIndex+1):nrow(chromosomeGenes)) {

      #cat(neigbouringGeneIndex)

      neigbouringGene=chromosomeGenes[neigbouringGeneIndex, ]$gene
      neighbouringGeneDistance=gene_distance_index(topGeneIndex, neigbouringGeneIndex, chromosomeGenes)
      chromosomeCopyNumbers$geneDistance <- ifelse(chromosomeCopyNumbers$gene == neigbouringGene, neighbouringGeneDistance, chromosomeCopyNumbers$geneDistance)

      chromosomeCopyNumbers$adjacent <- ifelse(chromosomeCopyNumbers$geneDistance == neighbouringGeneDistance, adjacent_to_deleted(neighbouringGeneDistance-1, chromosomeCopyNumbers), chromosomeCopyNumbers$adjacent)
      if (nrow(chromosomeCopyNumbers[geneDistance == neighbouringGeneDistance & adjacent == TRUE]) == 0) {
        break
      }
    }
  }

  if (topGeneIndex != 1) {
    for(neigbouringGeneIndex in (topGeneIndex-1):1) {

      neigbouringGene=chromosomeGenes[neigbouringGeneIndex, ]$gene
      neighbouringGeneDistance=gene_distance_index(neigbouringGeneIndex, topGeneIndex, chromosomeGenes)
      chromosomeCopyNumbers$geneDistance <- ifelse(chromosomeCopyNumbers$gene == neigbouringGene, -neighbouringGeneDistance, chromosomeCopyNumbers$geneDistance)

      chromosomeCopyNumbers$adjacent <- ifelse(chromosomeCopyNumbers$geneDistance == -neighbouringGeneDistance, adjacent_to_deleted(-neighbouringGeneDistance+1, chromosomeCopyNumbers), chromosomeCopyNumbers$adjacent)
      if (nrow(chromosomeCopyNumbers[geneDistance == -neighbouringGeneDistance & adjacent == TRUE]) == 0) {
        break
      }
    }
  }

  return (chromosomeCopyNumbers$adjacent %in% TRUE)
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

copy_number_deletions<-function(allGenes, allGeneCopyNumbers) {
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


  for (currentChromosome in c(1:22, 'X')) {
    #currentChromosome = 1
    cat("Processing chromosome:", currentChromosome, "\n")

    chromosomeGenes = allGenes[chromosome == currentChromosome][order(start)]
    chromosomeGeneCopyNumbers = allGeneCopyNumbers[chromosome == currentChromosome]
    chromosomeSummary = aggregate_gene_copy_numbers(chromosomeGeneCopyNumbers)

    localDeletes = chromosomeGeneCopyNumbers

    chromosomeResult = list()
    for (i in 1:20) {
      peakSummary = aggregate_gene_copy_numbers(localDeletes)
      if (nrow(peakSummary) == 0) {
        break
      }

      topGene = head(peakSummary[order(-score)], 1)
      genesWithSimilarScore = peakSummary[score >= topGene$score -1]
      cat("    Isolating gene:", topGene$gene, "\n")

      topGeneIndexLocalPeak = match(topGene$gene, peakSummary$gene)
      localPeak = peakSummary[max(1,topGeneIndexLocalPeak-7):min(topGeneIndexLocalPeak+7, nrow(peakSummary)),]
      #localPeakLeft = peakSummary[max(1, topGeneIndexLocalPeak-1)]
      #localPeakRight = peakSummary[min(nrow(peakSummary), topGeneIndexLocalPeak+1)]

      topGeneIndexTotalPeak = match(topGene$gene, chromosomeSummary$gene)
      totalPeak = chromosomeSummary[max(1,topGeneIndexTotalPeak-7):min(topGeneIndexTotalPeak+7, nrow(chromosomeSummary)),]

      adjacent = adjacent_to_gene(topGene$gene,  localDeletes, chromosomeGenes)

      removed = localDeletes[c(adjacent)]
      topGeneRemovals = rbind(topGeneRemovals, removed[gene==topGene$gene])

      removedSummary = aggregate_gene_copy_numbers_by_cancer_type(removed)
      localDeletes = localDeletes[!c(adjacent)]

      peak = list(topGene = topGene, similarScoredGenes = genesWithSimilarScore, localPeak = localPeak, totalPeak = totalPeak, removed = removedSummary)
      chromosomeResult[[i]] <- peak
    }

    copyNumberDeletions[[currentChromosome]] <- chromosomeResult
  }
  copyNumberDeletions[["summary"]] <- removed_summary(topGeneRemovals)[order(chromosome, -N)]
  return (copyNumberDeletions)
}



#load("~/hmf/geneCopyNumber.RData")
#allGenes = genes
#allGeneCopyNumbers = geneCopyNumberDeletes
#allGeneCopyNumbers = geneCopyNumberAmplifactions
#copyNumberDeletions = copy_number_deletions(allGenes, allGeneCopyNumbers)
#jon = copyNumberDeletions$summary

#detach("package:purple", unload=TRUE)
