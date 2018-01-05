
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

neighbours<-function(topGene, aggregatedCopyNumbers) {
  index = match(topGene$gene, aggregatedCopyNumbers$gene)
  neighbours = aggregatedCopyNumbers[max(1,index-7):min(index+7, nrow(aggregatedCopyNumbers)),]
  return(neighbours)
}

copy_number_drivers<-function(allGenes, allGeneCopyNumbers, maxDriversPerChromosome = 20, chromosomes = c(1:22, "X")) {

  # Clean input
  allGenes = data.table(allGenes)
  allGeneCopyNumbers = data.table(allGeneCopyNumbers)
  allGeneCopyNumbers$cancerType = ifelse(is.na(allGeneCopyNumbers$cancerType), "NA", allGeneCopyNumbers$cancerType)

  # Genome level variables
  genomeDrivers = list()
  genomeAdjacent <- data.frame()
  genomeCandidates = data.frame()

  for (currentChromosome in chromosomes) {
    #currentChromosome = 1
    genomeDrivers[[currentChromosome]] = list()
    cat("Processing chromosome:", currentChromosome, "\n")

    # Chromosome level variables
    chromosomeGenes = allGenes[chromosome == currentChromosome][order(start)]
    chromosomeGeneCopyNumbers = allGeneCopyNumbers[chromosome == currentChromosome]
    chromosomeSummary = aggregate_gene_copy_numbers(chromosomeGeneCopyNumbers)

    driverGeneCopyNumbers = chromosomeGeneCopyNumbers
    for (i in 1:maxDriversPerChromosome) {
      driverSummary = aggregate_gene_copy_numbers(driverGeneCopyNumbers)
      if (nrow(driverSummary) == 0) {
        break
      }

      topGene = head(driverSummary[order(-score)], 1)
      cat("       Isolating gene:", topGene$gene, "\n")

      driverNeighbours = neighbours(topGene, driverSummary)
      chromosomeNeighbours = neighbours(topGene, chromosomeSummary)

      # Determine adjacent
      isAdjacent = adjacent_to_gene(topGene$gene,  driverGeneCopyNumbers, chromosomeGenes)
      adjacent = driverGeneCopyNumbers[c(isAdjacent)]
      adjacentSummary = aggregate_gene_copy_numbers_by_cancer_type(adjacent)
      geneCandidates = candidates(topGene, adjacentSummary)

      # Update genome variables
      genomeAdjacent = rbind(genomeAdjacent, adjacent[gene==topGene$gene])
      genomeCandidates = rbind(genomeCandidates, geneCandidates)
      genomeDrivers[[currentChromosome]][[i]] <- list(topGene = topGene, candidates = geneCandidates$candidates, driverNeighbours = driverNeighbours, chromosomeNeighbours = chromosomeNeighbours, adjacent = adjacentSummary)

      # Remove adjacent copy numbers ready for next loop...
      driverGeneCopyNumbers = driverGeneCopyNumbers[!c(isAdjacent)]
    }
  }

  genomeAdjacentSummary <- aggregate_gene_copy_numbers_by_cancer_type(genomeAdjacent)[order(chromosome, -N)]
  genomeAdjacentSummary$candidates <- genomeCandidates[match(genomeAdjacentSummary$gene, genomeCandidates$gene), c("candidates")]
  genomeDrivers[["summary"]] <- genomeAdjacentSummary

  return (genomeDrivers)
}

#maxDriversPerChromosome = 3
#chromosomes = c(1:2)

#load("~/hmf/geneCopyNumber.RData")
#allGenes = genes
#allGeneCopyNumbers = geneCopyNumberDeletes
#allGeneCopyNumbers = geneCopyNumberAmplifactions
#driverAmplifications = copy_number_drivers(allGenes, geneCopyNumberAmplifactions)
#driverDeletions = copy_number_drivers(allGenes, geneCopyNumberDeletes)

#copyNumberDeletions[[1]][[1]]

#summary = copyNumberDeletions$summary
