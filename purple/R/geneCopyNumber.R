
adjacent_to_gene<-function(topGene, geneCopyNumbers, genes) {

  adjacent_to_deleted<-function(previous_distance, copyNumbers) {
    copyNumbers = data.table(copyNumbers)
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
      if (nrow(geneCopyNumbers[geneCopyNumbers$geneDistance == neighbouringGeneDistance & geneCopyNumbers$adjacent == TRUE, ]) == 0) {
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
      if (nrow(geneCopyNumbers[geneCopyNumbers$geneDistance == -neighbouringGeneDistance & geneCopyNumbers$adjacent == TRUE, ]) == 0) {
        break
      }
    }
  }

  return (geneCopyNumbers$adjacent %in% TRUE)
}

aggregate_gene_copy_numbers_by_cancer_type<-function(geneCopyNumbers) {
  aggregation = data.table(aggregate_gene_copy_numbers(geneCopyNumbers))[order(gene)]
  aggregationsByCancerType = data.table(dcast(geneCopyNumbers,  gene ~ cancerType, fun.aggregate = length, value.var = c("cancerType")))[order(gene)]

  return(cbind(aggregation, aggregationsByCancerType[, -c("gene")])[order(-score)])
}


aggregate_gene_copy_numbers<-function(geneCopyNumbers) {
  group_by(geneCopyNumbers, gene, chromosome, start, end, chromosomeBand) %>%
    summarise(score=sum(score), sd = sd(minCopyNumber), N = n()) %>%
    arrange(start)
}

candidates<-function(driverGene, adjacentAggregate) {
  candidates = c(adjacentAggregate[score >= driverGene$score -2 | score >= 0.95*driverGene$score]$gene)
  candidatesDF = data.frame(gene=driverGene$gene, candidates=I(list(candidates)), stringsAsFactors = FALSE)
  return (candidatesDF)
}

neighbours<-function(driverGene, aggregatedCopyNumbers) {
  index = match(driverGene$gene, aggregatedCopyNumbers$gene)
  neighbours = aggregatedCopyNumbers[max(1,index-7):min(index+7, nrow(aggregatedCopyNumbers)),]
  return(neighbours)
}

copy_number_drivers<-function(allGenes, allGeneCopyNumbers, maxDriversPerChromosome = 20, chromosomes = c(1:22, "X")) {

  # Clean input
  allGenes = data.table(allGenes)
  allGeneCopyNumbers = data.table(allGeneCopyNumbers)
  allGeneCopyNumbers$cancerType = ifelse(is.na(allGeneCopyNumbers$cancerType), "NA", allGeneCopyNumbers$cancerType)

  # Genome level variables
  allDrivers = list()
  allAdjacent <- data.frame()
  allGeneAttributes = data.frame()

  for (currentChromosome in chromosomes) {
    #currentChromosome = 17
    allDrivers[[currentChromosome]] = list()
    cat("Processing chromosome:", currentChromosome, "\n")

    # Chromosome level variables
    chromosomeGenes = filter(allGenes, chromosome == currentChromosome) %>% arrange(start)
    chromosomeGeneCopyNumbers = filter(allGeneCopyNumbers, chromosome == currentChromosome)
    chromosomeSummary = aggregate_gene_copy_numbers(chromosomeGeneCopyNumbers)

    driverGeneCopyNumbers = chromosomeGeneCopyNumbers
    for (i in 1:maxDriversPerChromosome) {

      driverSummary = aggregate_gene_copy_numbers(driverGeneCopyNumbers)
      if (nrow(driverSummary) == 0) {
        break
      }

      driverGene = head(arrange(driverSummary, -score), 1)
      if (driverGene$score == 0) {
        break
      }

      cat("       Isolating gene:", driverGene$gene,"\n")

      driverNeighbours = neighbours(driverGene, driverSummary)
      chromosomeNeighbours = neighbours(driverGene, chromosomeSummary)

      # Determine adjacent
      isAdjacent = adjacent_to_gene(driverGene$gene,  driverGeneCopyNumbers, chromosomeGenes)
      adjacent = driverGeneCopyNumbers[c(isAdjacent), ]
      adjacentSummary = aggregate_gene_copy_numbers_by_cancer_type(adjacent)

      # Gene Attributes
      geneAttributes = candidates(driverGene, adjacentSummary)
      geneAttributes$globalPeak <- chromosomeSummary[chromosomeSummary$gene == driverGene$gene, ]$score
      geneAttributes$localPeak <- driverSummary[driverSummary$gene == driverGene$gene, ]$score
      geneAttributes$neighbouringGenes <- sum(adjacentSummary$score)

      # Update genome variables
      allAdjacent = rbind(allAdjacent, adjacent[adjacent$gene==driverGene$gene, ])
      allGeneAttributes = rbind(allGeneAttributes, geneAttributes)
      allDrivers[[currentChromosome]][[i]] <- list(driverGene = driverGene, candidates = geneAttributes$candidates, driverNeighbours = driverNeighbours, chromosomeNeighbours = chromosomeNeighbours, adjacent = adjacentSummary)

      # Remove adjacent copy numbers ready for next loop...
      driverGeneCopyNumbers = driverGeneCopyNumbers[!isAdjacent, ]
    }
  }

  summary = aggregate_gene_copy_numbers_by_cancer_type(allAdjacent)[order(chromosome, -N)]
  summary$candidates <- allGeneAttributes[match(summary$gene, allGeneAttributes$gene), c("candidates")]
  summary$neighbouringGenes <- allGeneAttributes[match(summary$gene, allGeneAttributes$gene), c("neighbouringGenes")]
  summary$localPeak <- allGeneAttributes[match(summary$gene, allGeneAttributes$gene), c("localPeak")]
  summary$globalPeak <- allGeneAttributes[match(summary$gene, allGeneAttributes$gene), c("globalPeak")]

  allDrivers[["summary"]] <- summary

  return (allDrivers)
}

#maxDriversPerChromosome = 3
#chromosomes = c(1:2)

#load("~/hmf/geneCopyNumber.RData")
#allGenes = genes
#allGeneCopyNumbers = geneCopyNumberDeletes
#driverAmplifications = copy_number_drivers(genes, geneCopyNumberAmplifactions)
#driverDeletions = copy_number_drivers(genes, geneCopyNumberDeletes)
#driverDeletions = copy_number_drivers(allGenes, geneCopyNumberDeletes, maxDriversPerChromosome = 5, chromosomes = c(9:10))
#copyNumberDeletions = driverDeletions$summary
#copyNumberAmplifications = driverAmplifications$summary

#save(copyNumberDeletions, copyNumberAmplifications, file = "~/hmf/copyNumberSummaries.RData")
#save(driverAmplifications, driverDeletions, file = "~/hmf/copyNumberRaw.RData")
