
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
  return (paste(candidates, collapse = ","))
}

neighbours<-function(driverGene, aggregatedCopyNumbers) {
  index = match(driverGene$gene, aggregatedCopyNumbers$gene)
  neighbours = aggregatedCopyNumbers[max(1,index-7):min(index+7, nrow(aggregatedCopyNumbers)),]
  return(neighbours)
}


chromosome_copy_number_drivers <- function(currentChromosome, allGenes, allGeneCopyNumbers,  maxDriversPerChromosome = 20) {
  result = list()

  chromosomeGenes = filter(allGenes, chromosome == currentChromosome) %>% arrange(start)
  chromosomeGeneCopyNumbers = filter(allGeneCopyNumbers, chromosome == currentChromosome)
  chromosomeSummary = aggregate_gene_copy_numbers(chromosomeGeneCopyNumbers)

  for (i in 1:maxDriversPerChromosome) {

    summary = aggregate_gene_copy_numbers(chromosomeGeneCopyNumbers)
    if (nrow(summary) == 0) {
      break;
    }

    driverGene = head(arrange(summary, -score), 1)
    if (driverGene$score == 0) {
      break;
    }

    cat("Isolating gene", paste(driverGene$chromosome,driverGene$gene,sep=":"),"\n")

    chromosomeNeighbours = neighbours(driverGene, chromosomeSummary)
    driverNeighbours = neighbours(driverGene, summary)
    driverGeneCopyNumbers = chromosomeGeneCopyNumbers[chromosomeGeneCopyNumbers$gene == driverGene$gene, ]

    # Determine adjacent
    isAdjacent = adjacent_to_gene(driverGene$gene,  chromosomeGeneCopyNumbers, chromosomeGenes)
    adjacent = chromosomeGeneCopyNumbers[isAdjacent, ]
    adjacentSummary = aggregate_gene_copy_numbers_by_cancer_type(adjacent)

    # Update genome variables
    result[[i]] <- list(
      chromsome = driverGene$chromosome,
      gene = driverGene$gene,
      localPeak =  summary[summary$gene == driverGene$gene, ]$score,
      globalPeak = chromosomeSummary[chromosomeSummary$gene == driverGene$gene, ]$score,
      candidates = candidates(driverGene, adjacentSummary),
      driverNeighbours = driverNeighbours,
      chromosomeNeighbours = chromosomeNeighbours,
      driverGeneCopyNumbers = driverGeneCopyNumbers,
      adjacent = adjacentSummary,
      neighbouringGenes = sum(adjacentSummary$score))

    # Remove adjacent copy numbers ready for next loop...
    chromosomeGeneCopyNumbers = chromosomeGeneCopyNumbers[!isAdjacent, ]
  }

  return (result)
}


copy_number_drivers<-function(allGenes, allGeneCopyNumbers, maxDriversPerChromosome = 20, chromosomes = c(1:22, "X")) {

  # Clean input
  allGenes = data.table(allGenes)
  allGeneCopyNumbers = data.table(allGeneCopyNumbers)
  allGeneCopyNumbers$cancerType = ifelse(is.na(allGeneCopyNumbers$cancerType), "NA", allGeneCopyNumbers$cancerType)

  chromosomeDrivers = sapply(chromosomes, function(x) {chromosome_copy_number_drivers(x, allGenes, geneCopyNumberDeletes, maxDriversPerChromosome)})

  driverGeneCopyNumbers = lapply(chromosomeDrivers, function(x) {x$driverGeneCopyNumbers})
  driverGeneCopyNumbersSummary = aggregate_gene_copy_numbers_by_cancer_type(do.call(rbind, driverGeneCopyNumbers))
  remainderSummary = data.frame(t(sapply(chromosomeDrivers, function(x) {c(gene=x$gene, globalPeak=x$globalPeak, localPeak = x$localPeak, candidates = x$candidates, neighbouringGenes = x$neighbouringGenes)})))
  summary = left_join(driverGeneCopyNumbersSummary, remainderSummary, by = "gene")
  chromosomeDrivers[["summary"]] <- summary

  return (chromosomeDrivers)
}


