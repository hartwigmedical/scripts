
adjacent_to_gene_knockout<-function(topGene, geneCopyNumbers, genes) {

  adjacent_to_deleted<-function(targetGene, previousDistance, targetCopyNumbers) {
    result <- group_by(targetCopyNumbers, sampleId) %>%
      filter((geneDistance == previousDistance & adjacent) | gene == targetGene)  %>%
      summarise(priorDeleted = sum(ifelse(geneDistance == previousDistance, deleted, 0)) > 0,
                priorLOH = sum(ifelse(geneDistance == previousDistance, loh, 0)) > 0,
                currentDeleted = sum(ifelse(gene == targetGene, deleted, 0)) > 0,
                currentLOH = sum(ifelse(gene == targetGene, loh, 0)) > 0) %>%
      ungroup() %>%
      filter(priorDeleted | (priorLOH & !currentDeleted))
    return(result)
  }

  geneCopyNumbers$geneDistance <- ifelse(geneCopyNumbers$gene == topGene, 0, NA)
  geneCopyNumbers$adjacent <- ifelse(geneCopyNumbers$gene == topGene, geneCopyNumbers$deleted > 0 | geneCopyNumbers$biallelicVariant > 0, NA)

  topGeneIndex = match(topGene, genes$gene)
  if (topGeneIndex != nrow(genes)) {
    for(neigbouringGeneIndex in (topGeneIndex+1):nrow(genes)) {
      neigbouringGene=genes[neigbouringGeneIndex, ]$gene
      neighbouringGeneDistance=gene_distance_index(topGeneIndex, neigbouringGeneIndex, genes)
      geneCopyNumbers$geneDistance <- ifelse(geneCopyNumbers$gene == neigbouringGene, neighbouringGeneDistance, geneCopyNumbers$geneDistance)

      adjacent = adjacent_to_deleted(neigbouringGene, neighbouringGeneDistance - 1, geneCopyNumbers)
      geneCopyNumbers$adjacent <- ifelse(geneCopyNumbers$gene == neigbouringGene, geneCopyNumbers$sampleId %in% adjacent$sampleId, geneCopyNumbers$adjacent)
      if (nrow(adjacent) == 0) {
        break
      }
    }
  }

  if (topGeneIndex != 1) {
    for(neigbouringGeneIndex in (topGeneIndex-1):1) {

      neigbouringGene=genes[neigbouringGeneIndex, ]$gene
      neighbouringGeneDistance=gene_distance_index(neigbouringGeneIndex, topGeneIndex, genes)
      geneCopyNumbers$geneDistance <- ifelse(geneCopyNumbers$gene == neigbouringGene, -neighbouringGeneDistance, geneCopyNumbers$geneDistance)

      adjacent = adjacent_to_deleted(neigbouringGene, -neighbouringGeneDistance + 1, geneCopyNumbers)
      geneCopyNumbers$adjacent <- ifelse(geneCopyNumbers$gene == neigbouringGene, geneCopyNumbers$sampleId %in% adjacent$sampleId, geneCopyNumbers$adjacent)
      if (nrow(adjacent) == 0) {
        break
      }
    }
  }

  return (geneCopyNumbers$adjacent %in% TRUE)
}

aggregate_gene_copy_number_knockouts<-function(geneCopyNumbers) {
  group_by(geneCopyNumbers, gene, chromosome, start, end, chromosomeBand) %>%
    summarise(score=sum(score), sd = sd(minCopyNumber), N = n(), deleted = sum(deleted), loh = sum(loh), biallelicVariant = sum(biallelicVariant)) %>%
    arrange(start)
}

copy_number_knockouts<-function(allGenes, allGeneCopyNumbers, maxDriversPerChromosome = 20, chromosomes = c(1:22, "X")) {

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
    chromosomeSummary = aggregate_gene_copy_number_knockouts(chromosomeGeneCopyNumbers)

    driverGeneCopyNumbers = chromosomeGeneCopyNumbers
    for (i in 1:maxDriversPerChromosome) {

      driverSummary = aggregate_gene_copy_number_knockouts(driverGeneCopyNumbers)
      if (nrow(driverSummary) == 0) {
        break
      }

      driverGene = head(arrange(driverSummary, -score, -deleted), 1)
      if (driverGene$score == 0) {
        break
      }

      cat("       Isolating gene:", driverGene$gene,"\n")

      driverNeighbours = neighbours(driverGene, driverSummary)
      chromosomeNeighbours = neighbours(driverGene, chromosomeSummary)

      # Determine adjacent
      isAdjacent = adjacent_to_gene_knockout(driverGene$gene,  driverGeneCopyNumbers, chromosomeGenes)
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
