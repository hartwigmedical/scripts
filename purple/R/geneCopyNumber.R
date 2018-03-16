adjacent_to_gene<-function(topGene, geneCopyNumbers, genes) {

  geneCopyNumbers$geneDistance <- ifelse(geneCopyNumbers$gene == topGene, 0, NA)
  geneCopyNumbers$adjacent <- ifelse(geneCopyNumbers$gene == topGene, TRUE, NA)

  topGeneIndex = match(topGene, genes$gene)
  if (topGeneIndex != nrow(genes)) {
    for(neigbouringGeneIndex in (topGeneIndex+1):nrow(genes)) {

      neigbouringGene=genes[neigbouringGeneIndex, ]$gene
      neighbouringGeneDistance=gene_distance_index(topGeneIndex, neigbouringGeneIndex, genes)
      geneCopyNumbers$geneDistance <- ifelse(geneCopyNumbers$gene == neigbouringGene, neighbouringGeneDistance, geneCopyNumbers$geneDistance)

      previouslyAdjacent = geneCopyNumbers %>% filter(geneDistance == neighbouringGeneDistance-1 & adjacent) %>% group_by(sampleId) %>% count();
      if (nrow(previouslyAdjacent) == 0) {
        break
      }

      geneCopyNumbers$adjacent <- ifelse(geneCopyNumbers$geneDistance == neighbouringGeneDistance, geneCopyNumbers$sampleId %in% previouslyAdjacent$sampleId, geneCopyNumbers$adjacent)
    }
  }

  if (topGeneIndex != 1) {
    for(neigbouringGeneIndex in (topGeneIndex-1):1) {

      neigbouringGene=genes[neigbouringGeneIndex, ]$gene
      neighbouringGeneDistance=gene_distance_index(neigbouringGeneIndex, topGeneIndex, genes)
      geneCopyNumbers$geneDistance <- ifelse(geneCopyNumbers$gene == neigbouringGene, -neighbouringGeneDistance, geneCopyNumbers$geneDistance)

      previouslyAdjacent = geneCopyNumbers %>% filter(geneDistance == -neighbouringGeneDistance+1 & adjacent) %>% group_by(sampleId) %>% count();
      if (nrow(previouslyAdjacent) == 0) {
        break
      }
      geneCopyNumbers$adjacent <- ifelse(geneCopyNumbers$geneDistance == -neighbouringGeneDistance, geneCopyNumbers$sampleId %in% previouslyAdjacent$sampleId, geneCopyNumbers$adjacent)
    }
  }

  return (geneCopyNumbers$adjacent %in% TRUE)
}

aggregate_gene_copy_numbers<-function(geneCopyNumbers) {
  geneCopyNumbers %>%
    mutate(unsupported = minRegionStartSupport == "NONE" & minRegionEndSupport == "NONE") %>%
    mutate(telomereSupported = minRegionStartSupport == "TELOMERE" | minRegionEndSupport == "TELOMERE") %>%
    mutate(centromereSupported = minRegionStartSupport == "CENTROMERE" | minRegionEndSupport == "CENTROMERE") %>%
    group_by(gene, chromosome, start, end, chromosomeBand) %>%
    summarise(score=sum(score), sd = sd(minCopyNumber), N = n(), unsupported = sum(unsupported), telomereSupported = sum(telomereSupported), centromereSupported = sum(centromereSupported)) %>%
    arrange(start)
}

aggregate_gene_copy_numbers_by_cancer_type<-function(geneCopyNumbers) {
  aggregation = aggregate_gene_copy_numbers(geneCopyNumbers)
  aggregationByCancerType = geneCopyNumbers  %>%
    group_by(gene, cancerType) %>%
    summarise(count = n()) %>% spread(cancerType, count)

  return(left_join(aggregation, aggregationByCancerType, by = "gene"))
}

candidates<-function(driverGene, adjacentAggregate) {
  candidates = c(adjacentAggregate[adjacentAggregate$score >= driverGene$score -2 | adjacentAggregate$score >= 0.85 * driverGene$score, ]$gene)
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
      neighbouringGenes = sum(adjacentSummary$N))

    # Remove adjacent copy numbers ready for next loop...
    chromosomeGeneCopyNumbers = chromosomeGeneCopyNumbers[!isAdjacent, ]
  }

  return (result)
}


copy_number_drivers<-function(allGenes, allGeneCopyNumbers, maxDriversPerChromosome = 20, chromosomes = c(1:22, "X"), cl = NA) {

  # Clean input
  allGenes = data.table(allGenes)
  allGeneCopyNumbers = data.table(allGeneCopyNumbers)
  allGeneCopyNumbers$cancerType = ifelse(is.na(allGeneCopyNumbers$cancerType), "NA", allGeneCopyNumbers$cancerType)

  if (is.na(cl)) {
    chromosomeDrivers = sapply(chromosomes, function(x) {chromosome_copy_number_drivers(x, allGenes, allGeneCopyNumbers, maxDriversPerChromosome)})
  } else {
    cat("Going parallel")
    chromosomeDrivers = parSapply(cl, chromosomes, function(x) {chromosome_copy_number_drivers(x, allGenes, allGeneCopyNumbers, maxDriversPerChromosome)})
  }

  driverGeneCopyNumbers = lapply(chromosomeDrivers, function(x) {x$driverGeneCopyNumbers})
  driverGeneCopyNumbersSummary = aggregate_gene_copy_numbers_by_cancer_type(do.call(rbind, driverGeneCopyNumbers))
  remainderSummary = data.frame(t(sapply(chromosomeDrivers, function(x) {c(gene=x$gene, globalPeak=x$globalPeak, localPeak = x$localPeak, neighbouringGenes = x$neighbouringGenes, candidates = x$candidates)})), stringsAsFactors = F)
  remainderSummary[, 2:4] <- lapply(remainderSummary[, 2:4], as.numeric)

  summary = left_join(driverGeneCopyNumbersSummary, remainderSummary, by = "gene")
  summary$chromosome = factor(summary$chromosome, levels = c(1:22,"X","Y"))

  chromosomeDrivers[["summary"]] <- summary

  return (chromosomeDrivers)
}


#### RUBBISH

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

#currentChromosome = 9
#allGeneCopyNumbers = geneCopyNumberAmplifications

