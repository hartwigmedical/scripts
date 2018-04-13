adjacent_to_arm<-function(driverGene, geneCopyNumbers, genes) {
  driverArm = substr(driverGene$chromosomeBand, 1, 1)
  driverChromosome = driverGene$chromosome
  samplesWithDriverGeneCNV = geneCopyNumbers %>% filter(gene == driverGene$gene) %>% select(sampleId)
  result = geneCopyNumbers %>% mutate(
    arm = substr(geneCopyNumbers$chromosomeBand, 1, 1),
    adjacent = (chromosome == driverChromosome & arm == driverArm & sampleId %in% samplesWithDriverGeneCNV$sampleId)) %>% select(adjacent)
  return (result$adjacent)
}

adjacent_to_gene<-function(driverGene, geneCopyNumbers, genes) {
  topGene = driverGene$gene

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


chromosome_copy_number_drivers <- function(currentChromosome, allGenes, allGeneCopyNumbers,  maxDriversPerChromosome = 20, adjacent) {
  result = list()

  chromosomeGenes = filter(allGenes, chromosome == currentChromosome) %>% arrange(start)
  chromosomeGeneCopyNumbers = filter(allGeneCopyNumbers, chromosome == currentChromosome)
  chromosomeSummary = aggregate_gene_copy_numbers(chromosomeGeneCopyNumbers)

  if (maxDriversPerChromosome == 0) {
    maxDriversPerChromosome = 30000
  }

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
    if (adjacent == "gene") {
      isAdjacent = adjacent_to_gene(driverGene,  chromosomeGeneCopyNumbers, chromosomeGenes)
    } else {
      isAdjacent = adjacent_to_arm(driverGene,  chromosomeGeneCopyNumbers, chromosomeGenes)
    }

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


copy_number_drivers<-function(allGenes, allGeneCopyNumbers, maxDriversPerChromosome = 0, chromosomes = c(1:22, "X"), cl = NA, adjacent = "gene") {

  if (!adjacent %in% c("gene","arm")) {
    stop("Argument adjacent must be one of: 'gene','arm'")
  }

  # Clean input
  allGenes = data.table(allGenes)
  allGeneCopyNumbers = data.table(allGeneCopyNumbers)
  allGeneCopyNumbers$cancerType = ifelse(is.na(allGeneCopyNumbers$cancerType), "NA", allGeneCopyNumbers$cancerType)

  if (is.na(cl)) {
    chromosomeDrivers = sapply(chromosomes, function(x) {chromosome_copy_number_drivers(x, allGenes, allGeneCopyNumbers, maxDriversPerChromosome, adjacent)})
  } else {
    cat("Going parallel")
    chromosomeDrivers = parSapply(cl, chromosomes, function(x) {chromosome_copy_number_drivers(x, allGenes, allGeneCopyNumbers, maxDriversPerChromosome, adjacent)})
  }

  allChromosomeDrivers = unlist(chromosomeDrivers, recursive = F)
  driverGeneCopyNumbers = lapply(allChromosomeDrivers, function(x) {x$driverGeneCopyNumbers})
  driverGeneCopyNumbersSummary = aggregate_gene_copy_numbers_by_cancer_type(do.call(rbind, driverGeneCopyNumbers))
  remainderSummary = data.frame(t(sapply(allChromosomeDrivers, function(x) {c(gene=x$gene, globalPeak=x$globalPeak, localPeak = x$localPeak, neighbouringGenes = x$neighbouringGenes, candidates = x$candidates)})), stringsAsFactors = F)
  remainderSummary[, 2:4] <- lapply(remainderSummary[, 2:4], as.numeric)

  summary = left_join(driverGeneCopyNumbersSummary, remainderSummary, by = "gene")
  summary$chromosome = factor(summary$chromosome, levels = c(1:22,"X","Y"))

  chromosomeDrivers[["summary"]] <- summary

  return (chromosomeDrivers)
}


#### WORKING
#library(dplyr)
#library(tidyr)
#load("~/hmf/RData/geneCopyNumberAmplificationsData.RData")
#chromosomes = c(1:22, "X")
#chromosomes = c(8)
#allGeneCopyNumbers = geneCopyNumberAmplifications
#maxDriversPerChromosome = 5
#cl = NA
#jon = copy_number_drivers(allGenes, allGeneCopyNumbers, maxDriversPerChromosome = 0, chromosomes = c(22))
#currentChromosome = 8

