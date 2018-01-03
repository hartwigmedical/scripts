library(RMySQL)
library(data.table)

####### Queries #######
query_gene_copy_number_deletes<-function(dbConnect) {
  query = paste(
    "SELECT g.sampleId, g.chromosome, g.start, g.end, g.gene, g.minCopyNumber",
    "  FROM geneCopyNumber g, purity p",
    "WHERE g.sampleId = p.sampleId",
    "AND p.qcStatus = 'PASS'",
    "AND p.status != 'NO_TUMOR'",
    "AND g.germlineHetRegions = 0",
    "AND g.germlineHomRegions = 0",
    "AND g.minCopyNumber < 0.5",
    "AND g.chromosome <> 'Y'",
    sep = " ")
  return (dbGetQuery(dbConnect, query))
}

query_all_genes<-function(dbConnect, sample) {
  query = paste(
    "SELECT g.chromosome, g.start, g.end, g.gene",
    "  FROM geneCopyNumber g",
    " WHERE g.sampleId = '", sample, "'",
    "   AND g.chromosome <> 'Y'",
    sep = "")
  return (dbGetQuery(dbConnect, query))
}

query_clinical_data<-function(dbConnect) {
  query = paste(
    "SELECT c.sampleId, c.cancertype, c.birthYear, c.biopsyDate",
    " FROM clinical c",
    sep = " ")
  return ((dbGetQuery(dbConnect, query)))
}

####### Function #######
geneDistance<-function(topGeneName, geneCopyNumbers, chromosomeGenes) {
  topGeneIndex = match(topGeneName, chromosomeGenes$gene)
  distance = sapply(geneCopyNumbers$gene, function(x) {match(x, chromosomeGenes$gene) - topGeneIndex} )
  return (distance)
}


adjacent_to_gene<-function(topGene, chromosomeCopyNumbers, chromosomeGenes) {
  
  adjacent_to_deleted<-function(previous_distance, copyNumbers) {
    previouslyAdjacentSamples = copyNumbers[geneDistance == previous_distance & adjacent, .N, by = .(sampleId)]
    result = !is.na(match(copyNumbers$sampleId, previouslyAdjacentSamples$sampleId))
    return (result)
  }
  
  chromosomeCopyNumbers$geneDistance <- geneDistance(topGene, chromosomeCopyNumbers, chromosomeGenes)
  chromosomeCopyNumbers$adjacent <- ifelse(chromosomeCopyNumbers$geneDistance == 0, TRUE, FALSE)
  
  for(i in 1:30000) {
    chromosomeCopyNumbers$adjacent <- ifelse(chromosomeCopyNumbers$geneDistance == i, adjacent_to_deleted(i-1, chromosomeCopyNumbers), chromosomeCopyNumbers$adjacent)
    if (nrow(chromosomeCopyNumbers[geneDistance == i & adjacent == TRUE]) == 0) {
      break
    }
  }
  
  for(i in 1:30000) {
    chromosomeCopyNumbers$adjacent <- ifelse(chromosomeCopyNumbers$geneDistance == -i, adjacent_to_deleted(-i+1, chromosomeCopyNumbers), chromosomeCopyNumbers$adjacent)
    if (nrow(chromosomeCopyNumbers[geneDistance == -i & adjacent == TRUE]) == 0) {
      break
    }
  }
  
  return (chromosomeCopyNumbers$adjacent)
}

removed_summary<-function(removed) {
  cancerTypes = unique(removed$cancerType)
  removedByCancerType = dcast(removed, ... ~ cancerType, fun.aggregate = length, value.var = c("cancerType"))
  removedByCancerTypeSummary = removedByCancerType[, lapply(.SD, sum), by=list(gene, chromosome, start, end), .SDcols = cancerTypes]
  removedByCancerTypeSummary$total <- rowSums(removedByCancerTypeSummary[, cancerTypes, with=FALSE])
  return (removedByCancerTypeSummary)
}


copy_number_deletions<-function(allGenes, allDeletes, removeSample = FALSE) {
  allGenes = data.table(allGenes)
  DT = data.table(allDeletes[!is.na(allDeletes$cancerType), ])
  
  copyNumberDeletions = list()
  for (currentChromosome in c(1:22, 'X')) {
    cat("Processing chromosome:", currentChromosome, "\n")
    
    chromosomeGenes = allGenes[chromosome == currentChromosome][order(start)] 
    chromosomeCopyNumbers = DT[chromosome == currentChromosome]
    chromosomeSummary = chromosomeCopyNumbers[, .N, by=list(gene, chromosome, start, end)][order(start)]
    peakCopyNumber = chromosomeCopyNumbers
    
    chromosomeResult = list()
    for (i in 1:20) {
      peakSummary = peakCopyNumber[, .N, by=list(gene, chromosome, start, end)][order(start)]
      if (nrow(peakSummary) == 0) {
        break
      }

      topGene = head(peakSummary[order(-N)], 1)
      genesWithSameCount = peakSummary[N >= topGene$N -1]
      
      topGeneIndexLocalPeak = match(topGene$gene, peakSummary$gene)
      localPeak = peakSummary[max(1,topGeneIndexLocalPeak-7):min(topGeneIndexLocalPeak+7, nrow(peakSummary)),]
      
      topGeneIndexTotalPeak = match(topGene$gene, chromosomeSummary$gene)
      totalPeak = chromosomeSummary[max(1,topGeneIndexTotalPeak-7):min(topGeneIndexTotalPeak+7, nrow(chromosomeSummary)),]
      
      adjacent = adjacent_to_gene(topGene$gene,  peakCopyNumber, chromosomeGenes)
      removed = peakCopyNumber[c(adjacent)]
      removed[sampleId == 'CPCT02110002T']
      
      removedSummary = removed_summary(removed)
      peakCopyNumber = peakCopyNumber[!c(adjacent)]
      
      peak = list(topGene = topGene, similarSizedGenes = genesWithSameCount, localPeak = localPeak, totalPeak = totalPeak, removed = removedSummary)
      chromosomeResult[[i]] <- peak
    }
    
    copyNumberDeletions[[currentChromosome]] <- chromosomeResult
  }
  
  return (copyNumberDeletions)
}

####### LOAD DATA FROM DATABASE #######
#prodDB = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
#clinicalData = query_clinical_data(prodDB)
#allDeletes = query_gene_copy_number_deletes(prodDB)
#allDeletes$cancerType <- sapply(allDeletes$sampleId, function(x) {clinicalData[match(x, clinicalData$sampleId), c("cancerType")] })
#allGenes = query_all_genes(prodDB, allDeletes[1, c('sampleId')])
#save(allDeletes, allGenes, file = "~/hmf/copyNumberDeletions.RData")
#dbDisconnect(prodDB)
#rm(prodDB)

####### LOAD DATA FROM DATABASE #######
load("~/hmf/copyNumberDeletions.RData")
#allGenes = data.table(allGenes)
#DT = data.table(allDeletes[!is.na(allDeletes$cancerType), ])
#cancerTypes = unique(DT$cancerType)
#cancerTypes = cancerTypes[!is.na(cancerTypes)]

copyNumberDeletions = copy_number_deletions(allGenes, allDeletes)



chrom9 = copyNumberDeletions[["9"]]
chrom9[[7]]

sum(chrom9[[1]]$removed$total)
sum(chrom9[[2]]$removed$total)
sum(chrom9[[3]]$removed$total)
sum(chrom9[[4]]$removed$total)
sum(chrom9[[5]]$removed$total)
sum(removed$total)

chrom13 = copyNumberDeletions[[13]]
length(chrom13)
chrom13[[1]]
chrom13[[1]]$topGene
chrom13[[1]]$removed
