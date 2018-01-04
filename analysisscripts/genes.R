library(RMySQL)
library(data.table)
library(plyr)

####### Helpers #######
applySamples<-function(cohort, sampleFunction) {
  result = ddply(cohort, 1, sampleFunction, .drop = FALSE, .progress = "text", .parallel = FALSE)
  return (result)
}

column_as_names<-function(df, .column = "gene", .drop = TRUE) {
  row.names(df) <- df[, c(.column)]
  if (.drop) {
    df[, c(.column)] <- NULL
  }
  return (df)
}

####### Gene Panel #######
query_gene_panel<-function(dbConnect, panel = "HMF Paper") {
  query = paste(
    "SELECT gene ",
    "  FROM genePanel ",
    " WHERE panel='", panel, "'",
    sep = "")
  
  return (column_as_names(dbGetQuery(dbConnect, query), .drop = FALSE))
}

####### Gene Copy Number #######
query_sample_copynumber<-function(dbConnect, genes, sample) {
  geneString = paste("'", genes$gene, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT gene, minCopyNumber, maxCopyNumber ",
    "  FROM geneCopyNumber ",
    "WHERE sampleId = '", sample$sampleId, "'",
    "  AND gene in (", geneString, ")",
    sep = "")
 
  return((column_as_names(dbGetQuery(dbConnect, query), .drop = FALSE))) 
}

cohort_gene_copynumber<-function(dbConnect, genes, cohort) {
  return (applySamples(cohort, function(x) {query_sample_copynumber(dbConnect, genes, x)}))
}

####### Somatics #######
query_position_somatics<-function(dbConnect, genes, sample) {
  geneString = paste("'", genes$gene, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT sampleId, gene, chromosome, position, adjustedCopyNumber as copyNumber, ",
    "       ROUND(adjustedVAF * adjustedCopyNumber,2) as ploidy, ",
    "       clonality, loh", 
    "  FROM somaticVariant s ",
    " WHERE sampleId = '", sample$sampleId, "'",
    "   AND filter = 'PASS'",
    "   AND effect not in ('non coding exon variant', 'synonymous variant', 'UTR variant', 'sequence feature','intron variant')",
    "   AND gene in (", geneString, ")",
    sep = "")
  return (dbGetQuery(dbConnect, query)) 
}

sample_position_somatics<-function(dbConnect, sampleGeneCopyNumbers, sample) {
  sampleSomatics = query_position_somatics(dbConnect, genes, sample)
  result = merge(sampleSomatics, sampleGeneCopyNumbers, by="gene",all.x=TRUE)
  result$biallelic <- ifelse(result$ploidy + 0.5 > result$minCopyNumber | result$loh, TRUE, FALSE)
  result$location <- paste(result$chromosome, result$position, sep = ":")
  return (result)
}

cohort_position_somatics<-function(dbConnect, cohortGeneCopyNumbers, cohort) {
  return (applySamples(cohort, function(x) {sample_position_somatics(dbConnect, cohortGeneCopyNumbers[cohortGeneCopyNumbers$sampleId == x$sampleId, ], x)}))
}

cohort_gene_somatics<-function(cohortPositionSomatics) {
  DT = data.table(cohortPositionSomatics)
  cohortGeneSomatics = DT[, .(somatics=length(loh),
                                 minSomaticCopyNumber=min(copyNumber), maxSomaticCopyNumber = max(copyNumber),
                                 minSomaticPloidy=min(ploidy), maxSomaticPloidy=max(ploidy), sumSomaticPloidy=sum(ploidy),                                  
                                 somaticSingleHitBiallelic=any(biallelic)), 
                             by=list(gene, sampleId)]
  return (cohortGeneSomatics)
}


####### Structural Variants #######
query_sample_structural_variants<-function(dbConnect, genes, sample) {
  geneString = paste("'", genes$gene, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT sampleId, gene, min(copyNumber) as minCopyNumber, max(copyNumber) as maxCopyNumber, min(ploidy) as minPloidy, max(ploidy) as maxPloidy FROM (",
    "  SELECT s.id, sampleId, startChromosome as chromosome, startPosition as position, startOrientation as orientation, ploidy, gene, adjustedStartCopyNumber as copyNumber ",
    "    FROM structuralVariant s, structuralVariantBreakend b",
    "   WHERE b.structuralVariantId = s.id",
    "     AND isStartEnd",
    "     AND isCanonicalTranscript",
    "     AND sampleId = '", sample$sampleId, "'",
    "   UNION",
    "  SELECT s.id, sampleId, endChromosome as chromosome, endPosition as position, endOrientation as orientation, ploidy, gene, adjustedEndCopyNumber as copyNumber  ",
    "    FROM structuralVariant s, structuralVariantBreakend b ",
    "   WHERE b.structuralVariantId = s.id",
    "     AND !isStartEnd",
    "     AND isCanonicalTranscript",
    "     AND sampleId = '", sample$sampleId, "'",
    " ) legs ",
    " WHERE gene in (", geneString, ")",
    "GROUP BY id, gene, sampleId", sep = "")
  return (dbGetQuery(dbConnect, query)) 
}

sample_gene_structual_variants<-function(dbConnect, genes, sample) {
  queryResult = query_sample_structural_variants(dbConnect, genes, sample)
  DT = data.table(queryResult)
  sampleGene = DT[, .(structuralVariants=length(minCopyNumber), 
                      minSVCopyNumber=min(minCopyNumber), maxSVCopyNumber=max(maxCopyNumber), 
                      minSVPloidy=min(minPloidy), maxSVPloidy=max(maxPloidy), sumSVPloidy=sum(maxPloidy))
                  , by=list(gene, sampleId)]
  return (sampleGene)
}

cohort_gene_structual_variants<-function(dbConnect, genes, cohort) {
  return (applySamples(cohort, function(x) {sample_gene_structual_variants(dbConnect, genes, x)}))
}

############# COMBINE AND ADD FEATURES
gene_summary<-function(cohort, cohortGeneCopyNumbers, cohortGeneSomatics, cohortGeneStructuralVariants) {
  
  # Merge
  result = merge(cohortGeneSomatics, cohortGeneStructuralVariants, by=c('gene', 'sampleId'), all=TRUE)
  result = merge(result, cohortGeneCopyNumbers, by=c('gene', 'sampleId'), all=TRUE)
  
  #Add cohort data
  result$ploidy <- sapply(result$sampleId, function(x) {cohort[match(x, cohort$sampleId), c("ploidy")] })
  result$cancerType <- sapply(result$sampleId, function(x) {cohort[match(x, cohort$sampleId), c("cancerType")] })
  
  # Copy number features 
  result$copyNumberAmplification <- ifelse(result$minCopyNumber > 8 | result$minCopyNumber > 2.2 * result$ploidy, TRUE, FALSE)
  result$copyNumberDeletion <- ifelse(result$minCopyNumber < 0.5, TRUE, FALSE)
  result$ploidy <- NULL
  
  # Filter out uninteresting rows
  rowsToKeep = !is.na(result$somatics) | !is.na(result$structuralVariants) | result$copyNumberDeletion | result$copyNumberAmplification
  result <- result[rowsToKeep]
  
  # Structural Variant Features
  result$svBiallelic <- ifelse(result$sumSVPloidy >  result$minCopyNumber, TRUE, FALSE)
  
  # Somatic Features
  result$somaticMultiHitBiallelic<-ifelse(result$somatics > 1 & result$sumSomaticPloidy + 0.5 > result$minCopyNumber, TRUE, FALSE)
  result$somaticSingleHitNoneBiallelic<-ifelse(result$somatics == 1 & !result$somaticSingleHitBiallelic, TRUE, FALSE)
  result$somaticMultiHitNonBiallelic<-ifelse(result$somatics > 1 & !result$somaticMultiHitBiallelic, TRUE, FALSE)
  
  return (result)
}



############# EXECUTION
load("~/hmf/pilot.RData")
rm(allSamples)
rm(backupSamples)
#cohort = cohort[1:5, ]

pilotDB = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
genes = query_gene_panel(pilotDB)
cohortGeneStructuralVariants = cohort_gene_structual_variants(pilotDB, genes, cohort)
cohortGeneCopyNumbers = cohort_gene_copynumber(pilotDB, genes, cohort)
cohortPositionSomatics = cohort_position_somatics(pilotDB, cohortGeneCopyNumbers, cohort)
cohortGeneSomatics = cohort_gene_somatics(cohortPositionSomatics)
dbDisconnect(pilotDB)
rm(pilotDB)

cohortGeneComplete = gene_summary(cohort, cohortGeneCopyNumbers, cohortGeneSomatics, cohortGeneStructuralVariants)


# save(cohort, cohortPositionSomatics, cohortGeneSomatics, cohortGeneCopyNumbers, cohortGeneStructuralVariants, cohortGeneComplete, file = "~/hmf/pilotGenes.RData")

save(cohort, cohortGeneCopyNumbers, cohortGeneStructuralVariants, file = "~/hmf/prodGenes.RData")

