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
 
  return (column_as_names(dbGetQuery(dbConnect, query), .drop = FALSE)) 
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
    "   AND effect not in ('non coding exon variant', 'synonymous variant', 'UTR variant', 'sequence feature')",
    "   AND gene in (", geneString, ")",
    sep = "")
  
  return (dbGetQuery(dbConnect, query)) 
}

sample_position_somatics<-function(dbConnect, genes, sample) {
  sampleGeneCopyNumbers = query_sample_copynumber(dbConnect, genes, sample)
  sampleSomatics = query_position_somatics(dbConnect, genes, sample)
  result = merge(sampleSomatics, sampleGeneCopyNumbers, by="gene",all.x=TRUE)
  result$biallelic <- ifelse(result$ploidy + 0.5 > result$minCopyNumber | result$loh, 1, 0)
  result$location <- paste(result$chromosome, result$position, sep = ":")
  return (result)
}

cohort_position_somatics<-function(dbConnect, genes, cohort) {
  return (applySamples(cohort, function(x) {sample_position_somatics(dbConnect, genes, x)}))
}

cohort_gene_somatics<-function(cohort_position_somatics) {
  DT = data.table(cohort_position_somatics)
  cohortGeneSomatics = DT[, .(somatics=length(loh), loh=sum(loh), 
                                 minSomaticCopyNumber=min(copyNumber), maxSomaticCopyNumber = max(copyNumber),
                                 minSomaticPloidy=min(ploidy), maxSomaticPloidy=max(ploidy), sumSomaticPloidy=sum(ploidy)), 
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

############# EXECUTION
load("~/hmf/purple.RData")
rm(allSamples)
rm(backupSamples)
cohort = cohort[cohort$sampleId %in% c('CPCT02180005T', 'CPCT02080055T', 'CPCT02180008T'), ]

pilotDB = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
genes = query_gene_panel(pilotDB)
dbDisconnect(pilotDB)
rm(pilotDB)

prodDB = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")

cohortPositionSomatics = cohort_position_somatics(prodDB, genes, cohort)
cohortGeneSomatics = cohort_gene_somatics(cohortPositionSomatics)
cohortGeneCopyNumbers = cohort_gene_copynumber(prodDB, genes, cohort)
cohortGeneStructuralVariants = cohort_gene_structual_variants(prodDB, genes, cohort)
cohortGeneComplete = merge(cohortGeneSomatics, cohortGeneStructuralVariants, by=c('gene', 'sampleId'), all=TRUE)
cohortGeneComplete = merge(cohortGeneComplete, cohortGeneCopyNumbers, by=c('gene', 'sampleId'), all.x=TRUE)

dbDisconnect(prodDB)
rm(prodDB)



