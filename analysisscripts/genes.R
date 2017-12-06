
####### Helpers #######
applySamples<-function(cohort, sampleFunction) {
  result = ddply(cohort, 1, sampleFunction, .drop = FALSE, .progress = "text", .parallel = FALSE)
  return (result)
}

leftJoin<-function(left, right) {
  temp = merge(x = left, y = right, by="row.names", all.x=TRUE)
  return (column_as_names(temp, "Row.names"))
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
query_gene_copynumber<-function(dbConnect, genes, sample) {
  geneString = paste("'", genes$gene, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT gene, minCopyNumber, maxCopyNumber ",
    "  FROM geneCopyNumber ",
    "WHERE sampleId = '", sample$sampleId, "'",
    "  AND gene in (", geneString, ")",
    sep = "")
 
  return (column_as_names(dbGetQuery(dbConnect, query))) 
}

attach_gene_copynumber<-function(dbConnect, genes, sample) {
  copynumber = query_gene_copynumber(pilotDB, genes, sample)
  return (leftJoin(genes, copynumber))
}

####### Gene Somatics #######
query_gene_somatics<-function(dbConnect, genes, sample) {
  geneString = paste("'", genes$gene, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT gene, ",
    "       SUM(1) as total,",
    "       MIN(ROUND(adjustedVAF * adjustedCopyNumber,2)) as minPloidy,",
    "       MAX(ROUND(adjustedVAF * adjustedCopyNumber, 2)) as maxPloidy,",
    "       SUM(IF(clonality = 'CLONAL', 1, 0)) AS clonal,",
    "       SUM(IF(clonality = 'SUBCLONAL', 1, 0)) AS subclonal,",
    "       SUM(IF(clonality = 'INCONSISTENT', 1, 0)) AS inconsistentClonality,",
    "       SUM(IF(LOH, 1, 0)) as LOH",
    "  FROM somaticVariant ",
    " WHERE sampleId = '", sample$sampleId, "'",
    "   AND filter = 'PASS'",
    "   AND effect not in ('non coding exon variant', 'synonymous variant', 'UTR variant', 'sequence feature')",
    "   AND gene in (", geneString, ")",
    " GROUP BY gene",
    sep = "")
  
  return (column_as_names(dbGetQuery(dbConnect, query))) 
}

attach_gene_somatics<-function(dbConnect, genes, sample) {
  somatics = query_gene_somatics(pilotDB, genes, sample)
  return (leftJoin(genes, somatics))
}

####### Sample Drivers #######
sample_drivers<-function(dbConnect, genes, sample) {
  genes$sampleId <- sample$sampleId
  genes = attach_gene_copynumber(pilotDB, genes, sample)
  genes = attach_gene_somatics(pilotDB, genes, sample)
  genes = genes[!is.na(genes$total) & (genes$maxPloidy + 0.5 > genes$minCopyNumber | genes$LOH > 0) , ]
  genes$sampleId <- sample$sampleId
  return (genes)
}

####### Cohort Drivers #######
cohort_drivers<-function(dbConnect, genes, cohort) {
  return (applySamples(cohort, function(x) {sample_drivers(pilotDB, genes, x)}))
}

############# Load Purple Data
load("~/hmf/purple.RData")
rm(allSamples)
rm(backupSamples)

pilotDB = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
genes = query_gene_panel(pilotDB)
cohortDrivers = cohort_drivers(pilotDB, genes, cohort[1:5,])
dbDisconnect(pilotDB)
