library(RMySQL)

query_gene_copy_number_deletes<-function(dbConnect) {
  query = paste(
    "SELECT g.sampleId, g.chromosome, g.start, g.end, g.gene, g.chromosomeBand, g.minCopyNumber, 1 as score",
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

query_gene_copy_number_amplifications<-function(dbConnect) {
  query = paste(
    "SELECT g.sampleId, g.chromosome, g.start, g.end, g.gene, g.chromosomeBand, g.minCopyNumber, log2(g.minCopyNumber / p.ploidy / 2) as score",
    "  FROM geneCopyNumber g, purity p",
    "WHERE g.sampleId = p.sampleId",
    "AND p.qcStatus = 'PASS'",
    "AND p.status != 'NO_TUMOR'",
    "AND g.germlineHetRegions = 0",
    "AND g.germlineHomRegions = 0",
    "AND g.minCopyNumber / p.ploidy > 2",
    "AND g.chromosome <> 'Y'",
    sep = " ")
  return (dbGetQuery(dbConnect, query))
}

query_gene_copy_number<-function(dbConnect, sampleId) {
  query = paste(
    "SELECT g.*",
    "  FROM geneCopyNumber g",
    " WHERE g.sampleId = '", sampleId, "'",
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
