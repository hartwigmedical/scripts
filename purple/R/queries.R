query_gene_copy_number_deletes<-function(dbConnect) {
  query = paste(
    "SELECT g.sampleId, g.chromosome, g.start, g.end, g.gene, g.chromosomeBand, g.minCopyNumber, 1 as score",
    "  FROM geneCopyNumber g, purity p",
    " WHERE g.sampleId = p.sampleId",
    "   AND p.qcStatus = 'PASS'",
    "   AND p.status != 'NO_TUMOR'",
    "   AND p.purity > 0.15",
    "   AND g.germlineHetRegions = 0",
    "   AND g.germlineHomRegions = 0",
    "   AND g.minCopyNumber < 0.5",
    "   AND g.chromosome <> 'Y'",
    sep = " ")
  return (dbGetQuery(dbConnect, query))
}

query_gene_copy_number_amplifications<-function(dbConnect, cutoff = 3) {
  query = paste(
    "SELECT g.sampleId, g.chromosome, g.start, g.end, g.gene, g.chromosomeBand, g.minCopyNumber, log2(g.minCopyNumber / p.ploidy / ", cutoff, ") as score",
    "  FROM geneCopyNumber g, purity p",
    " WHERE g.sampleId = p.sampleId",
    "   AND p.qcStatus = 'PASS'",
    "   AND p.status != 'NO_TUMOR'",
    "   AND p.purity > 0.15",
    "   AND g.germlineHetRegions = 0",
    "   AND g.germlineHomRegions = 0",
    "   AND g.minCopyNumber / p.ploidy > ", cutoff,
    "   AND g.chromosome <> 'Y'",
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

query_purity<-function(dbConnect) {
  query = paste(
    "SELECT p.*",
    " FROM purity p",
    "WHERE qcStatus = 'PASS'",
    "  AND status <> 'NO_TUMOR'",
    sep = " ")
  return (dbGetQuery(dbConnect, query))
}

query_sample_data<-function(dbConnect) {
  query = paste(
    "SELECT sampleId, tumorPercentage AS pathologyPurity, arrivalDate, samplingDate",
    " FROM sample c",
    sep = " ")
  return (dbGetQuery(dbConnect, query))
}

query_somatic_variant_sample<-function(dbConnect, sampleId) {
  query = paste(
    "SELECT * ",
    " FROM somaticVariant ",
    "WHERE filter = 'PASS'",
    "   AND sampleId = '", sampleId, "'",
    sep = "")
  result = dbGetQuery(dbConnect, query)
  return (result)
}

query_msi_sample<-function(dbConnect, sampleId) {
  query = paste(
    "SELECT sampleId, count(*) as msiScore ",
    " FROM somaticVariant ",
    "WHERE filter = 'PASS'",
    "   AND type = 'INDEL'",
    "   AND repeatCount > 0 ",
    "   AND length(alt) <= 50",
    "   AND length(ref) <= 50",
    "   AND sampleId = '", sampleId, "'",
    " GROUP BY sampleId",
    sep = "")
  score = dbGetQuery(dbConnect, query)
  score$msiScore = score$msiScore/3095
  return (score)
}

