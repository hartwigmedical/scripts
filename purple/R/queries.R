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
    "SELECT *",
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
    "   AND length(alt) <= 50",
    "   AND length(ref) <= 50",
    "   AND repeatCount >= 4",
    "   AND sampleId = '", sampleId, "'",
    "   AND ((length(repeatSequence) BETWEEN 2 AND 4) OR (length(repeatSequence) = 1 AND repeatCount >= 5))",
    " GROUP BY sampleId",
    sep = "")
  score = dbGetQuery(dbConnect, query)
  score$msiScore = score$msiScore/3095
  score$msiStatus = ifelse(score$msiScore > 0.909, "MSI", "MSS")
  return (score)
}

query_msi<-function(dbConnect) {
  query = paste(
    "SELECT sampleId, count(*) as msiScore ",
    " FROM somaticVariant ",
    "WHERE filter = 'PASS'",
    "   AND type = 'INDEL'",
    "   AND length(alt) <= 50",
    "   AND length(ref) <= 50",
    "   AND repeatCount >= 4",
    "   AND ((length(repeatSequence) BETWEEN 2 AND 4) OR (length(repeatSequence) = 1 AND repeatCount >= 5))",
    " GROUP BY sampleId",
    sep = "")
  score = dbGetQuery(dbConnect, query)
  score$msiScore = score$msiScore/3095
  score$msiStatus = ifelse(score$msiScore > 0.909, "MSI", "MSS")
  return (score)
}

query_unsupported_segments<-function(dbConnect, sampleId) {
  query = paste(
    "SELECT sampleId, count(*) as unsupportedSegments ",
    " FROM copyNumber ",
    "WHERE sampleId = '",sampleId, "'",
    "  AND segmentStartSupport = 'NONE'",
    " GROUP BY sampleId",
    sep = "")
  return (dbGetQuery(dbConnect, query))
}

query_structural_variant_overview<-function(dbConnect, sampleId) {
  query = paste(
    "SELECT sampleId, type, count(*) as count",
    " FROM structuralVariant ",
    "WHERE sampleId = '", sampleId, "'",
    " GROUP BY sampleId, type",
    sep = "")
  result = data.table(dbGetQuery(dbConnect, query))
  structuralVariants = result[, list(structuralVariants=sum(count)), by=c("sampleId")]

  types = dcast(result, sampleId ~ type, value.var = "count")
  #colnames(types) <- paste("sv", colnames(types), sep="")
  #colnames(types)[1] <- "sampleId"
  structuralVariants = left_join(structuralVariants, types)

  return (structuralVariants)
}


query_somatic_overview_old<-function(dbConnect, sampleId) {
  query = paste(
    "SELECT sampleId, clonality, type, count(*) as count, SUM(IF (effect like '%missense%', 1, 0)) AS mutationalLoad",
    " FROM somaticVariant ",
    "WHERE filter = 'PASS' ",
    "   AND sampleId = '", sampleId, "'",
    " GROUP BY sampleId, clonality, type",
    sep = "")

  result = data.table(dbGetQuery(dbConnect, query))

  type = result[, list(count=sum(count)), by=c("sampleId","type")]
  typeLong = dcast(type, sampleId ~ type, value.var="count")

  clonality = result[, list(count=sum(count)), by=c("sampleId","clonality")]
  clonalityLong = dcast(clonality, sampleId ~ clonality, value.var="count")

  mutationalLoad = result[, list(mutationalLoad=sum(mutationalLoad)), by=c("sampleId")]

  somatics = result[, list(somaticVariants=sum(count)), by=c("sampleId")]
  somatics$mutationalLoad = mutationalLoad$mutationalLoad
  somatics = left_join(somatics, typeLong)
  somatics = left_join(somatics, clonalityLong)

  return (somatics)
}

query_somatic_overview<-function(dbConnect, sampleId) {
  query = paste(
    "SELECT sampleId, clonality, type, count(*) as count, SUM(IF (effect like '%missense%', 1, 0)) AS mutationalLoad",
    " FROM somaticVariant ",
    "WHERE filter = 'PASS' ",
    "   AND sampleId = '", sampleId, "'",
    " GROUP BY sampleId, clonality, type",
    sep = "")

  result = data.table(dbGetQuery(dbConnect, query))

  type = result[, list(count=sum(count)), by=c("sampleId","type", "clonality")]
  typeLong = dcast(type, sampleId ~ type+clonality, value.var="count")

  mutationalLoad = result[, list(mutationalLoad=sum(mutationalLoad)), by=c("sampleId")]

  somatics = result[, list(somaticVariants=sum(count)), by=c("sampleId")]
  somatics$mutationalLoad = mutationalLoad$mutationalLoad
  somatics = left_join(somatics, typeLong)

  return (somatics)
}

query_patient_id_lookup<-function(dbConnect) {
  query = paste(
    "SELECT CPCTCNT.patientId AS sampleId, concat('CPCT02',CPCTCNT.itemValue , IF(length(CPCTPN.itemValue) = 2, '00' + CPCTPN.itemValue, IF(left(CPCTPN.itemValue,1) = '-', right(CPCTPN.itemValue, 4), CPCTPN.itemValue ))) as patientId",
    "  FROM drupEcrf CTCT2YN,  drupEcrf CPCTCNT, drupEcrf CPCTPN ",
    " WHERE CPCTCNT.patientId = CPCTPN.patientId AND CTCT2YN.patientId = CPCTCNT.patientId ",
    "   AND CPCTCNT.item = 'FLD.REG.CPCTCNT' AND CPCTCNT.itemValue != ''",
    "   AND CPCTPN.item = 'FLD.REG.CPCTPN' AND CPCTPN.itemValue != ''",
    "   AND CTCT2YN.item = 'FLD.REG.CTCT2YN' AND CTCT2YN.itemValue = 'Yes'",
    sep = "")

  result = dbGetQuery(dbConnect, query)
  return (rbind(manual_patient_id(), result))
}

query_whole_genome_duplication<-function(dbConnect, sampleId) {
  query = paste(
    "SELECT sampleId, count(*) as duplicatedAutosomes",
    "  FROM ",
    " (select sampleId, chromosome, round(sum(bafCount*actualBaf*copyNumber)/sum(bafCount),1) as lwMajorAlleleAvg ",
    " from copyNumber ",
    " where sampleId = '", sampleId, "' and chromosome not in ('X', 'Y') ",
    " group by sampleId, chromosome ",
    " HAVING lwMajorAlleleAvg > 1.5) a ",
    " GROUP BY sampleId",
    sep = "")

  result = dbGetQuery(dbConnect, query)
  result$WGD <- ifelse(result$duplicatedAutosomes > 10, TRUE, FALSE)

  return (result)
}

query_snps_cohort<-function(dbConnect, cohort) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT s.sampleId, chromosome as chr, position as pos, ref, alt",
    "  FROM somaticVariant s ",
    " WHERE filter = 'PASS'",
    "   AND type = 'SNP'",
    "   AND gene <> ''",
    "   AND s.sampleId in (", sampleIdString, ")",
    sep = "")

  return (dbGetQuery(dbConnect, query))
}

query_snps_sample<-function(dbConnect, sample) {
  query = paste(
    "SELECT s.sampleId, chromosome as chr, position as pos, ref, alt",
    "  FROM somaticVariant s ",
    " WHERE filter = 'PASS'",
    "   AND type = 'SNP'",
    "   AND gene <> ''",
    "   AND s.sampleId = '", sample, "'",
    sep = "")

  return (dbGetQuery(dbConnect, query))
}

query_somatic_drivers<-function(dbConnect, cohort, genes) {
  geneString = paste("'", genes$gene, "'", collapse = ",", sep = "")
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT sampleId, gene, chromosome, position, adjustedCopyNumber as copyNumber, ",
    "       ROUND(adjustedVAF * adjustedCopyNumber,2) as ploidy, ",
    "       clonality, loh",
    "  FROM somaticVariant s ",
    " WHERE sampleId in (",sampleIdString, ")",
    "   AND filter = 'PASS'",
    "   AND effect not in ('non coding exon variant', 'synonymous variant', 'UTR variant', 'sequence feature','intron variant')",
    "   AND gene in (", geneString, ")",
    sep = "")
  return (dbGetQuery(dbConnect, query))
}


query_gene_panel<-function(dbConnect, panel = "HMF Paper") {
  query = paste(
    "SELECT gene ",
    "  FROM genePanel ",
    " WHERE panel='", panel, "'",
    sep = "")

  return (dbGetQuery(dbConnect, query))
}

