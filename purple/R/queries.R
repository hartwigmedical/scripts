query_highest_purity_cohort<-function(dbConnect) {

  cohort = purple::query_purity(dbConnect)

  # PatientIds
  patientIdLookups = query_patient_id_lookup(dbConnect)
  cohort$patientId <- sapply(cohort$sampleId, function(x) {purple::sample_to_patient_id(x, patientIdLookups)})

  #Clinical Data
  clinicalData = purple::query_clinical_data(dbConnect)
  cohort = left_join(cohort, clinicalData[, c("sampleId", "cancerType")])

  # Cohort
  highestPurityCohort = purple::highest_purity_patients(cohort)

  return (highestPurityCohort)
}

query_multiple_biopsy_cohort<-function(dbConnect) {

  cohort = purple::query_purity(dbConnect)

  # PatientIds
  patientIdLookups = query_patient_id_lookup(dbConnect)
  cohort$patientId <- sapply(cohort$sampleId, function(x) {purple::sample_to_patient_id(x, patientIdLookups)})


  #Clinical Data
  clinicalData = purple::query_clinical_data(dbConnect)
  cohort = left_join(cohort, clinicalData[, c("sampleId", "cancerType")])

  # Cohort
  multipleBiopsyCohort = purple::multiple_biopsy(cohort)

  return (multipleBiopsyCohort)
}

query_gene_copy_number_knockout<-function(dbConnect, cohort) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT g.sampleId, g.chromosome, g.start, g.end, g.gene, g.chromosomeBand, g.minCopyNumber, if(minCopyNumber < 0.5 || nonsenseBiallelicVariants + missenseBiallelicVariants + spliceBiallelicVariants > 0, 1, 0) as score, ",
    "       g.minCopyNumber < 0.5 as deleted, g.minMinorAllelePloidy < 0.5 as loh, nonsenseBiallelicVariants + missenseBiallelicVariants + spliceBiallelicVariants > 0 as biallelicVariant ",
    "  FROM geneCopyNumber g, purity p",
    " WHERE g.sampleId = p.sampleId",
    "   AND p.qcStatus = 'PASS'",
    "   AND p.status != 'NO_TUMOR'",
    "   AND p.purity >= 0.20",
    "   AND g.germlineHetRegions = 0",
    "   AND g.germlineHomRegions = 0",
    "   AND (g.minCopyNumber < 0.5 or g.minMinorAllelePloidy < 0.5 or nonsenseBiallelicVariants + missenseBiallelicVariants + spliceBiallelicVariants > 0)" ,
    "   AND g.chromosome <> 'Y'",
    "   AND p.sampleId in (",sampleIdString, ")",
    sep = " ")
  return (dbGetQuery(dbConnect, query))
}

query_gene_copy_number_loh<-function(dbConnect, cohort) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT g.sampleId, g.chromosome, g.start, g.end, g.gene, g.chromosomeBand, g.minCopyNumber, 1 as score",
    "  FROM geneCopyNumber g, purity p",
    " WHERE g.sampleId = p.sampleId",
    "   AND p.qcStatus = 'PASS'",
    "   AND p.status != 'NO_TUMOR'",
    "   AND p.purity >= 0.20",
    "   AND g.germlineHetRegions = 0",
    "   AND g.germlineHomRegions = 0",
    "   AND g.minMinorAllelePloidy < 0.5",
    "   AND g.chromosome <> 'Y'",
    "   AND p.sampleId in (",sampleIdString, ")",
    sep = " ")
  return (dbGetQuery(dbConnect, query))
}

query_gene_copy_number_deletes<-function(dbConnect, cohort) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT g.sampleId, g.chromosome, g.start, g.end, g.gene, g.chromosomeBand, g.minCopyNumber, 1 as score, g.minRegionStartSupport, g.minRegionEndSupport",
    "  FROM geneCopyNumber g, purity p",
    " WHERE g.sampleId = p.sampleId",
    "   AND p.qcStatus = 'PASS'",
    "   AND p.status != 'NO_TUMOR'",
    "   AND g.germlineHetRegions = 0",
    "   AND g.germlineHomRegions = 0",
    "   AND g.minCopyNumber < 0.5",
    "   AND g.chromosome <> 'Y'",
    "   AND p.sampleId in (",sampleIdString, ")",
    sep = " ")
  return (dbGetQuery(dbConnect, query))
}

query_gene_copy_number_amplifications<-function(dbConnect, cohort, cutoff = 3) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT g.sampleId, g.chromosome, g.start, g.end, g.gene, g.chromosomeBand, g.minCopyNumber, log2(2 * g.minCopyNumber / p.ploidy / ", cutoff, ") as score, g.minRegionStartSupport, g.minRegionEndSupport, p.ploidy",
    "  FROM geneCopyNumber g, purity p",
    " WHERE g.sampleId = p.sampleId",
    "   AND p.qcStatus = 'PASS'",
    "   AND p.status != 'NO_TUMOR'",
    "   AND g.germlineHetRegions = 0",
    "   AND g.germlineHomRegions = 0",
    "   AND g.minCopyNumber / p.ploidy > ", cutoff,
    "   AND p.sampleId in (",sampleIdString, ")",
    sep = " ")
  return (dbGetQuery(dbConnect, query))
}

query_gene_copy_number_drivers<-function(dbConnect, cohort, cutoff = 3) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT g.*, p.ploidy",
    "  FROM geneCopyNumber g, purity p",
    " WHERE g.sampleId = p.sampleId",
    "   AND p.qcStatus = 'PASS'",
    "   AND p.status != 'NO_TUMOR'",
    "   AND g.germlineHetRegions = 0",
    "   AND g.germlineHomRegions = 0",
    "   AND (g.minCopyNumber / p.ploidy > ", cutoff, " OR g.minCopyNumber < 0.5 OR g.minMinorAllelePloidy < 0.5)",
    "   AND p.sampleId in (",sampleIdString, ")",
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

query_purity<-function(dbConnect, purityCutoff = 0.2) {
  query = paste(
    "SELECT p.*",
    " FROM purity p",
    "WHERE qcStatus = 'PASS'",
    "  AND status <> 'NO_TUMOR'",
    "  AND p.purity >= ", purityCutoff,
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

### REMOVE??
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
    "SELECT CPCTCNT.patientId AS sampleId, concat('CPCT02', CPCTCNT.itemValue, LPAD(RIGHT(CPCTPN.itemValue,4), 4, '0')) as patientId",
    "  FROM drupEcrf CTCT2YN,  drupEcrf CPCTCNT, drupEcrf CPCTPN ",
    " WHERE CPCTCNT.patientId = CPCTPN.patientId AND CTCT2YN.patientId = CPCTCNT.patientId ",
    "   AND CPCTCNT.item = 'FLD.REG.CPCTCNT' AND CPCTCNT.itemValue != ''",
    "   AND CPCTPN.item = 'FLD.REG.CPCTPN' AND CPCTPN.itemValue != ''",
    "   AND CTCT2YN.item = 'FLD.REG.CTCT2YN' AND CTCT2YN.itemValue = 'Yes'",
    sep = "")

  result = dbGetQuery(dbConnect, query)
  return (result)
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

### NO LONGER USED???
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

query_somatic_variants <- function(dbConnect, cohort, filterEmptyGenes = FALSE, gene = NA) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT chromosome, position, sampleId, ref, alt, trinucleotideContext, adjustedVaf * adjustedCopyNumber as ploidy, clonality, type, biallelic, gene, canonicalEffect, canonicalCodingEffect",
    "FROM somaticVariant",
    "WHERE filter = 'PASS'",
    "  AND sampleId in (",sampleIdString, ")",
    sep = " ")

  if (filterEmptyGenes) {
    query = paste(query, " AND gene <> ''", sep = " ")
  }

  if (!is.na(gene)) {
    query = paste(query, " AND gene = '", gene ,"'",sep = "")
  }

  return (dbGetQuery(dbConnect, query))
}


query_structural_variants <- function(dbConnect, cohort) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT sampleId, startChromosome, endChromosome, startPosition, endPosition, startOrientation, endOrientation, type, ploidy",
    "FROM structuralVariant",
    "WHERE sampleId in (",sampleIdString, ")",
    sep = " ")

  return(dbGetQuery(dbConnect, query))
}

query_canonical_transcript <-function(dbConnect) {
  query = paste(
    "SELECT *",
    "FROM canonicalTranscript",
    sep = " ")

  return(dbGetQuery(dbConnect, query))
}

query_transcriptId_from_proteinId<-function(dbConnect, proteinIds) {
  sampleIdString = paste("'", proteinIds, "'", collapse = ",", sep = "")
  query = paste(
  " select t1.stable_id as transcriptId, t2.stable_id as proteinId from transcript t1, translation t2 where t1.transcript_id = t2.transcript_id  ",
  " AND t2.stable_id in (",sampleIdString, ")",
  sep = " ")
  return(dbGetQuery(dbConnect, query))
}

query_exons_from_ensembl<-function(dbConnect) {
  query = paste(
    "select",
    " seq_region.name as chromosome,",
    " gene.seq_region_start as gene_start,",
    " gene.seq_region_end as gene_end,",
    " gene.stable_id as gene_id,",
    " display_xref.display_label as gene_name,",
    " GROUP_CONCAT(DISTINCT entrez_xref.dbprimary_acc ORDER BY entrez_xref.dbprimary_acc SEPARATOR ',') as entrezId,",
    " GROUP_CONCAT(DISTINCT karyotype.band ORDER BY karyotype.band SEPARATOR '-') as chromosome_band,",
    " t.stable_id as transcript_id,",
    " t.version as transcript_version,",
    " t.seq_region_start as transcript_start,",
    "t.seq_region_end as transcript_end,",
    "e.stable_id as exon_id,",
    "e.seq_region_start as exon_start,",
    "e.seq_region_end as exon_end,",
    "t.seq_region_strand as strand,",
    "if(t.seq_region_strand = -1, ce.seq_region_end - tl.seq_end + 1, cs.seq_region_start + tl.seq_start - 1) as coding_start,",
    "if(t.seq_region_strand = -1, cs.seq_region_end - tl.seq_start + 1, ce.seq_region_start + tl.seq_end - 1) as coding_end",
    "from gene",
    "inner join object_xref on gene.gene_id=object_xref.ensembl_id and object_xref.ensembl_object_type = 'GENE'",
    "inner join xref as display_xref on (display_xref.xref_id=gene.display_xref_id)",
    "inner join karyotype on gene.seq_region_id=karyotype.seq_region_id",
    "inner join transcript t on gene.canonical_transcript_id = t.transcript_id",
    "inner join seq_region on t.seq_region_id = seq_region.seq_region_id",
    "inner join exon_transcript et on et.transcript_id = t.transcript_id",
    "inner join exon e on et.exon_id = e.exon_id",
    "left join xref as entrez_xref on (entrez_xref.xref_id=object_xref.xref_id and entrez_xref.external_db_id = 1300)",
    "left join translation tl on tl.transcript_id = t.transcript_id",
    "left join exon cs on cs.exon_id = tl.start_exon_id",
    "left join exon ce on ce.exon_id = tl.end_exon_id",
    "where",
    "seq_region.name in ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT') and",
    "((gene.seq_region_start >= karyotype.seq_region_start and gene.seq_region_start <= karyotype.seq_region_end)",
    "or (gene.seq_region_end >= karyotype.seq_region_start and gene.seq_region_end <= karyotype.seq_region_end))",
    "and gene.stable_id not in",
    "('ENSG00000250424','ENSG00000257028','ENSG00000184909','ENSG00000181464','ENSG00000244255','ENSG00000263203','ENSG00000268942','ENSG00000233280',",
    "'ENSG00000258465','ENSG00000251246','ENSG00000272414','ENSG00000258728','ENSG00000269783','ENSG00000273266','ENSG00000241489','ENSG00000269881',",
    "'ENSG00000231880','ENSG00000273291','ENSG00000269099','ENSG00000272781','ENSG00000249773','ENSG00000253117','ENSG00000263136','ENSG00000249967',",
    "'ENSG00000269846','ENSG00000259060','ENSG00000255154','ENSG00000270466','ENSG00000262304','ENSG00000268500','ENSG00000262660','ENSG00000258724',",
    "'ENSG00000250264','ENSG00000173366','ENSG00000254692','ENSG00000241690','ENSG00000198211','ENSG00000264668','ENSG00000232748','ENSG00000196826',",
    "'ENSG00000267179','ENSG00000188474','ENSG00000273045','ENSG00000255994','ENSG00000233050','ENSG00000256977','ENSG00000213906','ENSG00000273155',",
    "'ENSG00000228273','ENSG00000262621','ENSG00000233024','ENSG00000214967','ENSG00000272962','ENSG00000184040','ENSG00000173610','ENSG00000273439')",
    "and t.biotype = 'protein_coding'",
    "group by chromosome, gene_start, gene_end, gene_id, gene_name, transcript_id, transcript_version, transcript_start, transcript_end, exon_id, exon_start, exon_end, coding_start, coding_end, strand",
    "order by if(cast(chromosome as SIGNED) = 0, ascii(chromosome), cast(chromosome as SIGNED)), gene_start, gene_id, exon_start",
    sep = " "
  )

  return(dbGetQuery(dbConnect, query))
}

query_genes_from_ensembl<-function(dbConnect) {
  query = paste(
    "SELECT seq_region.name as chromosome, g_xref.display_label as gene_name, g.stable_id as gene_id, g.biotype as g_biotype,  t.stable_id as transcript_id, t.biotype as t_biotype, g.source,",
    "GROUP_CONCAT(DISTINCT entrez_xref.dbprimary_acc ORDER BY entrez_xref.dbprimary_acc SEPARATOR ',') as entrezId,",
    "GROUP_CONCAT(DISTINCT t_xref.dbprimary_acc ORDER BY t_xref.dbprimary_acc SEPARATOR ',') as ccdsId",
    "FROM gene g",
    "inner join xref as g_xref on g_xref.xref_id=g.display_xref_id",
    "inner join transcript t on g.canonical_transcript_id = t.transcript_id",
    "inner join object_xref on g.gene_id=object_xref.ensembl_id and object_xref.ensembl_object_type = 'GENE'",
    "inner join object_xref t_object_xref on t_object_xref.ensembl_id = t.transcript_id",
    "  inner join seq_region on t.seq_region_id = seq_region.seq_region_id",
    "    inner join karyotype on gene.seq_region_id=karyotype.seq_region_id ",
    "left join xref as entrez_xref on (entrez_xref.xref_id=object_xref.xref_id and entrez_xref.external_db_id = 1300)",
    "left join xref t_xref on t_xref.xref_id = t_object_xref.xref_id and t_xref.external_db_id = 3800 ",
    "WHERE seq_region.name in ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT')",
    "and ",
    "((gene.seq_region_start >= karyotype.seq_region_start and gene.seq_region_start <= karyotype.seq_region_end) ",
    "  or (gene.seq_region_end >= karyotype.seq_region_start and gene.seq_region_end <= karyotype.seq_region_end))	",
    "group by 1,2,3,4,5,6,7",
    sep = " ")

  return(dbGetQuery(dbConnect, query))
}



