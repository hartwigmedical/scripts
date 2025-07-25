#1 high tmb + high msi + somatic mutation
select purity.modified, purity.sampleId, tmbPerMb, tmbStatus, msIndelsPerMb, msStatus from purity
inner join driverCatalog on purity.sampleId = driverCatalog.sampleId
where tmbStatus = "HIGH" and msStatus = "MSI" and qcStatus = "PASS" and driver = "MUTATION" order by modified desc;

#2 hrd + germline mutation
select purity.modified, chord.sampleId, hrStatus, gene, canonicalTranscript, category, driver, likelihoodMethod from chord
inner join purity on purity.sampleId = chord.sampleId
inner join driverCatalog on purity.sampleId = driverCatalog.sampleId
where hrStatus = "HR_DEFICIENT" and driver = "GERMLINE_MUTATION" order by modified desc;

#3 mss + germline disruption
select purity.modified, purity.sampleId, msIndelsPerMb, msStatus, gene, canonicalTranscript, category, driver, likelihoodMethod  from purity
inner join driverCatalog on purity.sampleId = driverCatalog.sampleId
where msStatus = "MSS" and driver  = "GERMLINE_DISRUPTION" order by modified desc;

#4 del, amp, partial amp, disruption, hom_del_disruption
SELECT t.modified, t.sampleId, t.gene, t.canonicalTranscript, t.category, t.driver, t.likelihoodMethod
FROM driverCatalog t
WHERE t.driver IN ('AMP', 'DEL', 'DISRUPTION', 'HOM_DEL_DISRUPTION')
  AND t.sampleId IN (
    SELECT sampleId
    FROM driverCatalog
    WHERE driver IN ('AMP', 'DEL', 'DISRUPTION', 'HOM_DEL_DISRUPTION')
    GROUP BY sampleId
    HAVING COUNT(DISTINCT driver) = 4
  ) order by t.modified desc;

#5 germline deletion
select *  from  driverCatalog
inner join chord on driverCatalog.sampleId = chord.sampleId
where driver  = "GERMLINE_DELETION" and hrStatus = "HR_PROFICIENT" order by modified desc;

#6 hom dup disruption
select modified, sampleId,gene, canonicalTranscript, category, driver, likelihoodMethod from driverCatalog where driver = "HOM_DUP_DISRUPTION" order by modified desc;

#7 virus
select modified, sampleId, virusName, qcStatus, integrations, interpretation, likelihood, reported from virusAnnotation where reported and interpretation = "HPV" and likelihood = "HIGH" order by modified desc;

#8 fusion
SELECT *
FROM svFusion t
WHERE t.reportedType IN ('KNOWN_PAIR', 'PROMISCUOUS_3', 'PROMISCUOUS_5', 'EXON_DEL_DUP') and reported
  AND t.sampleId IN (
    SELECT sampleId
    FROM svFusion
    WHERE reportedType IN ('KNOWN_PAIR', 'PROMISCUOUS_3', 'PROMISCUOUS_5', 'EXON_DEL_DUP') and reported
    GROUP BY sampleId
    HAVING COUNT(DISTINCT reportedType) = 1
  ) order by modified desc;

#9 HRD cannot be determined + # low purity
select modified, chord.sampleId, hrStatus, qcStatus, purity from chord inner join purity on purity.sampleId = chord.sampleId where hrStatus = "CANNOT_BE_DETERMINED" and qcStatus = "WARN_LOW_PURITY" order by modified desc;

#10 germline variant without tumor support
select * from germlineVariant where reported and variantCopyNumber <= 0.5 order by modified desc;

#select XX samples for verification specific samples

#select 1 lab failure

#select 1 analyse failure
select modified, sampleId, qcStatus, purity from purity  where qcStatus = "FAIL_NO_TUMOR" order by modified desc;

#select 1 panel analyse

#select 1 panel lab failure