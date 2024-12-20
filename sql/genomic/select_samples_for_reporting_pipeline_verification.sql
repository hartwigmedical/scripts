# low purity + low TMB status
select modified, sampleId, tmbPerMb, tmbStatus  from purity where tmbStatus = "LOW" and qcStatus = "WARN_LOW_PURITY" order by modified desc;

# high tmb + high msi + somatic mutation
select purity.modified, purity.sampleId, tmbPerMb, tmbStatus, msIndelsPerMb, msStatus from purity
inner join driverCatalog on purity.sampleId = driverCatalog.sampleId
where tmbStatus = "HIGH" and msStatus = "MSI" and qcStatus = "PASS" and driver = "MUTATION" order by modified desc;

# hrd + germline mutation
select purity.modified, chord.sampleId, hrStatus, gene, canonicalTranscript, category, driver, likelihoodMethod from chord
inner join purity on purity.sampleId = chord.sampleId
inner join driverCatalog on purity.sampleId = driverCatalog.sampleId
where hrStatus = "HR_DEFICIENT" and driver = "GERMLINE_MUTATION" order by modified desc;

# mss + germline disruption
select purity.modified, purity.sampleId, msIndelsPerMb, msStatus, gene, canonicalTranscript, category, driver, likelihoodMethod  from purity
inner join driverCatalog on purity.sampleId = driverCatalog.sampleId
where msStatus = "MSS" and driver  = "GERMLINE_DISRUPTION" order by modified desc;

#del, amp, partial amp, disruption, hom_del_disruption
SELECT t.modified, t.sampleId, t.gene, t.canonicalTranscript, t.category, t.driver, t.likelihoodMethod
FROM driverCatalog t
WHERE t.driver IN ('AMP', 'DEL', 'DISRUPTION', 'HOM_DEL_DISRUPTION')
  AND t.sampleId IN (
    SELECT sampleId
    FROM driverCatalog
    WHERE driver IN ('AMP', 'DEL', 'DISRUPTION', 'HOM_DEL_DISRUPTION')
    GROUP BY sampleId
    HAVING COUNT(DISTINCT driver) = 4
  );

# germline deletion
select *  from  driverCatalog
inner join chord on driverCatalog.sampleId = chord.sampleId
where driver  = "GERMLINE_DELETION" and hrStatus = "HR_PROFICIENT" order by modified desc;

#hom dup disruption
select modified, sampleId,gene, canonicalTranscript, category, driver, likelihoodMethod from driverCatalog where driver = "HOM_DUP_DISRUPTION" order by modified desc;

# virus
select modified, sampleId, virusName, qcStatus, integrations, interpretation, likelihood, reported from virusAnnotation where reported and interpretation = "HPV" and likelihood = "HIGH" order by modified desc;

#fusion
SELECT *
FROM svFusion t
WHERE t.reportedType IN ('KNOWN_PAIR', 'PROMISCUOUS_3', 'PROMISCUOUS_5', 'EXON_DEL_DUP') and reported
  AND t.sampleId IN (
    SELECT sampleId
    FROM svFusion
    WHERE reportedType IN ('KNOWN_PAIR', 'PROMISCUOUS_3', 'PROMISCUOUS_5', 'EXON_DEL_DUP') and reported
    GROUP BY sampleId
    HAVING COUNT(DISTINCT reportedType) = 1
  );

# HRD cannot be determined
select modified, chord.sampleId, hrStatus from chord inner join purity on purity.sampleId = chord.sampleId where hrStatus = "CANNOT_BE_DETERMINED" order by modified desc;