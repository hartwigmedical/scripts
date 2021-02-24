CREATE OR REPLACE VIEW actionability AS (

SELECT
    profileName AS tumorProfile, therapyName AS treatment, name AS cancerType, termId AS cancerTypeId,
    evidenceType, responseType, ampCapAscoEvidenceLevel, approvalStatus, efficacyEvidence
FROM ckbEntry
INNER JOIN evidence ON evidence.ckbEntryId = ckbEntry.id
INNER JOIN therapyEvidence ON therapyEvidence.evidenceId = evidence.id
INNER JOIN therapy ON therapyEvidence.therapyId = therapy.id
INNER JOIN indicationEvidence ON indicationEvidence.evidenceId = evidence.id
INNER JOIN indication ON indicationEvidence.indicationId = indication.id);