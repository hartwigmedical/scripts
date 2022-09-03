CREATE OR REPLACE VIEW actionability AS (

SELECT
    ckbEntry.id as ckbEntryId, ckbEntry.ckbProfileId, profileName AS tumorProfile,
    therapyEvidence.therapyId, therapy.ckbTherapyId, therapy.createDate as therapyCreateDate, therapy.updateDate as therapyUpdateDate, therapyName AS treatment, description,
    indicationEvidence.indicationId, indication.ckbIndicationId, termId AS cancerTypeId, name AS cancerType,
    treatmentApproachEvidence.treatmentApproachEvidenceId, treatmentApproachDrugClass.treatmentApproachId, group_concat(drugClass) as drugClasses,
    evidence.ckbEvidenceId, evidenceType, responseType, ampCapAscoEvidenceLevel as level, approvalStatus, efficacyEvidence,
    evidenceReference.ckbReferenceId, group_concat(DISTINCT pubmedId) as pubmeds
FROM ckbEntry
INNER JOIN evidence ON evidence.ckbEntryId = ckbEntry.id
LEFT JOIN evidenceReference on evidenceReference.evidenceId = evidence.id
INNER JOIN therapyEvidence ON therapyEvidence.evidenceId = evidence.id
INNER JOIN therapy ON therapyEvidence.therapyId = therapy.id
INNER JOIN indicationEvidence ON indicationEvidence.evidenceId = evidence.id
INNER JOIN indication ON indicationEvidence.indicationId = indication.id
LEFT JOIN treatmentApproachEvidence ON treatmentApproachEvidence.evidenceId = evidence.id
LEFT JOIN treatmentApproachDrugClass ON treatmentApproachDrugClass.treatmentApproachId=treatmentApproachEvidence.treatmentApproachEvidenceId
group by ckbEntry.id, ckbEntry.ckbProfileId, therapyEvidence.therapyId, therapy.ckbTherapyId, indicationEvidence.indicationId,
indication.ckbIndicationId, evidence.ckbEvidenceId);