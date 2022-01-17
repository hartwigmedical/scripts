CREATE OR REPLACE VIEW actionability AS (

SELECT
    ckbEntry.id as ckbEntryId, therapy.createDate as therapyCreateDate, therapy.updateDate as therapyUpdateDate,
    profileName AS tumorProfile, therapyName AS treatment, group_concat(treatmentApproachDrugClass.drugClass) as drugClass,
    name AS cancerType, termId AS cancerTypeId,evidenceType, responseType, ampCapAscoEvidenceLevel, approvalStatus,
    efficacyEvidence, description,group_concat(DISTINCT pubmedId) as pubmeds
FROM ckbEntry
INNER JOIN evidence ON evidence.ckbEntryId = ckbEntry.id
LEFT JOIN evidenceReference on evidenceReference.evidenceId = evidence.id
INNER JOIN therapyEvidence ON therapyEvidence.evidenceId = evidence.id
INNER JOIN therapy ON therapyEvidence.therapyId = therapy.id
INNER JOIN indicationEvidence ON indicationEvidence.evidenceId = evidence.id
INNER JOIN indication ON indicationEvidence.indicationId = indication.id
INNER JOIN treatmentApproachEvidence ON treatmentApproachEvidence.treatmentApproachEvidenceId = evidence.id
INNER JOIN treatmentApproach ON treatmentApproachEvidence.treatmentApproachEvidenceId = treatmentApproach.id
INNER JOIN treatmentApproachDrugClass ON treatmentApproachDrugClass.treatmentApproachId = treatmentApproach.id
GROUP BY 1,2,3,4,5,6,7,8,9,10,11,12,13,14);