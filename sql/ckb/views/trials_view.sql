CREATE OR REPLACE VIEW trials AS (

SELECT
    ckbEntry.id as ckbEntryId, ckbEntry.ckbProfileId, profileName AS tumorProfile, requirementType,
	clinicalTrial.nctId, clinicalTrial.title, clinicalTrial.phase, clinicalTrial.recruitment, clinicalTrial.gender, clinicalTrial.sponsors,clinicalTrial.variantRequirement,
    group_concat(distinct ageGroup.ageGroup) as ageGroups, group_concat(distinct location.country) as countries,
    therapyName,indication.name as cancerType
--    treatmentApproachEvidence.treatmentApproachEvidenceId, treatmentApproachDrugClass.treatmentApproachId, group_concat(distinct drugClass) as drugClasses
FROM ckbEntry
LEFT JOIN clinicalTrial ON clinicalTrial.ckbEntryId = ckbEntry.id
LEFT JOIN variantRequirementDetail ON variantRequirementDetail.clinicalTrialId = clinicalTrial.id
LEFT JOIN ageGroup ON ageGroup.clinicalTrialId = clinicalTrial.id
LEFT JOIN location ON location.clinicalTrialId = clinicalTrial.id
LEFT JOIN therapyClinicalTrial ON therapyClinicalTrial.clinicalTrialId = clinicalTrial.id
LEFT JOIN therapy ON therapyClinicalTrial.therapyId = therapy.id
LEFT JOIN indicationClinicalTrial ON indicationClinicalTrial.clinicalTrialId = clinicalTrial.id
LEFT JOIN indication ON indicationClinicalTrial.indicationId = indication.id
--LEFT JOIN treatmentApproachEvidence ON treatmentApproachEvidence.evidenceId = clinicalTrial.id
--LEFT JOIN treatmentApproachDrugClass ON treatmentApproachDrugClass.treatmentApproachId=treatmentApproachEvidence.treatmentApproachEvidenceId
group by clinicalTrial.nctId, ckbEntry.ckbProfileId, therapyClinicalTrial.therapyId, indicationClinicalTrial.indicationId);