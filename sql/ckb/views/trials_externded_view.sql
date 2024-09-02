CREATE OR REPLACE VIEW trialsExtended AS (

SELECT DISTINCT
        clinicalTrial.nctId, clinicalTrial.title, clinicalTrial.phase, clinicalTrial.recruitment, clinicalTrial.gender, group_concat(distinct ageGroup.ageGroup) as ageGroups,
        clinicalTrial.variantRequirement, requirementType, profileName AS tumorProfile,
		therapyName,indication.name as cancerType
    FROM ckbEntry
    INNER JOIN clinicalTrial ON clinicalTrial.ckbEntryId = ckbEntry.id
    INNER JOIN variantRequirementDetail ON variantRequirementDetail.clinicalTrialId = clinicalTrial.id and variantRequirementDetail.ckbProfileId = ckbEntry.ckbProfileId
    INNER JOIN ageGroup ON ageGroup.clinicalTrialId = clinicalTrial.id
    INNER JOIN therapyClinicalTrial ON therapyClinicalTrial.clinicalTrialId = clinicalTrial.id
    INNER JOIN therapy ON therapyClinicalTrial.therapyId = therapy.id
    INNER JOIN indicationClinicalTrial ON indicationClinicalTrial.clinicalTrialId = clinicalTrial.id
    INNER JOIN indication ON indicationClinicalTrial.indicationId = indication.id
    group by clinicalTrial.nctId, ckbEntry.ckbProfileId, therapyClinicalTrial.therapyId, indicationClinicalTrial.indicationId