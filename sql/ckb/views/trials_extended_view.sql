CREATE OR REPLACE VIEW trialsExtended AS (

SELECT DISTINCT
        clinicalTrial.nctId, clinicalTrial.title, clinicalTrial.phase, clinicalTrial.recruitment, clinicalTrial.gender,
        clinicalTrial.variantRequirement, ckbEntry.profileName, variantRequirementDetail.requirementType,
        group_concat(DISTINCT ageGroup.ageGroup separator ',') AS ageGroups,
        group_concat(DISTINCT therapy.therapyName separator ',') AS therapyNames,
        group_concat(DISTINCT indication.name separator ',') AS cancerTypes
    FROM clinicalTrial
    INNER JOIN ckbEntry ON clinicalTrial.ckbEntryId = ckbEntry.id
    INNER JOIN variantRequirementDetail ON
        variantRequirementDetail.clinicalTrialId = clinicalTrial.id AND variantRequirementDetail.ckbProfileId = ckbEntry.ckbProfileId
    LEFT JOIN ageGroup ON ageGroup.clinicalTrialId = clinicalTrial.id
    LEFT JOIN therapyClinicalTrial ON therapyClinicalTrial.clinicalTrialId = clinicalTrial.id
    LEFT JOIN therapy ON therapyClinicalTrial.therapyId = therapy.id
    LEFT JOIN indicationClinicalTrial ON indicationClinicalTrial.clinicalTrialId = clinicalTrial.id
    LEFT JOIN indication ON indicationClinicalTrial.indicationId = indication.id
    GROUP BY 1,2,3,4,5,6,7,8
);