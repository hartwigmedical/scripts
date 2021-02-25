CREATE OR REPLACE VIEW trials AS (

SELECT DISTINCT
	profileName AS tumorProfile, clinicalTrial.updateDate AS trialUpdateDate, clinicalTrial.nctId, title, phase, recruitment, gender,
	group_concat(DISTINCT ageGroup) AS ageGroups,
	group_concat(DISTINCT country) AS countries,
    group_concat(DISTINCT name) as cancerType
FROM ckbEntry
INNER JOIN clinicalTrial ON clinicalTrial.ckbEntryId = ckbEntry.id
LEFT JOIN indicationClinicalTrial ON indicationClinicalTrial.clinicalTrialId=clinicalTrial.id
LEFT JOIN indication ON indication.id=indicationClinicalTrial.indicationId
LEFT JOIN ageGroup ON ageGroup.clinicalTrialId = clinicalTrial.id
LEFT JOIN location ON location.clinicalTrialId = clinicalTrial.id
GROUP BY 1,2,3,4,5,6,7);