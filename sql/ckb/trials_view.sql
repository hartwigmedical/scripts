CREATE OR REPLACE VIEW trials AS (

SELECT DISTINCT
	profileName, clinicalTrial.updateDate AS trialUpdateDate, clinicalTrial.nctId, title, phase, recruitment, gender,
	group_concat(DISTINCT ageGroup) AS ageGroups,
	group_concat(DISTINCT country) AS countries
FROM ckbEntry
INNER JOIN clinicalTrial ON clinicalTrial.ckbEntryId = ckbEntry.id
LEFT JOIN ageGroup ON ageGroup.clinicalTrialId = clinicalTrial.id
LEFT JOIN location ON location.clinicalTrialId = clinicalTrial.id
GROUP BY 1,2,3,4,5,6,7);