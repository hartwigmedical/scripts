CREATE OR REPLACE VIEW trials AS (

SELECT DISTINCT
	clinicalTrial.id, profileName AS tumorProfile, clinicalTrial.updateDate AS trialUpdateDate,
	clinicalTrial.nctId, title, phase, recruitment, gender
FROM ckbEntry INNER JOIN clinicalTrial ON clinicalTrial.ckbEntryId = ckbEntry.id);