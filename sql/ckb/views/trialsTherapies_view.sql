CREATE OR REPLACE VIEW trialsTherapies AS (

SELECT DISTINCT
    clinicalTrial.updateDate AS clinicalTrialUpdateDate, nctId, therapy.updateDate AS therapyUpdateDate, therapyName, description
FROM clinicalTrial
LEFT JOIN therapyClinicalTrial ON clinicalTrial.id = therapyClinicalTrial.clinicalTrialId
LEFT JOIN therapy ON therapyClinicalTrial.therapyId = therapy.id)