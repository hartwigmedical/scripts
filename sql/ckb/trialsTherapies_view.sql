CREATE OR REPLACE VIEW trialsTherapies AS (

select distinct clinicalTrial.updateDate as clinicalTrialUpdateDate, nctId, therapy.updateDate as therapyUpdateDate, therapyName, description  from clinicalTrial
left join therapyClinicalTrial on clinicalTrial.id=therapyClinicalTrial.clinicalTrialId
left join therapy on therapyClinicalTrial.therapyId=therapy.id)