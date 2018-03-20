# registratie datum
select distinct patientId, registrationDate from clinical
inner join patient on patient.cpctId=clinical.patientId order by 1;

# Informed consent is eerste moment in de tijd
select * from clinicalFindings where ecrfItem = 'FLD.INFORMEDCONSENT.ICDTC;FLD.BIOPS.BIOPTDT' and
message ='at least 1 biopsy taken more than 30 days prior to informed consent date' and formLocked = 'true';

select * from clinicalFindings where ecrfitem ='FLD.INFORMEDCONSENT.ICDTC' and message ='informed consent date empty or in wrong format' and formLocked='true';

# Ziekenhuis is bekend?
select distinct clinical.patientId, hospital, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from clinical
inner join ecrf on ecrf.patientId = clinical.patientId
where isnull(hospital) and ecrf.item = 'FLD.SELCRIT.NHOSPITAL' and ecrf.locked='true';

# Geslacht is bekend?
select * from clinicalFindings where ecrfItem ='FLD.DEMOGRAPHY.SEX' and message = 'gender empty' and formLocked = 'true';

# Geboortejaar is bekend?
select * from clinicalFindings where ecrfItem ='FLD.SELCRIT.NBIRTHYEAR;FLD.ELIGIBILITY.BIRTHYEAR;FLD.ELIGIBILITY.BIRTHDTCES' and
message = 'birth year could not be determined' and formLocked = 'true';

# Kan de tumor locatie cureren?
select * from clinicalFindings where ecrfitem ='FLD.CARCINOMA.PTUMLOC;FLD.CARCINOMA.PTUMLOCS' and message ='primary tumor location empty' and formLocked = 'true';

# Systemic pre-therapy is bekend?
select patient.cpctId, clinical.hasSystemicPretreatment, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.hasSystemicPretreatment is null and ecrf.item='FLD.PRETHERAPY.SYSTEMIC' and locked ='true';

# Radiotherapy pre-therapy is bekend?
select patient.cpctId, clinical.hasRadiotherapyPreTreatment, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.hasRadiotherapyPreTreatment is null and ecrf.item ='FLD.PRETHERAPY.RADIOTHER' and locked ='true';

# Kan alle pre therapies cureren?
select patient.cpctId, clinical.preTreatments, clinical.hasSystemicPretreatment, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.preTreatments is null and clinical.hasSystemicPretreatment = 'yes' and ecrf.item='FLD.PRETHERAPY.SYSTEMICREG' and locked ='true';

