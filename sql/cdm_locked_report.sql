# registratie datum patient
select cpctId, registrationDate from patient;

# biopten afgenomen na informed consent date?
select * from clinicalFindings where ecrfItem = 'FLD.INFORMEDCONSENT.ICDTC;FLD.BIOPS.BIOPTDT' and
message ='at least 1 biopsy taken more than 30 days prior to informed consent date' and formLocked = 'true';

# bevat een informed constent datum in een correct format?
select * from clinicalFindings where ecrfitem ='FLD.INFORMEDCONSENT.ICDTC' and message ='informed consent date empty or in wrong format' and formLocked='true';

# bevat naam van het ziekenhuis?
select distinct clinical.patientId, hospital, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from clinical
inner join ecrf on ecrf.patientId = clinical.patientId
where isnull(hospital) and ecrf.item = 'FLD.SELCRIT.NHOSPITAL' and ecrf.locked='true';

# heeft geslacht?
select * from clinicalFindings where ecrfItem ='FLD.DEMOGRAPHY.SEX' and message = 'gender empty' and formLocked = 'true';

# heeft geboortejaar?
select * from clinicalFindings where ecrfItem ='FLD.SELCRIT.NBIRTHYEAR;FLD.ELIGIBILITY.BIRTHYEAR;FLD.ELIGIBILITY.BIRTHDTCES' and
message = 'birth year could not be determined' and formLocked = 'true';

# heeft tumor locatie?
select * from clinicalFindings where ecrfitem ='FLD.CARCINOMA.PTUMLOC;FLD.CARCINOMA.PTUMLOCS' and message ='primary tumor location empty' and formLocked = 'true';

# Heeft een systemicPreTreatment?
select patient.cpctId, clinical.hasSystemicPretreatment, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.hasSystemicPretreatment is null and ecrf.item='FLD.PRETHERAPY.SYSTEMIC' and locked ='true';

# Heeft RadiotherapyPreTreatment?
select patient.cpctId, clinical.hasRadiotherapyPreTreatment, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.hasRadiotherapyPreTreatment is null and ecrf.item ='FLD.PRETHERAPY.RADIOTHER' and locked ='true';

# Bevat gecureerde pretreatments?
select patient.cpctId, clinical.preTreatments, clinical.hasSystemicPretreatment, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.preTreatments is null and clinical.hasSystemicPretreatment = 'yes' and ecrf.item='FLD.PRETHERAPY.SYSTEMICREG' and locked ='true';

# geen ingevulde ecrf formulier en 1 biopt gesequenced
select * from clinicalFindings where ecrfItem = 'FRM.BIOPS' and message ='less ecrf biopsies than biopsies sequenced.'
and details='ecrf biopsies: 0; sequenced: 1' and formLocked = 'true';

# 2+ biopten gesequenced met genoeg biopt formulieren?
select * from clinicalFindings where ecrfItem = 'FRM.BIOPS' and message ='less ecrf biopsies than biopsies sequenced.'
and details like '%ecrf biopsies: 1%' and formLocked = 'true';

# Heeft een unieke biopt match?
select * from clinicalFindings where ecrfitem ='FRM.BIOPS' and message ='more than 1 possible clinical biopsy match for sequenced sample.' and formLocked = 'true';

# bevat een match van clinical biopt met sequenced sample?
select * from clinicalFindings where ecrfitem ='FRM.BIOPS' and message ='could not match any clinical biopsy with sequenced sample.'
and details not like '%ecrf biopsies: [].%' and formLocked = 'true';

# heeft biopsySite?
select * from clinicalFindings where ecrfItem ='FLD.BIOPS.BILESSITE;FLD.BIOPS.BIOTHLESSITE' and message='biopsy site empty' and formLocked = 'true';

# heeft biopsyLocation?
select * from clinicalFindings where ecrfItem ='FLD.BIOPS.BILESLOC' and message='biopsy location empty' and formLocked = 'true';

# bevat een biopt datum in een correct format?
select * from clinicalFindings where ecrfitem ='FLD.BIOPS.BIOPTDT' and message ='biopsy date empty or in wrong format' and formLocked='true';

# heeft een gelijk sampling datum en biopt datum?
select * from clinicalFindings where ecrfitem ='FRM.BIOPS' and message ='sampling date does not equal biopsy date in matched biopsy' and formLocked = 'true';

# bevat treatmentGiven?
select patient.cpctId, clinical.treatmentGiven, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.treatmentGiven is null and ecrf.item='FLD.TRTAFTER.SYSTEMICST' and locked ='true';

# Bevat eind datum treatment?
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='end of at least 1 non-final treatment is missing' and formLocked = 'true';

# Bevat start datum treatment?
select patient.cpctId, clinical.treatmentStartDate, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.treatmentStartDate is null and ecrf.item='FLD.TRTAFTER.SYSSTDT' and locked='true';

# Bevat drug treatment?
select distinct patient.cpctId, clinical.treatment, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.treatment is null and ecrf.item='FLD.TRTAFTER.SYSREGPOST' and locked='true';

# Heeft treatment start dat na biopt datum?
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER;FRM.BIOPS' and message ='first treatment start is before first biopsy date' and formLocked = 'true';

# Treatment matched met 1 biopt?
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='multiple biopsy matches for treatment' and formLocked = 'true';

# heeft een treatment match met biopt?
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='no biopsy match for treatment' and formLocked = 'true';

# 2e treatment start is na eind datum 1e treatment
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='subsequent treatment starts before the end of previous treatment' and formLocked = 'true';

# bevat radiotherapyGiven?
select patient.cpctId, clinical.radiotherapyGiven, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.radiotherapyGiven is null and ecrf.item='FLD.TRTAFTER.RADIOTHERST' and locked ='true';

# heeft overlijdensdatum na einde treatment?
select * from clinicalFindings where ecrfitem ='FLD.DEATH.DDEATHDTC;FRM.TRTAFTER' and message ='death date before end of last treatment' and formLocked = 'true';

# Bevat response voor 1 treatment?
select * from clinicalFindings where ecrfitem ='FRM.TUMORMEASUREMENT' and message ='treatments are overlapping. Cannot match any response.' and formLocked = 'true';

# Bevat een ingevulde measurement?
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.TMYN' and message ='measurement done field empty' and formLocked = 'true';

# bevat een response of asessment datum in een correct format
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.ASSDTC;FLD.TUMORMEASUREMENT.RESPONSEDTC' and
 message ='response date and assessment date empty or in wrong format' and formLocked = 'true';

# Bevat een measurement en de uitkomst van response?
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.BESTRESPON' and message ='measurement done is yes, but response is empty (non-first response)'
and formLocked = 'true';

# Bevat treatment response voor treatment?
select * from clinicalFindings where ecrfitem ='FRM.TUMORMEASUREMENT' and message ='no treatment response for at least 1 treatment' and formLocked = 'true';

# Bevat treatment response data met een ingevulde treatment
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='treatment response filled in, but treatment data missing' and formLocked = 'true';

# Bevat eerste response meting voor treatment start datum?
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER;FRM.TUMORMEASUREMENT' and message ='first (baseline) measurement date is after first treatment start'
and formLocked = 'true';;

# Bevat baseline response meting voor start datum treatment
select * from clinicalFindings where ecrfitem ='FRM.TUMORMEASUREMENT' and message ='response after new baseline and before next treatment' and formLocked = 'true';

# bevat response of assesment datum als response uitkomst is ingevuld
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.ASSDTC;FLD.TUMORMEASUREMENT.RESPONSEDTC' and
message ='response filled in, but no assessment date and response date found' and formLocked = 'true';

# Response uitkomst veld is leeg als geen measurement is gedaan
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.TMYN' and message ='measurement done is no, but response filled in' and formLocked = 'true';

# Response of assessment datum is leeg als er geen measurement is gedaan
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.TMYN' and message ='measurement done is no, but assessment date or response date is filled in'
 and formLocked = 'true';