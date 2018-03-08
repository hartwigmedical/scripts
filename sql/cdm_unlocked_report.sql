# bevat naam van het ziekenhuis?
select distinct clinical.patientId, hospital, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from clinical
inner join ecrf on ecrf.patientId = clinical.patientId
where isnull(hospital) and ecrf.item = 'FLD.SELCRIT.NHOSPITAL' and ecrf.locked='false';

# heeft geslacht?
select * from clinicalFindings where ecrfItem ='FLD.DEMOGRAPHY.SEX' and message = 'gender empty' and formLocked = 'false';

# heeft geboortejaar?
select * from clinicalFindings where ecrfItem ='FLD.SELCRIT.NBIRTHYEAR;FLD.ELIGIBILITY.BIRTHYEAR;FLD.ELIGIBILITY.BIRTHDTCES' and
message = 'birth year could not be determined' and formLocked = 'false';

# Heeft een systemicPreTreatment?
select patient.cpctId, clinical.hasSystemicPretreatment, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.hasSystemicPretreatment is null and locked ='false';

# Heeft RadiotherapyPreTreatment?
select patient.cpctId, clinical.hasRadiotherapyPreTreatment, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.hasRadiotherapyPreTreatment is null and locked ='false';

# Bevat preTreatments?
select patient.cpctId, clinical.preTreatments, clinical.hasSystemicPretreatment, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.preTreatments is null and clinical.hasSystemicPretreatment = 'yes' and locked ='false';

# heeft een correct informed constent datum?
select * from clinicalFindings where ecrfitem ='FLD.INFORMEDCONSENT.ICDTC' and message ='informed consent date empty or in wrong format' and formLocked='false';

# heeft een matching informed constent?
select * from clinicalFindings where ecrfItem = 'FLD.INFORMEDCONSENT.ICDTC;FLD.BIOPS.BIOPTDT' and
message ='at least 1 biopsy taken more than 30 days prior to informed consent date' and formLocked = 'false';

# heeft een biopt formulier
select * from clinicalFindings where ecrfitem ='FRM.BIOPS' and message ='no biopsies found' and formLocked = 'false';

# biopt datum is correct ingevuld
select * from clinicalFindings where ecrfitem ='FLD.BIOPS.BIOPTDT' and message ='biopsy date empty or in wrong format' and formLocked='false';

# Bevat geen meerdere biopten voor sequenced biopt
select * from clinicalFindings where ecrfitem ='FRM.BIOPS' and message ='more than 1 possible clinical biopsy match for sequenced sample.' and formLocked = 'false';

# Bevat geen treatment dat matched met null datum biopt
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='treatment matched biopsy with null date.' and formLocked = 'false';

# geen ingevulde ecrf formulier en 1 biopt gesequenced
select * from clinicalFindings where ecrfItem = 'FRM.BIOPS' and message ='less ecrf biopsies than biopsies sequenced.' and formLocked = 'false'
and details='ecrf biopsies: 0; sequenced: 1';

# 1 ingevulde ecrf formulier en 2 biopten gesequenced
select * from clinicalFindings where ecrfItem = 'FRM.BIOPS' and message ='less ecrf biopsies than biopsies sequenced.' and formLocked = 'false'
and details='ecrf biopsies: 1; sequenced: 2';

# heeft verdacht bioptdatum?
select * from clinicalFindings where ecrfitem ='FRM.BIOPS' and message ='could not match any clinical biopsy with sequenced sample.'
and details not like '%ecrf biopsies: [].%' and formLocked = 'false';

# sampling date en biopt date zijn gelijk
select * from clinicalFindings where ecrfitem ='FRM.BIOPS' and message ='sampling date does not equal biopsy date in matched biopsy' and formLocked = 'false';

# biopt matched met treatment
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='no biopsy match for treatment' and formLocked = 'false';

# Bevat geen meerdere biopten dat match met treatment
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='multiple biopsy matches for treatment' and formLocked = 'false';

# treatment start is na biopt datum
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER;FRM.BIOPS' and message ='first treatment start is before first biopsy date' and formLocked = 'false';

# bevat geen overlappende treatments
select * from clinicalFindings where ecrfitem ='FRM.TUMORMEASUREMENT' and message ='treatments are overlapping. Cannot match any response.' and formLocked = 'false';

# Startdatum 2e treatment is na eind datum 1e treatment
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='subsequent treatment starts before the end of previous treatment' and formLocked = 'false';

# Heeft een ingevulde bioptdatum op biopt formulier for matching treatment
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='end of at least 1 non-final treatment is missing' and formLocked = 'false';

# heeft overlijdensdatum na einde treatment
select * from clinicalFindings where ecrfitem ='FLD.DEATH.DDEATHDTC;FRM.TRTAFTER' and message ='death date before end of last treatment' and formLocked = 'false';

# heeft tumor locatie?
select * from clinicalFindings where ecrfitem ='FLD.CARCINOMA.PTUMLOC;FLD.CARCINOMA.PTUMLOCS' and message ='primary tumor location empty' and formLocked = 'false';

# heeft biopsySite?
select * from clinicalFindings where ecrfItem ='FLD.BIOPS.BILESSITE;FLD.BIOPS.BIOTHLESSITE' and message='biopsy site empty' and formLocked = 'false';

# heeft biopsyLocation?
select * from clinicalFindings where ecrfItem ='FLD.BIOPS.BILESLOC' and message='biopsy location empty' and formLocked = 'false';

# bevat treatmentGiven
select patient.cpctId, clinical.treatmentGiven, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.treatmentGiven is null and locked ='false';

# bevat radiotherapyGiven
select patient.cpctId, clinical.radiotherapyGiven, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.radiotherapyGiven is null and locked ='false';

# Bevat een ingevulde measurement?
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.TMYN' and message ='measurement done field empty' and formLocked = 'false';

# Bevat treatment response en een ingevulde treatment
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='treatment response filled in, but treatment data missing' and formLocked = 'false';

# Bevat een treatment response en tenminste 1 treatment
select * from clinicalFindings where ecrfitem ='FRM.TUMORMEASUREMENT' and message ='no treatment response for at least 1 treatment' and formLocked = 'false';

# Measurement done is ja, en response is ingevuld
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.BESTRESPON' and message ='measurement done is yes, but response is empty (non-first response)'
and formLocked = 'false';

# Mesurement done is no, en bevat geen assesment of response datum
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.TMYN' and message ='measurement done is no, but assessment date or response date is filled in'
 and formLocked = 'false';

# Measurement done is yes, en bevat response
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.TMYN' and message ='measurement done is no, but response filled in' and formLocked = 'false';

# Bevat een response en een assment of response datum met measurement is done
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.ASSDTC;FLD.TUMORMEASUREMENT.RESPONSEDTC' and
message ='response filled in, but no assessment date and response date found' and formLocked = 'false';

# Response of assessment datum is correct ingevuld
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.ASSDTC;FLD.TUMORMEASUREMENT.RESPONSEDTC' and
 message ='response date and assessment date empty or in wrong format' and formLocked = 'false';

# De eerste response datum is voor treatment start datum
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER;FRM.TUMORMEASUREMENT' and message ='first (baseline) measurement date is after first treatment start'
and formLocked = 'false';

# Bevat response na nieuw baseline response en voor nieuw treatment
select * from clinicalFindings where ecrfitem ='FRM.TUMORMEASUREMENT' and message ='response after new baseline and before next treatment' and formLocked = 'false';