# count samples
select patientId, group_concat(sampleId separator ', ') as sampleNames, count(sampleId) as countSamples from clinical group by patientId;

# registratie datum
select distinct patientId, registrationDate from clinical
inner join patient on patient.cpctId=clinical.patientId order by 1;

# Informed consent is eerste moment in de tijd
select * from clinicalFindings where ecrfItem = 'FLD.INFORMEDCONSENT.ICDTC;FLD.BIOPS.BIOPTDT' and
message ='at least 1 biopsy taken more than 30 days prior to informed consent date';

select * from clinicalFindings where ecrfitem ='FLD.INFORMEDCONSENT.ICDTC' and message ='informed consent date empty or in wrong format';

# Ziekenhuis is bekend?
select distinct clinical.patientId, hospital, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from clinical
inner join ecrf on ecrf.patientId = clinical.patientId
where isnull(hospital) and ecrf.item = 'FLD.SELCRIT.NHOSPITAL';

# Geslacht is bekend?
select * from clinicalFindings where ecrfItem ='FLD.DEMOGRAPHY.SEX' and message = 'gender empty';

# Geboortejaar is bekend?
select * from clinicalFindings where ecrfItem ='FLD.SELCRIT.NBIRTHYEAR;FLD.ELIGIBILITY.BIRTHYEAR;FLD.ELIGIBILITY.BIRTHDTCES' and
message = 'birth year could not be determined';

# Kan de tumor locatie cureren?
select * from clinicalFindings where ecrfitem ='FLD.CARCINOMA.PTUMLOC;FLD.CARCINOMA.PTUMLOCS' and message ='primary tumor location empty';

# Systemic pre-therapy is bekend?
select patient.cpctId, clinical.hasSystemicPretreatment, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.hasSystemicPretreatment is null and ecrf.item='FLD.PRETHERAPY.SYSTEMIC';

# Radiotherapy pre-therapy is bekend?
select patient.cpctId, clinical.hasRadiotherapyPreTreatment, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.hasRadiotherapyPreTreatment is null and ecrf.item ='FLD.PRETHERAPY.RADIOTHER';

# Kan alle pre therapies cureren?
select patient.cpctId, clinical.preTreatments, clinical.hasSystemicPretreatment, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.preTreatments is null and clinical.hasSystemicPretreatment = 'yes' and ecrf.item='FLD.PRETHERAPY.SYSTEMICREG';

# Heeft genoeg biopt formulieren?
select * from clinicalFindings where ecrfItem = 'FRM.BIOPS' and message ='less ecrf biopsies than biopsies sequenced.'
and details='ecrf biopsies: 0; sequenced: 1';

select * from clinicalFindings where ecrfItem = 'FRM.BIOPS' and message ='less ecrf biopsies than biopsies sequenced.'
and details like '%ecrf biopsies: 1%';

# Kan elk sequenced biopt matchen met een biopt formulier?
select * from clinicalFindings where ecrfitem ='FLD.BIOPS.BIOPTDT' and message ='biopsy date empty or in wrong format';

select * from clinicalFindings where ecrfitem ='FRM.BIOPS' and message ='could not match any clinical biopsy with sequenced sample.'
and details not like '%ecrf biopsies: [].%';

select * from clinicalFindings where ecrfitem ='FRM.BIOPS' and message ='more than 1 possible clinical biopsy match for sequenced sample.';

# Biopt datum matched met  HMF sampling datum?
select * from clinicalFindings where ecrfitem ='FRM.BIOPS' and message ='sampling date does not equal biopsy date in matched biopsy';

# Alle biopt sites zijn bekend?
select * from clinicalFindings where ecrfItem ='FLD.BIOPS.BILESSITE;FLD.BIOPS.BIOTHLESSITE' and message='biopsy site empty';

# Alle biopt locaties zijn bekend?
select * from clinicalFindings where ecrfItem ='FLD.BIOPS.BILESLOC' and message='biopsy location empty';

# Heeft genoeg treatment formulieren?
select patient.cpctId, clinical.treatmentGiven, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.treatmentGiven is null and ecrf.item='FLD.TRTAFTER.SYSTEMICST';

select patient.cpctId, clinical.treatmentStartDate, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.treatmentStartDate is null and ecrf.item='FLD.TRTAFTER.SYSSTDT';

select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='end of at least 1 non-final treatment is missing';

# Kan elk sequenced biopt matchen met een treatment?
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='no biopsy match for treatment';

select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='multiple biopsy matches for treatment';

# Kan alle treatments cureren?
select distinct patient.cpctId, clinical.treatment, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.treatment is null and ecrf.item='FLD.TRTAFTER.SYSREGPOST';

# Alle treatments zijn opvolgend in tijd?
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER;FRM.BIOPS' and message ='first treatment start is before first biopsy date';

select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='subsequent treatment starts before the end of previous treatment';

# Voor alle treatments is bekend of radiotherapy is gegeven?
select patient.cpctId, clinical.radiotherapyGiven, status, locked from clinical
inner join sample on clinical.sampleId=sample.sampleId
inner join patient on sample.patientId=patient.id
inner join ecrf on ecrf.patientId=patient.cpctId
where clinical.radiotherapyGiven is null and ecrf.item='FLD.TRTAFTER.RADIOTHERST';

# Einde laatste treatment is voor overlijdensdsatum?
select * from clinicalFindings where ecrfitem ='FLD.DEATH.DDEATHDTC;FRM.TRTAFTER' and message ='death date before end of last treatment';

# Ellk treatment heeft minimaal 1 geldig response formulier
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.TMYN' and message ='measurement done field empty';

select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.ASSDTC;FLD.TUMORMEASUREMENT.RESPONSEDTC' and
 message ='response date and assessment date empty or in wrong format';

select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.BESTRESPON' and message ='measurement done is yes, but response is empty (non-first response)';

select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.ASSDTC;FLD.TUMORMEASUREMENT.RESPONSEDTC' and
message ='response filled in, but no assessment date and response date found';

# Elk respons is geldig en kan gematched worden aan een treatment
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.TMYN' and message ='measurement done is no, but response filled in';

select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.TMYN' and message ='measurement done is no, but assessment date or response date is filled in';

select * from clinicalFindings where ecrfitem ='FRM.TUMORMEASUREMENT' and message ='treatments are overlapping. Cannot match any response.';

select * from clinicalFindings where ecrfitem ='FRM.TUMORMEASUREMENT' and message ='response after new baseline and before next treatment';

select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER;FRM.TUMORMEASUREMENT' and message ='first (baseline) measurement date is after first treatment start';

select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='treatment response filled in, but treatment data missing';

select * from clinicalFindings where ecrfitem ='FRM.TUMORMEASUREMENT' and message ='no treatment response for at least 1 treatment';

