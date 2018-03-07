# Locked CDM report

# heeft een correct informed constent datum?
select * from clinicalFindings where ecrfitem ='FLD.INFORMEDCONSENT.ICDTC' and message ='informed consent date empty or in wrong format' and formLocked='true';

# heeft tumor locatie?
select * from clinicalFindings where ecrfitem ='FLD.CARCINOMA.PTUMLOC;FLD.CARCINOMA.PTUMLOCS' and message ='primary tumor location empty' and formLocked = 'true';

# heeft geslacht?
select * from clinicalFindings where ecrfItem ='FLD.DEMOGRAPHY.SEX' and message = 'gender empty' and formLocked = 'true';

# heeft geboortejaar?
select * from clinicalFindings where ecrfItem ='FLD.SELCRIT.NBIRTHYEAR;FLD.ELIGIBILITY.BIRTHYEAR;FLD.ELIGIBILITY.BIRTHDTCES' and
message = 'birth year could not be determined' and formLocked = 'true';

# heeft een biopt formulier
select * from clinicalFindings where ecrfitem ='FRM.BIOPS' and message ='no biopsies found' and formLocked = 'true';

# geen ingevulde ecrf formulier en 1 biopt gesequenced
select * from clinicalFindings where ecrfItem = 'FRM.BIOPS' and message ='less ecrf biopsies than biopsies sequenced.' and formLocked = 'true'
and details='ecrf biopsies: 0; sequenced: 1';

# 1 ingevulde ecrf formulier en 2 biopten gesequenced
select * from clinicalFindings where ecrfItem = 'FRM.BIOPS' and message ='less ecrf biopsies than biopsies sequenced.' and formLocked = 'true'
and details='ecrf biopsies: 1; sequenced: 2';

# bevat 1 ecrf formulier voor 1 sequenced biopt
select * from clinicalFindings where ecrfitem ='FRM.BIOPS' and message ='more than 1 possible clinical biopsy match for sequenced sample.' and formLocked = 'true';

# heeft biopsySite?
select * from clinicalFindings where ecrfItem ='FLD.BIOPS.BILESSITE;FLD.BIOPS.BIOTHLESSITE' and message='biopsy site empty' and formLocked = 'true';

# heeft biopsyLocation?
select * from clinicalFindings where ecrfItem ='FLD.BIOPS.BILESLOC' and message='biopsy location empty' and formLocked = 'true';

# biopt datum dat correct is ingevuld
select * from clinicalFindings where ecrfitem ='FLD.BIOPS.BIOPTDT' and message ='biopsy date empty or in wrong format' and formLocked='true';

# heeft een ingevulde bioptdatum en het biopt formulier is ingevuld
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='treatment matched biopsy with null date.' and formLocked = 'true';

# heeft verdacht bioptdatum?
select * from clinicalFindings where ecrfItem = 'FLD.INFORMEDCONSENT.ICDTC;FLD.BIOPS.BIOPTDT' and
message ='at least 1 biopsy taken more than 30 days prior to informed consent date' and formLocked = 'true';

select * from clinicalFindings where ecrfitem ='FRM.BIOPS' and message ='could not match any clinical biopsy with sequenced sample.'
and details not like '%ecrf biopsies: [].%' and formLocked = 'true';

# sampling date en biopt date zijn gelijk
select * from clinicalFindings where ecrfitem ='FRM.BIOPS' and message ='sampling date does not equal biopsy date in matched biopsy' and formLocked = 'true';

# biopt matched met treatment
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='no biopsy match for treatment' and formLocked = 'true';

# treatment start is na biopt datum
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER;FRM.BIOPS' and message ='first treatment start is before first biopsy date' and formLocked = 'true';

# bevat start datum 2e treatment na eind datum vorige behandeling
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='subsequent treatment starts before the end of previous treatment' and formLocked = 'true';

# bevat 1 biopt voor treatment
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='multiple biopsy matches for treatment' and formLocked = 'true';

# bevat geen overlappende treatments
select * from clinicalFindings where ecrfitem ='FRM.TUMORMEASUREMENT' and message ='treatments are overlapping. Cannot match any response.' and formLocked = 'true';

# response gedaan veld is ingevuld
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.TMYN' and message ='measurement done field empty' and formLocked = 'true';

# de eerste response datum is voor treatment start datum
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER;FRM.TUMORMEASUREMENT' and message ='first (baseline) measurement date is after first treatment start'
and formLocked = 'true';

# bevat een eind datum van de behandeling
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='end of at least 1 non-final treatment is missing' and formLocked = 'true';

# bevat treatment response en heeft een treatment
select * from clinicalFindings where ecrfitem ='FRM.TUMORMEASUREMENT' and message ='no treatment response for at least 1 treatment' and formLocked = 'true';

# bevat treatment response en bevat een ingevulde treatment
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='treatment response filled in, but treatment data missing' and formLocked = 'true';

# bevat ja in response en response is ingevuld
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.BESTRESPON' and message ='measurement done is yes, but response is empty (non-first response)'
and formLocked = 'true';

# bevat een response voor nieuw treatment
select * from clinicalFindings where ecrfitem ='FRM.TUMORMEASUREMENT' and message ='response after new baseline and before next treatment' and formLocked = 'true';

# bevat response of assessment datum en heeft een response
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.ASSDTC;FLD.TUMORMEASUREMENT.RESPONSEDTC' and
message ='response filled in, but no assessment date and response date found' and formLocked = 'true';

# response of assessment datum is correct ingevuld
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.ASSDTC;FLD.TUMORMEASUREMENT.RESPONSEDTC' and
 message ='response date and assessment date empty or in wrong format' and formLocked = 'true';

# bevat nee in response en bevat een response of assessment datum
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.TMYN' and message ='measurement done is no, but assessment date or response date is filled in'
 and formLocked = 'true';

# bevat nee in response en bevat response datum
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.TMYN' and message ='measurement done is no, but response filled in' and formLocked = 'true';

# heeft overlijdensdatum na einde treatment
select * from clinicalFindings where ecrfitem ='FLD.DEATH.DDEATHDTC;FRM.TRTAFTER' and message ='death date before end of last treatment' and formLocked = 'true';

# bevat naam van het ziekenhuis?
select distinct clinical.patientId, hospital, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from clinical
inner join ecrf on ecrf.patientId = clinical.patientId
where isnull(hospital) and ecrf.item = 'FLD.SELCRIT.NHOSPITAL' and ecrf.locked='true';
