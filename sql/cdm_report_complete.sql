# Number of CPCT patients
select patientId, count(sampleId) as countSamples from clinical where sampleId like '%CPCT%'group by patientId ;

# Registration date of patient
select patientIdentifier, registrationDate from baseline
inner join patient on patient.id=baseline.patientId order by 1;

# Timeline starts with informed consent?
select * from clinicalFindings where ecrfItem = 'FLD.INFORMEDCONSENT.ICDTC;FLD.BIOPS.BIOPTDT' and
message ='at least 1 biopsy taken before informed consent date';

select * from clinicalFindings where ecrfitem ='FLD.INFORMEDCONSENT.ICDTC' and message ='informed consent date empty or in wrong format';

# Hospital is known?
select * from clinicalFindings where message like '%no hospital%';

# Gender is known?
select * from clinicalFindings where ecrfItem ='FLD.DEMOGRAPHY.SEX' and message = 'gender empty';

# Birthyear is known?
select * from clinicalFindings where ecrfItem ='FLD.SELCRIT.NBIRTHYEAR;FLD.ELIGIBILITY.BIRTHYEAR;FLD.ELIGIBILITY.BIRTHDTCES' and
message = 'birth year could not be determined';

# Can curate the cancer type?
select * from clinicalFindings where ecrfitem ='FLD.CARCINOMA.PTUMLOC;FLD.CARCINOMA.PTUMLOCS' and message ='primary tumor location empty';
select * from clinicalFindings where ecrfitem ='FLD.CARCINOMA.PTUMLOC;FLD.CARCINOMA.PTUMLOCS' and message ='failed to curate primary tumor location.';

# Systemic pre-therapy is known?
select * from clinicalFindings where ecrfitem ='FLD.PRETHERAPY.SYSTEMIC' and message = 'pre systemic treatment given empty';

# Radiotherapy pre-therapy is known?
select * from clinicalFindings where ecrfitem ='FLD.PRETHERAPY.RADIOTHER' and message = 'pre radio treatment given empty';

# Can curate all pre-therapies?
select * from clinicalFindings where ecrfitem ='FLD.PRETHERAPY.SYSTEMICREG' and message = 'failed to curate ecrf drug. Curated list contained no matching entry, or match was ambiguous.';

# Has enough ECRF biopsy forms?
select * from clinicalFindings where ecrfItem ='FRM.BIOPS' and message = 'no biopsies found';

select * from clinicalFindings where ecrfItem = 'FRM.BIOPS' and message ='less ecrf biopsies than biopsies sequenced.'
and details='ecrf biopsies: 0; sequenced: 1';

select * from clinicalFindings where ecrfItem = 'FRM.BIOPS' and message ='less ecrf biopsies than biopsies sequenced.'
and details like '%ecrf biopsies: 1%';

# Can match every sequenced sample with an ECRF biopsy?
select * from clinicalFindings where ecrfitem ='FLD.BIOPS.BIOPTDT' and message ='biopsy date empty or in wrong format';

select * from clinicalFindings where ecrfitem ='FRM.BIOPS' and message ='could not match any clinical biopsy with sequenced sample.'
and details not like '%ecrf biopsies: [].%';

select * from clinicalFindings where ecrfitem ='FRM.BIOPS' and message ='more than 1 possible clinical biopsy match for sequenced sample.';

# ECRF biopsy date matches with HMF sampling date?
select * from clinicalFindings where ecrfitem ='FRM.BIOPS' and message ='sampling date does not equal biopsy date in matched biopsy';

# All biopsy sites are known?
select * from clinicalFindings where ecrfItem ='FLD.BIOPS.BILESSITE;FLD.BIOPS.BIOTHLESSITE' and message='biopsy site empty';

# All biopsy locations are known?
select * from clinicalFindings where ecrfItem ='FLD.BIOPS.BILESLOC' and message='biopsy location empty';

# Has enough ECRF treatment forms?
select * from clinicalFindings where ecrfItem ='FLD.TRTAFTER.SYSTEMICST' and message='treatment given field empty';

select * from clinicalFindings where message like '%treatment given is not yes/no%';

select * from clinicalFindings where ecrfItem ='FLD.TRTAFTER.SYSSTDT' and message='drug start date empty or in wrong format';

select * from clinicalFindings where ecrfItem ='FLD.TRTAFTER.PLANNEDTRT;FLD.TRTAFTER.SYSREGPOST' and message = 'drug name empty';

select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='end of at least 1 non-final treatment is missing';

select * from clinicalFindings where ecrfItem ='FLD.TRTAFTER.SYSTEMICST' and message='treatment given is no, but treatment data is filled in';

select * from clinicalFindings where ecrfItem ='FRM.TRTAFTER' and message = 'treatment given is yes, but no treatment data filled in';

# Can match every sequenced biopsy to a treatment?
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='no biopsy match for treatment';

select * from clinicalFindings where ecrfItem ='FLD.TRTAFTER.SYSSTDT;FLD.TRTAFTER.SYSENDT' and message='drug startDate is after drug endDate';

select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='multiple biopsy matches for treatment';

select * from clinicalFindings where message like '%treatment matched biopsy with null date%';

# Can curate all treatments?
select * from clinicalFindings where ecrfitem ='FLD.TRTAFTER.PLANNEDTRT;FLD.TRTAFTER.SYSREGPOST' and message = 'failed to curate ecrf drug. Curated list contained no matching entry, or match was ambiguous.';

# Are all treatments sequentially in time?
select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER;FRM.BIOPS' and message ='first treatment start is before first biopsy date';

select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='subsequent treatment starts before the end of previous treatment';

# Do we know whether radiotherapy was given for every treatment?
select * from clinicalFindings where ecrfitem ='FLD.TRTAFTER.RADIOTHERST' and message ='radio therapy given field empty';

# Does the final treatment end before death date?
select * from clinicalFindings where ecrfitem ='FLD.DEATH.DDEATHDTC;FRM.TRTAFTER' and message ='death date before end of last treatment';

# Every treatment has at least 1 valid response?
select * from cinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.TMYN' and message ='measurement done field empty';

select * from clinicalFindings where message like '%measurement done is not yes/no%';

select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.ASSDTC;FLD.TUMORMEASUREMENT.RESPONSEDTC' and
 message ='response date and assessment date empty or in wrong format';

select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.BESTRESPON' and message ='measurement done is yes, but response is empty (non-first response)';

select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.ASSDTC;FLD.TUMORMEASUREMENT.RESPONSEDTC' and
message ='response filled in, but no assessment date and response date found';

# Every response is valid and can be matched to a treatment?
select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.TMYN' and message ='measurement done is no, but response filled in';

select * from clinicalFindings where ecrfitem ='FLD.TUMORMEASUREMENT.TMYN' and message ='measurement done is no, but assessment date or response date is filled in';

select * from clinicalFindings where ecrfitem ='FRM.TUMORMEASUREMENT' and message ='treatments are overlapping. Cannot match any response.';

select * from clinicalFindings where ecrfitem ='FRM.TUMORMEASUREMENT' and message ='response after new baseline and before next treatment';

select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER;FRM.TUMORMEASUREMENT' and message ='first (baseline) measurement date is after first treatment start';

select * from clinicalFindings where ecrfitem ='FRM.TRTAFTER' and message ='treatment response filled in, but treatment data missing';

select * from clinicalFindings where ecrfitem ='FRM.TUMORMEASUREMENT' and message ='no treatment response for at least 1 treatment';

