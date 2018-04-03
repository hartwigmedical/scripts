# Number of CPCT patients
select patientId, count(sampleId) as countSamples from clinical where sampleId like '%CPCT%'group by patientId ;

# Registration date of patient
select patientIdentifier, registrationDate from baseline
inner join patient on patient.id=baseline.patientId order by 1;

# Timeline starts with informed consent?
select * from clinicalFindings where ecrfItem = 'FLD.INFORMEDCONSENT.ICDTC;FLD.BIOPS.BIOPTDT' and
message ='at least 1 biopsy taken before informed consent date' and formLocked='true';

select * from clinicalFindings where ecrfitem ='FLD.INFORMEDCONSENT.ICDTC' and message ='informed consent date empty or in wrong format' and formLocked='true';

# Hospital is known?
select * from clinicalFindings where message like '%no hospital%' and formLocked='true';

# Gender is known?
select * from clinicalFindings where ecrfItem ='FLD.DEMOGRAPHY.SEX' and message = 'gender empty' and formLocked='true';

# Birthyear is known?
select * from clinicalFindings where ecrfItem ='FLD.SELCRIT.NBIRTHYEAR;FLD.ELIGIBILITY.BIRTHYEAR;FLD.ELIGIBILITY.BIRTHDTCES' and
message = 'birth year could not be determined' and formLocked='true';

# Can curate the cancer type?
select * from clinicalFindings where ecrfitem ='FLD.CARCINOMA.PTUMLOC;FLD.CARCINOMA.PTUMLOCS' and message ='primary tumor location empty' and formLocked='true';
select * from clinicalFindings where ecrfitem ='FLD.CARCINOMA.PTUMLOC;FLD.CARCINOMA.PTUMLOCS' and message ='failed to curate primary tumor location.' and formLocked='true';

# Systemic pre-therapy is known?
select * from clinicalFindings where ecrfitem ='FLD.PRETHERAPY.SYSTEMIC' and message = 'pre systemic treatment given empty' and formLocked='true';

# Radiotherapy pre-therapy is known?
select * from clinicalFindings where ecrfitem ='FLD.PRETHERAPY.RADIOTHER' and message = 'pre radio treatment given empty' and formLocked='true';

# Can curate all pre-therapies?
select * from clinicalFindings where ecrfitem ='FLD.PRETHERAPY.SYSTEMICREG'
and message = 'failed to curate ecrf drug. Curated list contained no matching entry, or match was ambiguous.' and formLocked='true';
