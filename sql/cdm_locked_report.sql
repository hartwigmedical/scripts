#  LISC: Number of CPCT patients
select patientId, count(sampleId) as countSamples from clinical where sampleId like '%CPCT%' group by patientId ;

#  LISC: Registration date of patient
select distinct patientId, registrationDate, informedConsentDate from clinical where patientId like '%CPCT%' group by sampleId;

#  LISC: Timeline starts with informed consent?
select * from clinicalFindings where
message ='At least 1 biopsy taken before informed consent date' and patientId like '%CPCT%' and formLocked = 'true';

select * from clinicalFindings where message ='informed consent date empty or in wrong format' and patientId like '%CPCT%'  and formLocked = 'true';

#  LISC: Hospital is known?
select * from clinicalFindings where message like '%Hospital could not be determined%' and patientId like '%CPCT%'  and formLocked = 'true';

#  LISC: Gender is known?
select * from clinicalFindings where message = 'Gender empty' and patientId like '%CPCT%'  and formLocked = 'true';

#  LISC: Birthyear is known?
select * from clinicalFindings where
message = 'birth year could not be determined' and patientId like '%CPCT%'  and formLocked = 'true';

#  LISC: Can curate the cancer type?
select * from clinicalFindings where message ='primary tumor location empty' and patientId like '%CPCT%'  and formLocked = 'true';

#  LISC: Systemic pre-therapy is known?
select * from clinicalFindings where message = 'pre systemic treatment given empty' and patientId like '%CPCT%'  and formLocked = 'true';

#  LISC: Radiotherapy pre-therapy is known?
select * from clinicalFindings where message = 'pre radio treatment given empty' and patientId like '%CPCT%'  and formLocked = 'true';