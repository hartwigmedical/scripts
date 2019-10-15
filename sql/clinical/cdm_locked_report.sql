#  Baseline information
select patientId, count(sampleId), registrationDate, informedConsentDate from clinical where patientId like '%CPCT%' and not(isnull(registrationDate)) group by 1,3,4;

#  Timeline starts with informed consent?
select patientId, message, details from clinicalFindings where message ='At least 1 biopsy taken before informed consent date' and patientId like '%CPCT%' and formLocked = 'true'

union select patientId, message, details from clinicalFindings where message ='informed consent date empty or in wrong format' and patientId like '%CPCT%'  and formLocked = 'true';

#  Hospital is known?
select patientId, message, details from clinicalFindings where message like '%Hospital could not be determined%' and patientId like '%CPCT%'  and formLocked = 'true';

#  Gender is known?
select patientId, message, details from clinicalFindings where message = 'Gender empty' and patientId like '%CPCT%'  and formLocked = 'true';

#  Birthyear is known?
select patientId, message, details from clinicalFindings where message = 'birth year could not be determined' and patientId like '%CPCT%'  and formLocked = 'true';

#  Can curate the cancer type?
select patientId, message, details from clinicalFindings where message ='primary tumor location empty' and patientId like '%CPCT%'  and formLocked = 'true';

#  Systemic pre-therapy is known?
select patientId, message, details from clinicalFindings where message = 'pre systemic treatment given empty' and patientId like '%CPCT%'  and formLocked = 'true';


#  Radiotherapy pre-therapy is known?
select patientId, message, details from clinicalFindings where message = 'pre radio treatment given empty' and patientId like '%CPCT%'  and formLocked = 'true';