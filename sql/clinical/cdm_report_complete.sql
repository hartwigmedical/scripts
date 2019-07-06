#  Baseline information
select patientId, count(sampleId), registrationDate, informedConsentDate from clinical where patientId like '%CPCT%' and not(isnull(registrationDate)) group by 1,3,4;

#  Timeline starts with informed consent?
select patientId, message, details from clinicalFindings where message ='At least 1 biopsy taken before informed consent date' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message ='informed consent date empty or in wrong format' and patientId like '%CPCT%';

#  Hospital is known?
select patientId, message, details from clinicalFindings where message like '%Hospital could not be determined%' and patientId like '%CPCT%';

#  Gender is known?
select patientId, message, details from clinicalFindings where message = 'Gender empty' and patientId like '%CPCT%';

#  Birthyear is known?
select patientId, message, details from clinicalFindings where message = 'birth year could not be determined' and patientId like '%CPCT%';

#  Can curate the cancer type?
select patientId, message, details from clinicalFindings where message ='primary tumor location empty' and patientId like '%CPCT%';

#  Systemic pre-therapy is known?
select patientId, message, details from clinicalFindings where message = 'pre systemic treatment given empty' and patientId like '%CPCT%';

#  Radiotherapy pre-therapy is known?
select patientId, message, details from clinicalFindings where message = 'pre radio treatment given empty' and patientId like '%CPCT%';

#  Has enough ECRF biopsy forms?
select patientId, message, details from clinicalFindings where message = 'Not enough clinical biopsy forms to match for every sequenced sample' and patientId like '%CPCT%';

#  Can match every sequenced sample with an ECRF biopsy?
select patientId, message, details from clinicalFindings where message ='biopsy date empty or in wrong format' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message ='could not match any clinical biopsy with sequenced sample.' and details not like '%ecrf biopsies: [].%' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message ='more than 1 possible clinical biopsy match for sequenced sample.' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message = 'Undetermined match issue in biopsy matcher' and patientId like '%CPCT%';

#  ECRF biopsy date matches with HMF sampling date?
select patientId, message, details from clinicalFindings where message ='sampling date does not equal biopsy date in matched biopsy' and patientId like '%CPCT%';

#  All biopsy sites and biopsy location are empty?
select patientId, message, details from clinicalFindings where message='Biopsy site and biopsy location are empty' and patientId like '%CPCT%';

#  Can find treatment for every sequenced & matched biopsy?
select patientId, message, details from clinicalFindings where message='Could not match treatment to matched biopsy' and patientId like '%CPCT%';

#  Are all matched treatments correctly filled in?
select patientId, message, details from clinicalFindings where message='Treatment given is yes, but no drugs are filled in' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message='treatment given field empty' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message like '%treatment given is not yes/no%' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message='drug start date empty or in wrong format' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message = 'drug name empty' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message ='end of at least 1 non-final treatment is missing' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message='Treatment given is no, but drugs are filled in' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message='Drug start date is after drug end date' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message ='radio therapy given field empty' and patientId like '%CPCT%';

#  Are all treatments sequentially in time?
select patientId, message, details from clinicalFindings where message ='subsequent treatment starts before the end of previous treatment' and patientId like '%CPCT%';

#  Does the final treatment end before death date?
select patientId, message, details from clinicalFindings where message ='death date before end of last treatment' and patientId like '%CPCT%';

#  Does every treatment >16 weeks have at least 1 treatment response form?
select patientId, message, details from clinicalFindings where message ='No treatment response for at least 1 matched treatment' and patientId like '%CPCT%';

#  Every response is valid and can be matched to a treatment?
select patientId, message, details from clinicalFindings where message ='measurement done field empty' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message ='Response date and/or assessment date empty or in wrong format' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message ='measurement done is yes, but response is empty (non-first response)' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message ='response filled in, but no assessment date and response date found' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message like '%measurement done is not yes/no%'

union select patientId, message, details from clinicalFindings where message ='measurement done is no, but response filled in' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message ='measurement done is no, but assessment date or response date is filled in' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message ='treatments are overlapping. Cannot match any response.' and patientId like '%CPCT%'

union select patientId, message, details from clinicalFindings where message ='response after new baseline and before next treatment' and patientId like '%CPCT%';
