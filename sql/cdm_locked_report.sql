# heeft registratie datum voor biopt afname?
select * from clinicalFindings where ecrfitem ="FLD.ELIGIBILITY.REGDTC;FLD.SELCRIT.NREGDTC;FLD.BIOPS.BIOPTDT" and
 message ="at least 1 biopsy date prior to registration date" and formLocked = 'true';

# heeft tumor locatie?
select * from clinicalFindings where ecrfitem ="FLD.CARCINOMA.PTUMLOC;FLD.CARCINOMA.PTUMLOCS" and message ="primary tumor location empty" and formLocked = 'true';

# heeft geslacht?
select * from clinicalFindings where ecrfItem ="FLD.DEMOGRAPHY.SEX" and formLocked = 'true';

# heeft geboortejaar?
select * from clinicalFindings where ecrfItem ="FLD.SELCRIT.NBIRTHYEAR;FLD.ELIGIBILITY.BIRTHYEAR;FLD.ELIGIBILITY.BIRTHDTCES" and formLocked = 'true';

# geen ingevulde ecrf formulier en 1 biopt gesequenced
select * from clinicalFindings where ecrfItem = "FRM.BIOPS" and message ="less ecrf biopsies than biopsies sequenced." and formLocked = 'true' 
and details="ecrf biopsies: 0; sequenced: 1";

# 1 ingevulde ecrf formulier en 2 biopten gesequenced
select * from clinicalFindings where ecrfItem = "FRM.BIOPS" and message ="less ecrf biopsies than biopsies sequenced." and formLocked = 'true' 
and details="ecrf biopsies: 1; sequenced: 2";

# bevat 1 mogelijk biopt datum dat matched met sequenced biopt
select * from clinicalFindings where ecrfitem ="FRM.BIOPS" and message ="more than 1 possible clinical biopsy match for sequenced sample." and formLocked = 'true';

# heeft biopsySite?
select * from clinicalFindings where ecrfItem ="FLD.BIOPS.BILESSITE;FLD.BIOPS.BIOTHLESSITE" and formLocked = 'true';

# heeft biopsyLocation?
select * from clinicalFindings where ecrfItem ="FLD.BIOPS.BILESLOC" and formLocked = 'true';

# biopt datum is correct ingevuld
select * from clinicalFindings where ecrfitem ="FLD.BIOPS.BIOPTDT" and message ="biopsy date empty or in wrong format" and formLocked = 'true';

# heeft 'verdachte' biopt afname?
select clinical.sampleId, clinical.patientId,clinical.biopsyDate, sample.samplingDate, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from clinical
inner join sample on sample.sampleId = clinical.sampleId 
inner join ecrf on ecrf.patientId = clinical.patientId
where clinical.biopsyDate <> sample.samplingDate and ecrf.item ="FLD.BIOPS.BIOPTDT" and ecrf.Locked='true';

select * from clinicalFindings where ecrfItem = "FRM.BIOPS" and message ="could not match any clinical biopsy with sequenced sample." and formLocked = 'true'; 

# bevat startdatum van behandeling na afname biopsy
select * from clinicalFindings where ecrfitem ="FRM.TRTAFTER;FRM.BIOPS" and message ="first treatment start is before first biopsy date" and formLocked = 'true';

# bevat biopt dat matched met 1 behandeling
select * from clinicalFindings where ecrfitem ="FRM.TRTAFTER" and message ="multiple biopsy matches for treatment" and formLocked = 'true';

# bevat startdatum van behandeling en een biopt datum 
select * from clinicalFindings where ecrfitem ="FRM.TRTAFTER" and message ="treatment matched biopsy with null date." and formLocked = 'true';

# heeft biopt afname datum dat matched met behandeling
select * from clinicalFindings where ecrfitem ="FRM.TRTAFTER" and message ="no biopsy match for treatment" and formLocked = 'true';

# heeft behandeling?
select distinct clinical.patientId, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from sample
inner join treatment on sample.patientId = treatment.patientId 
inner join clinical on clinical.sampleId = sample.sampleId
inner join ecrf on ecrf.patientId = clinical.patientId
where treatment.treatmentGiven is null and ecrf.item = "FLD.TRTAFTER.SYSTEMICST" and ecrf.locked='true' and itemValue ='';

# bevat behandelde medicijn?
select distinct clinical.patientId, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from sample
inner join treatment on sample.patientId = treatment.patientId 
inner join clinical on clinical.sampleId = sample.sampleId
inner join ecrf on ecrf.patientId = clinical.patientId
where treatment.treatmentGiven ="Yes" and treatment.name is null and ecrf.item = "FLD.TRTAFTER.PLANNEDTRT" and ecrf.itemValue =" " and locked = 'true';

select distinct clinical.patientId, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from sample
inner join treatment on sample.patientId = treatment.patientId 
inner join clinical on clinical.sampleId = sample.sampleId
inner join ecrf on ecrf.patientId = clinical.patientId
where treatment.treatmentGiven ="Yes" and treatment.name is null and ecrf.item = "FLD.TRTAFTER.PLANNEDTRT" and ecrf.itemValue = "Other";

select distinct clinical.patientId, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from sample
inner join treatment on sample.patientId = treatment.patientId 
inner join clinical on clinical.sampleId = sample.sampleId
inner join ecrf on ecrf.patientId = clinical.patientId
where clinical.patientId in ('XXX') and ecrf.item like "FLD.TRTAFTER.SYSREGPOST" and ecrf.itemValue =''
 and locked = 'true';

# bevat start datum behandeling? 
select distinct clinical.patientId, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from sample
inner join treatment on sample.patientId = treatment.patientId 
inner join clinical on clinical.sampleId = sample.sampleId
inner join ecrf on ecrf.patientId = clinical.patientId
where treatment.treatmentGiven="Yes" and treatment.startDate is null and ecrf.item = "FLD.TRTAFTER.SYSSTDT" and ecrf.locked='true' and ecrf.itemValue ="";

# bevat eind datum behandeling?
select distinct clinical.patientId, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from sample
inner join treatment on sample.patientId = treatment.patientId 
inner join clinical on clinical.sampleId = sample.sampleId
inner join ecrf on ecrf.patientId = clinical.patientId
where treatment.treatmentGiven="Yes" and treatment.endDate is null and ecrf.item = "FLD.TRTAFTER.SYSENDT" and ecrf.locked='true' and itemValue=""
and treatment.startDate < '2017-05-29';

# bevat startdatum behandeling na eind datum vorige behandeling 
select * from clinicalFindings where ecrfitem ="FRM.TRTAFTER" and message ="subsequent treatment starts before the end of previous treatment" and formLocked = 'true';

# bevat repsonse datum behandeling, en data van behandeling is ingevuld 
select * from clinicalFindings where ecrfitem ="FRM.TRTAFTER" and message ="treatment response filled in, but treatment data missing" and formLocked = 'true';

# bevat een response van behandeling en matched met 1 behandeling
select * from clinicalFindings where ecrfitem ="FRM.TUMORMEASUREMENT" and message ="treatment response matches multiple treatments" and formLocked = 'true';

# heeft response datum?
select distinct clinical.patientId, ecrf.locked from treatmentResponse
inner join sample on sample.patientId = treatmentResponse.patientId 
inner join treatment on sample.patientId = treatment.patientId 
inner join clinical on clinical.sampleId = sample.sampleId
inner join ecrf on ecrf.patientId = clinical.patientId
where treatmentResponse.responseDate is null and ecrf.item in ("FLD.TUMORMEASUREMENT.BONEYN", "FLD.TUMORMEASUREMENT.RESPONSEDTC") and itemValue in ("Yes", "")
 and ecrf.locked = 'true';

# heeft meting van response?
select distinct clinical.patientId, ecrf.locked from treatmentResponse
inner join sample on sample.patientId = treatmentResponse.patientId 
inner join treatment on sample.patientId = treatment.patientId 
inner join clinical on clinical.sampleId = sample.sampleId
inner join ecrf on ecrf.patientId = clinical.patientId
where treatmentResponse.measurementDone is null and ecrf.item in("FLD.TUMORMEASUREMENT.BONEYN", "FLD.TUMORMEASUREMENT.TMYN") and ecrf.locked = 'true' 
and ecrf.itemValue in ("Yes", "");

# heeft response?
select distinct clinical.patientId, ecrf.locked from treatmentResponse
inner join sample on sample.patientId = treatmentResponse.patientId 
inner join treatment on sample.patientId = treatment.patientId 
inner join clinical on clinical.sampleId = sample.sampleId
inner join ecrf on ecrf.patientId = clinical.patientId
where treatmentResponse.response is null and ecrf.item in ("FLD.TUMORMEASUREMENT.BONEYN", "FLD.TUMORMEASUREMENT.BESTRESPON") and ecrf.locked = 'true' 
and ecrf.itemValue in ("Yes", "");

# bevat naam van het ziekenhuis?
select distinct clinical.patientId, hospital, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from clinical 
inner join ecrf on ecrf.patientId = clinical.patientId
where isnull(hospital) and ecrf.item = "FLD.SELCRIT.NHOSPITAL" and ecrf.locked='true';

# heeft overlijdensdatum na einde treatment
select * from clinicalFindings where ecrfitem ="FLD.DEATH.DDEATHDTC;FRM.TRTAFTER" and message ="death date before end of last treatment" and formLocked = 'true';


