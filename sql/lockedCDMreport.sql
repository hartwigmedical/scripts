# heeft tumor locatie?
select * from clinicalFindings where ecrfitem ="FLD.CARCINOMA.PTUMLOC;FLD.CARCINOMA.PTUMLOCS" and message ="primary tumor location empty" and formLocked = 'true';

# heeft geslacht?
select * from clinicalFindings where ecrfItem ="FLD.DEMOGRAPHY.SEX" and formLocked = 'true';

# heeft geboortejaar?
select * from clinicalFindings where ecrfItem ="FLD.SELCRIT.NBIRTHYEAR;FLD.ELIGIBILITY.BIRTHYEAR;FLD.ELIGIBILITY.BIRTHDTCES" and formLocked = 'true';


# 1 biopt gesequenced, minimaal 1 formulier en 2+ biopten gesequenced, minimaal 2+ biopt formulieren
select * from clinicalFindings where ecrfItem = "FRM.BIOPS" and message ="less ecrf biopsies than biopsies sequenced." and formLocked = 'true';

# heeft biopsySite?
select * from clinicalFindings where ecrfItem ="FLD.BIOPS.BILESSITE;FLD.BIOPS.BIOTHLESSITE" and formLocked = 'true';

# heeft biopsyLocation?
select * from clinicalFindings where ecrfItem ="FLD.BIOPS.BILESLOC" and formLocked = 'true';

# heeft betrouwbare datum biopt afname?
select clinical.sampleId, clinical.patientId,clinical.biopsyDate, sample.samplingDate, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from clinical
inner join sample on sample.sampleId = clinical.sampleId 
inner join ecrf on ecrf.patientId = clinical.patientId
where clinical.biopsyDate <> sample.samplingDate and ecrf.item ="FLD.BIOPS.BIOPTDT" and ecrf.Locked='true';

select * from clinicalFindings where ecrfItem = "FRM.BIOPS" and message ="could not match any clinical biopsy with sequenced sample." and formLocked = 'true'; # lijstje Immy

# , sample.sampleId, treatment.treatmentGiven, treatment.startDate, treatment.endDate, treatment.name, treatment.type

# heeft behandeling
select distinct clinical.patientId, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from sample
inner join treatment on sample.patientId = treatment.patientId 
inner join clinical on clinical.sampleId = sample.sampleId
inner join ecrf on ecrf.patientId = clinical.patientId
where treatment.treatmentGiven is null and ecrf.item = "FLD.TRTAFTER.SYSTEMICST" and ecrf.locked='true' and itemValue ='';


# heeft behandelde medicijn 
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

# heeft start datum 
select distinct clinical.patientId, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from sample
inner join treatment on sample.patientId = treatment.patientId 
inner join clinical on clinical.sampleId = sample.sampleId
inner join ecrf on ecrf.patientId = clinical.patientId
where treatment.treatmentGiven="Yes" and treatment.startDate is null and ecrf.item = "FLD.TRTAFTER.SYSSTDT" and ecrf.locked='true' and ecrf.itemValue ="";


# heeft eind datum 
select distinct clinical.patientId, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from sample
inner join treatment on sample.patientId = treatment.patientId 
inner join clinical on clinical.sampleId = sample.sampleId
inner join ecrf on ecrf.patientId = clinical.patientId
where treatment.treatmentGiven="Yes" and treatment.endDate is null and ecrf.item = "FLD.TRTAFTER.SYSENDT" and ecrf.locked='true' and itemValue=""
and treatment.startDate < '2017-05-29';

 #heeft  response date datum?
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

# heeft ziekenhuis
select distinct clinical.patientId, hospital, ecrf.status, ecrf.locked, ecrf.item, ecrf.itemValue from clinical 
inner join ecrf on ecrf.patientId = clinical.patientId
where isnull(hospital) and ecrf.item = "FLD.SELCRIT.NHOSPITAL" and ecrf.locked='true';

select * from ecrf where patientId = "CPCT02140007";
