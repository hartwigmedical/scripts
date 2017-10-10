select patients.patientId, baselinedtc, icdtc, adconsent, adbiopsy, adgen, adtissues, addna, adpid, adtum, 
adblood, amcons, amd17gr, tstam84,
tstam85, tstam86, tstam87,  tstam81, tstam82,
tstam83, tstam88, tstam89
 
from 
	(select distinct patientId from ecrf) patients
left join
   (select patientId, group_concat(itemValue separator ', ') as baselinedtc 
    from ecrf where fieldName = 'BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.BASELINEDTC' group by patientId) BASELINEDTC
on patients.patientId = BASELINEDTC.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as icdtc
     from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ICDTC' group by patientId) ICDTC
on BASELINEDTC.patientId = ICDTC.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as adconsent
     from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ADCONSENT' group by patientId) ADCONSENT
on ICDTC.patientId = ADCONSENT.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as adbiopsy
     from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ADBIOPSY' group by patientId) ADBIOPSY
on ADCONSENT.patientId = ADBIOPSY.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as adgen
     from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ADGEN' group by patientId) ADGEN
on ADBIOPSY.patientId = ADGEN.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as adtissues
     from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ADTISSUES' group by patientId) ADTISSUES
on ADGEN.patientId = ADTISSUES.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as addna
     from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ADDNA' group by patientId) ADDNA
on ADTISSUES.patientId = ADDNA.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as adpid
     from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ADPID' group by patientId) ADPID
on ADDNA.patientId = ADPID.patientId
left join
	(select patientId, group_concat(itemValue separator ', ') as adtum
	from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ADTUM' group by patientId) ADTUM
on ADPID.patientId = ADTUM.patientId
left join
	(select patientId, group_concat(itemValue separator ', ') as adblood
	from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ADBLOOD' group by patientId) ADBLOOD
on ADTUM.patientId = ADBLOOD.patientId
left join
	(select patientId, group_concat(itemValue separator ', ') as amcons
	from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.AMCONS' group by patientId) AMCONS
on ADBLOOD.patientId = AMCONS.patientId
left join
	(select patientId, group_concat(itemValue separator ', ') as amd17gr
	from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.AMD17GR' group by patientId) AMD17GR
on AMCONS.patientId = AMD17GR.patientId
left join
	(select patientId, group_concat(itemValue separator ', ') as tstam84
	from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM84' group by patientId) TSTAM84
on AMD17GR.patientId = TSTAM84.patientId
left join
	(select patientId, group_concat(itemValue separator ', ') as tstam85
	from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM85' group by patientId) TSTAM85
on TSTAM84.patientId = TSTAM85.patientId
left join
	(select patientId, group_concat(itemValue separator ', ') as tstam86
	from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM86' group by patientId) TSTAM86
on TSTAM85.patientId = TSTAM86.patientId
left join
	(select patientId, group_concat(itemValue separator ', ') as tstam87
	from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM87' group by patientId) TSTAM87
on TSTAM86.patientId = TSTAM87.patientId
left join
	(select patientId, group_concat(itemValue separator ', ') as tstam81
	from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM81' group by patientId) TSTAM81
on TSTAM87.patientId = TSTAM81.patientId
left join
	(select patientId, group_concat(itemValue separator ', ') as tstam82
	from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM82' group by patientId) TSTAM82
on TSTAM81.patientId = TSTAM82.patientId
left join
	(select patientId, group_concat(itemValue separator ', ') as tstam83
	from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM83' group by patientId) TSTAM83
on TSTAM82.patientId = TSTAM83.patientId
left join
	(select patientId, group_concat(itemValue separator ', ') as tstam88
	from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM88' group by patientId) TSTAM88
on TSTAM83.patientId = TSTAM88.patientId
left join
	(select patientId, group_concat(itemValue separator ', ') as tstam89
	from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM89' group by patientId) TSTAM89
on TSTAM88.patientId = TSTAM89.patientId

