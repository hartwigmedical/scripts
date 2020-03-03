select patients.patientId, icdtc, adconsent, adtissues, adpid, adtum, tstam81, tstam82, tstam83, tstam84, tstam85, tstam86, tstam87, tstam88, tstam89
from
    (select distinct patientId from cpctEcrf) patients
left join
   (select patientId, group_concat(itemValue separator ', ') as icdtc
    from cpctEcrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ICDTC' group by patientId) ICDTC
on patients.patientId = ICDTC.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as adconsent
    from cpctEcrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ADCONSENT' group by patientId) ADCONSENT
on ICDTC.patientId = ADCONSENT.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as adtissues
    from cpctEcrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ADTISSUES' group by patientId) ADTISSUES
on ADCONSENT.patientId = ADTISSUES.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as adpid
    from cpctEcrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ADPID' group by patientId) ADPID
on ADTISSUES.patientId = ADPID.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as adtum
    from cpctEcrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ADTUM' group by patientId) ADTUM
on ADPID.patientId = ADTUM.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as tstam81
    from cpctEcrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM81' group by patientId) TSTAM81
on ADTUM.patientId = TSTAM81.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as tstam82
    from cpctEcrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM82' group by patientId) TSTAM82
on ADTUM.patientId = TSTAM82.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as tstam83
    from cpctEcrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM83' group by patientId) TSTAM83
on ADTUM.patientId = TSTAM83.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as tstam84
    from cpctEcrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM84' group by patientId) TSTAM84
on ADTUM.patientId = TSTAM84.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as tstam85
    from cpctEcrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM85' group by patientId) TSTAM85
on ADTUM.patientId = TSTAM85.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as tstam86
    from cpctEcrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM86' group by patientId) TSTAM86
on ADTUM.patientId = TSTAM86.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as tstam87
    from cpctEcrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM87' group by patientId) TSTAM87
on ADTUM.patientId = TSTAM87.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as tstam88
    from cpctEcrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM88' group by patientId) TSTAM88
on ADTUM.patientId = TSTAM88.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as tstam89
    from cpctEcrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM89' group by patientId) TSTAM89
on TSTAM86.patientId = TSTAM89.patientId