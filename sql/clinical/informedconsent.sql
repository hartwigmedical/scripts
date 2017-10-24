select patients.patientId, icdtc, adconsent, adtissues, adpid, adtum, tstam86, tstam87

from
    (select distinct patientId from ecrf) patients
left join
   (select patientId, group_concat(itemValue separator ', ') as icdtc
    from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ICDTC' group by patientId) ICDTC
on patients.patientId = ICDTC.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as adconsent
    from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ADCONSENT' group by patientId) ADCONSENT
on ICDTC.patientId = ADCONSENT.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as adtissues
    from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ADTISSUES' group by patientId) ADTISSUES
on ADCONSENT.patientId = ADTISSUES.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as adpid
    from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ADPID' group by patientId) ADPID
on ADTISSUES.patientId = ADPID.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as adtum
    from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.ADTUM' group by patientId) ADTUM
on ADPID.patientId = ADTUM.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as tstam86
    from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM86' group by patientId) TSTAM86
on ADTUM.patientId = TSTAM86.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as tstam87
    from ecrf where fieldName ='BASELINE.INFORMEDCONSENT.INFORMEDCONSENT.TSTAM87' group by patientId) TSTAM87
on TSTAM86.patientId = TSTAM87.patientId