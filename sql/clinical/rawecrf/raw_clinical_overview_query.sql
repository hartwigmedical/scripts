select patients.patientId, 
informedConsentDate, regDate1, regDate2, hospital1, hospital2,
primaryTumorLocation, primaryTumorLocationOther, entryStage,
biopsyDate, biopsyLocation, biopsyLocationOther,
treatmentGiven, treatmentStart, treatmentEnd, treatment, treatmentOther,
measurement, responseAssessmentDate, responseDate, response
from
	(select distinct patientId from cpctEcrf) patients
left join
   (select patientId, group_concat(itemValue separator ', ') as informedConsentDate
    from cpctEcrf where item = 'FLD.INFORMEDCONSENT.ICDTC' group by patientId) informed
on patients.patientId = informed.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as regDate1
    from cpctEcrf where item = 'FLD.ELIGIBILITY.REGDTC' group by patientId) reg1
on patients.patientId = reg1.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as regDate2
    from cpctEcrf where item = 'FLD.SELCRIT.NREGDTC' group by patientId) reg2
on patients.patientId = reg2.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as hospital1
    from cpctEcrf where item = 'FLD.ELIGIBILITY.HOSPITAL' group by patientId) hosp1
on patients.patientId = hosp1.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as hospital2
    from cpctEcrf where item = 'FLD.SELCRIT.NHOSPITAL' group by patientId) hosp2
on patients.patientId = hosp2.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as primaryTumorLocation
    from cpctEcrf where item = 'FLD.CARCINOMA.PTUMLOC' group by patientId) ptumloc
on patients.patientId = ptumloc.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as primaryTumorLocationOther
    from cpctEcrf where item = 'FLD.CARCINOMA.PTUMLOCS' group by patientId) ptumlocother
on patients.patientId = ptumlocother.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as entryStage 
    from cpctEcrf where item = 'FLD.CARCINOMA.ENTRYSTAGE' group by patientId) entryst
on patients.patientId = entryst.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as biopsyDate
     from cpctEcrf where item ='FLD.BIOPS.BIOPTDT' group by patientId) bioptdt
on ptumloc.patientId = bioptdt.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as biopsyLocation
     from cpctEcrf where item ='FLD.BIOPS.BILESSITE' group by patientId) bilessite
on bioptdt.patientId = bilessite.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as biopsyLocationOther
     from cpctEcrf where item ='FLD.BIOPS.BIOTHLESSITE' group by patientId) biothlessite
on bilessite.patientId = biothlessite.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as treatmentGiven
     from cpctEcrf where item ='FLD.TRTAFTER.SYSTEMICST' group by patientId) systemicst
on biothlessite.patientId = systemicst.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as treatmentStart
     from cpctEcrf where item ='FLD.TRTAFTER.SYSSTDT' group by patientId) sysstdt
on systemicst.patientId = sysstdt.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as treatmentEnd
     from cpctEcrf where item ='FLD.TRTAFTER.SYSENDT' group by patientId) sysendt
on sysstdt.patientId = sysendt.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as treatment
     from cpctEcrf where item ='FLD.TRTAFTER.PLANNEDTRT' group by patientId) plannedtrt
on sysendt.patientId = plannedtrt.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as treatmentOther
    from cpctEcrf where item ='FLD.TRTAFTER.SYSREGPOST' group by patientId) sysregpost
on plannedtrt.patientId = sysregpost.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as measurement
    from cpctEcrf where item ='FLD.TUMORMEASUREMENT.TMYN' group by patientId) tmyn
on sysregpost.patientId = tmyn.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as responseAssessmentDate
    from cpctEcrf where item ='FLD.TUMORMEASUREMENT.ASSDTC' group by patientId) assdtc
on tmyn.patientId = assdtc.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as responseDate
    from cpctEcrf where item ='FLD.TUMORMEASUREMENT.RESPONSEDTC' group by patientId) responsedtc
on assdtc.patientId = responsedtc.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as response
    from cpctEcrf where item ='FLD.TUMORMEASUREMENT.BESTRESPON' group by patientId) bestrespon
on responsedtc.patientId = bestrespon.patientId