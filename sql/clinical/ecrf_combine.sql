select patients.patientId, primaryTumorLocation, biopsyDate, biopsyLocation, biopsyLocationOther
-- , treatmentGiven, treatmentStart, treatmentEnd, treatment
-- , treatmentOther,
-- measurement, responseAsseddmentDate, responseDate, response
from 
	(select distinct patientId from ecrf) patients
left join
   (select patientId, group_concat(itemValue separator ', ') as primaryTumorLocation 
    from ecrf where item = 'FLD.CARCINOMA.PTUMLOC' group by patientId) ptumloc
on patients.patientId = ptumloc.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as biopsyDate
     from ecrf where item ='FLD.BIOPS.BIOPTDT' group by patientId) bioptdt
on ptumloc.patientId = bioptdt.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as biopsyLocation
     from ecrf where item ='FLD.BIOPS.BILESSITE' group by patientId) bilessite
on bioptdt.patientId = bilessite.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as biopsyLocationOther
     from ecrf where item ='FLD.BIOPS.BIOTHLESSITE' group by patientId) biothlessite
on bilessite.patientId = biothlessite.patientId
-- inner join
--     (select patientId, group_concat(itemValue separator ', ') as treatmentGiven
--      from ecrf where item ='FLD.TRTAFTER.SYSTEMICST' group by patientId) systemicst
-- on biothlessite.patientId = systemicst.patientId
-- inner join
--     (select patientId, group_concat(itemValue separator ', ') as treatmentStart
--      from ecrf where item ='FLD.TRTAFTER.SYSSTDT' group by patientId) sysstdt
-- on systemicst.patientId = sysstdt.patientId
-- inner join
--     (select patientId, group_concat(itemValue separator ', ') as treatmentEnd
--      from ecrf where item ='FLD.TRTAFTER.SYSENDT' group by patientId) sysendt
-- on sysstdt.patientId = sysendt.patientId
-- inner join
--     (select patientId, group_concat(itemValue separator ', ') as treatment
--      from ecrf where item ='FLD.TRTAFTER.PLANNEDTRT' group by patientId) plannedtrt
-- on sysendt.patientId = plannedtrt.patientId
-- -- inner join
--     (select patientId, group_concat(itemValue separator ', ') as treatmentOther
--      from ecrf where item ='FLD.TRTAFTER.SYSREGPOST' group by patientId) sysregpost
-- on plannedtrt.patientId = sysregpost.patientId
-- inner join
--     (select patientId, group_concat(itemValue separator ', ') as measurement
--      from ecrf where item ='FLD.TUMORMEASUREMENT.TMYN' group by patientId) tmyn
-- on sysregpost.patientId = tmyn.patientId
-- inner join
--     (select patientId, group_concat(itemValue separator ', ') as responseAsseddmentDate
--      from ecrf where item ='FLD.TUMORMEASUREMENT.ASSDTC' group by patientId) assdtc
-- on tmyn.patientId = assdtc.patientId
-- inner join
--     (select patientId, group_concat(itemValue separator ', ') as responseDate
--      from ecrf where item ='FLD.TUMORMEASUREMENT.RESPONSEDTC' group by patientId) responsedtc
-- on assdtc.patientId = responsedtc.patientId
-- inner join
--     (select patientId, group_concat(itemValue separator ', ') as response
--      from ecrf where item ='FLD.TUMORMEASUREMENT.BESTRESPON' group by patientId) bestrespon
-- on responsedtc.patientId = bestrespon.patientId
-- where ptumloc.patientId in ("CPCT02010003", "CPCT02010007", "CPCT02010008")
