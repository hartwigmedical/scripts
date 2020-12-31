select patients.patientId,
location, locationOther, selectLocation
from
	(select distinct patientId from cpctEcrf) patients
left join
   (select patientId, group_concat(itemValue separator ', ') as location
    from cpctEcrf where item = 'FLD.CARCINOMA.PTUMLOC' group by patientId) location
on patients.patientId = location.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as locationOther
    from cpctEcrf where item = 'FLD.CARCINOMA.PTUMLOCS' group by patientId) locationOther
on patients.patientId = locationOther.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as selectLocation
    from cpctEcrf where item = 'FLD.SELCRIT.SELPTM' group by patientId) selectLocation
on patients.patientId = selectLocation.patientId
