select patients.patientId, 
CTCT2YN, CPCTCNT, CPCTPN
from 
	(select distinct patientId from drupEcrf) patients
left join
   (select patientId, group_concat(itemValue separator ', ') as CTCT2YN 
    from drupEcrf where item = 'FLD.CTCT2YN' group by patientId) yesno
on patients.patientId = yesno.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as CPCTCNT
     from drupEcrf where item ='FLD.CPCTCNT' group by patientId) center
on patients.patientId = center.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as CPCTPN
     from drupEcrf where item ='FLD.CPCTPN' group by patientId) id
on patients.patientId = id.patientId