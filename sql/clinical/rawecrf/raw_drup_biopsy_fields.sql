select patients.patientId,
TBDAT, TBSITE, TBSITEOTH, TBLOC
from
	(select distinct patientId from drupEcrf) patients
left join
   (select patientId, group_concat(itemValue separator ', ') as TBDAT
    from drupEcrf where item = 'FLD.TBDAT' group by patientId) tbdat
on patients.patientId = tbdat.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as TBSITE
     from drupEcrf where item ='FLD.TBSITE' group by patientId) tbsite
on patients.patientId = tbsite.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as TBSITEOTH
     from drupEcrf where item ='FLD.TBSITEOTH' group by patientId) tbsiteoth
on patients.patientId = tbsiteoth.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as TBLOC
     from drupEcrf where item ='FLD.TBLOC' group by patientId) tbloc
on patients.patientId = tbloc.patientId
where patients.patientId in (select distinct patientId from clinical) order by patients.patientId