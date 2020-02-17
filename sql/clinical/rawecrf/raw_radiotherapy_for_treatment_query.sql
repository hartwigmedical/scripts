select patients.patientId,
RADIOTHERST, RADIOSITE, RADIOAIM, RADIOTHERAIMST, RADIOIRRAD, RADIODOSE, RADIOSTDTC, RADIOENDTC, RADIORESP
from 
	(select distinct patientId from cpctEcrf) patients
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOTHERST 
    from cpctEcrf where item = 'FLD.TRTAFTER.RADIOTHERST' group by patientId) therst
on patients.patientId = therst.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOSITE 
    from cpctEcrf where item = 'FLD.TRTAFTER.RADIOSITE' group by patientId) site
on patients.patientId = site.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOAIM 
    from cpctEcrf where item = 'FLD.TRTAFTER.RADIOAIM' group by patientId) aim
on patients.patientId = aim.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOTHERAIMST 
    from cpctEcrf where item = 'FLD.TRTAFTER.RADIOTHERAIMST' group by patientId) aimst
on patients.patientId = aimst.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOIRRAD 
    from cpctEcrf where item = 'FLD.TRTAFTER.RADIOIRRAD' group by patientId) irrad
on patients.patientId = irrad.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIODOSE
    from cpctEcrf where item = 'FLD.TRTAFTER.RADIODOSE' group by patientId) dose
on patients.patientId = dose.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOSTDTC 
    from cpctEcrf where item = 'FLD.TRTAFTER.RADIOSTDTC' group by patientId) stdtc
on patients.patientId = stdtc.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOENDTC 
    from cpctEcrf where item = 'FLD.TRTAFTER.RADIOENDTC' group by patientId) endtc
on patients.patientId = endtc.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIORESP 
    from cpctEcrf where item = 'FLD.TRTAFTER.RADIORESP' group by patientId) resp
on patients.patientId = resp.patientId
