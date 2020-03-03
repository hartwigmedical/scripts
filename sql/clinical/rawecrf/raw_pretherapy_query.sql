select patients.patientId,
SYSTEMIC, SYSTEMICLNR, SYSTEMICREG, SYSTEMICLET, SYSTEMICTYPE, SYSTEMICCYCLES, SYSTEMICSTDTC, SYSTEMICENDTC, SYSTEMICRESP,
RADIOTHER,RADIOTHERLNR,RADIOTHERSITE,RADIOTHERAIMSP,RADIOTHERAIMSP,RADIOTHEREDTC,RADIOTHERDOSE,RADIOTHERENDTC
from 
	(select distinct patientId from cpctEcrf) patients
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMIC 
    from cpctEcrf where item = 'FLD.PRETHERAPY.SYSTEMIC' group by patientId) ic
on patients.patientId = ic.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMICLNR 
    from cpctEcrf where item = 'FLD.PRETHERAPY.SYSTEMICLNR' group by patientId) iclnr
on patients.patientId = iclnr.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMICREG 
    from cpctEcrf where item = 'FLD.PRETHERAPY.SYSTEMICREG' group by patientId) icreg
on patients.patientId = icreg.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMICLET 
    from cpctEcrf where item = 'FLD.PRETHERAPY.SYSTEMICLET' group by patientId) iclet
on patients.patientId = iclet.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMICTYPE
    from cpctEcrf where item = 'FLD.PRETHERAPY.SYSTEMICTYPE' group by patientId) ictype
on patients.patientId = ictype.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMICCYCLES 
    from cpctEcrf where item = 'FLD.PRETHERAPY.SYSTEMICCYCLES' group by patientId) iccycles
on patients.patientId = iccycles.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMICSTDTC 
    from cpctEcrf where item = 'FLD.PRETHERAPY.SYSTEMICSTDTC' group by patientId) icstdtc
on patients.patientId = icstdtc.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMICENDTC 
    from cpctEcrf where item = 'FLD.PRETHERAPY.SYSTEMICENDTC' group by patientId) icendtc
on patients.patientId = icendtc.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMICRESP 
    from cpctEcrf where item = 'FLD.PRETHERAPY.SYSTEMICRESP' group by patientId) icresp
on patients.patientId = icresp.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOTHER 
    from cpctEcrf where item = 'FLD.PRETHERAPY.RADIOTHER' group by patientId) radiother
on patients.patientId = radiother.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOTHERLNR 
    from cpctEcrf where item = 'FLD.PRETHERAPY.RADIOTHERLNR' group by patientId) radiolnr
on patients.patientId = radiolnr.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOTHERSITE
    from cpctEcrf where item = 'FLD.PRETHERAPY.RADIOTHERSITE' group by patientId) radiosite
on patients.patientId = radiosite.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOTHERAIMSP
    from cpctEcrf where item = 'FLD.PRETHERAPY.RADIOTHERAIMSP' group by patientId) radioaimsp
on patients.patientId = radioaimsp.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOTHEREDTC
    from cpctEcrf where item = 'FLD.PRETHERAPY.RADIOTHEREDTC' group by patientId) radioedtc
on patients.patientId = radioedtc.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOTHERDOSE 
    from cpctEcrf where item = 'FLD.PRETHERAPY.RADIOTHERDOSE' group by patientId) radiodose
on patients.patientId = radiodose.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOTHERENDTC 
    from cpctEcrf where item = 'FLD.PRETHERAPY.RADIOTHERENDTC' group by patientId) radioendtc
on patients.patientId = radioendtc.patientId