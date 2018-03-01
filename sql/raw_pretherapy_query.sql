select patients.patientId,
SYSTEMIC, SYSTEMICLNR, SYSTEMICREG, SYSTEMICLET, SYSTEMICTYPE, SYSTEMICCYCLES, SYSTEMICSTDTC, SYSTEMICENDTC, SYSTEMICRESP,
RADIOTHER,RADIOTHERLNR,RADIOTHERSITE,RADIOTHERAIMSP,RADIOTHERAIMSP,RADIOTHEREDTC,RADIOTHERDOSE,RADIOTHERENDTC
from 
	(select distinct patientId from ecrf) patients
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMIC 
    from ecrf where item = 'FLD.PRETHERAPY.SYSTEMIC' group by patientId) ic
on patients.patientId = ic.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMICLNR 
    from ecrf where item = 'FLD.PRETHERAPY.SYSTEMICLNR' group by patientId) iclnr
on patients.patientId = iclnr.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMICREG 
    from ecrf where item = 'FLD.PRETHERAPY.SYSTEMICREG' group by patientId) icreg
on patients.patientId = icreg.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMICLET 
    from ecrf where item = 'FLD.PRETHERAPY.SYSTEMICLET' group by patientId) iclet
on patients.patientId = iclet.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMICTYPE
    from ecrf where item = 'FLD.PRETHERAPY.SYSTEMICTYPE' group by patientId) ictype
on patients.patientId = ictype.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMICCYCLES 
    from ecrf where item = 'FLD.PRETHERAPY.SYSTEMICCYCLES' group by patientId) iccycles
on patients.patientId = iccycles.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMICSTDTC 
    from ecrf where item = 'FLD.PRETHERAPY.SYSTEMICSTDTC' group by patientId) icstdtc
on patients.patientId = icstdtc.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMICENDTC 
    from ecrf where item = 'FLD.PRETHERAPY.SYSTEMICENDTC' group by patientId) icendtc
on patients.patientId = icendtc.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as SYSTEMICRESP 
    from ecrf where item = 'FLD.PRETHERAPY.SYSTEMICRESP' group by patientId) icresp
on patients.patientId = icresp.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOTHER 
    from ecrf where item = 'FLD.PRETHERAPY.RADIOTHER' group by patientId) radiother
on patients.patientId = radiother.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOTHERLNR 
    from ecrf where item = 'FLD.PRETHERAPY.RADIOTHERLNR' group by patientId) radiolnr
on patients.patientId = radiolnr.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOTHERSITE
    from ecrf where item = 'FLD.PRETHERAPY.RADIOTHERSITE' group by patientId) radiosite
on patients.patientId = radiosite.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOTHERAIMSP
    from ecrf where item = 'FLD.PRETHERAPY.RADIOTHERAIMSP' group by patientId) radioaimsp
on patients.patientId = radioaimsp.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOTHEREDTC
    from ecrf where item = 'FLD.PRETHERAPY.RADIOTHEREDTC' group by patientId) radioedtc
on patients.patientId = radioedtc.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOTHERDOSE 
    from ecrf where item = 'FLD.PRETHERAPY.RADIOTHERDOSE' group by patientId) radiodose
on patients.patientId = radiodose.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as RADIOTHERENDTC 
    from ecrf where item = 'FLD.PRETHERAPY.RADIOTHERENDTC' group by patientId) radioendtc
on patients.patientId = radioendtc.patientId