select patients.patientId, 
BISLICENR, BISLICE, BISITERADI, BIOTHLESSITE, CPCT, BIOTHIMAGE, BIOPTDT, BILESTYPE, BILESSITE, BILESNR, BILESLOC, BIIMAGE
from 
	(select distinct patientId from ecrf) patients
left join
   (select patientId, group_concat(itemValue separator ', ') as BISLICENR 
    from ecrf where item = 'FLD.BIOPS.BISLICENR' group by patientId) slicenr
on patients.patientId = slicenr.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as BISLICE
     from ecrf where item ='FLD.BIOPS.BISLICE' group by patientId) slice
on patients.patientId = slice.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as BISITERADI
     from ecrf where item ='FLD.BIOPS.BISITERADI' group by patientId) radi
on patients.patientId = radi.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as BIOTHLESSITE
     from ecrf where item ='FLD.BIOPS.BIOTHLESSITE' group by patientId) thles
on patients.patientId = thles.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as BIOTHIMAGE
     from ecrf where item ='FLD.BIOPS.BIOTHIMAGE' group by patientId) thimage
on patients.patientId = thimage.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as CPCT
    from ecrf where item ='FLD.BIOPS.CPCT' group by patientId) bt
on patients.patientId = bt.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as BIOPTDT
     from ecrf where item ='FLD.BIOPS.BIOPTDT' group by patientId) tdt
on patients.patientId = tdt.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as BILESTYPE
     from ecrf where item ='FLD.BIOPS.BILESTYPE' group by patientId) btype
on patients.patientId = btype.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as BILESSITE
     from ecrf where item ='FLD.BIOPS.BILESSITE' group by patientId) site
on patients.patientId = site.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as BILESNR
     from ecrf where item ='FLD.BIOPS.BILESNR' group by patientId) bnr
on patients.patientId = bnr.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as BILESLOC
     from ecrf where item ='FLD.BIOPS.BILESLOC' group by patientId) loc
on patients.patientId = loc.patientId
left join
    (select patientId, group_concat(itemValue separator ', ') as BIIMAGE
     from ecrf where item ='FLD.BIOPS.BIIMAGE' group by patientId) imag
on patients.patientId = imag.patientId
where patients.patientId in (select distinct patientId from clinical) order by patients.patientId