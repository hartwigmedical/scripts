select patients.patientId,
ICDTC, INST, GEN, YOB, BASTTYP, BASTTOSP, TBTAKEN, TBDAT, TBSITE, TBLOC,DEATHDTC
from 
	(select distinct patientId from drupEcrf) patients
left join
   (select patientId, group_concat(itemValue separator ', ') as ICDTC 
    from drupEcrf where item = 'FLD.BAS.ICDTC' group by patientId) ic
on patients.patientId = ic.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as INST 
    from drupEcrf where item = 'FLD.REG.INST' group by patientId) hospital
on patients.patientId = hospital.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as GEN 
    from drupEcrf where item = 'FLD.CSF.GEN' group by patientId) gender
on patients.patientId = gender.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as YOB 
    from drupEcrf where item = 'FLD.CSF.YOB' group by patientId) yearofbirth
on patients.patientId = yearofbirth.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as BASTTYP 
    from drupEcrf where item = 'FLD.BAS.BASTTYP' group by patientId) tumloc
on patients.patientId = tumloc.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as BASTTOSP 
    from drupEcrf where item = 'FLD.BAS.BASTTOSP' group by patientId) tumlocother
on patients.patientId = tumlocother.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as TBTAKEN 
    from drupEcrf where item = 'FLD.BIOPT.TBTAKEN' group by patientId) biopsytaken
on patients.patientId = biopsytaken.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as TBDAT 
    from drupEcrf where item = 'FLD.BIOPT.TBDAT' group by patientId) biopsydate
on patients.patientId = biopsydate.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as TBSITE 
    from drupEcrf where item = 'FLD.BIOPT.TBSITE' group by patientId) biopsysite
on patients.patientId = biopsysite.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as TBLOC 
    from drupEcrf where item = 'FLD.BIOPT.TBLOC' group by patientId) biopsyloc
on patients.patientId = biopsyloc.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as DEATHDTC 
    from drupEcrf where item = 'FLD.EOT.DEATHDTC' group by patientId) death
on patients.patientId = death.patientId;