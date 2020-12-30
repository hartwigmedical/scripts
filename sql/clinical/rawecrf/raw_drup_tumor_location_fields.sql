select patients.patientId,
codeICD, desciptionICD, bastType, bastTosp, cohortName, ttype
from
	(select distinct patientId from drupEcrf) patients
left join
   (select patientId, group_concat(itemValue separator ', ') as codeICD
    from drupEcrf where item = 'FLD.ICD10Code' group by patientId) codeICD
on patients.patientId = codeICD.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as desciptionICD
    from drupEcrf where item = 'FLD.ICD10Descr' group by patientId) desciptionICD
on patients.patientId = desciptionICD.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as bastType
    from drupEcrf where item = 'FLD.BASTTYP' group by patientId) bastType
on patients.patientId = bastType.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as bastTosp
    from drupEcrf where item = 'FLD.BASTTOSP' group by patientId) bastTosp
on patients.patientId = bastTosp.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as cohortName
    from drupEcrf where item = 'FLD.COHORTNAME' group by patientId) cohortName
on patients.patientId = cohortName.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as ttype
    from drupEcrf where item = 'FLD.TTYPE' group by patientId) ttype
on patients.patientId = ttype.patientId