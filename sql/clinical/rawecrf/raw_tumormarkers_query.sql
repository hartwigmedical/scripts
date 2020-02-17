select patients.patientId,
LBSTAT_TUMORMARKERS,LBDTC_TUMORMARKERS, LBNAM_TUMORMARKERS, LBSEQ_TUMORMARKERS, LBTERM_TUMORMARKERS, LBTERMSP_TUMORMARKERS, 
LBORRES_TUMORMARKERS, LBORRESU_TUMORMARKERS,LBCAT_TUMORMARKERS
from 
	(select distinct patientId from cpctEcrf) patients
left join
   (select patientId, group_concat(itemValue separator ', ') as LBSTAT_TUMORMARKERS 
    from cpctEcrf where item = 'FLD.RESPONSE.LBSTAT_TUMORMARKERS' group by patientId) lbstat
on patients.patientId = lbstat.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as LBDTC_TUMORMARKERS 
    from cpctEcrf where item = 'FLD.RESPONSE.LBDTC_TUMORMARKERS' group by patientId) lbdtc
on patients.patientId = lbdtc.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as LBNAM_TUMORMARKERS 
    from cpctEcrf where item = 'FLD.RESPONSE.LBNAM_TUMORMARKERS' group by patientId) lbnam
on patients.patientId = lbnam.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as LBSEQ_TUMORMARKERS 
    from cpctEcrf where item = 'FLD.RESPONSE.LBSEQ_TUMORMARKERS' group by patientId) lbseq
on patients.patientId = lbseq.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as LBTERM_TUMORMARKERS 
    from cpctEcrf where item = 'FLD.RESPONSE.LBTERM_TUMORMARKERS' group by patientId) lbterm
on patients.patientId = lbterm.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as LBTERMSP_TUMORMARKERS 
    from cpctEcrf where item = 'FLD.RESPONSE.LBTERMSP_TUMORMARKERS' group by patientId) lbtermsp
on patients.patientId = lbtermsp.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as LBORRES_TUMORMARKERS
    from cpctEcrf where item = 'FLD.RESPONSE.LBORRES_TUMORMARKERS' group by patientId) lborres
on patients.patientId = lborres.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as LBORRESU_TUMORMARKERS 
    from cpctEcrf where item = 'FLD.RESPONSE.LBORRESU_TUMORMARKERS ' group by patientId) lborresu
on patients.patientId = lborresu.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as LBCAT_TUMORMARKERS 
    from cpctEcrf where item = 'FLD.RESPONSE.LBCAT_TUMORMARKERS' group by patientId) lbcat
on patients.patientId = lbcat.patientId
