select patients.patientId,
INFORMED_CONSENT_DATE, PLANNED_TREATMENT, PRIMARY_TUMOR_LOCATION, PRIMARY_TUMOR_LOCATION_OTHER,
BIOPT_TAKEN, BIOPT_DATE, BIOPSY_SITE, BIOPSY_SITE_OTHER, BIOPSY_LOCATION, BIOPSY_EVALUABLE
from
	(select distinct patientId from ecrf) patients
left join
   (select patientId, group_concat(itemValue separator ', ') as INFORMED_CONSENT_DATE
    from ecrf where item = 'FLD.INFORMEDCONSENT.ICDTC' group by patientId) icdtc
on patients.patientId = icdtc.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as PLANNED_TREATMENT
    from ecrf where item = 'FLD.ELIGIBILITY.PLANTRT' group by patientId) plantrt
on patients.patientId = plantrt.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as PRIMARY_TUMOR_LOCATION
    from ecrf where item = 'FLD.CARCINOMA.PTUMLOC' group by patientId) ptumloc
on patients.patientId = ptumloc.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as PRIMARY_TUMOR_LOCATION_OTHER
    from ecrf where item = 'FLD.CARCINOMA.PTUMLOCS' group by patientId) ptumlocs
on patients.patientId = ptumlocs.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as BIOPT_TAKEN
    from ecrf where item = 'FLD.BIOPS.CPCT' group by patientId) cpcttkn
on patients.patientId = cpcttkn.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as BIOPT_DATE
    from ecrf where item = 'FLD.BIOPS.BIOPTDT' group by patientId) biopttdt
on patients.patientId = biopttdt.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as BIOPSY_SITE
    from ecrf where item = 'FLD.BIOPS.BILESSITE' group by patientId) bilessite
on patients.patientId = bilessite.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as BIOPSY_SITE_OTHER
    from ecrf where item = 'FLD.BIOPS.BIOTHLESSITE' group by patientId) biothlessite
on patients.patientId = biothlessite.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as BIOPSY_LOCATION
    from ecrf where item = 'FLD.BIOPS.BILESLOC' group by patientId) bilesloc
on patients.patientId = bilesloc.patientId
left join
   (select patientId, group_concat(itemValue separator ', ') as BIOPSY_EVALUABLE
    from ecrf where item = 'FLD.BIOPS.BIOPEFS' group by patientId) biopefs
on patients.patientId = biopefs.patientId
where patients.patientId IN ("XXX") order by patients.patientId