select patientId, sampleId, informedConsentDate, primaryTumorLocation, cancerSubtype, gender, birthyear, biopsyDate, biopsySite, biopsyLocation
from clinical
where patientId like '%WIDE%' order by 1;