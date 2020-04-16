select patientIdentifier, treatmentGiven as treatmentInAvl, startDate, endDate, name, type, mechanism
from treatment inner join patient on patient.id = treatment.patientId
where patientIdentifier like '%WIDE%' order by 1,2,3,4;