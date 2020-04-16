select patientIdentifier, responseDate, response
from treatmentResponse inner join patient on patient.id = treatmentResponse.patientId
where patientIdentifier like '%WIDE%' order by 1,2;