USE hmfpatients;

	SELECT 
    count(*) AS samples, 'all' AS category
    FROM sample INNER JOIN patient ON sample.patientId = patient.id
    
UNION

	SELECT count(*), 'with known cancer type' 
    FROM sample INNER JOIN patient ON sample.patientId = patient.id
    WHERE NOT(isnull(primaryTumorCategory))
    
UNION

    SELECT count(*), 'with matched biopsy' 
    FROM sample INNER JOIN patient ON sample.patientId = patient.id
    LEFT JOIN biopsy ON biopsy.sampleId = sample.sampleId
    WHERE NOT(isnull(biopsy.id))

UNION

    SELECT count(*), 'with matched biopsy and known biopsy site' 
    FROM sample INNER JOIN patient ON sample.patientId = patient.id
    LEFT JOIN biopsy ON biopsy.sampleId = sample.sampleId
    WHERE NOT(isnull(biopsy.id)) and NOT(isnull(biopsy.biopsySite))
    
UNION

    SELECT count(*), 'with matched biopsy and treatment' 
    FROM sample INNER JOIN patient ON sample.patientId = patient.id
    LEFT JOIN biopsy ON biopsy.sampleId = sample.sampleId
    LEFT JOIN treatment on treatment.biopsyId = biopsy.id
    WHERE NOT(isnull(treatment.id))

UNION

	SELECT count(*), 'with matched biopsy, treatment and response' 
    FROM sample INNER JOIN patient ON sample.patientId = patient.id
    LEFT JOIN biopsy ON biopsy.sampleId = sample.sampleId
    LEFT JOIN treatment on treatment.biopsyId = biopsy.id
    LEFT JOIN firstMatchedTreatmentResponse on treatment.id = firstMatchedTreatmentResponse.treatmentId
    WHERE NOT(isnull(firstMatchedTreatmentResponse.id))