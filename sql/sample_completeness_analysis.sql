	SELECT
    count(*) AS samples, 'all' AS category
    FROM sample
	WHERE sampleId like '%CPCT%'

UNION

	SELECT count(*), 'with known cancer type'
    FROM baseline INNER JOIN patient ON baseline.patientId = patient.id
    WHERE NOT(isnull(cancerType))
    
UNION

    SELECT count(*), 'with matched biopsy' 
    FROM sample LEFT JOIN biopsy ON biopsy.sampleId = sample.sampleId
    WHERE NOT(isnull(biopsy.id))

UNION

SELECT count(*), 'with matched biopsy and known biopsy site' 
    FROM sample LEFT JOIN biopsy ON biopsy.sampleId = sample.sampleId
    WHERE NOT(isnull(biopsy.id)) and NOT(isnull(biopsy.biopsySite))
    
UNION

    SELECT count(*), 'with matched biopsy and treatment' 
    FROM sample
    LEFT JOIN biopsy ON biopsy.sampleId = sample.sampleId
    LEFT JOIN treatment on treatment.biopsyId = biopsy.id
    WHERE NOT(isnull(treatment.id))

UNION

	SELECT count(*), 'with matched biopsy, treatment and response' 
    FROM sample
    LEFT JOIN biopsy ON biopsy.sampleId = sample.sampleId
    LEFT JOIN treatment on treatment.biopsyId = biopsy.id
    LEFT JOIN firstMatchedTreatmentResponse on treatment.id = firstMatchedTreatmentResponse.treatmentId
    WHERE NOT(isnull(firstMatchedTreatmentResponse.id))