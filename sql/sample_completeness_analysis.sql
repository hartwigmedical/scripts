	SELECT
    count(*) AS samples, 'all' AS category
    FROM sample
	WHERE sample.sampleId like '%CPCT%'

UNION

	SELECT count(*), 'with known cancer type'
	FROM sample
	INNER JOIN patient ON sample.patientId = patient.id
	INNER JOIN baseline ON baseline.patientId = patient.id
	WHERE sample.sampleId like '%CPCT%' AND NOT(isnull(cancerType))
    
UNION

    SELECT count(*), 'with matched biopsy' 
    FROM sample LEFT JOIN biopsy ON biopsy.sampleId = sample.sampleId
    WHERE sample.sampleId like '%CPCT%' AND NOT(isnull(biopsy.id))

UNION

SELECT count(*), 'with matched biopsy and known biopsy site' 
    FROM sample LEFT JOIN biopsy ON biopsy.sampleId = sample.sampleId
    WHERE sample.sampleId like '%CPCT%' AND NOT(isnull(biopsy.id)) and NOT(isnull(biopsy.biopsySite))
    
UNION

    SELECT count(*), 'with matched biopsy and treatment' 
    FROM sample
    LEFT JOIN biopsy ON biopsy.sampleId = sample.sampleId
    LEFT JOIN treatment on treatment.biopsyId = biopsy.id
    WHERE sample.sampleId like '%CPCT%' AND NOT(isnull(treatment.id))

UNION

	SELECT count(*), 'with matched biopsy, treatment and response' 
    FROM sample
    LEFT JOIN biopsy ON biopsy.sampleId = sample.sampleId
    LEFT JOIN treatment on treatment.biopsyId = biopsy.id
    LEFT JOIN firstMatchedTreatmentResponse on treatment.id = firstMatchedTreatmentResponse.treatmentId
    WHERE sample.sampleId like '%CPCT%' AND NOT(isnull(firstMatchedTreatmentResponse.id))