	SELECT count(*) AS 'form count', concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND locked='TRUE') / count(*)) , '%') AS 'locked percentage', 'all' AS form
	FROM formsMetadata
	WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug'

UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'informedConsent' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'informedConsent'

UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'eligibility' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'eligibility'

UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'selectionCriteria' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'selectionCriteria'

UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'demography' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'demography'
    
UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'primaryTumor' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'primaryTumor'

UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'pretreatment' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'pretreatment'
    
UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'biopsy' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'biopsy'

UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'treatment' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'treatment'
            
UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'treatmentResponse' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'treatmentResponse'
    
UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'tumorMarker' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'tumorMarker'
    
UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'death' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' AND tableName != 'preTreatmentDrug' AND form = 'death'

