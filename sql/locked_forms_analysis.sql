	SELECT count(*) AS 'form count', concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND locked='TRUE') / count(*)) , '%') AS 'locked percentage', 'all' AS form
	FROM formsMetadata
	WHERE tableName != 'drug'

UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND form = 'biopsy' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' and form = 'biopsy'

UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND form = 'treatment' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' and form = 'treatment'
    
UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND form = 'death' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' and form = 'death'
    
UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND form = 'demography' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' and form = 'demography'
    
UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND form = 'eligibility' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' and form = 'eligibility'
    
UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND form = 'primaryTumor' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' and form = 'primaryTumor'
    
UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND form = 'selectionCriteria' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' and form = 'selectionCriteria'
    
UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) FROM formsMetadata WHERE tableName != 'drug' AND form = 'treatmentResponse' AND locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' and form = 'treatmentResponse'

