USE hmfpatients;

	SELECT count(*) as 'form count', concat(round(100 * (SELECT count(*) from formsMetadata where tableName != 'drug' and locked='TRUE') / count(*)) , '%') as 'locked percentage', 'all' as form
	FROM formsMetadata
	WHERE tableName != 'drug'

UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) from formsMetadata where tableName != 'drug' and form = 'biopsy' and locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' and form = 'biopsy'

UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) from formsMetadata where tableName != 'drug' and form = 'treatment' and locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' and form = 'treatment'
    
UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) from formsMetadata where tableName != 'drug' and form = 'death' and locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' and form = 'death'
    
UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) from formsMetadata where tableName != 'drug' and form = 'demography' and locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' and form = 'demography'
    
UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) from formsMetadata where tableName != 'drug' and form = 'eligibility' and locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' and form = 'eligibility'
    
UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) from formsMetadata where tableName != 'drug' and form = 'primaryTumor' and locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' and form = 'primaryTumor'
    
UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) from formsMetadata where tableName != 'drug' and form = 'selectionCriteria' and locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' and form = 'selectionCriteria'
    
UNION

	SELECT count(*), concat(round(100 * (SELECT count(*) from formsMetadata where tableName != 'drug' and form = 'treatmentResponse' and locked='TRUE') / count(*)), '%'), form
	FROM formsMetadata
	WHERE tableName != 'drug' and form = 'treatmentResponse'

