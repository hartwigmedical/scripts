USE hmfpatients;

	SELECT count(*) as findings, 'any' as level, 'any' as item, 'any' as locked FROM clinicalFindings

UNION

	SELECT count(*), 'any', 'any', 'yes'
	FROM clinicalFindings
	WHERE formLocked = 'TRUE'

UNION

	SELECT count(*), level, 'any', 'any'
	FROM clinicalFindings
	GROUP by level

UNION

	SELECT count(*), 'any', ecrfItem, 'any'
	FROM clinicalFindings
	GROUP BY level, ecrfItem