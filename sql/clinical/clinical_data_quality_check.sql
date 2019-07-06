	SELECT count(*) AS findings, 'any' AS level, 'any' AS message, 'any' AS locked FROM clinicalFindings

UNION

	SELECT count(*), 'any', 'any', 'yes'
	FROM clinicalFindings
	WHERE formLocked = 'TRUE'

UNION

	SELECT count(*), level, 'any', 'any'
	FROM clinicalFindings
	GROUP BY level

UNION

	SELECT count(*), 'any', message, 'any'
	FROM clinicalFindings
	GROUP BY level, message