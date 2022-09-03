SELECT
    MIN(pipeline) AS pipeline,
    gene,
    transcript,
    isCanonical,
    event,
    eventIsHighDriver,
    germline,
    reported,
    treatment,
    onLabel,
    level,
    direction,
    source,
    sourceEvent,
    sourceUrls,
    evidenceType,
    rangeRank,
    evidenceUrls,
    COUNT(*)
FROM
    (SELECT
        id,
            'OnlyInTruth' AS pipeline,
            gene,
			transcript,
			isCanonical,
			event,
			eventIsHighDriver,
			germline,
			reported,
			treatment,
			onLabel,
			level,
			direction,
			source,
			sourceEvent,
			sourceUrls,
			evidenceType,
			rangeRank,
			evidenceUrls
    FROM
        VARIABLE_TRUTH_DB_SCHEMA.protect
    WHERE
        sampleId = 'VARIABLE_TRUTH_SAMPLE_ID'
            AND reported = 1 UNION SELECT
        id,
            'OnlyInNew' AS pipeline,
            gene,
			transcript,
			isCanonical,
			event,
			eventIsHighDriver,
			germline,
			reported,
			treatment,
			onLabel,
			level,
			direction,
			source,
			sourceEvent,
			sourceUrls,
			evidenceType,
			rangeRank,
			evidenceUrls
    FROM
        VARIABLE_NEW_DB_SCHEMA.protect
    WHERE
        sampleId = 'VARIABLE_NEW_SAMPLE_ID'
            AND reported = 1) AS a
GROUP BY gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction, source, sourceEvent, sourceUrls, evidenceType, rangeRank, evidenceUrls
HAVING COUNT(*) != 2
ORDER BY event;