SELECT
    MIN(pipeline) AS pipeline,
    gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction,
    source, sourceEvent, sourceUrls, evidenceType, rangeRank, evidenceUrls, COUNT(*)
FROM
    (SELECT
        id,
            'OnlyInTruth' AS pipeline,
    gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction,
    source, sourceEvent, sourceUrls, evidenceType, rangeRank, evidenceUrls
    FROM
        protect
    WHERE
        sampleId = 'VARIABLE_TRUTH_SAMPLE_ID'
        UNION SELECT
        id,
            'OnlyInNew' AS pipeline,
            gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction,
    source, sourceEvent, sourceUrls, evidenceType, rangeRank, evidenceUrls
    FROM
        protect
    WHERE
        sampleId = 'VARIABLE_NEW_SAMPLE_ID') AS a
GROUP BY gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction,
    source, sourceEvent, sourceUrls, evidenceType, rangeRank, evidenceUrls
    HAVING COUNT(*) != 2
ORDER BY event;