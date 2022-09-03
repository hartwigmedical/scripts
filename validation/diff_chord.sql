SELECT
    MIN(pipeline) AS pipeline, hrd, hrStatus, hrdType, COUNT(*)
FROM
    (SELECT
        'OnlyInTruth' AS pipeline, hrd, hrStatus, hrdType
    FROM
        VARIABLE_TRUTH_DB_SCHEMA.chord
    WHERE
        sampleId = 'VARIABLE_TRUTH_SAMPLE_ID' UNION SELECT
        'OnlyInNew' AS pipeline, hrd, hrStatus, hrdType
    FROM
        VARIABLE_NEW_DB_SCHEMA.chord
    WHERE
        sampleId = 'VARIABLE_NEW_SAMPLE_ID') AS a
GROUP BY hrd , hrStatus , hrdType
HAVING COUNT(*) != 2;