SELECT
    MIN(pipeline) AS pipeline, event, germline, reported, treatment, onLabel, level, direction, sources, COUNT(*)
FROM
        (SELECT
            id, 'OnlyInTruth' AS pipeline, event, germline, reported, treatment, onLabel, level, direction, sources
        FROM
            VARIABLE_TRUTH_DB_SCHEMA.protect
        WHERE
            sampleId = 'VARIABLE_TRUTH_SAMPLE_ID' AND reported = 1
    UNION
        SELECT
            id, 'OnlyInNew' AS pipeline, event, germline, reported, treatment, onLabel, level, direction, sources
        FROM
            VARIABLE_NEW_DB_SCHEMA.protect
        WHERE
            sampleId = 'VARIABLE_NEW_SAMPLE_ID' AND reported = 1)
    AS a
GROUP BY 2 , 3 , 4, 5, 6, 7, 8, 9
HAVING COUNT(*) != 2 ORDER BY pipeline;