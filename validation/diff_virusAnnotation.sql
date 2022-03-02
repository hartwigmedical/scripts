SELECT
    MIN(pipeline) AS pipeline, taxid, virusName, interpretation, integrations, qcStatus, reported, COUNT(*)
FROM
    (SELECT
        'OnlyInTruth' AS pipeline, taxid, virusName, interpretation, integrations, qcStatus, reported
    FROM
        VARIABLE_TRUTH_DB_SCHEMA.virusAnnotation
    WHERE
        sampleId = 'VARIABLE_TRUTH_SAMPLE_ID' UNION SELECT
        'OnlyInNew' AS pipeline, taxid, virusName, interpretation, integrations, qcStatus, reported
    FROM
        VARIABLE_NEW_DB_SCHEMA.virusAnnotation
    WHERE
        sampleId = 'VARIABLE_NEW_SAMPLE_ID') AS a
GROUP BY taxid, virusName, interpretation, integrations, qcStatus, reported
HAVING COUNT(*) != 2;