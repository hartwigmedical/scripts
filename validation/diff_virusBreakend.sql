SELECT
    MIN(pipeline) AS pipeline, taxidGenus, nameGenus, nameSpecies, qcStatus, COUNT(*)
FROM
    (SELECT
        'OnlyInTruth' AS pipeline, taxidGenus, nameGenus, nameSpecies, qcStatus
    FROM
        VARIABLE_TRUTH_DB_SCHEMA.virusBreakend
    WHERE
        sampleId = 'VARIABLE_TRUTH_SAMPLE_ID' UNION SELECT
        'OnlyInNew' AS pipeline, taxidGenus, nameGenus, nameSpecies, qcStatus
    FROM
        VARIABLE_NEW_DB_SCHEMA.virusBreakend
    WHERE
        sampleId = 'VARIABLE_NEW_SAMPLE_ID') AS a
GROUP BY taxidGenus, nameGenus, nameSpecies, qcStatus
HAVING COUNT(*) != 2;