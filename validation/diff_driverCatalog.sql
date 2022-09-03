SELECT
    MIN(pipeline) AS pipeline,
    chromosome,
    chromosomeBand,
    gene,
    transcriptId,
    category,
    driver,
    biallelic,
    COUNT(*),
    MIN(driverLikelihood)
FROM
    (SELECT
        'OnlyInTruth' AS pipeline,
            chromosome,
            chromosomeBand,
            gene,
            transcriptId,
            category,
            driver,
            driverLikelihood,
            biallelic
    FROM
        VARIABLE_TRUTH_DB_SCHEMA.driverCatalog
    WHERE
        sampleId = 'VARIABLE_TRUTH_SAMPLE_ID' UNION SELECT
        'OnlyInNew' AS pipeline,
            chromosome,
            chromosomeBand,
            gene,
            transcriptId,
            category,
            driver,
            driverLikelihood,
            biallelic
    FROM
        VARIABLE_NEW_DB_SCHEMA.driverCatalog
    WHERE
        sampleId = 'VARIABLE_NEW_SAMPLE_ID') AS a
GROUP BY chromosome , chromosomeBand , gene , transcriptId , category , driver , biallelic
HAVING COUNT(*) != 2;