SELECT
    'Truth' AS pipeline, sampleId, refMeanCoverage, refCoverage10xPercentage, refCoverage20xPercentage, tumorMeanCoverage,
    tumorCoverage30xPercentage, tumorCoverage60xPercentage
FROM
    VARIABLE_TRUTH_DB_SCHEMA.metric
WHERE
    sampleId = 'VARIABLE_TRUTH_SAMPLE_ID'
UNION SELECT
    'New' AS pipeline, sampleId, refMeanCoverage, refCoverage10xPercentage, refCoverage20xPercentage, tumorMeanCoverage,
    tumorCoverage30xPercentage, tumorCoverage60xPercentage
FROM
    VARIABLE_NEW_DB_SCHEMA.metric
WHERE
    sampleId = 'VARIABLE_NEW_SAMPLE_ID';
