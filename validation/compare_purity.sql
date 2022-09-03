SELECT
    'Truth' AS pipeline, version, sampleId, gender, fitMethod, qcStatus, purity, normFactor, score, somaticPenalty, ploidy,
    diploidProportion, polyclonalProportion, wholeGenomeDuplication, minPurity, maxPurity, minDiploidProportion, maxDiploidProportion,
    msIndelsPerMb, msStatus
FROM
    VARIABLE_TRUTH_DB_SCHEMA.purity
WHERE
    sampleId = 'VARIABLE_TRUTH_SAMPLE_ID'
UNION SELECT 
    'New' AS pipeline, version, sampleId, gender, fitMethod, qcStatus, purity, normFactor,score, somaticPenalty, ploidy,
    diploidProportion, polyclonalProportion, wholeGenomeDuplication,minPurity, maxPurity, minDiploidProportion, maxDiploidProportion,
    msIndelsPerMb, msStatus
FROM
    VARIABLE_NEW_DB_SCHEMA.purity
WHERE
    sampleId = 'VARIABLE_NEW_SAMPLE_ID';
