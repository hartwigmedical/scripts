SELECT 
    'Truth' AS pipeline,
    version,
    sampleId,
    gender,
    status,
    qcStatus,
    purity,
    normFactor,
    score,
    somaticPenalty,
    ploidy,
    diploidProportion,
    polyclonalProportion,
    wholeGenomeDuplication,
    minPurity,
    maxPurity,
    minDiploidProportion,
    maxDiploidProportion,
    msIndelsPerMb,
    msStatus
FROM
    reference_validation_sets.purity
WHERE
    sampleId = 'COLO829v003T' 
UNION SELECT 
    'New' AS pipeline,
    version,
    sampleId,
    gender,
    status,
    qcStatus,
    purity,
    normFactor,
    score,
    somaticPenalty,
    ploidy,
    diploidProportion,
    polyclonalProportion,
    wholeGenomeDuplication,
    minPurity,
    maxPurity,
    minDiploidProportion,
    maxDiploidProportion,
    msIndelsPerMb,
    msStatus
FROM
    pipeline_v5_validation.purity
WHERE
    sampleId = 'COLO829v003T';
