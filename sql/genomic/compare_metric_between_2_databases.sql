SELECT "Truth" AS pipeline, sampleId, refMeanCoverage, refCoverage10xPercentage, refCoverage20xPercentage, tumorMeanCoverage, tumorCoverage30xPercentage, tumorCoverage60xPercentage
FROM reference_validation_sets.metric WHERE sampleId = 'XXX'
UNION
SELECT "New" AS pipeline,  sampleId, refMeanCoverage, refCoverage10xPercentage, refCoverage20xPercentage, tumorMeanCoverage, tumorCoverage30xPercentage, tumorCoverage60xPercentage
FROM pipeline_v5_validation.metric WHERE sampleId = 'XXX';