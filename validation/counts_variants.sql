## Check variant counts
SELECT 
(SELECT COUNT(*) FROM pipeline_v5_validation.somaticVariant WHERE sampleId = "COLO829v003T" AND filter = "PASS") AS NewSom,
(SELECT COUNT(*) FROM reference_validation_sets.somaticVariant WHERE sampleId = "COLO829v003T" AND filter = "PASS") AS TruthSom,
(SELECT COUNT(*) FROM pipeline_v5_validation.structuralVariant WHERE sampleId = "COLO829v003T" and filter <> "PON") AS NewSv,
(SELECT COUNT(*) FROM reference_validation_sets.structuralVariant WHERE sampleId = "COLO829v003T" and filter <> "PON") AS TruthSv,
(SELECT COUNT(*) FROM pipeline_v5_validation.svBreakend WHERE sampleId = "COLO829v003T") AS NewSvbe,
(SELECT COUNT(*) FROM reference_validation_sets.svBreakend WHERE sampleId = "COLO829v003T") AS TruthSvbe;
