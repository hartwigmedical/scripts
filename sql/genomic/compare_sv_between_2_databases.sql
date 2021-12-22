SELECT min(pipeline), startChromosome, startPosition, endChromosome, endPosition,TYPE, endPosition-startPosition AS len, count(*), max(junctionCopyNumber), max(qualScore)
FROM
	(SELECT id, 'Truth' AS pipeline, startChromosome, startPosition, endChromosome, endPosition, TYPE, junctionCopyNumber, qualScore
	FROM reference_validation_sets.structuralVariant WHERE sampleId = 'XXX' AND filter <> "PON"
	UNION
	SELECT id, 'New' AS pipeline, startChromosome, startPosition, endChromosome, endPosition, TYPE, junctionCopyNumber, qualScore
	FROM pipeline_v5_validation.structuralVariant WHERE sampleId = 'XXX' AND filter <> "PON") AS a
GROUP BY 2,3,4,5,6,7 HAVING count(*) = 1;