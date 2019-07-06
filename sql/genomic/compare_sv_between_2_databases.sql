SELECT min(pipeline), startChromosome, startPosition, endChromosome, endPosition,TYPE,endPosition-startPosition AS len, count(*),max(ploidy),max(qualScore) FROM
(SELECT id, startChromosome, startPosition, endChromosome, endPosition,TYPE, 'Pv4' AS pipeline,ploidy,qualScore FROM reference_validation_sets.structuralVariant WHERE sampleId = 'XXX' AND filter <> "PON"
UNION
SELECT id, startChromosome, startPosition, endChromosome, endPosition,TYPE, 'Pv5' AS pipeline,ploidy,qualScore FROM pipeline_v5.structuralVariant WHERE sampleId = 'XXX' AND filter <> "PON") AS a
GROUP BY 2,3,4,5,6,7 HAVING count(*) =1;