SELECT min(pipeline) as pipeline, sampleId, gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction, source, sourceEvent, evidenceType, rangeRank, count(pipeline)
FROM (
( SELECT 'OnlyInTruth' AS pipeline, SUBSTRING(sampleId,1,13) as sampleId, gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction, source, sourceEvent, evidenceType, rangeRank FROM hmfpatients_pilot.protect where (sampleId in (
'_XXX'
)) and reported = 1 group by sampleId, gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction, source, sourceEvent, evidenceType, rangeRank )
UNION
( SELECT 'OnlyInNew' AS pipeline, SUBSTRING(sampleId,1,13) as sampleId, gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction, source, sourceEvent, evidenceType, rangeRank FROM hmfpatients_pilot.protect where (sampleId in (
'XXX')) and reported = 1 group by sampleId, gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction, source, sourceEvent, evidenceType, rangeRank ) ) as a
group by sampleId, gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction, source, sourceEvent, evidenceType, rangeRank
having count(pipeline) = 1
order by sampleId, gene, event, treatment, pipeline desc;