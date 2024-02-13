SELECT pipeline, sampleId, gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction, source, sourceEvent, evidenceType, rangeRank
FROM (
( SELECT 'OnlyInTruth' AS pipeline, sampleId, gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction, source, sourceEvent, evidenceType, rangeRank FROM hmfpatients_pilot.protect where (sampleId in (
'XXX'
)) and reported = 1 group by sampleId, gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction, source, sourceEvent, evidenceType, rangeRank )
UNION
( SELECT 'OnlyInNew' AS pipeline, SUBSTRING(sampleId,1,13) as sampleId, gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction, source, sourceEvent, evidenceType, rangeRank FROM hmfpatients_pilot.protect where (sampleId in (
'XXX'
)) and reported = 1 group by sampleId, gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction, source, sourceEvent, evidenceType, rangeRank ) ) as a
order by sampleId, gene, event, treatment, pipeline desc;