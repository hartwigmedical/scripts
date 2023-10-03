SELECT MIN(pipeline) AS pipeline, sampleId, gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment,
onLabel, level, direction, source, sourceEvent, evidenceType, rangeRank, COUNT(pipeline) AS count
FROM (
(SELECT 'OnlyInTruth' AS pipeline, sampleId, gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment,
onLabel, level, direction, source, sourceEvent, evidenceType, rangeRank
FROM VARIABLE_TRUTH_DB_SCHEMA.protect WHERE (sampleId IN ('VARIABLE_TRUTH_SAMPLE_ID')) AND reported = 1
GROUP BY sampleId, gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment,
onLabel , level , direction , source , sourceEvent , evidenceType , rangeRank)
UNION (SELECT 'OnlyInNew' AS pipeline, SUBSTRING(sampleId, 1, 13) AS sampleId, gene, transcript, isCanonical, event, eventIsHighDriver,
germline, reported, treatment, onLabel, level, direction, source, sourceEvent, evidenceType, rangeRank
FROM VARIABLE_NEW_DB_SCHEMA.protect WHERE (sampleId IN ('VARIABLE_NEW_SAMPLE_ID')) AND reported = 1
GROUP BY sampleId, gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction,
source, sourceEvent, evidenceType, rangeRank)
) AS a
GROUP BY sampleId, gene, transcript, isCanonical, event, eventIsHighDriver, germline, reported, treatment, onLabel, level, direction,
source, sourceEvent, evidenceType, rangeRank
HAVING count != 2
ORDER BY sampleId, gene, event, treatment, pipeline DESC;