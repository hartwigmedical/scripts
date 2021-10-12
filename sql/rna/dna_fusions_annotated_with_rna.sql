SELECT * FROM svFusion s
LEFT JOIN rnaFusion f ON s.sampleId=f.sampleId AND s.name=f.name
WHERE s.sampleId IN ('XXX') AND reported;