SELECT geneId, transcriptId, a.* FROM
(SELECT * FROM novelSpliceJunction
WHERE type IN ('NOVEL_INTRON','NOVEL_EXON') AND gene IN (SELECT gene FROM driverGenePanel) AND fragmentCount>5 AND cohortFrequency<10
AND sampleId IN ('XXX')) AS a
INNER JOIN canonicalTranscript
ON a.gene=canonicalTranscript.gene;