SELECT geneId, transcriptId, a.* FROM
(SELECT * FROM novelSpliceJunction
WHERE type = 'NOVEL_EXON' AND gene IN (SELECT gene FROM driverGenePanel) AND fragmentCount>5 AND cohortFrequency<5
AND sampleId IN ('XXX')) AS a
INNER JOIN canonicalTranscript
ON a.gene=canonicalTranscript.gene;