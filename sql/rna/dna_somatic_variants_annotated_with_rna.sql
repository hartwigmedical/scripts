SELECT a.*, driversInGene, purity
FROM
	(SELECT filter,sampleId,chromosome,position,gene,ref,alt,canonicalCodingEffect,alleleReadCount/totalReadCount AS DNAVaf,rnaAlleleReadCount/rnaTotalReadCount AS
	RNAVaf,alleleReadCount,totalReadCount,rnaAlleleReadCount,rnaTotalReadCount,round(variantCopyNumber,1) AS variantCopyNumberDna,copyNumber,reported
	FROM somaticVariant)
    AS a
    LEFT JOIN (SELECT sampleId, gene, group_concat(likelihoodMethod) AS driversInGene FROM driverCatalog d GROUP BY 1,2) AS b
		ON a.sampleId=b.sampleId AND a.gene=b.gene
    LEFT JOIN purity p
		ON a.sampleId=p.sampleId
WHERE a.sampleId IN ('XXX') AND (reported OR (canonicalCodingEffect in ('SPLICE','MISSENSE','NONSENSE_OR_FRAMESHIFT') AND a.gene IN (SELECT gene FROM driverGenePanel) AND a.filter = 'PASS'))
ORDER BY a.sampleId, a.gene;