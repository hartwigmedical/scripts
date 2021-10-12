SELECT a.*, group_concat(likelihoodMethod) AS driversInGene, purity
FROM
	(SELECT sampleId,chromosome,position,gene,ref,alt,canonicalCodingEffect,alleleReadCount/totalReadCount AS DNAVaf,rnaAlleleReadCount/rnaTotalReadCount AS
	RNAVaf,alleleReadCount,totalReadCount,rnaAlleleReadCount,rnaTotalReadCount,copyNumber,reported
	FROM hmfpatients_pilot.somaticVariant)
    AS a
	LEFT JOIN driverCatalog d ON a.sampleId=d.sampleId AND a.gene=d.gene
    LEFT JOIN purity p ON a.sampleId=p.sampleId
WHERE a.sampleId IN ('XXX') AND reported
GROUP BY a.sampleId, a.gene
ORDER BY a.sampleId, a.gene;