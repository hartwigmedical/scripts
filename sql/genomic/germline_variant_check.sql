SELECT * FROM germlineVariant2
WHERE (
(gene IN (SELECT gene FROM driverGenePanel WHERE reportGermlineVariant=1 OR reportGermlineHotspot=1) AND pathogenic AND filter <> 'PASS')
OR (gene IN (SELECT gene FROM driverGenePanel WHERE reportGermlineVariant=1 OR reportGermlineHotspot=1) AND pathogenic AND filter = 'PASS' AND NOT reported)
OR (gene IN (SELECT gene FROM driverGenePanel WHERE reportGermlineVariant=1 OR reportGermlineHotspot=1) AND pathogenicity = 'CONFLICTING' AND filter = 'PASS')
OR (sampleId IN (SELECT sampleId FROM germlineVariant2 WHERE gene = 'MUTYH' AND pathogenic GROUP BY sampleId HAVING count(sampleId)>1))
OR (sampleId IN (SELECT sampleId FROM germlineVariant2 WHERE gene = 'NTHL1' AND pathogenic GROUP BY sampleId HAVING count(sampleId)>1)))
AND sampleId IN ('XXX');
