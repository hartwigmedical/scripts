SELECT sampleId, count(*)/2859 AS TMB, if (count(*)/2859 > 10, "High", "Low") AS status
FROM somaticVariant WHERE filter = 'PASS' AND sampleId IN ('XXX')
GROUP BY 1;
