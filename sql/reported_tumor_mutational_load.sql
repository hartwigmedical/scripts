SELECT sampleId, count(*) AS TML, if (count(*) > 140, "High", "Low") AS status
FROM somaticVariant WHERE filter = "PASS" AND worstCodingEffect = "MISSENSE" AND sampleId IN ('XXX')
GROUP BY 1;