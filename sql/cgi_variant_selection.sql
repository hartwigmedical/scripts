SELECT chromosome, position, ref, alt
FROM somaticVariant
WHERE filter = 'PASS' AND (canonicalEffect LIKE '%protein%' OR canonicalEffect LIKE '%structural%' OR canonicalEffect LIKE '%missense%' OR canonicalEffect
LIKE '%frame%' OR canonicalEffect LIKE '%stop%')
AND sampleId = 'XXX';