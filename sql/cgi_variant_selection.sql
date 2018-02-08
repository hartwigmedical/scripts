SELECT chromosome, position, ref, alt
FROM somaticVariant
WHERE filter = 'PASS' AND (effect LIKE '%protein%' OR effect LIKE '%structural%' OR effect LIKE '%missense%' OR effect LIKE '%frame%' OR effect LIKE '%stop%')
AND sampleId = 'XXX';