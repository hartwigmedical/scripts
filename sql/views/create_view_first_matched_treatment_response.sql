CREATE OR REPLACE view firstMatchedTreatmentResponse AS
SELECT *
FROM treatmentResponse tr
WHERE tr.treatmentId IS NOT NULL
  AND NOT (exists(SELECT 1
                  FROM treatmentResponse tr1
                  WHERE tr1.treatmentId = tr.treatmentId
                    AND tr1.responseDate <= tr.responseDate
                    AND tr1.id <> tr.id));