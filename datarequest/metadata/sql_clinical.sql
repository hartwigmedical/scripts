USE hmfpatients;

SELECT
   purity.sampleId AS 'puritySampleId',
   purity.qcStatus AS 'purityQC',
   purity.status AS 'purityStatus',
   ROUND(purity.purity, 2) AS 'tumorPurity',
   c.*,
   allPostTreatments,
   allPostTreatmentTypes
FROM
   purity
   INNER JOIN
      c
      ON purity.sampleId = c.sampleId
   LEFT JOIN
      patient p
      ON c.patientId = p.patientIdentifier
   LEFT JOIN
      (
         SELECT
            patientid,
            Group_concat(t.name SEPARATOR '/') AS allPostTreatments,
            Group_concat(t.type SEPARATOR '/') AS allPostTreatmentTypes
         FROM
            treatment t
         GROUP BY
            patientId
      )
      treatments
      ON p.id = treatments.patientId
;
