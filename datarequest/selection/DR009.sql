SELECT
DISTINCT(dr.patientId) AS '#patientId'
FROM
datarequest as dr
  INNER JOIN
  patient AS p on dr.patientId = p.patientIdentifier
  INNER JOIN
  baseline AS b ON b.patientId = p.id
  LEFT JOIN
  drug AS d ON d.patientId = p.id
  INNER JOIN
  treatment AS t ON d.treatmentId = t.id
WHERE
d.type LIKE '%immuno%'
AND
t.treatmentGiven = "Yes"
AND NOT
isnull(t.startDate)
ORDER BY 1;