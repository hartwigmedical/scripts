SELECT DISTINCT
    c.sampleId AS '#sampleId',
    d.startDate,
    d.endDate,
    d.name,
    d.type,
    d.mechanism
FROM
    drug AS d
        LEFT JOIN
    treatment AS t ON d.treatmentId = t.id
        LEFT JOIN
    biopsy AS b ON t.biopsyId = b.id
        LEFT JOIN
    clinical AS c ON b.sampleId = c.sampleId
WHERE
    c.informedConsentDate > '2016-04-20'
AND
    t.treatmentGiven = "Yes"
AND
    not isnull(d.startDate)
ORDER BY
    c.sampleId
;
