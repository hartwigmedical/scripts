SELECT DISTINCT
    c.sampleId AS '#sampleId',
    tr.responseDate,
    tr.response,
    t.startDate,
    t.endDate,
    t.name,
    t.type,
    t.mechanism
FROM
    treatmentResponse tr
        LEFT JOIN
    treatment AS t ON tr.treatmentId = t.id
        LEFT JOIN
    biopsy AS b ON t.biopsyId = b.id
        LEFT JOIN
    clinical AS c ON c.sampleId = b.sampleId
WHERE
    c.informedConsentDate > '2016-04-20'
AND
    t.treatmentGiven = "Yes"
AND
    tr.measurementDone = "Yes"
ORDER BY
    c.sampleId, tr.responseDate
;
