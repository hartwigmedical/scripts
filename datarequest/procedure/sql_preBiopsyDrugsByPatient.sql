SELECT DISTINCT
    p.patientIdentifier AS '#patientId',
    d.startDate,
    d.endDate,
    d.name,
    d.type,
    d.mechanism
FROM
    preTreatmentDrug AS d
        LEFT JOIN
    patient AS p ON d.patientId = p.id
        LEFT JOIN
    clinical AS c ON c.patientId = p.patientIdentifier
WHERE
    c.informedConsentDate > '2016-04-20'
AND
    not isnull(d.startDate)
ORDER BY
    p.patientIdentifier
;