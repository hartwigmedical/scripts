SELECT patients.patientId,
CTCT2YN, CPCTCNT, CPCTPN, WIDEPNR
FROM
	(SELECT DISTINCT patientId FROM drupEcrf) patients
LEFT JOIN
   (SELECT patientId, group_concat(itemValue separator ', ') AS CTCT2YN
    FROM drupEcrf WHERE item = 'FLD.CTCT2YN' GROUP BY patientId) yesno
ON patients.patientId = yesno.patientId
LEFT JOIN
    (SELECT patientId, group_concat(itemValue separator ', ') AS CPCTCNT
     FROM drupEcrf WHERE item ='FLD.CPCTCNT' GROUP BY patientId) center
ON patients.patientId = center.patientId
LEFT JOIN
    (SELECT patientId, group_concat(itemValue separator ', ') AS CPCTPN
     FROM drupEcrf WHERE item ='FLD.CPCTPN' GROUP BY patientId) id
ON patients.patientId = id.patientId
LEFT JOIN
    (SELECT patientId, group_concat(itemValue separator ', ') AS WIDEPNR
     FROM drupEcrf WHERE item ='FLD.WIDEPNR' GROUP BY patientId) wideid
ON patients.patientId = wideid.patientId