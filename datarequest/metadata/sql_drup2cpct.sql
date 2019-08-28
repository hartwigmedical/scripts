USE hmfpatients;
SELECT patients.patientId,
CPCT2YN, CPCTCNT, CPCTPN
FROM
  (SELECT DISTINCT patientId FROM drupEcrf) patients
LEFT JOIN
  (SELECT patientId, group_concat(itemValue separator ', ') as CPCT2YN
  FROM drupEcrf WHERE item = 'FLD.CTCT2YN' GROUP BY patientId) yesno
on patients.patientId = yesno.patientId
LEFT JOIN
  (SELECT patientId, group_concat(itemValue separator ', ') as CPCTCNT
  FROM drupEcrf WHERE item ='FLD.CPCTCNT' GROUP BY patientId) center
on patients.patientId = center.patientId
LEFT JOIN
  (SELECT patientId, group_concat(itemValue separator ', ') as CPCTPN
  FROM drupEcrf WHERE item ='FLD.CPCTPN' GROUP BY patientId) id
ON patients.patientId = id.patientId
WHERE CPCT2YN = "Yes" AND CHAR_LENGTH(CPCTCNT) = 2 AND CHAR_LENGTH(CPCTPN) = 4
;
