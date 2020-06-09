SELECT DISTINCT patientId as '#patientId'
FROM (
  SELECT DISTINCT sampleId
  FROM datarequest
  WHERE hasRNA=1
  UNION
  select distinct(svBreakend.sampleId) from svBreakend
  where gene = 'FGFR2') as a
  inner join datarequest on a.sampleId = datarequest.sampleId;