select distinct therapyName, pubMedId, title, url, authors, journal, year from therapy
inner join therapyReference on therapyReference.therapyId = therapy.id
where therapyName = 'XXXX' and not isnull(pubMedId);