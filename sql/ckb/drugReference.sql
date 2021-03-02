select distinct drugName, tradeName, casregistryNum, ncitId, pubMedId, title, url, authors, journal, year from drug
inner join drugReference on drugReference.drugId = drug.id
where drugName = 'XXX' and not isnull(pubMedId);