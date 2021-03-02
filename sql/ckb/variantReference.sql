select distinct fullName, impact, pubMedId, title, url, authors, journal, year from variant
inner join variantReference on variantReference.variantId = variant.id
where fullName = 'XXX' and not isnull(pubMedId);