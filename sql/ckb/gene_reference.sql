select distinct geneSymbol, geneRole, entrezId, canonicalTranscript, title, url, authors, journal, year from gene
inner join geneReference on geneReference.geneId = gene.id
where geneSymbol = 'XXX' and not isnull(pubMedId);