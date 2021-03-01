select * from germlineVariant
where sampleId in("XXX") and reported and (sampleId, gene, hgvsCoding)
not in (select sampleId, gene, canonicalHgvsCodingImpact from germlineVariant2 where sampleId in("XXX") and reported);

select * from germlineVariant2
where sampleId in("XXX") and reported and (sampleId, gene, canonicalHgvsCodingImpact)
not in (select sampleId, gene, hgvsCoding from germlineVariant where sampleId in("XXX") and reported);