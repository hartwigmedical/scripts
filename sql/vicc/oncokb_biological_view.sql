CREATE OR REPLACE VIEW oncokbBiologicalView AS

SELECT
    viccEntryId, gene,
    oncokbVariantBiological.name AS variantName, alteration, proteinStart, proteinEnd, refResidues, variantResidues,
    term, description, isGenerallyTruncating,
    oncogenic, mutationEffect, mutationEffectPmids, mutationEffectAbstracts,
    hugoSymbol, GROUP_CONCAT(geneAlias SEPARATOR ","),
    oncokbGeneBiological.name AS geneName, oncokbGeneBiological.entrezGeneId AS geneEntrezGeneId,
    oncokbBiological.entrezGeneId AS variantEntrezGeneId, isoform, refseq,
    curatedIsoform, curatedRefSeq, oncogene, tsg
FROM oncokbBiological
INNER JOIN oncokb ON oncokbBiological.oncokbId = oncokb.id
INNER JOIN oncokbVariantBiological ON oncokbBiological.id = oncokbVariantBiological.oncokbBiologicalId
INNER JOIN oncokbConsequenceBiological ON oncokbVariantBiological.id = oncokbConsequenceBiological.oncokbVariantBiologicalId
INNER JOIN oncokbGeneBiological ON oncokbVariantBiological.id = oncokbGeneBiological.oncokbVariantBiologicalId
LEFT JOIN oncokbGeneAliasBiological ON oncokbGeneBiological.id = oncokbGeneAliasBiological.oncokbGeneBiologicalId
GROUP BY 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24,25,26;