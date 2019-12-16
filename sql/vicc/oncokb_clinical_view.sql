CREATE OR REPLACE VIEW oncokbClinicalView AS

SELECT
    viccEntryId, gene,
    oncokbVariantClinical.name AS variantName, alteration, proteinStart, proteinEnd, refResidues, variantResidues,
    term, description, isGenerallyTruncating,
    cancerType, drug, drugPmids, level, levelLabel,
    GROUP_CONCAT(text SEPARATOR ",") AS drugTexts, GROUP_CONCAT(link SEPARATOR ",") AS drugLinks,
    hugoSymbol, GROUP_CONCAT(geneAlias SEPARATOR ",") AS geneAliases,
    oncokbGeneClinical.name AS geneName, oncokbGeneClinical.entrezGeneId AS geneEntrezGeneId,
    oncokbClinical.entrezGeneId AS variantEntrezGeneId, isoform, refseq,
    curatedIsoform, curatedRefSeq, oncogene, tsg
FROM oncokbClinical
INNER JOIN oncokb ON oncokbClinical.oncokbId = oncokb.id
INNER JOIN oncokbVariantClinical ON oncokbClinical.id = oncokbVariantClinical.oncokbClinicalId
INNER JOIN oncokbConsequenceClinical ON oncokbVariantClinical.id = oncokbConsequenceClinical.oncokbVariantClinicalId
INNER JOIN oncokbGeneClinical ON oncokbVariantClinical.id = oncokbGeneClinical.oncokbVariantClinicalId
LEFT JOIN oncokbGeneAliasClinical ON oncokbGeneClinical.id = oncokbGeneAliasClinical.oncokbGeneClinicalId
LEFT JOIN oncokbDrugAbstractClinical ON oncokbClinical.id = oncokbDrugAbstractClinical.oncokbClinicalId
GROUP BY 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,19,21,22,23,24,25,26,27,28,29;