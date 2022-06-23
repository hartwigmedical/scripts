CREATE OR REPLACE VIEW studyMutationConditions AS

select idDB, acronym, title, eudra, nct, ipn, ccmo, gene, mutation
from study
inner join mutationConditions
on study.id= mutationConditions.mutationConditionId;