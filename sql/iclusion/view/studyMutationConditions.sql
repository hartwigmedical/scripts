CREATE OR REPLACE VIEW studyMutationCondition AS

select idDB, acronym, title, eudra, nct, ipn, ccmo, gene, mutation
from study
left join mutationCondition
on study.id= mutationCondition.mutationConditionId;