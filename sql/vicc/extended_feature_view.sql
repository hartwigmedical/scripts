select feature.*, sequenceOntology.*, hierarchyNames, provenanceNames, synonymNames, linkNames from feature
left join sequenceOntology on sequenceOntology.featureId = feature.id
left join
	(select sequenceOntologyId, GROUP_CONCAT(hierarchyName SEPARATOR ",") as hierarchyNames from hierarchy group by 1)
	hierarchies on hierarchies.sequenceOntologyId = sequenceOntology.id
left join
	(select featureId, GROUP_CONCAT(provenanceName SEPARATOR ",") as provenanceNames from provenance group by 1)
	provenances on provenances.featureId = feature.id
left join
	(select featureId, GROUP_CONCAT(synonymName SEPARATOR ",") as synonymNames from synonym group by 1)
	synonyms on synonyms.featureId = feature.id
left join
	(select featureId, GROUP_CONCAT(linkName SEPARATOR ",") as linkNames from link group by 1)
	links on links.featureId = feature.id;