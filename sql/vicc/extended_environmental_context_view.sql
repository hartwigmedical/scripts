select environmentalContext.*, taxonomy.*, approvedCountries from environmentalContext
left join taxonomy on environmentalContext.id = taxonomy.environmentContextId
left join
	(select environmentContextId, GROUP_CONCAT(approvedCountryName SEPARATOR ",") as approvedCountries from approvedCountry group by 1)
	approvedCountries on approvedCountries.environmentContextId = environmentalContext.id;