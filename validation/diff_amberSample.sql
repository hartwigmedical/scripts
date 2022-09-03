SELECT
    MIN(pipeline) AS pipeline,
    site1, site2, site3, site4, site5, site6, site7, site8, site9, site10, site11, site12, site13, site14, site15, site16, site17, site18, site19, site20, site21, site22, site23, site24, site25, site26, site27, site28, site29, site30, site31, site32, site33, site34, site35, site36, site37, site38, site39, site40, site41, site42, site43, site44, site45, site46, site47, site48, site49, site50, site51, site52, site53, site54, site55, site56, site57, site58, site59, site60, site61, site62, site63, site64, site65, site66, site67, site68, site69, site70, site71, site72, site73, site74, site75, site76, site77, site78, site79, site80, site81, site82, site83, site84, site85, site86, site87, site88, site89, site90, site91, site92, site93, site94, site95, site96, site97, site98, site99, site100
FROM
    (SELECT
        'OnlyInTruth' AS pipeline,
        site1, site2, site3, site4, site5, site6, site7, site8, site9, site10, site11, site12, site13, site14, site15, site16, site17, site18, site19, site20, site21, site22, site23, site24, site25, site26, site27, site28, site29, site30, site31, site32, site33, site34, site35, site36, site37, site38, site39, site40, site41, site42, site43, site44, site45, site46, site47, site48, site49, site50, site51, site52, site53, site54, site55, site56, site57, site58, site59, site60, site61, site62, site63, site64, site65, site66, site67, site68, site69, site70, site71, site72, site73, site74, site75, site76, site77, site78, site79, site80, site81, site82, site83, site84, site85, site86, site87, site88, site89, site90, site91, site92, site93, site94, site95, site96, site97, site98, site99, site100
    FROM
        VARIABLE_TRUTH_DB_SCHEMA.amberSample
    WHERE
        sampleId = 'VARIABLE_TRUTH_SAMPLE_ID' UNION SELECT
        'OnlyInNew' AS pipeline,
        site1, site2, site3, site4, site5, site6, site7, site8, site9, site10, site11, site12, site13, site14, site15, site16, site17, site18, site19, site20, site21, site22, site23, site24, site25, site26, site27, site28, site29, site30, site31, site32, site33, site34, site35, site36, site37, site38, site39, site40, site41, site42, site43, site44, site45, site46, site47, site48, site49, site50, site51, site52, site53, site54, site55, site56, site57, site58, site59, site60, site61, site62, site63, site64, site65, site66, site67, site68, site69, site70, site71, site72, site73, site74, site75, site76, site77, site78, site79, site80, site81, site82, site83, site84, site85, site86, site87, site88, site89, site90, site91, site92, site93, site94, site95, site96, site97, site98, site99, site100
    FROM
        VARIABLE_NEW_DB_SCHEMA.amberSample
    WHERE
        sampleId = 'VARIABLE_NEW_SAMPLE_ID') AS a
GROUP BY site1, site1, site2, site3, site4, site5, site6, site7, site8, site9, site10, site11, site12, site13, site14, site15, site16, site17, site18, site19, site20, site21, site22, site23, site24, site25, site26, site27, site28, site29, site30, site31, site32, site33, site34, site35, site36, site37, site38, site39, site40, site41, site42, site43, site44, site45, site46, site47, site48, site49, site50, site51, site52, site53, site54, site55, site56, site57, site58, site59, site60, site61, site62, site63, site64, site65, site66, site67, site68, site69, site70, site71, site72, site73, site74, site75, site76, site77, site78, site79, site80, site81, site82, site83, site84, site85, site86, site87, site88, site89, site90, site91, site92, site93, site94, site95, site96, site97, site98, site99, site100
HAVING COUNT(*) != 2;