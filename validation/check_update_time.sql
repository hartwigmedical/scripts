SELECT TABLE_SCHEMA, UPDATE_TIME, TABLE_NAME, TABLE_ROWS
FROM   information_schema.tables
WHERE  TABLE_SCHEMA = 'VARIABLE_NEW_DB_SCHEMA'
ORDER by UPDATE_TIME DESC;
