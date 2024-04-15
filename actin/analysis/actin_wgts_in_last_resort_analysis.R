library(dplyr)
library(tibble)

# Retrieve data ------------------------------------------------------------------
dbActin <- dbConnect(MySQL(), dbname='actin_pilot', groups="RAnalysis")

query_sample <-"select * from actinLastResortPaperSamples;"
query_molecular <-"select * from molecular where sampleId in (select sampleId from actinLastResortPaperSamples);"
query_tumor <-"select * from tumor where patientId in (select patientId from actinLastResortPaperSamples);"

sample <- dbGetQuery(dbActin, query_sample)
molecular <- dbGetQuery(dbActin, query_molecular)
tumor <- dbGetQuery(dbActin, query_tumor)

dbDisconnect(dbActin)

