library(RMySQL)
library(dndscv)
library(IRanges)

query_clinical_data<-function(dbConnect) {
  query = paste(
    "SELECT c.sampleId, c.cancertype, c.birthYear, c.biopsyDate",
    " FROM clinical c",
    sep = " ")
  return ((dbGetQuery(dbConnect, query)))
}

query_somatics<-function(dbConnect) {
  query = paste(
    "SELECT s.sampleId, chromosome as chr, position as pos, ref, alt ",
    "  FROM somaticVariant s, purity p ",
    " WHERE s.sampleId = p.sampleId",
    "   AND status <> 'NO_TUMOR' AND qcstatus = 'PASS' ",
    "   AND right(s.sampleId,1) ='T'",
    "   AND p.modified > '2017-12-21'",
    "   AND filter = 'PASS'",
    "   AND type = 'SNP'",
    "   AND gene <> ''",
    # PETE - COMMENT OUT NEXT LINE FOR LOTS OF FUN :)
    "   AND s.sampleId in ('CPCT02020213T','CPCT02030255T', 'CPCT02070012T')",
    sep = "")
  
  return (dbGetQuery(dbConnect, query)) 
}

#Select clinical data
dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
clinicalData = query_clinical_data(dbProd)
dbDisconnect(dbProd)
rm(dbProd)

#Select somatics
dbPilot = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
somatics = query_somatics(dbPilot)
dbDisconnect(dbPilot)
rm(dbPilot)

# Attach cancertype to somatics
somatics$cancerType <- sapply(somatics$sampleId, function(x) {clinicalData[match(x, clinicalData$sampleId), c("cancerType")] })

# Which cancer types to we have?
unique(somatics$cancerType)

allInput = somatics[, c("sampleId", "chr", "pos", "ref", "alt")]
allOutput<-dndscv(allInput)

sarcomaInput = somatics[somatics$cancerType == 'Sarcoma', c("sampleId", "chr", "pos", "ref", "alt")];
sarcomaOutput<-dndscv(sarcomaInput)

prostateInput = somatics[somatics$cancerType == 'Prostate', c("sampleId", "chr", "pos", "ref", "alt")];
prostateOutput<-dndscv(prostateInput)

str(allOutput)
allOutput$globaldnds
allOutput$sel_cv
