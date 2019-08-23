#!/usr/bin/Rscript
library(getopt)
suppressPackageStartupMessages(library(RMySQL))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(jsonlite))

DB_NAME = 'hmfpatients'
DB_GROUP = 'RAnalysis'
SAMPLE_TO_SET_MAPPING_FILE = '/data/data_archive/datarequests/procedure/files/sampleid2setname.tsv'
LOG_FILE = "./log.txt"

main <- function(configFile) {
    
    dateTime <- Sys.time()
    write( sprintf("Date: %s", dateTime), file=LOG_FILE, append=FALSE)

    catInfo( sprintf("Using config (%s) and database (%s)", configFile, DB_NAME) )
    config <- parseConfig(configFile)
    id2set <- parseMappingTsv(SAMPLE_TO_SET_MAPPING_FILE)
    
    db = dbConnect(MySQL(), dbname = DB_NAME, groups = DB_GROUP)
    samplesRaw = queryDatarequestBase(db, logFile)
    samplesFlt = queryDatarequestFiltered(db, config, logFile)
    
    catInfo( sprintf("Adding Set names by sampleId") )
    samplesRaw <- addSetNames(samplesRaw, id2set)
    samplesFlt <- addSetNames(samplesFlt, id2set)    
    sIds <- samplesFlt$sampleId
    pIds <- samplesFlt$patientId

    catInfo( sprintf("Retrieving PostBiopsyTreatment responses") )
    postResp <- queryPostBiopsyResponsesBySample(db, sIds, logFile)
    postDrug <- queryPostBiopsyDrugsBySample(db, sIds, logFile)

    catInfo( sprintf("Retrieving PreBiopsyTreatment drugs") )
    preDrug <- queryPreBiopsyDrugsByPatient(db, pIds, logFile)

    ## cleanup
    dbDisconnect(db)
    
    ## output to files
    writeOutputFile(samplesRaw, "./metadata_raw.json", "json")
    writeOutputFile(samplesRaw, "./metadata_raw.tsv", "tsv")
    writeOutputFile(samplesFlt, "./metadata_filtered.json", "json")
    writeOutputFile(samplesFlt, "./metadata_filtered.tsv", "tsv")
    writeOutputFile(  postResp, "./postBiopsyTreatmentResponses.json", "json")
    writeOutputFile(  postResp, "./postBiopsyTreatmentResponses.tsv", "tsv")
    writeOutputFile(  postDrug, "./postBiopsyTreatmentDrugs.json", "json")
    writeOutputFile(  postDrug, "./postBiopsyTreatmentDrugs.tsv", "tsv")
    writeOutputFile(   preDrug, "./preBiopsyTreatmentDrugs.json", "json")
    writeOutputFile(   preDrug, "./preBiopsyTreatmentDrugs.tsv", "tsv")
    catInfo(sprintf("Done. See log for exact queries used (%s)", LOG_FILE))
    
}

addSetNames <- function( metadataTable, id2set ) {
    setkey(metadataTable,sampleId)
    setkey(id2set,sampleId)
    out = merge( metadataTable, id2set, all.x = TRUE )
    return(out)
}

parseMappingTsv <- function( filePath ) {
    catInfo( sprintf("Parsing file (%s)", filePath) )
    data <- fread(filePath, header=FALSE)
    colnames(data) <- c("sampleId", "setName")
    return(data)
}

parseConfig <- function( filePath ) {
    catInfo( sprintf("Parsing config (%s)", filePath) )
    config <- read_json(filePath, simplifyVector = TRUE)
    return(config)
}

createSqlConditionsStringFromFilters <- function( filtersList ) {
    catInfo( "Parsing Filters" )
    conditions <- character()
    for (field in names(filtersList)){
        content <- filtersList[[field]]
        count <- length(content)
        if (count < 1) {
            catError(sprintf("Filter field found (%s) but no content?", field))
        } else if (count == 1) { 
            if (field == "study") {
                ## special field 
                statement <- sprintf( "%s LIKE \"%s%s\"", "sampleId", content, "%" )
            } else {
                statement <- sprintf( "%s = \"%s\"", field, content )
            }
            catInfo( sprintf("  Filter found: %s", statement))
            conditions <- append(conditions, statement)
        } else {
            list_string <- paste( '"', content, '"', sep="", collapse=",")
            statement <- sprintf( "%s IN (%s)", field, list_string)
            catInfo( sprintf("  Filter found: %s", statement))
            conditions <- append(conditions, statement)
        }
    }
    final_conditions_string <- ""
    if ( length(conditions) > 0 ){
        final_conditions_string <- paste0( " WHERE ", paste( conditions, collapse=" AND "))
    }
    return(final_conditions_string)
}

writeOutputFile <- function( obj, filePath, fileType ) {
    objNrow <- nrow(obj)
    if ( fileType == "tsv" ){
        originalColname = colnames(obj)[1]
        colnames(obj)[1] <- paste0("#",originalColname)
        write.table(obj, file=filePath, quote=FALSE, sep='\t', col.names=T, row.names=F)
        colnames(obj)[1] <- originalColname
    } else if ( fileType == "json" ){
        write_json(obj, path=filePath, pretty=TRUE, na="null")
    } else {
        stop( sprintf("[EXIT] Unknown output file type (%s)", fileType) )
    }
    cat( sprintf("[INFO] Written output file with %s objects of type %s (%s)\n", objNrow, fileType, filePath) )
}

queryDatabase <- function(db, query) {
    raw_data = dbGetQuery(db, query)
    first_column = colnames(raw_data)[1]
    DT = data.table(raw_data, key=first_column)
    return(DT)
}

queryDatarequestBase <- function(db, logFile) {
    query = paste0(
        "SELECT datarequestBase.*, metric.refMeanCoverage, metric.tumorMeanCoverage",
        " FROM datarequestBase LEFT JOIN metric ON datarequestBase.sampleId = metric.sampleId"
    )
    log(title="queryDatarequestBase", msg=query)
    DT = queryDatabase(db, query)
    return(DT)
}

queryDatarequestFiltered <- function(db, config, logFile) {
    conditional_string <- createSqlConditionsStringFromFilters(config$FILTERS)
    ## sampleId column name present multiple times due to JOIN
    conditional_string <- gsub( " sampleId IN", " datarequestFiltered.sampleId IN", conditional_string)
    conditional_string <- gsub( " sampleId =", " datarequestFiltered.sampleId =", conditional_string)
    conditional_string <- gsub( " sampleId LIKE", " datarequestFiltered.sampleId LIKE", conditional_string)
    query = paste0(
        "SELECT datarequestFiltered.*, metric.refMeanCoverage, metric.tumorMeanCoverage",
        " FROM datarequestFiltered LEFT JOIN metric ON datarequestFiltered.sampleId = metric.sampleId",
        conditional_string
    )
    log(title="queryDatarequestFiltered", msg=query)
        
    DT = queryDatabase(db, query)
    return(DT)
}

queryPostBiopsyResponsesBySample <- function(dbConnect, sampleIds, logFile) {
    inClauseString <- paste0( paste0( '"', sampleIds, '"'), collapse=',' )
    query = paste(
    'SELECT DISTINCT',
        'c.sampleId AS "sampleId",',
        'tr.responseDate, tr.response,',
        't.startDate, t.endDate, t.name, t.type, t.mechanism',
    'FROM',
    'treatmentResponse AS tr',
        'LEFT JOIN treatment AS t ON tr.treatmentId = t.id',
        'LEFT JOIN biopsy    AS b ON t.biopsyId = b.id',
        'LEFT JOIN clinical  AS c ON c.sampleId = b.sampleId',
    'WHERE',
        'c.informedConsentDate > "2016-04-20"',
        'AND t.treatmentGiven = "Yes"',
        'AND tr.measurementDone = "Yes"',
        'AND c.sampleId IN (', inClauseString, ')',
    'ORDER BY',
        'c.sampleId, tr.responseDate',
    sep = ' ')
    log( title="FinalPostResponsesBySampleQuery", msg=query)
    return(queryDatabase(dbConnect, query))
}

queryPostBiopsyDrugsBySample <- function(dbConnect, sampleIds, logFile) {
    inClauseString <- paste0( paste0( '"', sampleIds, '"'), collapse=',' )
    query = paste(
    'SELECT DISTINCT c.sampleId AS "sampleId", d.startDate, d.endDate, d.name, d.type, d.mechanism',
    'FROM drug AS d',
        'LEFT JOIN treatment AS t ON d.treatmentId = t.id',
        'LEFT JOIN biopsy    AS b ON t.biopsyId    = b.id',
        'LEFT JOIN clinical  AS c ON c.sampleId    = b.sampleId',
    'WHERE c.informedConsentDate > "2016-04-20"',
        'AND t.treatmentGiven = "Yes"',
        'AND not isnull(d.startDate)',
        'AND c.sampleId IN (', inClauseString, ')',
    'ORDER BY c.sampleId',
    sep = ' ')
    log( title="FinalPostDrugsByPatientQuery", msg=query)
    return(queryDatabase(dbConnect, query))
}

queryPreBiopsyDrugsByPatient <- function(dbConnect, patientIds, logFile) {
    inClauseString <- paste0( paste0( '"', patientIds, '"'), collapse=',' )
    query = paste( 'SELECT DISTINCT p.patientIdentifier AS "patientId", d.startDate, d.endDate, d.name, d.type, d.mechanism',
    'FROM preTreatmentDrug AS d',
        'LEFT JOIN patient  AS p ON d.patientId = p.id',
        'LEFT JOIN clinical AS c ON c.patientId = p.patientIdentifier',
    'WHERE c.informedConsentDate > "2016-04-20"',
        'AND not isnull(d.startDate)',
        'AND c.patientId IN (', inClauseString, ')',
    'ORDER BY p.patientIdentifier',
    sep = ' ')
    log( title="FinalPreDrugsByPatientQuery", msg=query)
    return(queryDatabase(dbConnect, query))
}

catInfo <- function( msg ) {
    cat(sprintf("[INFO] %s\n", msg))
}
catError <- function( msg ) {
    stop(sprintf("[EXIT] %s\n", msg))
}

log <- function( title, msg ){
    write( 
        sprintf( "-----\n%s\n%s", title, msg), 
        LOG_FILE, 
        append=TRUE 
    )
}

spec = matrix(c(
    'help',    'h', 0, "logical",
    'config',  'c', 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

if ( length(commandArgs(TRUE)) < 1 || !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE))
    cat( "---\n" )
    cat( "  Description:\n" )
    cat( "    Creates data request metadata tsv/json from database.\n" )
    cat( "  Usage: --config /path/to/configuration.json\n" )
    cat( "    (see exampleConfig*.json for examples)\n" )
    cat( "---\n" )
    q(status=1)
}

main(opt$config)



