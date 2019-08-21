#!/usr/bin/Rscript
library(getopt)
suppressPackageStartupMessages(library(RMySQL))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(jsonlite))

DB_NAME = 'hmfpatients'
DB_GROUP = 'RAnalysis'
SAMPLE_TO_SET_MAPPING_FILE = '/data/data_archive/datarequests/procedure/sampleid2setname.tsv'

main <- function(configFile) {
    
    cat( sprintf("[INFO] Using config (%s) and database (%s)\n", configFile, DB_NAME) )
    config <- parseConfig(configFile)
    id2set <- parseMappingTsv(SAMPLE_TO_SET_MAPPING_FILE)
    
    db = dbConnect(MySQL(), dbname = DB_NAME, groups = DB_GROUP)
    samplesRaw = queryDatarequestBase(db)
    samplesFlt = queryDatarequestFiltered(db, config)
    dbDisconnect(db)
    
    catInfo( sprintf("Adding Set names by sampleId") )
    samplesRaw <- addSetNames(samplesRaw, id2set)
    samplesFlt <- addSetNames(samplesFlt, id2set)

    outPathJsnRaw = "./metadata_raw.json"
    outPathTsvRaw = "./metadata_raw.tsv"
    outPathJsnFlt = "./metadata_filtered.json"
    outPathTsvFlt = "./metadata_filtered.tsv"

    writeOutputFile(samplesRaw, outPathJsnRaw, "json")
    writeOutputFile(samplesRaw, outPathTsvRaw, "tsv")
    
    writeOutputFile(samplesFlt, outPathJsnFlt, "json")
    writeOutputFile(samplesFlt, outPathTsvFlt, "tsv")    
    
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
    catInfo( "Data inspection:" )
    print(head(data))
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
    cat( sprintf("[INFO] Written output file of type %s with %s objects (%s)\n", fileType, objNrow, filePath) )
}

queryDatabase <- function(db, query) {
    raw_data = dbGetQuery(db, query)
    DT = data.table(raw_data, key="sampleId")
    return(DT)
}

queryDatarequestBase <- function(db) {
    query = paste0(
        "SELECT datarequestBase.*, metric.refMeanCoverage, metric.tumorMeanCoverage",
        " FROM datarequestBase LEFT JOIN metric ON datarequestBase.sampleId = metric.sampleId"
    )
    DT = queryDatabase(db, query)
    return(DT)
}

queryDatarequestFiltered <- function(db, config) {
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
    catInfo( sprintf("Final query (filtered): %s", query) )
        
    DT = queryDatabase(db, query)
    return(DT)
}

catInfo <- function( msg ) {
    cat(sprintf("[INFO] %s\n", msg))
}
catError <- function( msg ) {
    stop(sprintf("[EXIT] %s\n", msg))
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


