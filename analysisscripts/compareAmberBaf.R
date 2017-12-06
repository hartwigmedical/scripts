#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table) # used to speed up merging of dataframes
options(scipen = 999) # avoids scientific notation

SCRIPT_NAME="AMBERBAFCHECK"
OUTPUT_HEADER = c( "Input1", "Input2", "Sample1", "Sample2", "Count1", "Count2", "OverlapCount", "OverlapProp" )
OUTPUT_DELIM=","
OUTPUT_DECIMALS=3
OUTPUT_FILE="amber_baf_comparisons_table.csv"

## ----------
## MAIN
main <- function(){
    
    ## general input checks
    amberBafFilePaths <- args
    if (length(amberBafFilePaths)<2) {
        stop( "Provide two or more *.amber.baf files: eg <file-1.baf> <file-2.baf> [<file-n.baf>]\n", call.=FALSE)
    }
    if ( file.exists(OUTPUT_FILE) ){ 
        stop( paste( "[EXIT] Output file already exists:", OUTPUT_FILE, sep=" ") ) 
    }
    
    ## preparation
    write.table(data.frame(), file=OUTPUT_FILE, col.names=FALSE)
    
    ## input parsing
    catMsg( c("Parsing input of", length(amberBafFilePaths), "files") )
    bafObjects <- getObjectsFromBafInputFiles( amberBafFilePaths )
    objCount <- length(bafObjects)
    pairCount <- choose(objCount,2)
    pairs <- list()
    
    ## pairwise comparisons
    atComparison=1
    catMsg( c("Performing a total", pairCount, "pairwise comparisons for", objCount, "objects") )
    catMsg( c("For progress check the line count of ouptut file:", OUTPUT_FILE) )
    sink( OUTPUT_FILE, append=FALSE, split=FALSE )
    cat( "#", paste( c( OUTPUT_HEADER), collapse=OUTPUT_DELIM ), "\n", sep="" )
    
    for (i in 1:(objCount)){
        for (j in i:objCount){
            if ( i != j ){
                if ( atComparison %% 10000 == 0 ){ catMsg( c( "At comparison round", atComparison )) }
                fname1 <- bafObjects[[ i ]]$name
                fname2 <- bafObjects[[ j ]]$name
                sname1 <- substr( fname1, 1, 12 )
                sname2 <- substr( fname2, 1, 12 )
                count1 <- nrow( bafObjects[[ i ]]$data )
                count2 <- nrow( bafObjects[[ j ]]$data )
                countB <- nrow( bafObjects[[ i ]]$data[ bafObjects[[ j ]]$data, nomatch=0L, on = c("CHR", "POS") ] )
                 propB <- sprintf( "%.3f", round( countB / min( count1, count2), 3 ) )
                comparison <- list( 
                    "file1"   = fname1,
                    "file2"   = fname2,
                    "sample1" = sname1,
                    "sample2" = sname2,
                    "countFile1"  = count1,
                    "countFile2"  = count2,
                    "OverlapCount"= countB,
                    "OverlapProp" = propB
                )
                cat( paste( c(fname1, fname2, sname1, sname2, count1, count2, countB, propB), collapse=OUTPUT_DELIM), "\n", sep="" )
                atComparison <- atComparison+1
            }
        }
    }
    
    sink()
    catMsg( c("Finished: output in", OUTPUT_FILE) )
    
}
## /MAIN
## ----------

catMsg <- function( msg ){    
   cat( SCRIPT_NAME, "INFO", date(), msg, "\n", sep=" " ) 
}

getObjectsFromBafInputFiles <- function( filePaths ){
    
    fileCount <- length(filePaths)
    objects <- list()
    for (idx in 1:fileCount){
        filePath <- filePaths[idx]
        fileBase <- basename( filePath )
        datatable <- read.csv( 
            file=filePath, 
            sep="\t", 
            colClasses=c(NA, NA, "NULL", "NULL") # this ignores column 3 and 4
        )
        colnames(datatable) <- c( "CHR", "POS" )
        object <- list( name = fileBase, data = as.data.table(datatable) )
        objects[[ idx ]] <- object
    }
    
    return(objects)
}

main()



