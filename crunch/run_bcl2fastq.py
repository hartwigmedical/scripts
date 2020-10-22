#!/usr/bin/env python

import fnmatch
import os

DATA_DIRS=[
    '/data1/illumina_data/',
    '/data1/illumina_data/TestRuns/'
]

NEXTSEQ_IDS=[
    'NB500901',
    'NB500902'
]

ISEQ_IDS=[
    'FS10001173'
]

BCL2FASTQ_TOOL='/data/fastqconversion/bcl2fastq_v2.20.0.422/bin/bcl2fastq'
METRICS_TOOL='/data/repos/scripts/crunch/check_bcl2fastq_conversion.pl'
TMP_RUNNING_FILE='/tmp/bcl2fastq_running'

## DEFS
def ConvertBcl(run_dir):

    LOG_FILE  = run_dir + '/conversionLog.txt'
    ERR_FILE  = run_dir + '/conversionError.txt'
    STRT_FILE = run_dir + '/conversionStarted.txt'
    DONE_FILE = run_dir + '/conversionDone.txt'
    FQDEL_FILE = run_dir + '/conversionFastqDeleted.txt'

    out_dir = run_dir + '/Fastq'
    METRICS_FILE = out_dir + '/Stats/conversion_metrics_table.txt'

    FLOWCELL_ID = os.path.split(run_dir)[-1].split("_")[-1]
    print("[INFO]   Flowcell ID: " + FLOWCELL_ID)
    print("[INFO]   Output directory: " + out_dir)
 
    os.system( 'touch ' + STRT_FILE )
    os.system( 'echo "[INFO] Started conversion of ' + FLOWCELL_ID + ' ($(date))" >>' + LOG_FILE )
    os.system( BCL2FASTQ_TOOL + " --runfolder-dir %s --output-dir %s -l WARNING 1>>%s 2>>%s" % (run_dir, out_dir, LOG_FILE, ERR_FILE) )

    ## Add flowcell id to fastq file names
    os.system( 'echo "[INFO] Output directory set to ' + out_dir + '" >>' + LOG_FILE )
    os.system( 'echo "[INFO] Adding flowcell ID to output FASTQ filenames" >>' + LOG_FILE )
    for root, dirs, files in os.walk(run_dir+"/Data/Intensities/BaseCalls"):
        for name in files:
            if fnmatch.fnmatch(name, '*.fastq.gz'):
                name_parts = name.split('_')
                if name_parts[1] != FLOWCELL_ID:
                    newname = name.split('_')[0]+'_'+FLOWCELL_ID+'_'+'_'.join(name.split('_')[1:])
                    os.rename( os.path.join(root, name), os.path.join(root, newname) )

    ## Clean up fastq (unless test run or deliberate choice to keep)
    ## TODO: keep fastq for RNA flowcell
    testRunString = 'TestRuns'
    keepFastqFileName = 'KEEP_FASTQ'
    
    isTestRun = testRunString in run_dir
    doKeepFastq = os.path.isfile( run_dir + '/' + keepFastqFileName )

    if( doKeepFastq ):
        os.system( 'echo "[INFO] Skipped FASTQ cleanup because keep file exists (' + keepFastqFileName + ')" >>' + LOG_FILE )
    elif( isTestRun ):
        os.system( 'echo "[INFO] Skipped FASTQ cleanup because path contains test string (' + testRunString + ')" >>' + LOG_FILE )
    else:
        os.system( 'echo "[INFO] Cleaning up all FASTQ files:" >>' + LOG_FILE )
        os.system( 'find ' + run_dir + ' -name "*.fastq.gz" -delete -print >>' + LOG_FILE )
        os.system( 'touch ' + FQDEL_FILE )

    ## Collect metrics (creates table format of yield/q30 per object)
    os.system( "%s -run_dir %s 1>>%s 2>&1" % (METRICS_TOOL, run_dir, METRICS_FILE) )

    ## Finishing touch
    os.system( 'echo "[INFO] Finished conversion ($(date))" >>' + LOG_FILE )
    os.system( 'touch ' + DONE_FILE )

def SubDirPath (d):
    return filter(os.path.isdir, [os.path.join(d,f) for f in os.listdir(d)])


## MAIN 
for searchDir in DATA_DIRS:
    
    ## Not all machines have all search dirs
    if( not os.path.exists(searchDir) ):
        continue

    ## now visit each subdir and run conversion if all conditions apply
    for runDirPath in SubDirPath(searchDir):

        runDirName = os.path.basename(runDirPath)
        hasSampleSheet = os.path.isfile(runDirPath + '/SampleSheet.csv')
        sequencingDone = os.path.isfile(runDirPath + '/RTAComplete.txt')
        conversionStrd = os.path.isfile(runDirPath + '/conversionStarted.txt')
        conversionDone = os.path.isfile(runDirPath + '/conversionDone.txt')

        ## Check if novaseq conversion is already running
        if(not hasSampleSheet):
            print("[INFO] Skipping: no SampleSheet.csv in " + runDirName)
            continue
        elif(not sequencingDone):
            print("[INFO] Skipping: no RTAComplete.txt in " + runDirName)
            continue
        elif(conversionDone):
            print("[INFO] Skipping: conversion already done for " + runDirName)
            continue
        elif(conversionStrd):
            print("[INFO] Skipping: conversion already started for " + runDirName)
            continue
        else:
            conversionIsRunning = os.path.isfile( TMP_RUNNING_FILE )
            sequencerId = runDirName.split("_")[1]
            platform = "NA"
            if(sequencerId in ISEQ_IDS or sequencerId in NEXTSEQ_IDS):
                platform = "IseqOrNextseq"

            print("[INFO] Considering conversion for: " + runDirName + " (sequencer=" + sequencerId +", platform=" + platform +")")
            
            if(conversionIsRunning and not platform == "IseqOrNextseq"):
                print("[INFO]   Skipping because another conversion is already running")
                continue
            else:
                if(platform == "IseqOrNextseq"):
                    print("[INFO]   Starting Iseq/Nextseq conversion")
                    ConvertBcl(runDirPath)
                else:
                    print("[INFO]   Starting non Iseq/Nextseq conversion")
                    print("[INFO]   Creating " + TMP_RUNNING_FILE)
                    os.system('touch ' + TMP_RUNNING_FILE)
                    ConvertBcl(runDirPath)
                    ## cleanup tmp file to allow next round of conversions to start
                    print("[INFO] Removing " + TMP_RUNNING_FILE)
                    os.system( 'rm ' + TMP_RUNNING_FILE )    
