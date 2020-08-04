#!/usr/bin/env python

import fnmatch
import os

DATA_DIRS=[
    '/data1/illumina_data/',
    '/data1/illumina_data/TestRuns/'
]

BCL2FASTQ_TOOL='/data/fastqconversion/bcl2fastq_v2.20.0.422/bin/bcl2fastq'
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
    os.system( "check_bcl2fastq_conversion.pl -run_dir %s 1>>%s 2>&1" % (run_dir, METRICS_FILE) )

    ## Finishing touch
    os.system( 'echo "[INFO] Finished conversion ($(date))" >>' + LOG_FILE )
    os.system( 'touch ' + DONE_FILE )

def SubDirPath (d):
    return filter(os.path.isdir, [os.path.join(d,f) for f in os.listdir(d)])


## MAIN 
for search_dir in DATA_DIRS:
    
    ## Not all machines have all search dirs
    if( not os.path.exists(search_dir) ):
        continue
    
    ## Check if conversion is already running
    conversionRunning = os.path.isfile( TMP_RUNNING_FILE )
    if( conversionRunning ):
        print("[INFO] Skipping conversions because conversion is running (" + TMP_RUNNING_FILE + ")")
        continue
    else:
        os.system( 'touch ' + TMP_RUNNING_FILE )

    ## now visit each subdir and run conversion if all conditions apply
    for run_dir in SubDirPath(search_dir):
        hasSampleSheet = os.path.isfile( run_dir + '/SampleSheet.csv' )
        sequencingDone = os.path.isfile( run_dir + '/RTAComplete.txt' )
        conversionStrd = os.path.isfile( run_dir + '/conversionStarted.txt' )
        conversionDone = os.path.isfile( run_dir + '/conversionDone.txt' )

        if( hasSampleSheet and sequencingDone and not conversionStrd and not conversionDone):
            print("[INFO] Starting conversion for flowcell: " + run_dir)
            ConvertBcl(run_dir)
        else:
            print("[INFO] Skipping conversion for flowcell: " + run_dir)

    ## cleanup tmp file to allow next round of conversions to start later
    os.system( 'rm ' + TMP_RUNNING_FILE )
    
