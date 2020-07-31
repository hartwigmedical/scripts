#!/usr/bin/env python

import fnmatch
import os

DATA_DIRS=[
    '/data1/illumina_data/',
    '/data1/illumina_data/TestRuns/'
]

BCL2FASTQ_TOOL='/data/fastqconversion/bcl2fastq_v2.20.0.422/bin/bcl2fastq'

## DEFS
def ConvertBcl(run_dir):

    LOG_FILE  = run_dir + '/conversionLog.txt';
    ERR_FILE  = run_dir + '/conversionError.txt';
    STRT_FILE = run_dir + '/conversionStarted.txt';
    DONE_FILE = run_dir + '/conversionDone.txt';
    FQDEL_FILE = run_dir + '/conversionFastqDeleted.txt';

    out_dir = run_dir + '/Fastq'

    FLOWCELL_ID = os.path.split(run_dir)[-1].split("_")[-1]
    print "[INFO]   Flowcell ID: " + FLOWCELL_ID
    print "[INFO]   Output directory: " + out_dir

    os.system( 'touch '  + STRT_FILE )
    os.system( 'echo "[INFO] Started at $(date)" >>' + LOG_FILE )
    os.system( BCL2FASTQ_TOOL + " --runfolder-dir %s --output-dir %s -r 18 -p 18 -w 4 -l WARNING 1>>%s 2>>%s" % (run_dir, out_dir, LOG_FILE, ERR_FILE) )

    ## Add flowcell id to fastq file names
    os.system( 'echo "[INFO] Adding flowcell ID to FASTQ files" >>' + LOG_FILE )
    for root, dirs, files in os.walk(run_dir+"/Data/Intensities/BaseCalls"):
        for name in files:
            if fnmatch.fnmatch(name, '*.fastq.gz'):
                name_parts = name.split('_')
                if name_parts[1] != FLOWCELL_ID:
                    newname = name.split('_')[0]+'_'+FLOWCELL_ID+'_'+'_'.join(name.split('_')[1:])
                    os.rename( os.path.join(root, name), os.path.join(root, newname) )

    ## Clean up fastq
    keepFastqFileName = 'KEEP_FASTQ'
    hasKeepFastqFile = os.path.isfile( run_dir + '/' + keepFastqFileName )
    if( not hasKeepFastqFile ):
        os.system( 'echo "[INFO] Cleaning up all FASTQ files:" >>' + LOG_FILE )
        os.system( 'find ' + run_dir + ' -name "*.fastq.gz" -delete -print >>' + LOG_FILE )
        os.system( 'touch ' + FQDEL_FILE )
    else:
        os.system( 'echo "[INFO] Skipping FASTQ cleanup because keep file exists (' + keepFastqFileName + ')" >>' + LOG_FILE )
    
    ## finish up
    os.system( 'echo "[INFO] Finished at  $(date)" >>' + LOG_FILE )
    os.system( 'touch ' + DONE_FILE )

def SubDirPath (d):
    return filter(os.path.isdir, [os.path.join(d,f) for f in os.listdir(d)])

## MAIN 
for search_dir in DATA_DIRS:
    
    ## not all machines have all search dirs
    if( not os.path.exists(search_dir) ):
        continue

    ## now visit each subdir and run conversion if all conditions apply
    for run_dir in SubDirPath(search_dir):
        hasSampleSheet = os.path.isfile( run_dir + '/SampleSheet.csv' )
        sequencingDone = os.path.isfile( run_dir + '/RTAComplete.txt' )
        conversionStrd = os.path.isfile( run_dir + '/conversionStarted.txt' )
        conversionDone = os.path.isfile( run_dir + '/conversionDone.txt' )

        if( hasSampleSheet and sequencingDone and not conversionStrd and not conversionDone):
            print "[INFO] Starting conversion for flowcell: " + run_dir
            ConvertBcl(run_dir)
        else:
            print "[INFO] Skipping conversion for flowcell: " + run_dir
