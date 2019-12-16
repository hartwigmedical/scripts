#!/usr/bin/env python

import fnmatch
import os

#from email.mime.text import MIMEText
#from xml.dom.minidom import parse
#import smtplib

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

    FLOWCELL_ID = os.path.split(run_dir)[-1].split("_")[-1]
    print "[INFO] Found flowcell: " + FLOWCELL_ID

    os.system( "touch "  + STRT_FILE )
    os.system( "date >>" + LOG_FILE )
    os.system( BCL2FASTQ_TOOL + " --runfolder-dir %s -r 18 -p 18 -w 4 -l WARNING 1>>%s 2>>%s" % (run_dir, LOG_FILE, ERR_FILE) )

    ## Add flowcell id to fastq file names
    for root, dirs, files in os.walk(run_dir+"/Data/Intensities/BaseCalls"):
        for name in files:
            if fnmatch.fnmatch(name, '*fastq.gz'):
                name_parts = name.split('_')
                if name_parts[1] != FLOWCELL_ID:
                    newname = name.split('_')[0]+'_'+FLOWCELL_ID+'_'+'_'.join(name.split('_')[1:])
                    os.rename( os.path.join(root, name), os.path.join(root, newname) )

    ## finish up
    os.system( "date >>" + LOG_FILE )
    os.system( "touch "  + DONE_FILE )

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
            print "[INFO] Start conversion of: " + run_dir
            ConvertBcl(run_dir)


