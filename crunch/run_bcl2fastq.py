#!/usr/bin/env python

import fnmatch
import os

DATA_DIRS = ['/data1/illumina_data/', '/data1/illumina_data/TestRuns/']
NEXTSEQ_IDS = ['NB500901', 'NB500902']
ISEQ_IDS = ['FS10001173']

BCL2FASTQ_TOOL = '/data/fastqconversion/bcl2fastq_v2.20.0.422/bin/bcl2fastq'
METRICS_TOOL = '/data/repos/scripts/crunch/check_bcl2fastq_conversion.pl'
TMP_RUNNING_FILE = '/tmp/bcl2fastq_running'


def convert_bcl(run_dir):
    log_file = run_dir + '/conversionLog.txt'
    err_file = run_dir + '/conversionError.txt'
    started_file = run_dir + '/conversionStarted.txt'
    done_file = run_dir + '/conversionDone.txt'
    fqdel_file = run_dir + '/conversionFastqDeleted.txt'

    out_dir = run_dir + '/Fastq'
    metrics_file = out_dir + '/Stats/conversion_metrics_table.txt'

    flowcell_id = os.path.split(run_dir)[-1].split("_")[-1]
    print("[INFO]   Flowcell ID: " + flowcell_id)
    print("[INFO]   Output directory: " + out_dir)

    os.system('touch ' + started_file)
    os.system('echo "[INFO] Started conversion of ' + flowcell_id + ' ($(date))" >>' + log_file)
    os.system(BCL2FASTQ_TOOL + " --runfolder-dir %s --output-dir %s -l WARNING 1>>%s 2>>%s" % (
        run_dir, out_dir, log_file, err_file))

    # Add flowcell id to fastq file names
    os.system('echo "[INFO] Output directory set to ' + out_dir + '" >>' + log_file)
    os.system('echo "[INFO] Adding flowcell ID to output FASTQ filenames" >>' + log_file)
    for root, dirs, files in os.walk(run_dir + "/Data/Intensities/BaseCalls"):
        for name in files:
            if fnmatch.fnmatch(name, '*.fastq.gz'):
                name_parts = name.split('_')
                if name_parts[1] != flowcell_id:
                    new_name = name.split('_')[0] + '_' + flowcell_id + '_' + '_'.join(name.split('_')[1:])
                    os.rename(os.path.join(root, name), os.path.join(root, new_name))

    # Clean up fastq (unless test run or deliberate choice to keep)
    # TODO: keep fastq for RNA flowcell
    test_run_string = 'TestRuns'
    keep_fastq_file = 'KEEP_FASTQ'

    is_test_run = test_run_string in run_dir
    do_keep_fastq = os.path.isfile(run_dir + '/' + keep_fastq_file)

    if do_keep_fastq:
        os.system(
            'echo "[INFO] Skipping FASTQ cleanup because keep file exists (' + keep_fastq_file + ')" >>' + log_file)
    elif is_test_run:
        os.system(
            'echo "[INFO] Skipping FASTQ cleanup because is testRuns flowcell (' + test_run_string + ')" >>' + log_file)
    else:
        os.system('echo "[INFO] Cleaning up all FASTQ files:" >>' + log_file)
        os.system('find ' + run_dir + ' -name "*.fastq.gz" -delete -print >>' + log_file)
        os.system('touch ' + fqdel_file)

    # Collect metrics (creates table format of yield/q30 per object)
    os.system("%s -run_dir %s 1>>%s 2>&1" % (METRICS_TOOL, run_dir, metrics_file))

    # Finishing touch
    os.system('echo "[INFO] Finished conversion ($(date))" >>' + log_file)
    os.system('touch ' + done_file)


def sub_dir_path(d):
    return filter(os.path.isdir, [os.path.join(d, f) for f in os.listdir(d)])


# MAIN
for search_dir in DATA_DIRS:

    # Not all machines have all search dirs
    if not os.path.exists(search_dir):
        continue

    # Now visit each subdir and run conversion if all conditions apply
    for run_dir_path in sub_dir_path(search_dir):

        run_dir_name = os.path.basename(run_dir_path)
        has_samplesheet = os.path.isfile(run_dir_path + '/SampleSheet.csv')
        sequencing_done = os.path.isfile(run_dir_path + '/RTAComplete.txt')
        conversion_started = os.path.isfile(run_dir_path + '/conversionStarted.txt')
        conversion_done = os.path.isfile(run_dir_path + '/conversionDone.txt')

        # Check if novaseq conversion is already running
        if not has_samplesheet:
            print("[INFO] Skipping: no SampleSheet.csv in " + run_dir_name)
            continue
        elif not sequencing_done:
            print("[INFO] Skipping: no RTAComplete.txt in " + run_dir_name)
            continue
        elif conversion_done:
            print("[INFO] Skipping: conversion already done for " + run_dir_name)
            continue
        elif conversion_started:
            print("[INFO] Skipping: conversion already started for " + run_dir_name)
            continue
        else:
            conversion_is_running = os.path.isfile(TMP_RUNNING_FILE)
            sequencer_id = run_dir_name.split("_")[1]
            platform = "NA"
            if sequencer_id in ISEQ_IDS or sequencer_id in NEXTSEQ_IDS:
                platform = "IseqOrNextseq"

            print('[INFO] Considering conversion for: ' + run_dir_name)

            if conversion_is_running and not platform == "IseqOrNextseq":
                print("[INFO]   Skipping because another conversion is already running")
                continue
            else:
                if platform == "IseqOrNextseq":
                    print("[INFO]   Starting Iseq/Nextseq conversion")
                    convert_bcl(run_dir_path)
                else:
                    print("[INFO]   Starting non Iseq/Nextseq conversion")
                    print("[INFO]   Creating " + TMP_RUNNING_FILE)
                    os.system('touch ' + TMP_RUNNING_FILE)
                    convert_bcl(run_dir_path)

                    # Cleanup tmp file to allow next round of conversions to start
                    print("[INFO] Removing " + TMP_RUNNING_FILE)
                    os.system('rm ' + TMP_RUNNING_FILE)
