import subprocess
import os
import shutil
import argparse
import fnmatch
from gsutil import get_bucket_and_blob_from_gs_path, get_file_name_from_blob_name
from google.cloud.storage import Bucket, Client


def main():
    parser = argparse.ArgumentParser(description='Uploads the relevant report files to NC')
    parser.add_argument('sample_barcode')

    args = parser.parse_args()

    reports_to_nc(sample_barcode=args.sample_barcode)

def reports_to_nc(sample_barcode):
    """
    This script can be used to upload the reporting artifacts to Next Cloud.

    When running this script, it will prompt the user for a set of different files they might want to upload. After the
    user has described which artifacts they want on Next Cloud, the program will automatically perform the uploading.

    :param sample_barcode: the sample barcode of the report to upload the artifacts for.
    """
    temp_dir_path = f'{os.path.expanduser("~")}/temp'

    # temp folder clean
    if os.path.isdir(temp_dir_path):
        shutil.rmtree(temp_dir_path)
    os.mkdir(temp_dir_path)

    # alle patterns combineren in een enkele regex
    gsutil_pattern = f"gs://patient-reporter-final-prod-1/**/{sample_barcode}/patient-reporter/*_oncoact_wgs_report.*"
    corr_pattern = f"gs://patient-reporter-final-prod-1/**/corr-{sample_barcode}/patient-reporter/*_oncoact_wgs_report.*"

    cmd = [
        "bash", "-c",
        f"gsutil ls -l '{gsutil_pattern}' '{corr_pattern}' "
        r"| grep -v '^Total:' "
        r"| sort -k2,2nr "
        r"| head -n 3 "
        r"| awk '{ print substr($0, index($0, $3)) }'"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    files = [f for f in result.stdout.strip().split("\n") if f]

    print("Files to download:", files)

    for filepath in files:
        subprocess.run(["gsutil", "cp", filepath, temp_dir_path], check=True)
        print("Downloaded:", filepath)

    upload_to_nextcloud(temp_dir_path)
    shutil.rmtree(temp_dir_path)

    print("All done!")
    exit(0)



def upload_to_nextcloud(path):
    """
    Uploads the contents located at the path to nextcloud.
    :param path: the path to the content. If it is a single file, it will upload that file. If it is a directory,
    it will upload all its content (recursively).
    """
    if os.path.isdir(path):
        files = os.listdir(path)
        for file in files:
            upload_to_nextcloud(f'{path}/{file}')
    else:
        _upload_file_to_nextcloud(path)


def _upload_file_to_nextcloud(filepath):
    subprocess.check_output(['upload_file_to_nc_for_sharing', filepath])


if __name__ == '__main__':
    main()
