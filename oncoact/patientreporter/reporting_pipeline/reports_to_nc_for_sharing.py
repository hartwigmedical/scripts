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
    patterns = [
        f"{sample_barcode}/patient-reporter/*_oncoact_wgs_report.xml",
        f"{sample_barcode}/patient-reporter/*_oncoact_wgs_report.json",
        f"{sample_barcode}/patient-reporter/*_oncoact_wgs_report.pdf",
        f"corr-{sample_barcode}/patient-reporter/*_oncoact_wgs_report.xml",
        f"corr-{sample_barcode}/patient-reporter/*_oncoact_wgs_report.json",
        f"corr-{sample_barcode}/patient-reporter/*_oncoact_wgs_report.pdf"
    ]

    temp_dir_path = f'{os.path.expanduser("~")}/temp'

    try:
        os.mkdir(temp_dir_path)
    except FileExistsError as e:
        print(e)
        delete_and_try_again = input('Do you want to delete this folder and try again? (y/n)\n').lower() == 'y'
        if not delete_and_try_again:
            exit(1)
        shutil.rmtree(temp_dir_path)
        os.mkdir(temp_dir_path)

    client = Client()
    bucket = client.bucket("patient-reporter-final-prod-1")

    for pattern in patterns:
        prefix = pattern.split("*")[0]

        blobs = client.list_blobs(bucket, prefix=prefix)

        for blob in blobs:
            if fnmatch.fnmatch(blob.name, pattern):
                destination = os.path.join(temp_dir_path, os.path.basename(blob.name))
                with open(destination, "wb") as f:
                    blob.download_to_file(f)
                print("Downloaded:", blob.name)

    upload_to_nextcloud(temp_dir_path)
    shutil.rmtree(temp_dir_path)
    print('All done!')
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
