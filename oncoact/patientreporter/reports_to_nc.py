import subprocess
import os
import shutil
import argparse
from rest_util import RestClient
from gsutil import get_bucket_and_blob_from_gs_path, get_file_name_from_blob
from google.cloud.storage import Bucket, Client


def main():
    parser = argparse.ArgumentParser(description='Uploads the relevant report files to NC')
    parser.add_argument('sample_barcode')

    args = parser.parse_args()

    pipeline_output_bucket = 'diagnostic-pipeline-output-prod-1'
    reports_to_nc(sample_barcode=args.sample_barcode,
                  pipeline_output_bucket=pipeline_output_bucket)


def reports_to_nc(sample_barcode, pipeline_output_bucket):
    """
    This script can be used to upload the reporting artifacts to Next Cloud.

    When running this script, it will prompt the user for a set of different files they might want to upload. After the
    user has described which artifacts they want on Next Cloud, the program will automatically perform the uploading.

    :param sample_barcode: the sample barcode of the report to upload the artifacts for.
    :param pipeline_output_bucket: the output bucket of the pipeline. This is where the artifacts are located.
    """
    api_util = RestClient(profile='prod')
    report_created = api_util.get_report_created(sample_barcode)
    sample_name = report_created['sample_name']
    sample_set = api_util.get_sample_set_by_sample_name(sample_name)
    set_name = sample_set['name']

    reports = [report_file for report_file in report_created['report_files'] if
               report_file['datatype'] in {'report_pdf', 'report_json'}]

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

    upload_report_json = input('Do you want to upload the final OncoAct report (PDF and json)? y or n\n')
    if upload_report_json.lower() == 'y':
        for report in reports:
            path = report['path']
            (bucket_name, blob) = get_bucket_and_blob_from_gs_path(path)
            blob_file = get_file_name_from_blob(blob)
            remote_bucket: Bucket = client.bucket(bucket_name)
            destination_file_name = f'{temp_dir_path}/{blob_file}'
            with open(destination_file_name, 'xb') as file:
                remote_bucket.blob(blob).download_to_file(file)

    upload_orange_report = input("Do you want to upload the ORANGE report? y or n\n")
    if upload_orange_report.lower() == 'y':
        bucket: Bucket = client.bucket(pipeline_output_bucket)
        blob_name = f'{set_name}/orange/{sample_name}.orange.pdf'
        blob_file = get_file_name_from_blob(blob_name)
        destination_file_name = f'{temp_dir_path}/{blob_file}'
        with open(destination_file_name, 'xb'):
            bucket.blob(blob_name).download_to_filename(destination_file_name)

    upload_cuppa = input('Do you want to upload the CUPPA RUO report? y or n\n')
    if upload_cuppa.lower() == 'y':
        bucket: Bucket = client.bucket(pipeline_output_bucket)
        blob_name = f'{set_name}/cuppa/{sample_name}_cup_report.pdf'
        blob_file = get_file_name_from_blob(blob_name)
        destination_file_name = f'{temp_dir_path}/{blob_file}'
        with open(destination_file_name, 'xb'):
            bucket.blob(blob_name).download_to_filename(destination_file_name)

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
        upload_to_nextcloud(path)


def _upload_file_to_nextcloud(filepath):
    subprocess.check_output(['upload_file_to_nc_for_viewing', filepath])


if __name__ == '__main__':
    main()
