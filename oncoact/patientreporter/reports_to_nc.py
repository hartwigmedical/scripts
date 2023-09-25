import subprocess
import os
import argparse
from api_util import get_report_created, get_sample_set
from gsutil import get_bucket_and_blob_from_gs_path
from google.cloud.storage import Bucket, Client

PIPELINE_OUTPUT_BUCKET = 'diagnostic-pipeline-output-prod-1'


def main(sample_barcode):
    report_created = get_report_created(sample_barcode)
    sample_name = report_created['sample_name']
    sample_set = get_sample_set(sample_name)
    set_name = sample_set['name']

    reports = [report_file for report_file in report_created['report_files'] if
               report_file['datatype'] in {'report_pdf', 'report_json'}]

    temp_dir_path = '~/temp'
    os.mkdir(temp_dir_path)

    client = Client()

    upload_report_json = input('Do you want to upload the final OncoAct report (PDF and json)? y or n')
    if upload_report_json.lower() == 'y':
        for report in reports:
            path = report['path']
            (bucket_name, blob) = get_bucket_and_blob_from_gs_path(path)
            remote_bucket: Bucket = client.bucket(bucket_name)
            destination_file_name = f'{temp_dir_path}/{blob}'
            remote_bucket.blob(blob).download_to_file(destination_file_name)

    upload_orange_report = input("Do you want to upload the ORANGE report? y or n")
    if upload_orange_report == 'y':
        bucket: Bucket = client.bucket(PIPELINE_OUTPUT_BUCKET)
        blob_name = f'{set_name}/orange/{sample_name}.orange.pdf'
        destination_file_name = f'{temp_dir_path}/{blob_name}'
        bucket.blob(blob_name).download_to_file(destination_file_name)

    upload_cuppa = input('Do you want to upload the CUPPA RUO report? y or n')
    if upload_cuppa == 'y':
        bucket: Bucket = client.bucket(PIPELINE_OUTPUT_BUCKET)
        blob_name = f'{set_name}/cuppa/{sample_name}_cup_report.pdf'
        destination_file_name = f'{temp_dir_path}/{blob_name}'
        bucket.blob(blob_name).download_to_file(destination_file_name)

    upload_to_nextcloud(temp_dir_path)
    os.removedirs(temp_dir_path)
    print('All done!')
    exit(0)


def upload_file_to_nextcloud(filepath):
    subprocess.check_output(['upload_file_to_nc_for_viewing', filepath])


def upload_to_nextcloud(path):
    if os.path.isdir(path):
        file = os.listdir(path)
        upload_file_to_nextcloud(f'{path}/{file}')
    else:
        upload_to_nextcloud(path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Uploads the relevant report files to NC')
    parser.add_argument('sample_barcode')
    args = parser.parse_args()
    main(sample_barcode=args.sample_barcode)
