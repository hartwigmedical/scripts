import json
import requests
import subprocess
import os
import argparse


def main(sample_barcode):
    if not sample_barcode:
        raise ValueError('No sample barcode provided')
    run = get_run(sample_barcode)
    set_name = get_set_name_from_run(run)
    sample_name = set['tumor-sample']
    report_files: list[json] = [report_file for report_file in run['report_files'] if
                                report_file['stage'] == 'reporting-pipeline']

    temp_dir_path = '~/temp'
    os.mkdir(temp_dir_path)
    upload_report_json = input('Do you want to upload the final OncoAct report (PDF and json)? y or n')
    if upload_report_json == 'y':
        for file in report_files:
            path = file['path']
            if path.endswith('.pdf') or path.endswith('.json'):
                subprocess.run(['gsutil', 'cp', path, temp_dir_path])

    upload_orange_report = input("Do you want to upload the ORANGE report? y or n")
    if upload_orange_report == 'y':
        subprocess.run(
            ['gsutil', 'cp', f'gs://diagnostic-pipeline-output-prod-1/{set_name}/orange/{sample_name}*.orange.pdf',
             temp_dir_path])

    upload_cuppa = input('Do you want to upload the CUPPA RUO report? y or n')
    if upload_cuppa == 'y':
        subprocess.run(
            ['gsutil', 'cp', f'gs://diagnostic-pipeline-output-prod-1/{set_name}/cuppa/{sample_name}*_cup_report.pdf',
             temp_dir_path])

    upload_to_nextcloud(temp_dir_path)
    os.removedirs(temp_dir_path)
    print('All done!')
    exit(0)


def get_run(sample_barcode):
    response = requests.get('http://api.prod-1/hmf/v1/reports/2/created', params={'sample_barcode': sample_barcode})
    if not response.ok:
        raise ValueError(f'Response status {response.status_code}')

    return response.json()[-1]


def get_set_name_from_run(run) -> str:
    return run['set']['name']


def upload_file_to_nextcloud(filepath):
    subprocess.run(['upload_file_to_nc_for_viewing', filepath])
    print(f"Uploaded '{filepath}' to NC")


def upload_to_nextcloud(path):
    if os.path.isdir(path):
        file = os.listdir(path)
        upload_file_to_nextcloud(f'{path}/{file}')
    else:
        upload_to_nextcloud(path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Uploads the relevant report files to NC")
    parser.add_argument('sample_barcode')
    args = parser.parse_args()
    main(sample_barcode=args.sample_barcode)
