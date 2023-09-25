import subprocess
import os
import argparse
from api_util import get_report_created, get_set


def main(sample_barcode):
    report_created = get_report_created(sample_barcode)
    sample_name = report_created['sample_name']
    sample_set = get_set(sample_name)
    set_name = sample_set['name']

    reports = [report_file for report_file in report_created['report_files'] if
               report_file['datatype'] in {'report_pdf', 'report_json'}]

    temp_dir_path = '~/temp'
    os.mkdir(temp_dir_path)
    upload_report_json = input('Do you want to upload the final OncoAct report (PDF and json)? y or n')
    if upload_report_json.lower() == 'y':
        for report in reports:
            path = report['path']
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


def upload_file_to_nextcloud(filepath):
    subprocess.check_output(['upload_file_to_nc_for_viewing', filepath])


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
