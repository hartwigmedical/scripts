import requests
import re
import argparse
from api_util import get_report_created, API_BASE_URL
from google.cloud import storage

# Constants

FINAL_BUCKET_NAME = "patient-reporter-final-prod-1"
PORTAL_BUCKET_NAME = ""

# Regex for gs file locations: gs:://<bucket-name>/<blob-name>
GS_PATH_REGEX = re.compile(r'^gs://([a-z0-9._-]+)((?:/[a-zA-Z0-9_.-]+)*)$')


def main(sample_barcode):
    """
    This program does the following:

    - Copy the relevant report files from the run bucket to the final bucket.
    - Copy the relevant report files from the run bucket to the portal bucket.
    - Update shared status of the report in the HMF API.

    :param sample_barcode: The sample_barcode whose run to process.
    """
    report_created = get_report_created(sample_barcode)
    report_files = report_created["report_files"]

    reports = [file for file in report_files if file['datatype'] in {'report_pdf', 'report_xml', 'report_json'}]
    # driver_catalog = next(file for file in report_files if file['path'].endswith('Driver.catalog.somatic.tsv'))
    # purple_somatic = next(file for file in report_files if file['path'].endswith('purple.somatic.vcf'))
    # purple_sv = next(file for file in report_files if file['path'].endswith('purple.sv.vcf'))
    # orange = next(file for file in report_files if file['path'].endswith('Orange.pdf'))
    # TODO we also need Linx fusion output. Unfortunately, I do not know the file name

    storage_client = storage.Client()
    target_bucket_portal: storage.Bucket = storage_client.bucket(PORTAL_BUCKET_NAME)
    target_bucket_final: storage.Bucket = storage_client.bucket(FINAL_BUCKET_NAME)

    for file in reports:
        (bucket, blob) = get_bucket_and_blob_from_gs_path(file['path'])
        bucket_instance: storage.Bucket = storage_client.bucket(bucket)
        print(f"Copying '{file}' to '{target_bucket_portal}'...")
        bucket_instance.copy_blob(blob, target_bucket_portal)
        print(f"Copying '{file}' to '{target_bucket_final}'...")
        bucket_instance.copy_blob(blob, target_bucket_final)

    print('Updating report shared status in the API...')
    params = {
        'report_created_id': report_created['id'],
        'notify_users': False,
        'publish_to_portal': True
    }
    requests.post(f'{API_BASE_URL}/hmf/v1/reports/2/shared', params=params)
    print('API updated!')
    print('all Done (ﾉ◕ヮ◕)ﾉ*:・ﾟ✧ !')
    exit(0)


def get_bucket_and_blob_from_gs_path(gs_path: str) -> (str, str):
    """
    Returns the bucket name and the blob name from a given google storage path.

    :param gs_path: the path in the following format: (gs://<bucket-name>/<blob-name>).
    :return: a tuple containing both the bucket name and blob name: (bucket, blob).
    """
    match = GS_PATH_REGEX.match(gs_path)
    if not match:
        raise ValueError("No match!")

    groups = match.groups()
    bucket = groups[0]
    blob = groups[1][1:] if groups[1] else None  # the [1:] is to remove the leading slash from the blob name

    return bucket, blob


if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument('sample_barcode')
    args = argument_parser.parse_args()
    main(args.sample_barcode)
