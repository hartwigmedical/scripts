import requests
import re
import argparse
from api_util import get_report_created, API_BASE_URL, get_set
from google.cloud import storage

# Constants
PIPELINE_OUTPUT_BUCKET = 'diagnostic-pipeline-output-prod-1'
FINAL_BUCKET_NAME = "patient-reporter-final-prod-1"
PORTAL_BUCKET_NAME = "hmf-customer-portal-report-shared-prod"

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

    storage_client = storage.Client()
    target_bucket_portal: storage.Bucket = storage_client.bucket(PORTAL_BUCKET_NAME)
    target_bucket_final: storage.Bucket = storage_client.bucket(FINAL_BUCKET_NAME)

    for report in reports:
        (bucket, blob) = get_bucket_and_blob_from_gs_path(report['path'])
        bucket_instance: storage.Bucket = storage_client.bucket(bucket)
        copy_and_log(bucket_instance, target_bucket_portal, blob)
        copy_and_log(bucket_instance, target_bucket_final, blob)

    sample_name = report_created['sample_name']
    sample_set = get_set(sample_name)
    set_name = sample_set['name']

    orange_pdf = f'{set_name}/orange/{sample_name}.orange.pdf'
    purple_sv_vcf = f'{set_name}/purple/{sample_name}.purple.sv.vcf.gz'
    purple_somatic_vcf = f'{set_name}/purple/{sample_name}.purple.somatic.vcf.gz'
    purple_catalog = f'{set_name}/purple/{sample_name}.driver.catalog.somatic.tsv'
    linx_fusion = f'{set_name}/linx/{sample_name}.linx.fusion.tsv'
    linx_catalog = f'{set_name}/linx/{sample_name}.linx.driver.catalog.tsv'

    pipline_output_bucket: storage.Bucket = storage_client.bucket(PIPELINE_OUTPUT_BUCKET)
    for blob in [orange_pdf, purple_sv_vcf, purple_somatic_vcf, purple_catalog, linx_fusion, linx_catalog]:
        copy_and_log(pipline_output_bucket, target_bucket_portal, blob)

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


def copy_and_log(source_bucket: storage.Bucket, target_bucket: storage.Bucket, blob_name: str):
    blob = storage.Blob(blob_name, source_bucket)
    print(f"Copying '{blob_name}' to '{target_bucket}'...")
    source_bucket.copy_blob(blob=blob, destination_bucket=target_bucket)


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
