import json
import requests
from google.cloud import storage
import re

# Constants
API_BASE_URL = "http://api.prod-1"
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
    created_reports = get_created_reports_from_api(sample_barcode)
    report_files = created_reports["report_files"]

    reports = next(file for file in report_files if file['datatype'] in ['report_pdf', 'report_xml', 'report_json'])
    driver_catalog = next(file for file in report_files if file['path'].endswith('Driver.catalog.somatic.tsv'))
    purple_somatic = next(file for file in report_files if file['path'].endswith('purple.somatic.vcf'))
    purple_sv = next(file for file in report_files if file['path'].endswith('purple.sv.vcf'))
    orange = next(file for file in report_files if file['path'].endswith('Orange.pdf'))
    # TODO we also need Linx fusion output. Unfortunately, I do not know the file name

    storage_client = storage.Client()
    target_bucket_portal: storage.Bucket = storage_client.bucket(PORTAL_BUCKET_NAME)
    target_bucket_final: storage.Bucket = storage_client.bucket(FINAL_BUCKET_NAME)

    for file in [*reports, driver_catalog, purple_somatic, purple_sv, orange]:
        if not file:
            continue
        (bucket, blob) = get_bucket_and_blob_from_gs_path(file['path'])
        bucket_instance: storage.Bucket = storage_client.bucket(bucket)
        print(f"Copying '{file}' to '{target_bucket_portal}'...")
        bucket_instance.copy_blob(blob, target_bucket_portal)
        print(f"Copying '{file}' to '{target_bucket_final}'...")
        bucket_instance.copy_blob(blob, target_bucket_final)

    print('Updating report shared status in the API...')

    requests.post(f'{API_BASE_URL}/hmf/v1/reports/2/shared')  # TODO

    print('all Done (ﾉ◕ヮ◕)ﾉ*:・ﾟ✧ !')
    exit(0)


def get_args_as_dictionary(args):
    """
    Transforms the given command line arguments to a dictionary, skipping optional arguments that were not set.

    :param args: the command line arguments.
    :return: a dictionary the key is the argument name and the value is the argument value.
    """
    return {k: (v[0] if len(v) == 1 else v) for (k, v) in args.__dict__.items() if v is not None}


def get_created_reports_from_api(sample_barcode: str) -> json:
    """
    Queries the 'reports/2/created' endpoint with the given sample_barcode.

    :param sample_barcode: the sample_barcode to query for.
    :return: The query result.
    :raises ValueError: If the response returned a non 2XX status code or if the query returned more than 1 sample.
    """
    response = requests.get(url=f"{API_BASE_URL}/hmf/v1/reports/2/created", params={'sample_barcode': sample_barcode})
    if not response.ok:
        raise ValueError(f"Response was not ok: {response.status_code}")
    response_json = response.json()
    if len(response_json) > 1:
        raise ValueError(f"Query returned more than one sample: '{len(response_json)}'")
    return response_json[0]


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
    print(get_bucket_and_blob_from_gs_path("gs://jwfijewf/wjfwepijf/wjfipewjfwep/wfjiwe"))
