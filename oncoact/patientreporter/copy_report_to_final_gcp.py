import requests
import argparse
from api_util import ApiUtil
from google.cloud.storage import Bucket, Blob, Client
from gsutil import get_bucket_and_blob_from_gs_path, get_file_name_from_blob
from cli_util import perform_prod_test


def main():
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument('sample_barcode')
    argument_parser.add_argument('--profile', choices=['pilot', 'prod'], default='pilot')
    args = argument_parser.parse_args()

    profile = args.profile
    perform_prod_test(profile)

    # pipeline_output_bucket = f'diagnostic-pipeline-output-{profile}-1'
    pipeline_output_bucket = f'diagnostic-pipeline-output-prod-1'
    final_bucket = f"patient-reporter-final-{profile}-1"
    portal_bucket = f'hmf-customer-portal-report-shared-{profile}'

    copy_report_to_final_gcp(args.sample_barcode,
                             profile,
                             pipeline_output_bucket=pipeline_output_bucket,
                             final_bucket=final_bucket,
                             portal_bucket=portal_bucket)


def copy_report_to_final_gcp(sample_barcode, profile, portal_bucket, final_bucket, pipeline_output_bucket):
    """
    This program does the following:

    - Copy the relevant report files from the run bucket to the final bucket.
    - Copy the relevant report files from the run bucket to the portal bucket.
    - Update shared status of the report in the HMF API.

    :param profile: the profile to run this program in (pilot, prod, etc.).
    :param portal_bucket: the name of the portal bucket.
    :param final_bucket: the name of the final bucket for internal quality purposes.
    :param pipeline_output_bucket: the bucket where the pipeline output is stored.
    :param sample_barcode: The sample_barcode whose run to process.
    """
    api_util = ApiUtil(profile)
    report_created = api_util.get_report_created(sample_barcode)
    report_files = report_created["report_files"]
    run_id = report_created['run_id']
    run = api_util.get_run(run_id)
    if not run:
        cont = input(
            f"No associate run found for tumor barcode '{sample_barcode}'. Are you sure you want to continue? (y/n)\n")
        if cont.lower() != 'y':
            exit(1)
    if run['status'] != 'Validated':
        cont = input(f"Run status for tumor barcode '{sample_barcode}' is not yet 'Validated' "
                     f"(actual status: {run['status']}). Are you sure you want to continue? (y/n)\n")
        if cont.lower() != 'y':
            exit(1)
    storage_client = Client()
    target_bucket_portal: Bucket = storage_client.bucket(portal_bucket,
                                                         user_project='hmf-customer-portal')
    target_bucket_final: Bucket = storage_client.bucket(final_bucket)

    delete_old_report(target_bucket_portal, sample_barcode)

    reports = [file for file in report_files if file['datatype'] in {'report_pdf', 'report_xml', 'report_json'}]

    for report in reports:
        (bucket, blob) = get_bucket_and_blob_from_gs_path(report['path'])
        file_name = get_file_name_from_blob(blob)
        bucket_instance: Bucket = storage_client.bucket(bucket)
        copy_and_print(bucket_instance, target_bucket_portal, blob, f'{sample_barcode}/{file_name}')
        # copy_and_print(bucket_instance, target_bucket_final, blob, blob)

    sample_name = report_created['sample_name']
    sample_set = api_util.get_sample_set_by_sample_name(sample_name)
    set_name = sample_set['name']

    orange_pdf = f'{set_name}/orange/{sample_name}.orange.pdf'
    purple_sv_vcf = f'{set_name}/purple/{sample_name}.purple.sv.vcf.gz'
    purple_somatic_vcf = f'{set_name}/purple/{sample_name}.purple.somatic.vcf.gz'
    purple_catalog = f'{set_name}/purple/{sample_name}.driver.catalog.somatic.tsv'
    linx_fusion = f'{set_name}/linx/{sample_name}.linx.fusion.tsv'
    linx_catalog = f'{set_name}/linx/{sample_name}.linx.driver.catalog.tsv'

    pipline_output_bucket: Bucket = storage_client.bucket(pipeline_output_bucket)
    for blob in [orange_pdf, purple_sv_vcf, purple_somatic_vcf, purple_catalog, linx_fusion, linx_catalog]:
        # Replace the set_name prefix with sample_barcode
        new_blob_name = f'{sample_barcode}{blob[len(set_name):]}'
        copy_and_print(pipline_output_bucket, target_bucket_portal, blob, new_blob_name)

    print('Updating report shared status in the API...')
    body = {
        'report_created_id': report_created['id'],
        'notify_users': 'false',
        'publish_to_portal': 'true'
    }
    response = requests.post(f'{api_util.api_base_url()}/hmf/v1/reports/2/shared', json=body)
    if not response.ok:
        print(f"API update failed: '{response.status_code}' reason: '{response.reason}")
        exit(1)
    print('API updated!')
    print('All done・ﾟ✧ !')
    exit(0)


def copy_and_print(source_bucket: Bucket, target_bucket: Bucket, blob_name: str, new_blob_name: str):
    """
    Copies the blob from the source bucket to the target bucket and prints this in the stdout.

    :param source_bucket: the bucket to copy from.
    :param target_bucket: the bucket to copy to.
    :param blob_name: the name of the blob.
    :param new_blob_name: the new name of the blob within the target bucket.
    """
    print(f"Copying '{blob_name}' from '{source_bucket}' to '{target_bucket}'...")
    blob = Blob(blob_name, source_bucket)
    source_bucket.copy_blob(blob=blob, destination_bucket=target_bucket, new_name=new_blob_name)


def delete_old_report(portal_bucket: Bucket, sample_barcode):
    """
    Prompts the user if they want to delete any old report artifacts for the same sample_barcode.

    :param portal_bucket: the portal bucket where all the report artifacts are copied to.
    :param sample_barcode: the sample barcode for this report.
    """
    blobs_old_run = list(portal_bucket.list_blobs(prefix=sample_barcode))
    if len(blobs_old_run) > 0:
        print(f"Old report artifacts found for report '{sample_barcode}':", [blob.name for blob in blobs_old_run])
        delete = input(f"Do you want to delete these? If you choose 'n' the program will exit now (y/n)\n")
        if delete.lower() != 'y':
            exit(1)
        portal_bucket.delete_blobs(blobs=blobs_old_run)


if __name__ == "__main__":
    main()
