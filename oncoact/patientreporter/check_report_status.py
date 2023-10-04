import argparse
import pandas as pd
import subprocess

import requests

from rest_util import RestClient
from gsutil import get_bucket_and_blob_from_gs_path
from google.cloud.storage import Bucket, Blob, Client


class ReportWithStatus:
    def __init__(self, report, status):
        self.report = report
        self.status = status


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', choices=['pilot', 'prod'], default='pilot')
    args = parser.parse_args()

    storage_client = Client(project='hmf-pipeline-prod')
    check_report_status(storage_client, args.profile)


def check_report_status(storage_client: Client, profile: str):
    rest_client = RestClient(profile)

    all_reports = pd.DataFrame(rest_client.get_all_reports_created())
    # these if checks are to ensure that if the dataframe is empty,
    # it still has the required columns to prevent KeyErrors.
    if len(all_reports) == 0:
        all_reports = pd.DataFrame(columns=['sample_barcode', 'run_id'])
    all_shared_reports = pd.DataFrame([shared['report_created'] for shared in rest_client.get_all_reports_shared()])
    if len(all_shared_reports) == 0:
        all_shared_reports = pd.DataFrame(columns=['sample_barcode', 'run_id'])
    all_runs = pd.DataFrame(rest_client.get_runs())
    if len(all_runs) == 0:
        all_runs = pd.DataFrame(columns=['status', 'id'])

    validated_runs = all_runs[all_runs['status'] == 'Validated']
    finished_runs = all_runs[all_runs['status'] == 'Finished']
    failed_runs = all_runs[all_runs['status'] == 'Failed']

    reports_not_shared = all_reports[~all_reports['sample_barcode'].isin(all_shared_reports['sample_barcode'])]

    finished_reports_but_shared = all_shared_reports[all_shared_reports['run_id'].isin(finished_runs)]
    failed_report_shared = all_shared_reports[all_shared_reports['run_id'].isin(failed_runs)]

    validated_reports_and_shared = all_shared_reports[all_shared_reports['run_id'].isin(validated_runs['id'])].iloc[:5]
    validated_reports_but_not_shared = reports_not_shared[reports_not_shared['run_id'].isin(validated_runs['id'])].iloc[
                                       :5]
    validated_warnings_and_errors = check_errors_for_validated_reports(storage_client, pd.concat(
        [validated_reports_and_shared, validated_reports_but_not_shared], ignore_index=True),
                                                                       validated_runs, profile)

    finished_reports_not_shared = reports_not_shared[reports_not_shared['run_id'].isin(finished_runs)]

    failed_report_but_not_shared = reports_not_shared[reports_not_shared['run_id'].isin(failed_runs)]

    runs_without_report = all_runs[~all_runs['id'].isin(all_reports['run_id'])]
    reports_without_runs = all_reports[~all_reports['run_id'].isin(all_runs['id'])]

    failed_executions = rest_client.get_failed_executions()
    failed_sample_barcodes_with_reason = {rest_client.get_tumor_sample_barcode_from_run_id(run_id): errors for
                                          (run_id, errors) in failed_executions.items()}

    print('----------')
    print("Failed reports that have not yet been shared:\n")
    for i, report in failed_report_but_not_shared.iterrows():
        # fail_type = failed_runs[failed_runs['id'] == report['sample_barcode']['run_id']]['failure']['type']
        print(f"{i + 1} - {report['sample_barcode']}")
        print(f"\tSet name: {get_set_name_from_report(report, failed_runs)}")

    print('----------')
    print("Failed reports that have been shared:\n")
    for i, report in failed_report_shared.iterrows():
        # fail_type = failed_runs[failed_runs['id'] == report['sample_barcode']['run_id']]['failure']['type']
        print(f"{i + 1} - {report['sample_barcode']}")
        print(f"\tset name: {get_set_name_from_report(report, failed_runs)}")

    print('----------')
    print("Finished reports that have not yet been shared:\n")
    for i, report in finished_reports_not_shared.iterrows():
        print(f"{i + 1} - {report['sample_barcode']}")
        print(f"\tset name: {get_set_name_from_report(report, finished_runs)}")

    print('----------')
    print("Finished reports that have already been shared (BUT NOT VALIDATED YET!):\n")
    for i, report in finished_reports_but_shared.iterrows():
        print(f"{i + 1} - {report['sample_barcode']}")
        print(f"\tset name: {get_set_name_from_report(report, finished_runs)}")

    print('----------')
    print("Validated reports that have not yet been shared:\n")
    for i, report in validated_reports_but_not_shared.iterrows():
        print(f"{i + 1} - {report['sample_barcode']}")
        print(f"\tset name: {get_set_name_from_report(report, validated_runs)}")
        report_id = report['id']
        if validated_warnings_and_errors[report_id]:
            print(f"\tWarnings: {validated_warnings_and_errors[report_id]}")

    print('----------')
    print("Reports without a run:\n")
    for i, report in enumerate(reports_without_runs['sample_barcode']):
        print(f"{i + 1} - {report}")

    # Runs that failed in the patient reporter:
    print('----------')
    print("Runs whose execution failed in the reporting pipeline")
    for i, (tumor_barcode, errors) in enumerate(failed_sample_barcodes_with_reason):
        print(f"{i+1} - {tumor_barcode} : {errors}")


def check_errors_for_validated_reports(client: Client, validated_reports: pd.DataFrame,
                                       validated_runs: pd.DataFrame, profile):
    res = {}
    for i, report in validated_reports.iterrows():
        print('Checking errors and warnings for validated runs...', i, '/', len(validated_reports))
        run = validated_runs[validated_runs['id'] == report['run_id']].iloc[:1]
        if run.empty:
            continue

        report_files = report['report_files']
        log = next((file for file in report_files if file['datatype'] == 'report_log'), None)
        if not log:
            continue
        (bucket_name, blob_name) = get_bucket_and_blob_from_gs_path(log['path'])

        bucket: Bucket = client.bucket(bucket_name)
        blob: Blob = bucket.get_blob(blob_name=blob_name)
        log_content = blob.download_as_string().decode()

        set_name = run['set'].values[0]['name']
        health_error = get_health_error_validated_run(set_name)

        errors = []

        if has_rose_error(client, set_name, report['sample_name'], profile):
            errors.append("Summary was incorrect (rose error)")
        if "WARN" in log_content:
            errors.append("Warn was present in log file")
        if "Consent" in log_content or "Mismatching ref sample name" in log_content or "do not match" in log_content or "Missing or invalid hospital" in log_content:
            errors.append("Lims error")
        if "WARN" in health_error:
            errors.append("Health check error")

        res[report['id']] = errors
    return res


def get_health_error_validated_run(set_name: str) -> str:
    return subprocess.check_output(['health_check_validated_run', set_name, '2>&1']).decode()


def has_rose_error(client: Client, set_name: str, sample_name: str, profile: str):
    """
    Finds whether there is an error with the 'rose' tool for a given set and sample.

    :param client: the gcp storage client to access the rose log.
    :param set_name: the set to check for.
    :param sample_name: the sample to check for within the set.
    :param profile: the profile to run this script in (pilot/prod/etc).
    :return: true if a 'rose' error was found, else false.
    """
    bucket_name = 'diagnostic-pipeline-output-prod-1' if profile == 'prod' else 'diagnostic-pipeline-output-pilot-1'

    bucket: Bucket = client.bucket(bucket_name)
    blob: Blob = bucket.get_blob(f'{set_name}/purple/{sample_name}.driver.catalog.somatic.tsv')

    content = blob.download_as_string().decode()
    return content.find('AMP') > 0


def get_set_name_from_report(report, runs):
    return runs[runs['id'] == report['run_id']].iloc[:1]['set'].values[0]['name']


if __name__ == '__main__':
    main()
