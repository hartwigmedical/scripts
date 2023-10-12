import argparse
import pandas as pd
import subprocess

from rest_util import RestClient
from gsutil import get_bucket_and_blob_names_from_gs_path
from google.cloud.storage import Bucket, Blob, Client


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', choices=['pilot', 'prod'], default='pilot')
    args = parser.parse_args()

    storage_client = Client(project='hmf-pipeline-prod')
    check_report_and_run_status(storage_client, args.profile)


def check_report_and_run_status(storage_client: Client, profile: str):
    rest_client = RestClient(profile)

    print("Fetching reports and runs...")
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

    not_validated_reports_but_shared = all_shared_reports[all_shared_reports['run_id'].isin(finished_runs)]
    failed_report_shared = all_shared_reports[all_shared_reports['run_id'].isin(failed_runs)]

    validated_reports_and_shared = all_shared_reports[all_shared_reports['run_id'].isin(validated_runs['id'])].iloc[:5]
    validated_reports_but_not_shared = reports_not_shared[reports_not_shared['run_id'].isin(validated_runs['id'])].iloc[
                                       :5]
    validated_warnings_and_errors = check_warnings_for_validated_reports(storage_client, pd.concat(
        [validated_reports_and_shared, validated_reports_but_not_shared], ignore_index=True),
                                                                         validated_runs, profile)

    finished_reports_not_shared = reports_not_shared[reports_not_shared['run_id'].isin(finished_runs)]
    failed_report_but_not_shared = reports_not_shared[reports_not_shared['run_id'].isin(failed_runs)]

    runs_without_report = all_runs[~all_runs['id'].isin(all_reports['run_id'])]
    reports_without_runs = all_reports[~all_reports['run_id'].isin(all_runs['id'])]

    failed_executions = rest_client.get_failed_executions()

    failed_reports_with_errors = {}
    for run_id, errors in failed_executions:
        sample_barcode = rest_client.get_tumor_sample_barcode_from_run_id(run_id)
        failed_reports_with_errors[sample_barcode] = errors

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
    for i, report in not_validated_reports_but_shared.iterrows():
        print(f"{i + 1} - {report['sample_barcode']}")
        print(f"\tset name: {get_set_name_from_report(report, finished_runs)}")

    print('----------')
    print("Validated reports that have not yet been shared:\n")
    for i, report in validated_reports_but_not_shared.iterrows():
        print(f"{i + 1} - {report['sample_barcode']}")
        print(f"\tset name: {get_set_name_from_report(report, validated_runs)}")
        report_id = report['id']
        print("\tWarnings:")
        if validated_warnings_and_errors[report_id]:
            for warning in validated_warnings_and_errors[report_id]:
                print(f'\t{warning}')
        else:
            print("\t\tNone, Ready for sharing :)!")

    print('----------')
    print("Reports without a run:\n")
    for i, report in enumerate(reports_without_runs['sample_barcode']):
        print(f"{i + 1} - {report}")

    # Runs that failed in the patient reporter:
    print('----------')
    print("Runs whose execution failed in the reporting pipeline")
    for i, (tumor_barcode, errors) in enumerate(failed_reports_with_errors):
        print(f"{i + 1} - {tumor_barcode} : {errors}")

    print('----------')
    print('Runs without a report')
    for i, run_id in enumerate(runs_without_report):
        sample_barcode = rest_client.get_tumor_sample_barcode_from_run_id(run_id)
        print(f'{i + 1} - {sample_barcode}')
        print(f'\trun_id: {run_id}')


class Warning:
    def __init__(self, warning, action):
        self._warning = warning
        self._action = action

    def __str__(self):
        return f"""
        -----
        Warning: {self._warning}
        Recommended action: 
            - {self._action}
        -----
        """


def check_warnings_for_validated_reports(client: Client, validated_reports: pd.DataFrame,
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
        (bucket_name, blob_name) = get_bucket_and_blob_names_from_gs_path(log['path'])

        bucket: Bucket = client.bucket(bucket_name)
        blob: Blob = bucket.get_blob(blob_name=blob_name)
        log_content = blob.download_as_string().decode()

        set_name = run['set'].values[0]['name']
        health_error = get_health_error_validated_run(set_name)

        warnings = []

        if has_rose_error(client, set_name, report['sample_name'], profile):
            rose_location = get_rose_location_as_gs_string(profile, set_name, report['sample_name'])
            warnings.append(
                Warning('Summary was incorrect (rose error).', f"Check for errors 'gsutil cat {rose_location}'"))
        if "WARN" in log_content:
            warnings.append(
                Warning('Warn was present in log file.', f"Check for warns 'gsutil cat {log['path']} | grep WARN'"))
        if (
                "Consent" in log_content or "Mismatching ref sample name" in log_content or "do not match" in log_content or "Missing or invalid hospital" in log_content):
            warnings.append(Warning('Lims error.', f"Check for lims problems 'gsutil cat {log['path']}'"))
        if "WARN" in health_error:
            warnings.append(Warning(f'Health check error.', f"Hint: [{health_error}]"))

        res[report['id']] = warnings
    return res


def get_health_error_validated_run(set_name: str) -> str:
    return subprocess.check_output(['health_check_validated_run', set_name, '2>&1']).decode()


def has_rose_error(storage_client: Client, set_name: str, sample_name: str, profile: str):
    """
    Finds whether there is an error with the 'rose' tool for a given set and sample.

    :param storage_client: the gcp storage client to access the rose log.
    :param set_name: the set to check for.
    :param sample_name: the sample to check for within the set.
    :param profile: the profile to run this script in (pilot/prod/etc).
    :return: true if a 'rose' error was found, else false.
    """
    bucket_name = 'diagnostic-pipeline-output-prod-1' if profile == 'prod' else 'diagnostic-pipeline-output-pilot-1'

    bucket: Bucket = storage_client.bucket(bucket_name)
    blob: Blob = bucket.get_blob(f'{set_name}/purple/{sample_name}.driver.catalog.somatic.tsv')

    content = blob.download_as_string().decode()
    return content.find('AMP') > 0


def get_rose_location_as_gs_string(profile, set_name, sample_name):
    bucket_name = 'diagnostic-pipeline-output-prod-1' if profile == 'prod' else 'diagnostic-pipeline-output-pilot-1'
    return f'gs://{bucket_name}/{set_name}/purple/{sample_name}.driver.catalog.somatic.tsv'


def get_set_name_from_report(report, all_runs):
    """
    Given a report (as found in the report created endpoint), this method tries to find the
    sample set name associated with the report.

    In order to do this, it needs a list of runs to get the SOP string from.

    :param report: the report
    :param all_runs: all the runs
    :return: the sample_set name
    """
    return all_runs[all_runs['id'] == report['run_id']].iloc[:1]['set'].values[0]['name']


if __name__ == '__main__':
    main()
