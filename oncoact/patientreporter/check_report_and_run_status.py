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

    StatusChecker(profile=args.profile).print_report_and_run_summary()


class StatusChecker:

    def __init__(self, profile):
        self.rest_client = RestClient(profile)
        self.storage_client = Client()

        self.all_reports = pd.DataFrame(self.rest_client.get_all_reports_created())
        # these if checks are to ensure that if the dataframe is empty,
        # it still has the required columns to prevent KeyErrors.
        if len(self.all_reports) == 0:
            self.all_reports = pd.DataFrame(columns=['sample_barcode', 'run_id'])

        self.shared_reports = pd.DataFrame([shared['report_created']
                                            for shared in self.rest_client.get_all_reports_shared()])
        if len(self.shared_reports) == 0:
            self.shared_reports = pd.DataFrame(columns=['id', 'sample_barcode', 'run_id'])
        self.not_shared_reports = self.all_reports[~self.all_reports['id'].isin(self.shared_reports['id'])]

        self.all_runs = pd.DataFrame(self.rest_client.get_diagnostic_somatic_cpct_runs())
        if len(self.all_runs) == 0:
            self.all_runs = pd.DataFrame(columns=['status', 'id'])

        self.failed_runs = self.all_runs[self.all_runs['status'] == 'Failed']
        self.finished_runs = self.all_runs[self.all_runs['status'] == 'Finished']
        self.validated_runs = self.all_runs[self.all_runs['status'] == 'Validated']
        self.runs_without_report = self.all_runs[~self.all_runs['id'].isin(self.all_reports['run_id'])]

    def print_report_and_run_summary(self):
        print("REPORT AND RUN SUMMARY")
        self._handle_failed_runs()
        self._handle_finished_runs()
        self._handle_validated_runs()

        self._handle_reporting_pipeline_failures()

    def _handle_failed_runs(self):
        print('\n---- Failed runs ----\n')

        self._handle_failed_runs_without_report()
        self._handle_failed_runs_with_report_not_shared()
        self._handle_failed_runs_with_report_shared()

    def _handle_failed_runs_without_report(self):
        failed_runs_without_report = self.runs_without_report[self.runs_without_report['status'] == 'Failed']
        if len(failed_runs_without_report) > 0:
            failed_runs_without_report = failed_runs_without_report.reset_index(drop=True)
            print('** FAILED RUNS WITHOUT A REPORT **')
            print('These are runs that have failed, and have no report generated for them.',
                  'Consider investigating the failure and/or generating a QC failed report.', sep='\n')
            for i, run_record in failed_runs_without_report.iterrows():
                output_data = {
                    'run_id': run_record['id'],
                    'context': run_record['context'],
                    'fail_reason': _get_fail_reason(run_record)
                }
                StatusChecker._print_entry(i, output_data)

    def _handle_failed_runs_with_report_not_shared(self):
        not_shared_failed_reports = self.not_shared_reports[
            self.not_shared_reports['run_id'].isin(self.failed_runs['id'])]
        if len(not_shared_failed_reports) > 0:
            not_shared_failed_reports.reset_index(drop=True)
            print("** FAILED RUNS WHOSE REPORT HAS NOT BEEN SHARED YET **")
            print("These failed runs have an generated report, but this report has not been shared yet.")
            for i, report_record in not_shared_failed_reports.iterrows():
                self._print_failed_report(index=i, report_record=report_record)

    def _handle_failed_runs_with_report_shared(self):
        shared_failed_reports = self.shared_reports[self.shared_reports['run_id'].isin(self.failed_runs['id'])]
        if len(shared_failed_reports) > 0:
            shared_failed_reports = shared_failed_reports.reset_index(drop=True)
            print("** FAILED RUNS WHOSE REPORTS HAVE BEEN SHARED **")
            print(
                "Typically, the API status for these kind of runs is set to validated. Here, this is not yet the case",
                "Consider setting the API to status 'validated'", sep='\n')
            for i, report_record in shared_failed_reports.iterrows():
                self._print_failed_report(index=i, report_record=report_record)

    def _print_failed_report(self, index, report_record):
        associated_run = self.failed_runs[self.failed_runs['id'] == report_record['run_id']]
        if len(associated_run) == 0:
            fail_reason = 'Unknown (run not present in API)'
        else:
            fail_reason = _get_fail_reason(associated_run.iloc[0])

        StatusChecker._print_entry(index, {'run_id': int(report_record['run_id']),
                                           'sample_barcode': report_record['sample_barcode'],
                                           'isolation_barcode': report_record['barcode'],
                                           'fail_reason': fail_reason})

    def _handle_finished_runs(self):
        print('\n---- Finished runs ----\n')
        self._handle_finished_runs_without_report()
        self._handle_finished_runs_with_report_not_shared()
        self._handle_finished_runs_with_report_shared()

    def _handle_finished_runs_without_report(self):
        finished_runs_without_report = self.runs_without_report[self.runs_without_report['status'] == 'Finished']
        if len(finished_runs_without_report) > 0:
            finished_runs_without_report = finished_runs_without_report.reset_index(drop=True)
            print("** FINISHED RUNS WITHOUT REPORT **")
            print("These pipeline runs are finished, but no report exists yet.",
                  "This could be because the patient reporter is still processing this run.",
                  sep='\n')
            for i, run_record in finished_runs_without_report.iterrows():
                print_body = {
                    'run_id': run_record['id'],
                    'context': run_record['context']
                }
                StatusChecker._print_entry(i, print_body)

    def _handle_finished_runs_with_report_not_shared(self):
        not_shared_finished_reports = self.not_shared_reports[
            self.not_shared_reports['run_id'].isin(self.finished_runs['id'])]
        if len(not_shared_finished_reports) > 0:
            not_shared_finished_reports = not_shared_finished_reports.reset_index(drop=True)
            print('** FINISHED RUNS WHOSE REPORTS HAVE NOT BEEN SHARED **')
            print('Typically, these runs have not been validated yet by the snp-check. '
                  'Hence no action has to be performed.')
            for i, report_record in not_shared_finished_reports.iterrows():
                self._print_finished_report(index=i, report_record=report_record)

    def _handle_finished_runs_with_report_shared(self):
        shared_finished_reports = self.shared_reports[self.shared_reports['run_id'].isin(self.finished_runs['id'])]
        if len(shared_finished_reports) > 0:
            shared_finished_reports = shared_finished_reports.reset_index(drop=True)
            print('** FINISHED RUNS WHOSE REPORTS HAVE ALREADY BEEN SHARED **')
            print('WARNING: these reports have not been validated yet, but have already been shared!')
            for i, report_record in shared_finished_reports.iterrows():
                self._print_finished_report(index=i, report_record=report_record)

    @staticmethod
    def _print_finished_report(index, report_record):
        StatusChecker._print_entry(index, {'sample barcode': report_record["sample_barcode"],
                                           'isolation barcode': report_record["barcode"]})

    def _handle_validated_runs(self):
        print('\n---- Validated runs ----\n')
        self._handle_validated_runs_without_report()
        self._handle_validated_runs_whith_report_not_shared()

    def _handle_validated_runs_without_report(self):
        validated_runs_without_report = self.runs_without_report[self.runs_without_report['status'] == 'Validated']
        if len(validated_runs_without_report) > 0:
            validated_runs_without_report = validated_runs_without_report.reset_index(drop=True)
            print('** VALIDATED RUNS WITHOUT A REPORT **')
            print('These runs are already validated, but have no report yet.',
                  'It could be that the reporting pipeline is still processing them.',
                  sep='\n')
            for i, run_record in validated_runs_without_report.iterrows():
                print_data = {
                    'run_id': run_record['id'],
                    'context': run_record['context']
                }
                StatusChecker._print_entry(i, print_data)

    def _handle_validated_runs_whith_report_not_shared(self):
        not_shared_validated_reports = self.not_shared_reports[
            self.not_shared_reports['run_id'].isin(self.validated_runs['id'])]
        if len(not_shared_validated_reports) > 0:
            not_shared_validated_reports = not_shared_validated_reports.reset_index(drop=True)
            print('** VALIDATED RUNS WHOSE REPORT IS HAS NOT YET BEEN SHARED **')
            print('These reports are (almost) ready to be shared, any pending warnings are displayed below.')
            for i, report_record in not_shared_validated_reports.iterrows():
                self._print_validated_report(index=i, report_record=report_record)

    def _print_validated_report(self, index, report_record):
        warnings = self._get_all_report_associated_warnings(report_record)
        StatusChecker._print_entry(index, {'sample barcode': report_record['sample_barcode'],
                                           'isolation barcode': report_record['barcode'],
                                           'warnings': warnings})

    def _get_all_report_associated_warnings(self, report_record):
        return (self._get_patient_reporter_log_related_warnings(report_record) +
                # self._get_health_checker_related_warnings(report_record) +
                self._get_rose_related_warnings(report_record))

    def _get_patient_reporter_log_related_warnings(self, report_record):
        warnings = []

        report_files = report_record['report_files']
        log_file = next((file for file in report_files if file['datatype'] == 'report_log'), None)
        if not log_file:
            return warnings
        bucket_name, blob_name = get_bucket_and_blob_names_from_gs_path(log_file['path'])
        log_content = self.storage_client.get_bucket(bucket_name).get_blob(blob_name).download_as_string().decode()
        if 'WARN' in log_content:
            warnings.append('A warning was found in the patient-reporter log. '
                            f"Try running 'gsutil cat {log_file}' | grep WARN' to find out more.")

        if ("Consent" in log_content or
                "Mismatching ref sample name" in log_content or
                "do not match" in log_content or
                "Missing or invalid hospital" in log_content):
            warnings.append('A lims error was found in the patient-reporter log. '
                            f"Check for lims problems 'gsutil cat {log_file}'.")

        return warnings

    def _get_health_checker_related_warnings(self, report_record):
        warnings = []
        associated_set_name = self._get_set_name_from_report(report_record)
        health_checker_log = _get_health_error_validated_run(associated_set_name)
        if 'WARN' in health_checker_log:
            warnings.append('A warning was found in the health checker log. '
                            f'See: {health_checker_log}')
        return warnings

    def _get_rose_related_warnings(self, report_record):
        warnings = []
        associated_run = self._get_associated_run(report_record=report_record)
        run_files = self.rest_client.get_run_files(run_id=associated_run['id'])

        purple_driver_catalog_file = next(
            (file for file in run_files if file['datatype'] == 'purple_somatic_driver_catalog'), None)
        if not purple_driver_catalog_file:
            return warnings
        path = purple_driver_catalog_file['filepath']
        bucket_name, blob_name = get_bucket_and_blob_names_from_gs_path(path)
        content = self.storage_client.bucket(bucket_name).get_blob(blob_name).download_as_string().decode()
        if content.find('AMP'):
            warnings.append(f"A rose related warning was found. Try running 'gsutil cat {path}'")
        return warnings

    def _get_associated_run(self, report_record):
        return self.all_runs[self.all_runs['id'] == report_record['run_id']]

    def _get_set_name_from_report(self, report_record):
        associated_run = self.all_runs[self.all_runs['id'] == report_record['run_id']]
        if len(associated_run) == 0:
            return None
        return associated_run.iloc[:1]['set'].values[0]['name']

    def _handle_reporting_pipeline_failures(self):
        failed_executions = self.rest_client.get_failed_executions()
        if len(failed_executions) > 0:
            print("** RUNS WITH A PATIENT REPORTER FAIL **")
            print("For these runs there is no report due to a patient report failure.")
            for i, (run_id, failure) in enumerate(failed_executions):
                if run_id in self.all_reports['run_id']:
                    continue
                sample_barcode = self.rest_client.get_tumor_sample_barcode_from_run_id(run_id)
                StatusChecker._print_entry(i, {'sample_barcode': sample_barcode,
                                               'failure': failure})

    @staticmethod
    def _print_entry(index, body: dict):
        print('',
              f'\t{index + 1}.',
              *[f'\t{k}: {v}' for (k, v) in body.items()],
              sep='\n')


def _get_fail_reason(run_record):
    return run_record['failure']  # .array[0]


def _get_health_error_validated_run(set_name: str) -> str:
    return subprocess.check_output(['health_check_validated_run', set_name, '2>&1']).decode()


class Chapter:

    def __init__(self, name):
        self.name = name
        self.sections = []

    def add_section(self, section):
        self.sections.append(section)


class Section:
    def __init__(self, name, content):
        self.name = name
        self.content = content


if __name__ == '__main__':
    main()
