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

    StatusChecker(profile=args.profile).generate_and_print_report_summary()


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

    def generate_and_print_report_summary(self):
        chapters = [self._failed_runs_chapter(),
                    self._finished_runs_chapter(),
                    self._validated_runs_chapter(),
                    self._reporting_pipeline_failures_chapter()]

        _print_chapters(chapters)

    def _failed_runs_chapter(self):
        failed_runs_chapter = Chapter(name="Failed runs")
        failed_runs_chapter.add_section(self._failed_runs_without_report_section())
        failed_runs_chapter.add_section(self._failed_runs_with_report_not_shared_section())
        failed_runs_chapter.add_section(self._failed_runs_with_report_shared_section())

        return failed_runs_chapter

    def _finished_runs_chapter(self):
        chapter = Chapter(name="Finished runs")

        chapter.add_section(self._finished_runs_without_report_section())
        chapter.add_section(self._finished_runs_with_report_not_shared_section())
        chapter.add_section(self._finished_runs_with_report_shared_section())
        return chapter

    def _validated_runs_chapter(self):
        chapter = Chapter(name='Validated runs')
        chapter.add_section(self._validated_runs_without_report_section())
        chapter.add_section(self._validated_runs_with_report_not_shared_section())
        return chapter

    def _reporting_pipeline_failures_chapter(self):
        chapter = Chapter(name='Reporting pipeline fails')
        section = Section(name='Runs with a patient reporter fail',
                          description='For these runs there is no report due to a patient report failure.')
        failed_executions = self.rest_client.get_failed_executions()

        for _, (run_id, failure) in enumerate(failed_executions):
            sample_barcode = self.rest_client.get_tumor_sample_barcode_from_run_id(run_id)
            section.add_content({'sample_barcode': sample_barcode,
                                 'failure': failure})
        chapter.add_section(section)
        return chapter

    def _failed_runs_without_report_section(self):
        section = Section(name='Failed runs without a report',
                          description='These are runs that have failed, and have no report generated for them. '
                                      'Consider investigating the failure and/or generating a QC failed report.')

        failed_runs_without_report = self.runs_without_report[
            self.runs_without_report['status'] == 'Failed']
        for _, run_record in failed_runs_without_report.iterrows():
            section.add_content({
                'run_id': run_record['id'],
                'context': run_record['context'],
                'fail_reason': _get_fail_reason(run_record)
            })
        return section

    def _failed_runs_with_report_not_shared_section(self):
        section = Section(name="Failed runs whose report has not been shared",
                          description="These failed runs have an generated report, "
                                      "but this report has not been shared yet.")

        not_shared_failed_reports = self.not_shared_reports[
            self.not_shared_reports['run_id'].isin(self.failed_runs['id'])]
        for _, report_record in not_shared_failed_reports.iterrows():
            content = self._generate_failed_report_content(report_record=report_record)
            section.add_content(content)
        return section

    def _failed_runs_with_report_shared_section(self):
        section = Section(name="Failed runs whose reports have been shared",
                          description="Typically, the API status for these kind of runs is set to validated. "
                                      "Here, this is not yet the case Consider setting the API to status 'validated'")

        shared_failed_reports = self.shared_reports[self.shared_reports['run_id'].isin(self.failed_runs['id'])]
        for _, report_record in shared_failed_reports.iterrows():
            content = self._generate_failed_report_content(report_record=report_record)
            section.add_content(content)
        return section

    def _generate_failed_report_content(self, report_record):
        associated_run = self.failed_runs[self.failed_runs['id'] == report_record['run_id']]
        if len(associated_run) == 0:
            fail_reason = 'Unknown (run not present in API)'
        else:
            fail_reason = _get_fail_reason(associated_run.iloc[0])

        return {'run_id': int(report_record['run_id']),
                'sample_barcode': report_record['sample_barcode'],
                'isolation_barcode': report_record['barcode'],
                'fail_reason': fail_reason}

    def _finished_runs_without_report_section(self):
        section = Section(name="Finished runs without a report",
                          description="These pipeline runs are finished, but no report exists yet. "
                                      "This could be because the patient reporter is still processing this run.")

        finished_runs_without_report = self.runs_without_report[self.runs_without_report['status'] == 'Finished']
        for _, run_record in finished_runs_without_report.iterrows():
            content = {
                'run_id': run_record['id'],
                'context': run_record['context']
            }
            section.add_content(content)
        return section

    def _finished_runs_with_report_not_shared_section(self):
        section = Section(name='Finished runs whose reports have not been shared',
                          description='Typically, these runs have not been validated yet by the snp-check. '
                                      'Hence no action has to be performed.')

        not_shared_finished_reports = self.not_shared_reports[
            self.not_shared_reports['run_id'].isin(self.finished_runs['id'])]

        for _, report_record in not_shared_finished_reports.iterrows():
            content = {'sample barcode': report_record["sample_barcode"],
                       'isolation barcode': report_record["barcode"]}
            section.add_content(content)
        return section

    def _finished_runs_with_report_shared_section(self):
        section = Section(name='Finished runs whose reports have been shared already',
                          description='WARNING: these reports have not been validated yet, '
                                      'but have already been shared!')
        shared_finished_reports = self.shared_reports[self.shared_reports['run_id'].isin(self.finished_runs['id'])]

        for _, report_record in shared_finished_reports.iterrows():
            content = {'sample barcode': report_record["sample_barcode"],
                       'isolation barcode': report_record["barcode"]}
            section.add_content(content)
        return section

    def _validated_runs_without_report_section(self):
        section = Section(name='Validated runs without a report',
                          description='These runs are already validated, but have no report yet. '
                                      'It could be that the reporting pipeline is still processing them.')

        validated_runs_without_report = self.runs_without_report[self.runs_without_report['status'] == 'Validated']

        for _, run_record in validated_runs_without_report.iterrows():
            content = {
                'run_id': run_record['id'],
                'context': run_record['context']
            }
            section.add_content(content)
        return section

    def _validated_runs_with_report_not_shared_section(self):
        section = Section(name="Validated runs whose report has not been shared yet",
                          description='These reports are (almost) ready to be shared,'
                                      ' any pending warnings are displayed below.')

        not_shared_validated_reports = self.not_shared_reports[
            self.not_shared_reports['run_id'].isin(self.validated_runs['id'])]
        for _, report_record in not_shared_validated_reports.iterrows():
            section.add_content(self._print_validated_report(report_record))
        return section

    def _print_validated_report(self, report_record):
        warnings = self._get_all_report_associated_warnings(report_record)
        return {'sample barcode': report_record['sample_barcode'],
                'isolation barcode': report_record['barcode'],
                'warnings': warnings}

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


def _get_fail_reason(run_record):
    return run_record['failure']


def _get_health_error_validated_run(set_name: str) -> str:
    return subprocess.check_output(['health_check_validated_run', set_name, '2>&1']).decode()


class Chapter:
    def __init__(self, name):
        self.name = name
        self.sections = []

    def add_section(self, section):
        self.sections.append(section)


class Section:
    def __init__(self, name, description=None):
        self.name = name
        self.description = description
        self.contents = []

    def add_content(self, to_add):
        self.contents.append(to_add)


def _print_chapters(chapters):
    for chapter in chapters:
        print(f'\n\n--- {chapter.name} ---')
        for section in chapter.sections:
            print(f'\n** {section.name.upper()} **')
            print(section.description)
            for i, content in enumerate(section.contents):
                print('',
                      f'\t{i + 1}.',
                      *[f'\t{k}: {v}' for (k, v) in content.items()],
                      sep='\n')


if __name__ == '__main__':
    main()
