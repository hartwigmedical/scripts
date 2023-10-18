import argparse
import pandas as pd
import subprocess

from rest_util import RestClient
from gsutil import get_bucket_and_blob_names_from_gs_path
from google.cloud.storage import Client


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', choices=['pilot', 'prod'], default='pilot')
    args = parser.parse_args()

    StatusChecker(profile=args.profile).generate_and_print_summary()


class StatusChecker:

    def __init__(self, profile):
        self.rest_client = RestClient(profile)
        self.storage_client = Client()

        self.all_reports = ((pd.DataFrame(self.rest_client.get_all_reports_created()))
                            .drop_duplicates(subset='id', keep='last'))

        # these if checks are to ensure that if the dataframe is empty,
        # it still has the required columns to prevent KeyErrors.
        if len(self.all_reports) == 0:
            self.all_reports = pd.DataFrame(columns=['sample_barcode', 'run_id'])

        self.shared_reports = (pd.DataFrame([shared['report_created']
                                             for shared in self.rest_client.get_all_reports_shared()])
                               .drop_duplicates(subset='id', keep='last'))
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

    def generate_and_print_summary(self):
        chapters = [self._failed_runs_chapter(),
                    self._finished_runs_chapter(),
                    self._validated_runs_chapter(),
                    self._reporting_pipeline_failures_chapter(),
                    self._waiting_pending_processing_runs_chapter()]

        _print_chapters(chapters)

    def _failed_runs_chapter(self):
        failed_runs_chapter = Chapter(name="Failed runs")
        failed_runs_chapter.add_section(self._failed_runs_without_report_section())
        failed_runs_chapter.add_section(self._failed_runs_with_report_not_shared_section())
        failed_runs_chapter.add_section(self._failed_runs_with_report_shared_section())

        return failed_runs_chapter

    def _failed_runs_without_report_section(self):
        section = Section(name='Failed runs without a report',
                          description='These are runs that have failed, and have no report generated for them. '
                                      'Consider investigating the failure and/or generating a QC failed report.')

        failed_runs_without_report = self.runs_without_report[
            self.runs_without_report['status'] == 'Failed']
        for _, run_record in failed_runs_without_report.iterrows():
            run_content = _get_default_run_content(run_record)
            run_content['fail_reason'] = run_record['failure']
            section.add_content(run_content)
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
        run_record = self._get_run_from_report_record(report_record)
        fail_reason = run_record['failure']

        run_content = _get_default_run_content(run_record)
        run_content['fail_reason'] = fail_reason
        report_content = _get_default_report_content(report_record)
        return {**run_content,
                **report_content}

    def _finished_runs_chapter(self):
        chapter = Chapter(name="Finished runs")

        chapter.add_section(self._finished_runs_without_report_section())
        chapter.add_section(self._finished_runs_with_report_not_shared_section())
        chapter.add_section(self._finished_runs_with_report_shared_section())
        return chapter

    def _finished_runs_without_report_section(self):
        section = Section(name="Finished runs without a report",
                          description="These molecular pipeline runs are finished, but no report exists yet. "
                                      "This could be because the patient reporter is still processing this run.")

        finished_runs_without_report = self.runs_without_report[self.runs_without_report['status'] == 'Finished']
        for _, run_record in finished_runs_without_report.iterrows():
            content = _get_default_run_content(run_record)
            section.add_content(content)
        return section

    def _finished_runs_with_report_not_shared_section(self):
        section = Section(name='Finished runs whose reports have not been shared',
                          description='Typically, these runs have not been validated yet by the snp-check. '
                                      'Hence no action has to be performed.')

        not_shared_finished_reports = self.not_shared_reports[
            self.not_shared_reports['run_id'].isin(self.finished_runs['id'])]

        for _, report_record in not_shared_finished_reports.iterrows():
            section.add_content(self._get_default_report_and_run_content(report_record))
        return section

    def _finished_runs_with_report_shared_section(self):
        section = Section(name='Finished runs whose reports have been shared already',
                          description='WARNING: these reports have not been validated yet, '
                                      'but have already been shared!')
        shared_finished_reports = self.shared_reports[self.shared_reports['run_id'].isin(self.finished_runs['id'])]

        for _, report_record in shared_finished_reports.iterrows():
            section.add_content(self._get_default_report_and_run_content(report_record))
        return section

    def _validated_runs_chapter(self):
        chapter = Chapter(name='Validated runs')
        chapter.add_section(self._validated_runs_without_report_section())
        chapter.add_section(self._validated_runs_with_report_not_shared_section())
        return chapter

    def _validated_runs_without_report_section(self):
        section = Section(name='Validated runs without a report',
                          description='These runs are already validated, but have no report yet. '
                                      'It could be that the reporting pipeline is still processing them.')

        validated_runs_without_report = self.runs_without_report[self.runs_without_report['status'] == 'Validated']

        for _, run_record in validated_runs_without_report.iterrows():
            content = _get_default_run_content(run_record)
            section.add_content(content)
        return section

    def _validated_runs_with_report_not_shared_section(self):
        section = Section(name="Validated runs whose report has not been shared yet",
                          description='These reports are (almost) ready to be shared,'
                                      ' any pending warnings are displayed below.')

        not_shared_validated_reports = self.not_shared_reports[
            self.not_shared_reports['run_id'].isin(self.validated_runs['id'])]
        for _, report_record in not_shared_validated_reports.iterrows():
            warnings = self._get_all_report_associated_warnings(report_record)
            content = {**self._get_default_report_and_run_content(report_record),
                       'warnings': warnings}
            section.add_content(content)
        return section

    def _get_all_report_associated_warnings(self, report_record):
        return (self._get_patient_reporter_log_related_warnings(report_record) +
                self._get_health_checker_related_warnings(report_record) +
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
        run_record = self._get_run_from_report_record(report_record)
        set_name = run_record['set']['name']

        health_checker_log = subprocess.check_output(['health_check_validated_run', set_name, '2>&1']).decode()
        if 'WARN' in health_checker_log:
            warnings.append('A warning was found in the health checker log. '
                            f'See: {health_checker_log}')
        return warnings

    def _get_rose_related_warnings(self, report_record):
        warnings = []
        run_record = self._get_run_from_report_record(report_record=report_record)
        run_files = self.rest_client.get_run_files(run_id=run_record['id'])

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

    def _waiting_pending_processing_runs_chapter(self):
        chapter = Chapter(name='Waiting, pending and processing runs')
        for status in {'Waiting', 'Pending', 'Processing'}:
            chapter.add_section(self._simple_run_status_section(status))
        return chapter

    def _simple_run_status_section(self, status):
        section = Section(name=f"{status} runs")
        runs = self.all_runs[self.all_runs['status'] == status]
        for _, run_record in runs.iterrows():
            section.add_content(_get_default_run_content(run_record))
        return section

    def _get_run_from_report_record(self, report_record):
        filtered = self.all_runs[self.all_runs['id'] == report_record['run_id']]
        if len(filtered) == 0:
            return None
        return self.all_runs[self.all_runs['id'] == report_record['run_id']].iloc[0]

    def _get_default_report_and_run_content(self, report_record):
        run_record = self._get_run_from_report_record(report_record)
        if run_record is not None:
            return {**_get_default_run_content(run_record),
                    **_get_default_report_content(report_record)}
        return {'run_id': report_record['run_id'],
                **_get_default_report_content(report_record)}


def _get_default_report_content(report_created_record):
    return {
        'report_created_id': report_created_record['id'],
        'sample_barcode': report_created_record['sample_barcode'],
        'isolation_barcode': report_created_record['barcode'],
        'report_type': report_created_record['report_type']
    }


def _get_default_run_content(run_record):
    return {
        'run_id': run_record['id'],
        'run_finished_date': run_record['endTime']
    }


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

    def is_empty(self):
        return len(self.contents) == 0


def _print_chapters(chapters):
    for chapter in chapters:
        print(f'\n\n--- {chapter.name} ---')
        for section in chapter.sections:
            print(f'\n** {section.name.upper()} **')
            if section.description:
                print(section.description)
            if section.is_empty():
                print('\n\t- This section is empty -')
                continue
            for i, content in enumerate(section.contents):
                print('',
                      f'\t{i + 1}.',
                      *[f'\t{k}: {v}' for (k, v) in content.items()],
                      sep='\n')


if __name__ == '__main__':
    main()
