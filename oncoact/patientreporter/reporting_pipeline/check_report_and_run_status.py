import argparse
import pandas as pd
import subprocess

from rest_util import RestClient, _get_as_json
from gsutil import get_bucket_and_blob_from_gs_path
from google.cloud.storage import Client, Bucket
import json

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', choices=['pilot', 'prod'], default='pilot')
    args = parser.parse_args()

    StatusChecker(profile=args.profile).generate_and_print_summary()


class StatusChecker:

    def __init__(self, profile):
        self.rest_client = RestClient(profile)
        self.storage_client = Client()

        self.oncoact_bucket: Bucket = self.storage_client.bucket(bucket_name="patient-reporter-final-prod-1")

        print("Gathering data")
        self.all_reports_with_null = ((pd.DataFrame(self.rest_client.get_all_reports_created()))
                                      .drop_duplicates(subset='id', keep='last'))
        self.all_reports = self.all_reports_with_null[self.all_reports_with_null['sample_barcode'].notnull()]

        if len(self.all_reports) == 0:
            self.all_reports = pd.DataFrame(columns=['sample_barcode', 'run_id'])
        self.shared_reports_with_null = (pd.DataFrame([shared['report_created']
                                                       for shared in self.rest_client.get_all_reports_shared()])
                                         .drop_duplicates(subset='id', keep='last'))
        self.shared_reports = self.shared_reports_with_null[self.shared_reports_with_null['sample_barcode'].notnull()]

        if len(self.shared_reports) == 0:
            self.shared_reports = pd.DataFrame(columns=['id', 'sample_barcode', 'run_id'])

        not_shared_reports = self.all_reports[~self.all_reports['id'].isin(self.shared_reports['id'])]
        self.to_be_shared_reports = not_shared_reports[
            ~not_shared_reports['run_id'].isin(self.shared_reports['run_id'])]
        self.to_be_shared_reports = self.to_be_shared_reports.sort_values(by='id', ascending=True).drop_duplicates(
            subset=['sample_barcode'], keep='last')

        self.all_runs = pd.DataFrame(self.rest_client.get_all_relevant_runs())

        if len(self.all_runs) == 0:
            self.all_runs = pd.DataFrame(columns=['status', 'id'])

        recent_raw = []
        for ini_type in ["Somatic", "Targeted"]:
            recent_raw += _get_as_json(
                self.rest_client._runs_url,
                params={
                    "ini": f"{ini_type}.ini",
                    "context": "DIAGNOSTIC"
                }
            )

        extra_df = pd.DataFrame(recent_raw)

        if not extra_df.empty:
            missing_processing = extra_df[extra_df["status"] == "Processing"]
            missing_processing = missing_processing[
                ~missing_processing["id"].isin(self.all_runs["id"])
            ]

            missing_waiting = extra_df[extra_df["status"] == "Waiting"]
            missing_waiting = missing_waiting[
                ~missing_waiting["id"].isin(self.all_runs["id"])
            ]

            to_add = pd.concat(
                [missing_processing, missing_waiting],
                ignore_index=True
            )

            if not to_add.empty:
                self.all_runs = pd.concat([self.all_runs, to_add], ignore_index=True)


        self.processing_runs = self.all_runs[self.all_runs['status'] == 'Processing']
        self.waiting_runs = self.all_runs[self.all_runs['status'] == 'Waiting']
        self.failed_runs = self.all_runs[self.all_runs['status'] == 'Failed']
        self.finished_runs = self.all_runs[self.all_runs['status'] == 'Finished']
        self.validated_runs = self.all_runs[self.all_runs['status'] == 'Validated']
        self.runs_without_report = self.all_runs[~self.all_runs['id'].isin(self.all_reports_with_null['run_id'])]

    def generate_and_print_summary(self):
        print("Generating report summary")
        chapters = [self._running_pipeline_chapter(),
                    self._failed_runs_chapter(),
                    self._finished_runs_chapter(),
                    self._validated_runs_chapter(),
                    ##self._reporting_pipeline_chapter()
                    ]

        _print_chapters(chapters)

    def _running_pipeline_chapter(self):
        print("Processing running pipeline chapter")
        running_pipeline_chapter = Chapter(name="Waiting & Processing runs")

        running_pipeline_chapter.add_section(self._waiting_runs_section())
        running_pipeline_chapter.add_section(self._processing_runs_section())

        return running_pipeline_chapter


    def _waiting_runs_section(self):
        section = Section(name='Waiting runs',
                          description='These runs are still waiting to start. '
                                      'Consider investigating these if they have not started on production days.')

        for _, run_record in self.waiting_runs.iterrows():
            run_content = _get_default_run_content(run_record)
            section.add_content(run_content)

        return section


    def _processing_runs_section(self):
        section = Section(name='Processing runs',
                          description='These runs are still processing. '
                                      'Consider investigating these is this is taking longer than usual.')

        for _, run_record in self.processing_runs.iterrows():
            run_content = _get_default_run_content(run_record)
            section.add_content(run_content)

        return section

    def _failed_runs_chapter(self):
        print("Processing failed runs chapter")
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

        not_shared_failed_reports = self.to_be_shared_reports[
            self.to_be_shared_reports['run_id'].isin(self.failed_runs['id'])]
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
        print("Processing finished runs chapter")
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

        not_shared_finished_reports = self.to_be_shared_reports[
            self.to_be_shared_reports['run_id'].isin(self.finished_runs['id'])]

        for _, report_record in not_shared_finished_reports.iterrows():
            section.add_content(self._get_default_report_and_run_content(report_record))
        return section

    def _finished_runs_with_report_shared_section(self):
        section = Section(name='Finished runs whose reports have been shared already',
                          description='WARNING: these reports have not been validated yet, '
                                      'but have already been shared!')
        finished_non_ini_runs = self.finished_runs[self.finished_runs['ini'] != 'Targeted.ini']

        shared_finished_reports = self.shared_reports[self.shared_reports['run_id'].isin(finished_non_ini_runs['id'])]

        for _, report_record in shared_finished_reports.iterrows():
            section.add_content(self._get_default_report_and_run_content(report_record))
        return section

    def _validated_runs_chapter(self):
        print("Processing validated runs chapter")
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
        not_shared_validated_reports = self.to_be_shared_reports[
            self.to_be_shared_reports['run_id'].isin(self.validated_runs['id'])]
        for _, report_record in not_shared_validated_reports.iterrows():

            warnings = self._get_all_report_associated_warnings(report_record)
            content = {**self._get_default_report_and_run_content(report_record),
                       'warnings': warnings}
            section.add_content(content)
        return section

    def _get_all_report_associated_warnings(self, report_record):
        return (self._get_patient_reporter_log_related_warnings(report_record) +
                self._get_health_checker_related_warnings(report_record) +
                self._get_virus_names(report_record) +
                self._get_doid_warnings(report_record) +
                self._get_rose_warnings(report_record) +
                self._get_protect_warnings(report_record) +
                self._get_virus_warnings(report_record))

    def _get_patient_reporter_log_related_warnings(self, report_record):
        warnings = []
        sample_barcode = report_record["sample_barcode"]
        path = f"gs://{self.oncoact_bucket.name}/{sample_barcode}/patient-reporter.log"
        if not path:
            return warnings
        _, log_blob = get_bucket_and_blob_from_gs_path(self.storage_client, path)

        if log_blob is None:
            warnings.append('The patient reporter log file was not found!')
            return warnings
        log_content = log_blob.download_as_string().decode()
        if 'WARN' in log_content:
            warnings.append('A warning was found in the patient-reporter log. '
                            f"Try running 'gsutil cat {path} | grep WARN' to find out more.")

        if ("Consent" in log_content or
                "Mismatching ref sample name" in log_content or
                "do not match" in log_content or
                "Missing or invalid hospital" in log_content):
            warnings.append('A lims error was found in the patient-reporter log. '
                            f"Check for lims problems 'gsutil cat {path}'.")

        return warnings

    def _get_virus_names(self, report_record):
        warnings = []
        #get reporting id
        used_lama_data = self._get_lama_data_used_for_report(report_record)
        if used_lama_data is None:
            warnings.append(f"No Lama data was found in the patient reporter for {report_record['barcode']}.")
            return warnings

        patient_id = "reportingId"
        patient_id_value = used_lama_data[patient_id] if patient_id in used_lama_data else None

        pathology_id = "hospitalSampleLabel"
        pathology_id_value = used_lama_data[pathology_id] if pathology_id in used_lama_data else None

        if pathology_id_value and patient_id_value:
            reporting_id = f"{patient_id_value}-{pathology_id_value}"
        elif patient_id_value:
            reporting_id = patient_id_value
        else:
            warnings.append(f"Both reportingId and hospitalSampleLabel are missing in lama data for {report_record['barcode']}.")
            return warnings

        #Get credentials for sql
        command = 'bash -i -c "source database_functions && get_secret_from_secret_manager mysql-diagnostic-patients-sql-prod-1-reader"'
        creds_result = subprocess.run(command, shell=True, text=True, capture_output=True)
        creds = creds_result.stdout.strip()
        #Construct SQL query
        sql_query = (
            "SELECT sampleId, virusName, qcStatus, integrations, interpretation, reported, likelihood "
            "FROM virusAnnotation "
            "WHERE reported =0  "
            f"AND sampleId = '{reporting_id}'"
        )
        query_command = f"do_execute_sql_on_database \"{sql_query}\" hmfpatients '{creds}'"
        output = subprocess.run(query_command, shell=True, text=True, capture_output=True)
        # Process the output
        if output.stdout.strip():
            lines = output.stdout.strip().split('\n')
            if len(lines) > 1:  # Ensure there are both headers and data
                headers = lines[0].split('\t')
                values = lines[1].split('\t')

                # Combine headers and values into key-value pairs
                key_value_pairs = [f"{header}: {value}" for header, value in zip(headers, values)]

                # Create the readable output
                readable_output = ", ".join(key_value_pairs)
                warnings.append(f"The unreported virus information info is as follows: {readable_output}")
            else:
                warnings.append("The unreported virus information info could not be retrieved or is incomplete.")

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

    def _get_doid_warnings(self, report_record):
        warnings = []

        used_lama_data = self._get_lama_data_used_for_report(report_record)
        current_lama_data = self._get_current_lama_data(report_record)

        if used_lama_data is None or current_lama_data is None:
            if used_lama_data is None:
                warnings.append("No Lama data was found in the patient reporter for this report. "
                                "This is most certainly a bug in the scripts or in the patient reporter.")
            if current_lama_data is None:
                warnings.append("No lama data was found in the API for this report. "
                                "This is most certainly a bug in the scripts or in LAMA.")
            return warnings
        tumor_type_key = "primaryTumorType"
        used_primary_tumor_type = used_lama_data[tumor_type_key] if tumor_type_key in used_lama_data else None
        current_primary_tumor_type = current_lama_data[tumor_type_key] if tumor_type_key in current_lama_data else None

        if used_primary_tumor_type is None or current_primary_tumor_type is None:
            if used_primary_tumor_type is None:
                warnings.append("The primary tumor typed used for generating this report was 'None', "
                                "No Doid information was available during the generation of this report.")
            if current_primary_tumor_type is None:
                warnings.append("The current primary tumor type in LAMA for this report is 'None', "
                                "No doid information is available in LAMA as of right now.")
            return warnings

        if used_primary_tumor_type["doids"] != current_primary_tumor_type["doids"]:
            warnings.append("There was a difference detected between the doid data used for creating this report "
                            "and the doid data stored in LAMA. "
                            f"Used doids: {used_primary_tumor_type['doids']} "
                            f"Doids in LAMA: {current_primary_tumor_type['doids']}")

        if used_primary_tumor_type["location"].lower() == "unknown":
            warnings.append("The primary tumor type value was 'unknown' "
                            "at the time the patient reporter generated this report.")
        if current_primary_tumor_type["location"].lower() == "unknown":
            warnings.append("The primary tumor type value is 'unknown' at this time in LAMA for this sample.")

        return warnings

    def _get_lama_data_used_for_report(self, report_record):
        sample_barcode = report_record["sample_barcode"]
        path = f"gs://{self.oncoact_bucket.name}/{sample_barcode}/lama/patient-reporter.json"

        _, blob = get_bucket_and_blob_from_gs_path(self.storage_client, path)

        if blob is None:
            return None
        return json.loads(blob.download_as_string().decode())

    def _get_current_lama_data(self, report_record):
        isolation_barcode = report_record['barcode']
        return self.rest_client.get_lama_patient_reporter_data(isolation_barcode)

    def _get_rose_warnings(self, report_record):
        warnings = []
        sample_barcode = report_record["sample_barcode"]
        path = f"gs://{self.oncoact_bucket.name}/{sample_barcode}/rose.log"
        if not path:
            return warnings
        _, log_blob = get_bucket_and_blob_from_gs_path(self.storage_client, path)

        if log_blob is None:
            return ["Rose log not found."]
        rose_log = log_blob.download_as_string().decode()
        if 'WARN' in rose_log:
            warnings.append(
                f"A warning was found in the rose log. Use 'gsutil cat {path}'")
        return warnings

    def _get_protect_warnings(self, report_record):
        warnings = []
        sample_barcode = report_record["sample_barcode"]
        path = f"gs://{self.oncoact_bucket.name}/{sample_barcode}/protect.log"

        if not path:
            return warnings
        _, log_blob = get_bucket_and_blob_from_gs_path(self.storage_client, path)

        if log_blob is None:
            return ["Protect log not found"]
        protect_log = log_blob.download_as_string().decode()
        if 'WARN' in protect_log:
            warnings.append(
                f"A warning was found in the protect log. Use 'gsutil cat ${path}'")
        return warnings

    def _get_virus_warnings(self, report_record):
        warnings = []
        #get reporting id
        used_lama_data = self._get_lama_data_used_for_report(report_record)

        if used_lama_data is None:
            warnings.append(f"No Lama data was found in the patient reporter for {report_record['barcode']}.")
            return warnings

        patient_id = "reportingId"
        patient_id_value = used_lama_data[patient_id] if patient_id in used_lama_data else None

        pathology_id = "hospitalSampleLabel"
        pathology_id_value = used_lama_data[pathology_id] if pathology_id in used_lama_data else None

        if pathology_id_value and patient_id_value:
            reporting_id = f"{patient_id_value}-{pathology_id_value}"
        elif patient_id_value:
            reporting_id = patient_id_value
        else:
            warnings.append(f"Both reportingId and hospitalSampleLabel are missing in lama data for {report_record['barcode']}.")
            return warnings

        #Get credentials for sql
        command = 'bash -i -c "source database_functions && get_secret_from_secret_manager mysql-diagnostic-patients-sql-prod-1-reader"'
        creds_result = subprocess.run(command, shell=True, text=True, capture_output=True)
        creds = creds_result.stdout.strip()
        #Construct SQL query
        sql_query = (
            "SELECT breakend.* "
            "FROM virusAnnotation as annotation"
            "INNER JOIN virusBreakend as breakend on breakend.sampleId = annotation.sampleId"
            "WHERE annotation.reported  "
            "AND nameSpecies = nameAssigned and nameSpecies = 'Alphapapillomavirus 7'"
            f"AND sampleId = '{reporting_id}'"
        )
        query_command = f"do_execute_sql_on_database \"{sql_query}\" hmfpatients '{creds}'"
        output = subprocess.run(query_command, shell=True, text=True, capture_output=True)
        # Process the output
        if output.stdout.strip():
            lines = output.stdout.strip().split('\n')
            if len(lines) > 1:  # Ensure there are both headers and data
                headers = lines[0].split('\t')
                values = lines[1].split('\t')

                # Combine headers and values into key-value pairs
                key_value_pairs = [f"{header}: {value}" for header, value in zip(headers, values)]

                # Create the readable output
                readable_output = ", ".join(key_value_pairs)
                warnings.append(f"The virus Alphapapillomavirus 7 information info is as follows: {readable_output}")
            else:
                warnings.append("The virus Alphapapillomavirus 7 information info could not be retrieved or is incomplete.")

        return warnings

    # def _reporting_pipeline_chapter(self):
    #     print("Processing reporting pipeline failures chapter")
    #     chapter = Chapter(name='Reporting pipeline fails')
    #     chapter.add_section(self._reporting_pipeline_fail_section())
    #     chapter.add_section(self._reporting_pipeline_in_progress_section())
    #
    #     return chapter

    def _reporting_pipeline_fail_section(self):
        section = Section(name='Runs with a patient reporter fail',
                          description='For these entries there is no report due to a patient report failure.')
        failed_executions = self.rest_client.get_failed_executions()
        for _, execution in enumerate(failed_executions):
            content = {
                **execution
            }
            section.add_content(content)
        return section

    def _reporting_pipeline_in_progress_section(self):
        section = Section(name="Reporting pipeline in progress",
                          description="For these reports, the reporting pipeline is running")
        running_executions = self.rest_client.get_running_executions()
        for _, execution in enumerate(running_executions):
            content = {
                **execution
            }
            section.add_content(content)
        return section


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
        'reportingId': report_created_record['reportingId'],
        'hospitalSampleLabel': report_created_record['hospitalSampleLabel'],
        'report_type': report_created_record['report_type']
    }

def _get_default_run_content(run_record):
    return {
        'run_id': run_record['id'],
        'run_ini': run_record['ini'],
        'run_finished_date': run_record['endTime'],
        'set_name': run_record['set']['name']
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
    # Kleuren per chapter
    color_map = {
        "Failed runs": "\033[91m",       # red
        "Finished runs": "\033[94m",     # blue
        "Validated runs": "\033[92m",    # green
        "Waiting & Processing runs": "\033[33m",  # yello
        "Reporting pipeline fails": "\033[38;5;208m",   # orange
    }

    reset = "\033[0m"
    bold = "\033[1m"

    for chapter in chapters:
        color = color_map.get(chapter.name, "\033[97m")

        # Chapter header
        print("\n" + "=" * 70)
        print(f"{color}{bold}{chapter.name.upper()}{reset}")
        print("=" * 70)

        # Sections
        for section in chapter.sections:
            print(f"\n{bold}{section.name.upper()}{reset}")
            if section.description:
                print(section.description)

            if section.is_empty():
                print("\n\t- This section is empty -")
                continue

            for i, content in enumerate(section.contents):
                print('',
                      f'\t{i + 1}.',
                      *[f'\t{k}: {v}' for (k, v) in content.items()],
                      sep='\n')


if __name__ == '__main__':
    main()