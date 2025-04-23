import json
import requests
from datetime import datetime, timedelta


class RestClient:

    def __init__(self, profile):
        if profile == 'pilot':
            self._api_base_url = 'http://api.pilot-1'
            self._reporting_pipeline_url = "http://reporting-pipeline-launcher.pilot-1"
            self._lama_url = "http://lama.pilot-1"
        if profile == 'prod' or profile == 'preview':
            self._api_base_url = 'http://api.prod-1'
            self._reporting_pipeline_url = "http://reporting-pipeline-launcher.prod-1"
            self._lama_url = "http://lama.prod-1"

        self._report_created_url = f'{self._api_base_url}/hmf/v1/reports/created'
        self._report_shared_url = f'{self._api_base_url}/hmf/v1/reports/shared'
        self._sets_url = f'{self._api_base_url}/hmf/v1/sets'
        self._runs_url = f'{self._api_base_url}/hmf/v1/runs'
        self._files_url = f'{self._api_base_url}/hmf/v1/files'
        self._samples_url = f'{self._api_base_url}/hmf/v1/samples'

    def get_all_reports_created(self, lookback_days=90):
        """
        Queries the 'reports/created' endpoint and returns all results.

        :param lookback_days the amount of days to look back for the report creation time.
        """
        from_date = (datetime.today().date() - timedelta(days=lookback_days)).strftime('%Y-%m-%dT%H:%M:%S')
        return _get_as_json(url=self._report_created_url, params={'from_report_create_time': from_date})

    def get_report_created(self, sample_barcode):
        """
        Queries the 'reports/created' endpoint for the given sample_barcode.

        :param sample_barcode: the sample_barcode to query for.
        """
        response_json = _get_as_json(self._report_created_url, params={'sample_barcode': sample_barcode})
        if len(response_json) == 0:
            raise ValueError(f"No report created records found for sample_barcode: '{sample_barcode}'")
        if len(response_json) > 1:
            print(f"Report created endpoint returned more than one result")
            print('#', 'summary', sep='\t')
            for i, created in enumerate(response_json):
                print(i + 1, {
                    'id': created['id'],
                    'type': created['report_type'],
                    'create_time': created['create_time']
                }, sep='\t')
            to_return = input("Which one do you want to return? Please enter the #\n")
            return response_json[int(to_return) - 1]
        return response_json[0]

    def get_all_reports_shared(self, lookback_days=90):
        """
        Queries the 'reports/shared' endpoint and returns all results.

        :param lookback_days the amount of days to look back for the report creation time.
        """
        from_date = (datetime.today().date() - timedelta(days=lookback_days)).strftime('%Y-%m-%dT%H:%M:%S')
        return _get_as_json(url=self._report_shared_url, params={'from_report_create_time': from_date})

    def get_all_relevant_runs(self, lookback_days=90):
        """
        Queries the 'runs' endpoint for all diagnostic runs and returns somatic, CPCT and targeted inis.
        :param lookback_days the amount of days to look back for.
        """
        from_end_date = (datetime.today().date() - timedelta(days=lookback_days)).strftime('%Y-%m-%dT%H:%M:%S')

        res = []
        for ini_type in ['Somatic', 'Targeted']:
            res += _get_as_json(self._runs_url, params={'ini': f'{ini_type}.ini',
                                                        'context': 'DIAGNOSTIC',
                                                        'from_enddate': from_end_date})
        return res

    def get_run(self, run_id):
        """
        Gets the run by id.
        """
        return _get_as_json(f'{self._runs_url}/{run_id}')

    def get_yield(self, barcode):
        """
        Gets the sample information at sample barcode
        """
        sample = _get_as_json(self._samples_url, params={'barcode': barcode})
        print(sample)
        return sample['yld']

    def get_run_files(self, run_id):
        return _get_as_json(self._files_url, params={'run_id': run_id})

    def get_sample_barcode_from_isolation_barcode(self, tumor_isolation_barcode):
        """
        Gets the tumor sample barcode from the tumor isolation barcode.

        This method does not rely on the existence of a report. It uses Lama to retrieve the tumor sample barcode.
        """
        print(_get_as_json(f'{self._lama_url}/api/statuses/sample-barcode/{tumor_isolation_barcode}', params={'sampleBarcode': tumor_isolation_barcode}))
        return _get_as_json(f'{self._lama_url}/api/statuses/sample-barcode/{tumor_isolation_barcode}')['sampleBarcode']

    def get_tumor_sample_barcode_from_run_id(self, run_id):
        """
        For a given run_id, return the tumor sample barcode.

        """
        run = self.get_run(run_id)
        print(run)
        sample_set = self.get_sample_set_by_id(run['set']['id'])
        samples = sample_set['samples']
        tumor_samples = [s for s in samples if s['type'] == 'tumor']
        if len(tumor_samples) == 0:
            raise ValueError(f"No tumor sample found for run '{run_id}'")
        tumor_sample = tumor_samples[0]
        tumor_isolation_barcode = tumor_sample['barcode']
        print(tumor_isolation_barcode)
        return self.get_sample_barcode_from_isolation_barcode(tumor_isolation_barcode)

    def get_sample_set_by_sample_name(self, sample_name):
        """
        Queries the 'sets' endpoint with the given tumor sample name and returns the last entry.
        """
        response_json = _get_as_json(self._sets_url, params={'tumor_sample': sample_name})
        if len(response_json) == 0:
            raise ValueError(f"No sets found for sample_name '{sample_name}'")
        return response_json[-1]

    def get_sample_set_by_id(self, set_id):
        return _get_as_json(f'{self._sets_url}/{set_id}')

    def post_report_shared(self, report_created_id, notify_users, publish_to_portal):
        """
        Set a report as shared in the API.

        :param report_created_id: the reports created ID found at the created endpoint of the api.
        :param notify_users: boolean value that determines whether to notify the users that the report has been shared.
        :param publish_to_portal: boolean value that determines whether the report should be published to the portal.
        """

        body = {
            'report_created_id': report_created_id,
            'notify_users': notify_users,
            'publish_to_portal': publish_to_portal
        }
        response = requests.post(self._report_shared_url, json=body)
        response.raise_for_status()
        return response.json()

    def get_failed_executions(self):
        """
        Gets the failed executions from the reporting-pipeline.
        """
        json_response = _get_as_json(f'{self._reporting_pipeline_url}/executions', params={'success': 'false'})
        res = [execution for execution in
               json_response]
        return res

    def get_running_executions(self):
        json_response = _get_as_json(f'{self._reporting_pipeline_url}/executions')
        res = [execution for execution in json_response
               if execution['isRunning']]
        return res

    def delete_finished_execution(self, sample_barcode):
        return _delete_as_json(f'{self._reporting_pipeline_url}/report-executions', params={'sample_barcode': sample_barcode})

    def get_lama_patient_reporter_data(self, isolation_barcode):
        return _get_as_json(f"{self._lama_url}/api/queries/patient-reporter/isolation-barcode/{isolation_barcode}")


def _get_as_json(url, params=None):
    if params:
        response = requests.get(url, params=params)
    else:
        response = requests.get(url)
    response.raise_for_status()
    return response.json()


def _delete_as_json(url, params=None):
    if params:
        response = requests.delete(url, params=params)
    else:
        response = requests.delete(url)
    response.raise_for_status()
    return response.json()