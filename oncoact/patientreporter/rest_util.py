import json
from typing import List

import requests


class RestClient:

    def __init__(self, profile):
        self._api_base_url = 'https://api.prod-1'
        self._reporting_pipeline_url = ""
        self._lama_url = "https://lama.prod-1"

    def get_all_reports_created(self):
        """
        Queries the 'reports/2/created' endpoint and returns all results.

        :return: the query result.
        :raises ValueError: If the response returned a non 2XX status code.
        """
        response = requests.get(f'{self._api_base_url}/hmf/v1/reports/2/created')
        response.raise_for_status()
        json_response = response.json()
        return json_response

    def get_report_created(self, sample_barcode):
        """
        Queries the 'reports/2/created' endpoint with the given sample_barcode.

        :param sample_barcode: the sample_barcode to query for.
        :return: The query result.
        :raises ValueError: If the response returned a non 2XX status code or if the query returned more than 1 result.
        """
        response = requests.get(url=f'{self._api_base_url}/hmf/v1/reports/2/created',
                                params={'sample_barcode': sample_barcode})
        response.raise_for_status()
        response_json = response.json()
        if len(response_json) > 1:
            raise ValueError(f"Query returned more than one result: '{len(response_json)}'")
        return response_json[0]

    def get_all_reports_shared(self):
        """
        Queries the 'reports/2/shared' endpoint and returns all results.

        :param profile which profile to run this in (pilot/prod etc.)
        :return: the query result.
        :raises ValueError: If the response returned a non 2XX status code.
        """
        response = requests.get(f'{self._api_base_url}/hmf/v1/reports/2/shared')
        response.raise_for_status()
        json_response = response.json()
        return json_response

    def get_sample_set_by_sample_name(self, sample_name):
        """
        Queries the 'sets' endpoint with the given tumor sample name and returns the first entry.

        :param sample_name: the sample_name to query for
        :return: the first entry of the query result.
        :raises ValueError: If the response returned a non 2XX status code.
        """
        response = requests.get(url=f'{self._api_base_url}/hmf/v1/sets', params={'tumor_sample': sample_name})
        response.raise_for_status()
        response_json = response.json()
        return response_json[0]

    def get_sample_set_by_id(self, set_id):
        response = requests.get(url=f'{self._api_base_url}/hmf/v1/sets/{set_id}')
        response.raise_for_status()
        return response.json()

    def get_runs(self):
        """
        Queries the 'runs' endpoint and returns somatic and CPCT inis

        :return The query result
        :raises ValueError: If the response returned a non 2XX status code.
        """
        res = []
        for ini_type in ['CPCT', 'Somatic']:
            response = requests.get(url=f'{self._api_base_url}/hmf/v1/runs', params={'ini': f'{ini_type}.ini'})
            response.raise_for_status()
            res += response.json()
        return res

    def get_run(self, run_id):
        response = requests.get(url=f'{self._api_base_url}/hmf/v1/runs/{run_id}')
        response.raise_for_status()
        return response.json()

    def get_failed_executions(self):
        response = requests.get(f'{self._reporting_pipeline_url}/executions', params={'success': 'false'})
        response.raise_for_status()

        json_response: [json] = response.json()
        res = []

        for execution in json_response:
            res.append({
                "run_id": execution.runName,
                "stage_states": execution.stageStates
            })

        return res

    def get_tumor_sample_barcode(self, tumor_isolation_barcode):
        response = requests.get(f'{self._lama_url}/api/statuses/sample-barcode/{tumor_isolation_barcode}')
        response.raise_for_status()
        return response.json()['sampleBarcode']

    def get_tumor_sample_barcode_from_run_id(self, run_id):
        run = self.get_run(run_id)
        sample_set = self.get_sample_set_by_id(run['set']['id'])
        samples = sample_set['samples']
        tumor_samples = [s for s in samples if s['type'] == 'tumor']
        if len(tumor_samples) == 0:
            raise ValueError(f"No tumor sample found for run '{run_id}'")
        tumor_sample = tumor_samples[0]
        tumor_isolation_barcode = tumor_sample['barcode']

        self.get_tumor_sample_barcode(tumor_isolation_barcode)
