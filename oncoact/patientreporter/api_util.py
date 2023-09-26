import json
from typing import List

import requests


class ApiUtil:

    def __init__(self, profile):
        self._profile = profile

    def api_base_url(self):
        if self._profile == 'pilot':
            return "http://api.pilot-1"
        if self._profile == 'prod':
            "http://api.prod-1"
        else:
            raise ValueError(f"Unknown profile: '{self._profile}'")

    def get_all_reports_created(self):
        """
        Queries the 'reports/2/created' endpoint and returns all results.

        :return: the query result.
        :raises ValueError: If the response returned a non 2XX status code.
        """
        response = requests.get(f'{self.api_base_url()}/hmf/v1/reports/2/created')
        if not response.ok:
            raise ValueError(f"Response was not ok: '{response.status_code}'")
        json_response = response.json()
        return json_response

    def get_report_created(self, sample_barcode):
        """
        Queries the 'reports/2/created' endpoint with the given sample_barcode.

        :param sample_barcode: the sample_barcode to query for.
        :return: The query result.
        :raises ValueError: If the response returned a non 2XX status code or if the query returned more than 1 result.
        """
        response = requests.get(url=f'{self.api_base_url()}/hmf/v1/reports/2/created',
                                params={'sample_barcode': sample_barcode})
        if not response.ok:
            raise ValueError(f"Response was not ok: '{response.status_code}'")
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
        response = requests.get(f'{self.api_base_url()}/hmf/v1/reports/2/shared')
        if not response.ok:
            raise ValueError(f"Response was not ok: '{response.status_code}'")
        json_response = response.json()
        return json_response

    def get_sample_set(self, sample_name):
        """
        Queries the 'sets' endpoint with the given sample id and returns the first entry.

        :param sample_name: the sample_name to query for
        :return: the query result.
        :raises ValueError: If the response returned a non 2XX status code.
        """
        response = requests.get(url=f'{self.api_base_url()}/hmf/v1/sets', params={'tumor_sample': sample_name})
        if not response.ok:
            raise ValueError(f"Response was not ok: '{response.status_code}'")
        response_json = response.json()
        return response_json[0]
