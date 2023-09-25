import json
import requests

API_BASE_URL = "http://api.prod-1"


def get_all_reports_created() -> list[json]:
    """
    Queries the 'reports/2/created' endpoint and returns all results.

    :return: the query result.
    :raises ValueError: If the response returned a non 2XX status code.
    """
    response = requests.get(f'{API_BASE_URL}/hmf/v1/reports/2/created')
    if not response.ok:
        raise ValueError(f"Response was not ok: {response.status_code}")
    json_response = response.json()
    return json_response


def get_report_created(sample_barcode: str) -> json:
    """
    Queries the 'reports/2/created' endpoint with the given sample_barcode.

    :param sample_barcode: the sample_barcode to query for.
    :return: The query result.
    :raises ValueError: If the response returned a non 2XX status code or if the query returned more than 1 result.
    """
    response = requests.get(url=f"{API_BASE_URL}/hmf/v1/reports/2/created", params={'sample_barcode': sample_barcode})
    if not response.ok:
        raise ValueError(f"Response was not ok: {response.status_code}")
    response_json = response.json()
    if len(response_json) > 1:
        raise ValueError(f"Query returned more than one result: '{len(response_json)}'")
    return response_json[0]


def get_all_reports_shared() -> list[json]:
    """
    Queries the 'reports/2/shared' endpoint and returns all results.

    :return: the query result.
    :raises ValueError: If the response returned a non 2XX status code.
    """
    response = requests.get(f'{API_BASE_URL}/hmf/v1/reports/2/shared')
    if not response.ok:
        raise ValueError(f"Response was not ok: {response.status_code}")
    json_response = response.json()
    return json_response


def get_set(sample_id: str) -> json:
    """
    Queries the 'sets' endpoint with the given sample id.

    :param sample_id: the sample_id to query for
    :return: the query result.
    :raises ValueError: If the response returned a non 2XX status code or if the query returned more than 1 result.
    """
    response = requests.get(url=f"{API_BASE_URL}/hmf/v1/sets", params={'sample_id': sample_id})
    if not response.ok:
        raise ValueError(f"Response was not ok: {response.status_code}")
    response_json = response.json()
    if len(response_json) > 1:
        raise ValueError(f"Query returned more than one result: '{len(response_json)}'")
    return response_json[0]
