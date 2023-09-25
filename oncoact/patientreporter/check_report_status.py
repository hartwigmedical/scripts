import requests
import json

API_BASE_URL = 'http://api.prod-1/hmf/v1'


def main():
    all_reports = get_all_reports()
    all_shared_reports = {report['run_id'] for report in get_all_shared_reports()}
    reports_not_shared = [report for report in all_reports if report['run_id'] not in all_shared_reports]

    print("The following reports are created but not shared:")
    print(json.dumps(reports_not_shared, indent=2))


def get_all_reports():
    response = requests.get(f'{API_BASE_URL}/reports/2/created')
    if not response.ok:
        raise ValueError(f'An error has occurred while querying the API: {response.status_code}')

    json_response = response.json()
    return json_response


def get_all_shared_reports():
    response = requests.get(f'{API_BASE_URL}/reports/2/shared')
    if not response.ok:
        raise ValueError(f'An error has occurred while querying the API: {response.status_code}')

    json_response = response.json()
    return json_response


if __name__ == '__main__':
    main()
