import json
import argparse
from api_util import ApiUtil
from cli_util import perform_prod_test


def main(profile: str):
    api_util = ApiUtil(profile)
    all_reports = api_util.get_all_reports_created()
    all_shared_reports = {report['run_id'] for report in api_util.get_all_reports_shared()}
    reports_not_shared = [report for report in all_reports if report['run_id'] not in all_shared_reports]

    print('The following reports are created but not shared:')
    print(json.dumps(reports_not_shared, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', choices=['pilot', 'prod'], default='pilot')
    args = parser.parse_args()

    perform_prod_test(args.profile)
    main(args.profile)
