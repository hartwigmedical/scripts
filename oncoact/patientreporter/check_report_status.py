import json
import argparse
from api_util import ApiUtil


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

    if args.profile == 'prod':
        prod_warn = input("Warning: you are running in prod. Type 'y' to continue.")
        if prod_warn.lower() != 'y':
            print('Program aborted')
            exit(1)

    main(args.profile)
