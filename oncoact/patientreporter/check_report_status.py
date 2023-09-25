import json
from api_util import get_all_reports_created, get_all_reports_shared


def main():
    all_reports = get_all_reports_created()
    all_shared_reports = {report['run_id'] for report in get_all_reports_shared()}
    reports_not_shared = [report for report in all_reports if report['run_id'] not in all_shared_reports]

    print("The following reports are created but not shared:")
    print(json.dumps(reports_not_shared, indent=2))


if __name__ == '__main__':
    main()
