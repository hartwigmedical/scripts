import argparse
import pandas as pd

from api_util import ApiUtil


class ReportWithStatus:
    def __init__(self, report, status):
        self.report = report
        self.status = status


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', choices=['pilot', 'prod'], default='pilot')
    args = parser.parse_args()
    check_report_status(args.profile)


def check_report_status(profile: str):
    api_util = ApiUtil(profile)

    all_reports = pd.DataFrame(api_util.get_all_reports_created())
    # these if checks are to ensure that if the dataframe is empty,
    # it still has the required columns to prevent KeyErrors.
    if len(all_reports) == 0:
        all_reports = pd.DataFrame(columns=['sample_barcode', 'run_id'])
    all_shared_reports = pd.DataFrame([shared['report_created'] for shared in api_util.get_all_reports_shared()])
    if len(all_shared_reports) == 0:
        all_shared_reports = pd.DataFrame(columns=['sample_barcode', 'run_id'])
    all_runs = pd.DataFrame(api_util.get_runs())
    if len(all_runs) == 0:
        all_runs = pd.DataFrame(columns=['status', 'id'])

    validated_runs = all_runs[all_runs['status'] == 'Validated']
    finished_runs = all_runs[all_runs['status'] == 'Finished']
    failed_runs = all_runs[all_runs['status'] == 'Failed']

    reports_not_shared = all_reports[~all_reports['sample_barcode'].isin(all_shared_reports['sample_barcode'])]

    finished_reports_but_shared = all_shared_reports[all_shared_reports['run_id'].isin(finished_runs)]
    failed_report_shared = all_shared_reports[all_shared_reports['run_id'].isin(failed_runs)]

    validated_reports_but_not_shared = reports_not_shared[reports_not_shared['run_id'].isin(validated_runs['id'])]

    finished_reports_not_shared = reports_not_shared[reports_not_shared['run_id'].isin(finished_runs)]

    failed_report_but_not_shared = reports_not_shared[reports_not_shared['run_id'].isin(failed_runs)]

    runs_without_report = all_runs[~all_runs['id'].isin(all_reports['run_id'])]
    reports_without_runs = all_reports[~all_reports['run_id'].isin(all_runs['id'])]

    print('----------')
    print("Failed reports that have not yet been shared:\n")
    for i, report in enumerate(failed_report_but_not_shared['sample_barcode']):
        print(f"{i} - {report}")

    print('----------')
    print("Failed reports that have been shared:\n")
    for i, report in enumerate(failed_report_shared['sample_barcode']):
        print(f"{i} - {report}")

    print('----------')
    print("Validated reports that have not yet been shared:\n")
    for i, report in enumerate(validated_reports_but_not_shared['sample_barcode']):
        print(f"{i} - {report}")

    print('----------')
    print("Finished reports that have not yet been shared:\n")
    for i, report in enumerate(finished_reports_not_shared['sample_barcode']):
        print(f"{i} - {report}")

    print('----------')
    print("Finished reports that have already been shared (NOT VALIDATED YET!):\n")
    for i, report in enumerate(finished_reports_but_shared['sample_barcode']):
        print(f"{i} - {report}")

    print('----------')
    print("Reports without a run:\n")
    for i, report in enumerate(reports_without_runs['sample_barcode']):
        print(f"{i} - {report}")


if __name__ == '__main__':
    main()
