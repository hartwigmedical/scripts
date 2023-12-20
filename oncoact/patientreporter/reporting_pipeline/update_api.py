import argparse
from rest_util import RestClient
from google.cloud.storage import Bucket, Client

def main():
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument('sample_barcode')
    argument_parser.add_argument('--profile', choices=['pilot', 'preview', 'prod'], default='pilot')
    argument_parser.add_argument('--publish', default=False, action='store_true',
                                 help='whether to publish to portal or not')
    argument_parser.add_argument('--notify-users', default=False, action='store_true',
                                 help='whether to email the users of the share event')
    args = argument_parser.parse_args()

    ReportSharer(sample_barcode=args.sample_barcode, profile=args.profile).share_report(publish_to_portal=args.publish,
                                                                                        notify_users=args.notify_users)

class ReportSharer:

    def __init__(self, sample_barcode, profile='pilot'):
        self.sample_barcode = sample_barcode
        self.rest_client: RestClient = RestClient(profile)
        self.storage_client: Client = Client()

        self.pipeline_output_bucket: Bucket = self.storage_client.bucket(f"diagnostic-pipeline-output-{profile}-1")

        self.archive_bucket: Bucket = self.storage_client.bucket(f'patient-reporter-final-{profile}-1')
        self.portal_bucket: Bucket = self.storage_client.bucket(f'hmf-customer-portal-report-shared-{profile}')
        self.panel_pipeline_output_bucket: Bucket = self.storage_client.bucket(f'targeted-pipeline-output-{profile}-1')
        self.panel_share_bucket: Bucket = self.storage_client.bucket(f'oncoact-panel-files-nki')

        self.report_created_record = self.rest_client.get_report_created(self.sample_barcode)
        self.run_id = self.report_created_record['run_id']
        self.run = self.rest_client.get_run(self.run_id) if self.run_id else None

    def share_report(self, publish_to_portal, notify_users):
        """
        Shares the report by updating the remote portal bucket and the remote archive bucket. It also updates the API.

        :param publish_to_portal: whether to publish a pub/sub message to portal.
        :param notify_users: setting this to true will email the users regarding the portal update.
        """

        print('Updating api')
        response = self.rest_client.post_report_shared(report_created_id=self.report_created_record['id'],
                                                       publish_to_portal=publish_to_portal,
                                                       notify_users=notify_users)
        print("API response:", response)
        print("Done!")

if __name__ == "__main__":
    main()