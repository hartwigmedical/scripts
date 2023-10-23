import argparse
from rest_util import RestClient
from google.cloud.storage import Bucket, Client
from gsutil import get_bucket_and_blob_from_gs_path, get_file_name_from_blob_name


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

        self.archive_bucket: Bucket = self.storage_client.bucket(f'patient-reporter-final-{profile}-1')
        self.portal_bucket: Bucket = self.storage_client.bucket(f'hmf-customer-portal-report-shared-{profile}')
        self.panel_pipeline_output_bucket: Bucket = self.storage_client.bucket(f'targeted-pipeline-output-{profile}-1')

        self.report_created_record = self.rest_client.get_report_created(self.sample_barcode)
        self.run_id = self.report_created_record['run_id']
        self.run = self.rest_client.get_run(self.run_id) if self.run_id else None

    def share_report(self, publish_to_portal, notify_users):
        """
        Shares the report by updating the remote portal bucket and the remote archive bucket. It also updates the API.

        :param publish_to_portal: whether to publish a pub/sub message to portal.
        :param notify_users: setting this to true will email the users regarding the portal update.
        """
        if not self.run:
            self._prompt_user_no_run()
        elif self.run['status'] != 'Validated':
            self._prompt_user_non_validated_run()

        self._delete_old_artifacts_in_portal_bucket()
        if publish_to_portal:
            self._copy_files_to_remote_buckets_publish()
        else:
            self._copy_files_to_remote_buckets_no_publish()

        print('Updating api')
        response = self.rest_client.post_report_shared(report_created_id=self.report_created_record['id'],
                                                       publish_to_portal=publish_to_portal,
                                                       notify_users=notify_users)
        print("API response:", response)
        print("Done!")

    def _prompt_user_no_run(self):
        cont = input(
            f"No associated run found for tumor barcode '{self.sample_barcode}'. "
            f"Are you sure you want to continue? (y/n)\n")
        if cont.lower() != 'y':
            exit(1)

    def _prompt_user_non_validated_run(self):
        cont = input(f"Run status for tumor barcode '{self.sample_barcode}' is not yet 'Validated' "
                     f"(actual status: {self.run['status']}). Are you sure you want to continue? (y/n)\n")
        if cont.lower() != 'y':
            exit(1)

    def _delete_old_artifacts_in_portal_bucket(self):
        print(f"deleting old artifacts from bucket '{self.portal_bucket.name}'")
        blobs_old_run = list(self.portal_bucket.list_blobs(prefix=self.sample_barcode))
        print(*[blob.name for blob in blobs_old_run], sep='\n')
        self.portal_bucket.delete_blobs(blobs=blobs_old_run)

    def _copy_files_to_remote_buckets_publish(self):
        run_blobs = self._get_run_files_as_blobs()
        report_blobs = self._get_report_files_as_blobs()

        print(f"Copying a total of '{len(run_blobs)}' run files to remote buckets")
        for blob in run_blobs:
            self._copy_blob_to_portal_bucket(blob=blob, target_sub_folder='RUO')

        print(f"Copying a total of '{len(report_blobs)}' report files to remote buckets")
        for blob in report_blobs:
            self._copy_blob_to_portal_bucket(blob=blob, target_sub_folder='')
            self._copy_blob_to_archive_bucket(blob=blob)

    def _copy_files_to_remote_buckets_no_publish(self):
        report_blobs = self._get_report_files_as_blobs()
        print(f"Copying a total of '{len(report_blobs)}' report files to remote buckets")
        for blob in report_blobs:
            self._copy_blob_to_archive_bucket(blob=blob)

    def _get_run_files_as_blobs(self):
        if self.run is None:  # if there is no run there are also no run files to return.
            return []
        if self.report_created_record['report_type'] == 'oncopanel_result_report':
            return self._get_targeted_run_files_as_blobs()
        else:
            return self._get_wgs_run_files_as_blobs()

    def _get_wgs_run_files_as_blobs(self):
        all_run_files = self.rest_client.get_run_files(self.run_id)
        run_file_types = {'purple_somatic_driver_catalog',
                          'linx_driver_catalog',
                          'somatic_variants_purple',
                          'structural_variants_purple',
                          'linx_fusions',
                          'orange_output_pdf'}
        run_files_to_upload = [file for file in all_run_files if file['datatype'] in run_file_types]
        result = []
        for file in run_files_to_upload:
            _, blob = get_bucket_and_blob_from_gs_path(storage_client=self.storage_client, gs_path=file['filepath'])
            result.append(blob)
        return result

    def _get_targeted_run_files_as_blobs(self):
        """
        Targeted run files currently are NOT stored in the API because the archiver doesn't run.
        Hence, we have to hardcode almost all panel blob-names to retrieve them.
        """
        panel_file_suffixes = {
            "driver.catalog.somatic.tsv",
            "purple.somatic.vcf.gz"
            "purple.purity.tsv",
            "purple.cnv.gene.tsv",
            "reported.somatic.vcf",
            "sampleQcReport.pdf",
            "sampleQcReport.png",
            "sampleQcDeamination.pdf",
            "sampleQcDeamination.png",
            "reported.somatic.vcf"
        }
        result = []
        set_name = self.run['set']['name']
        blobs = list(self.panel_pipeline_output_bucket.list_blobs(prefix=set_name))
        for suffix in panel_file_suffixes:
            for blob in blobs:
                if blob.name[-len(suffix):] == suffix:  # this checks if the blob name ends with the suffix.
                    result.append(blob)
                    break

        return result

    def _get_report_files_as_blobs(self):
        all_report_files = self.report_created_record["report_files"]
        report_files_to_upload = [file for file in all_report_files if
                                  file['datatype'] in {'report_pdf', 'report_xml', 'report_json'}]
        result = []
        for file in report_files_to_upload:
            _, blob = get_bucket_and_blob_from_gs_path(storage_client=self.storage_client, gs_path=file['path'])
            result.append(blob)
        return result

    def _copy_blob_to_portal_bucket(self, blob, target_sub_folder):
        if target_sub_folder != '' and target_sub_folder[-1] != '/':
            target_sub_folder += '/'
        file_name = get_file_name_from_blob_name(blob.name)
        source_bucket = blob.bucket
        print(f"{blob.name} ---> {self.portal_bucket.name}")
        source_bucket.copy_blob(blob=blob,
                                destination_bucket=self.portal_bucket,
                                new_name=f'{self.sample_barcode}/{target_sub_folder}{file_name}')

    def _copy_blob_to_archive_bucket(self, blob):
        file_name = get_file_name_from_blob_name(blob.name)
        source_bucket = blob.bucket
        print(f"{blob.name} ---> {self.archive_bucket.name}")
        source_bucket.copy_blob(blob=blob,
                                destination_bucket=self.archive_bucket,
                                new_name=file_name)


if __name__ == "__main__":
    main()
