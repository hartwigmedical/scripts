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
        if not self.run:
            self._prompt_user_no_run()
        elif self.run['status'] != 'Validated':
            self._prompt_user_non_validated_run()

        self._delete_old_artifacts_in_portal_bucket()
        self._delete_old_artifacts_in_gcp_share_bucket()

        self._copy_files_to_remote_buckets(publish_to_portal=publish_to_portal)

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

    def _delete_old_artifacts_in_gcp_share_bucket(self):
        print(f"deleting old artifacts from bucket '{self.panel_share_bucket.name}'")
        blobs_old_run = list(self.panel_share_bucket.list_blobs(prefix=self.sample_barcode))
        print(*[blob.name for blob in blobs_old_run], sep='\n')
        self.panel_share_bucket.delete_blobs(blobs=blobs_old_run)

    def _copy_files_to_remote_buckets(self, publish_to_portal):
        report_blobs = self._get_report_blobs()
        self._archive_blobs(report_blobs)
        if publish_to_portal:
            if self._is_failure():
                self._share_failure_report(report_blobs)
            elif self._is_panel_failure():
                self._share_panel_failure_report(report_blobs)
            elif self._is_panel():
                self._share_panel_report(report_blobs)
            else:
                self._share_wgs_report(report_blobs)

    def _archive_blobs(self, blobs):
        print(f"Copying ${len(blobs)} blobs to the archive bucket")
        for blob in blobs:
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.archive_bucket)

    def _share_failure_report(self, report_blobs):
        print(f"Sharing ${len(report_blobs)} report files with the portal")
        for blob in report_blobs:
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.portal_bucket)

    def _share_panel_failure_report(self, report_blobs):
        print(f"Sharing ${len(report_blobs)} report files with the portal")
        for blob in report_blobs:
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.portal_bucket)
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.panel_share_bucket)

    def _share_panel_report(self, report_blobs):
        panel_blobs = self._get_blobs_from_bucket(bucket=self.panel_pipeline_output_bucket,
                                                  file_names=self._panel_files())

        print(f"Sharing ${len(report_blobs)} report files with the portal")
        for blob in report_blobs:
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.portal_bucket)
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.panel_share_bucket)
        print(f"Sharing a total of '{len(panel_blobs)}' panel files with the portal")
        for blob in panel_blobs:
            self._copy_blob_to_bucket(blob, self.portal_bucket, target_sub_folder='RUO')
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.panel_share_bucket, target_sub_folder='RUO')

    def _share_wgs_report(self, report_blobs):
        molecular_blobs = self._get_blobs_from_bucket(bucket=self.pipeline_output_bucket,
                                                      file_names=self._molecular_files())
        germline_blobs = self._get_blobs_from_bucket(bucket=self.pipeline_output_bucket,
                                                     file_names=self._germline_files())
        print(f"Sharing ${len(report_blobs)} report files with the portal")
        for blob in report_blobs:
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.portal_bucket)
        print(f"Sharing a total of '{len(molecular_blobs)}' molecular files with the portal")
        for blob in molecular_blobs:
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.portal_bucket, target_sub_folder='RUO')
        print(f"Sharing a total of '{len(germline_blobs)}' germline files with the portal")
        for blob in germline_blobs:
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.portal_bucket,
                                      target_sub_folder='RUO_germline')

    def _get_blobs_from_bucket(self, bucket, file_names):
        result = []
        bucket_contents = list(bucket.list_blobs(prefix=self._set_name()))
        for blob in bucket_contents:
            for file_name in file_names:
                if blob.name[-len(file_name):] == file_name:  # this checks if the blob name ends with the file name.
                    result.append(blob)
                    break
        return result

    def _panel_files(self):
        return {
            "driver.catalog.somatic.tsv",
            "purple.somatic.vcf.gz",
            "purple.purity.tsv",
            "purple.cnv.gene.tsv",
            "purple.sv.vcf.gz",
            "reported.somatic.vcf",
            "sampleQcReport.pdf",
            "sampleQcReport.png",
            "sampleQcDeamination.pdf",
            "sampleQcDeamination.png",
            "cobalt.ratio.tsv.gz",
            "cobalt.ratio.pcf",
            "amber.baf.tsv.gz",
            "purple.cnv.somatic.tsv",
        }

    def _molecular_files(self):
        set_name = self._set_name()
        reporting_id = self.rest_client.get_lama_patient_reporter_data(self.report_created_record['barcode'])[
            'reportingId']

        return {
            "purple.driver.catalog.somatic.tsv",
            "linx.driver.catalog.tsv",
            "purple.somatic.vcf.gz",
            "purple.sv.vcf.gz",
            "linx.fusion.tsv",
            f"{set_name}/orange_no_germline/{reporting_id}.orange.pdf"
        }

    def _germline_files(self):
        set_name = self._set_name()
        reporting_id = self.rest_client.get_lama_patient_reporter_data(self.report_created_record['barcode'])[
            'reportingId']

        return {
            "purple.germline.vcf.gz",
            f"{set_name}/purple/{reporting_id}.purple.germline.vcf.gz.tbi",
            "purple.germline.deletion.tsv",
            "purple.driver.catalog.germline.tsv",
            "gripss.filtered.germline.vcf.gz",
            f"{set_name}/gripss_germline/{reporting_id}.gripss.filtered.germline.vcf.gz.tbi",
            "linx.germline.disruption.tsv",
            "linx.germline.driver.catalog.tsv"
        }

    def _get_report_blobs(self):
        all_report_files = self.report_created_record["report_files"]
        report_files_to_upload = [file for file in all_report_files if
                                  file['datatype'] in {'report_pdf', 'report_xml', 'report_json'}]
        result = []
        for file in report_files_to_upload:
            _, blob = get_bucket_and_blob_from_gs_path(storage_client=self.storage_client, gs_path=file['path'])
            result.append(blob)
        return result

    def _copy_blob_to_bucket(self, blob, destination_bucket, target_sub_folder=''):
        if target_sub_folder != '' and target_sub_folder[-1] != '/':
            target_sub_folder += '/'
        file_name = get_file_name_from_blob_name(blob.name)
        source_bucket = blob.bucket
        print(f"{blob.name} ---> {destination_bucket.name}")
        source_bucket.copy_blob(blob=blob,
                                destination_bucket=destination_bucket,
                                new_name=f'{self.sample_barcode}/{target_sub_folder}{file_name}')

    def _set_name(self):
        return self.run['set']['name']

    def _is_panel(self):
        return self.report_created_record['report_type'] in ['panel_result_report', 'panel_result_report_fail']

    def _is_failure(self):
        return self.report_created_record['report_type'] in ['wgs_processing_issue', 'wgs_isolation_fail',
                                                             'wgs_tcp_shallow_fail', 'wgs_preparation_fail',
                                                             'wgs_tumor_processing_issue', 'wgs_pipeline_fail',
                                                             'wgs_tcp_fail']

    def _is_panel_failure(self):
        return self.report_created_record['report_type'] in ['panel_result_report_fail']


if __name__ == "__main__":
    main()
