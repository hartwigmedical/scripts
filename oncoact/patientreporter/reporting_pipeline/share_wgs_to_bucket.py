import argparse
import sys
from datetime import date
from rest_util import RestClient
from google.cloud.storage import Bucket, Client
from gsutil import get_bucket_and_blob_from_gs_path, get_file_name_from_blob_name
import subprocess

def main():
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument('sample_barcode')
    argument_parser.add_argument('--profile', choices=['pilot', 'preview', 'prod'], default='pilot')
    args = argument_parser.parse_args()

    ReportSharer(sample_barcode=args.sample_barcode, profile=args.profile).share_report()

class ReportSharer:

    def __init__(self, sample_barcode, profile='pilot'):
        self.sample_barcode = sample_barcode
        self.rest_client: RestClient = RestClient(profile)
        self.storage_client: Client = Client()

        self.pipeline_output_bucket: Bucket = self.storage_client.bucket(f'diagnostic-pipeline-output-{profile}-1')
        self.portal_bucket: Bucket = self.storage_client.bucket(f'hmf-customer-portal-report-shared-{profile}')
        self.wgs_share_bucket: Bucket = self.storage_client.bucket(f'oncoact-wgs-files-nki')


    def share_report(self):
        """
        Shares the report by updating the remote portal bucket and the remote archive bucket. It also updates the API.
        """
        print("Deleting existing wgs files of GCP shared bucket")
        self._delete_old_artifacts_in_gcp_share_bucket()

        print("Copy wgs files to GCP shared bucket")
        self._copy_files_to_remote_bucket()

        print("Creating completion marker file")
        self._create_complete_file()

        print("Done!")

    def _delete_old_artifacts_in_gcp_share_bucket(self):
        print(f"deleting old artifacts from bucket '{self.wgs_share_bucket.name}'")
        blobs_old_run = list(self.wgs_share_bucket.list_blobs(prefix=self.sample_barcode))
        print(*[blob.name for blob in blobs_old_run], sep='\n')
        self.wgs_share_bucket.delete_blobs(blobs=blobs_old_run)

    def _copy_files_to_remote_bucket(self):
        report_blobs = self._get_report_blobs()
        if not report_blobs:
            print("This is not a NKI sample so skipping copy files")
            sys.exit(1)
        else:
            self._share_report(report_blobs)

    def _get_report_blobs(self):
        sample_barcode = self.sample_barcode
        result = []
        cmd_pdf = f"gsutil ls gs://{self.portal_bucket.name}/{sample_barcode}/*_NKI-AVL_*_oncoact_wgs_report.pdf"
        pathPdf = subprocess.check_output(cmd_pdf, shell=True, text=True).strip()
        _, blobPdf = get_bucket_and_blob_from_gs_path(storage_client=self.storage_client, gs_path=pathPdf)

        cmd_json = f"gsutil ls gs://{self.portal_bucket.name}/{sample_barcode}/*_NKI-AVL_*_oncoact_wgs_report.json"
        pathJson = subprocess.check_output(cmd_json, shell=True, text=True).strip()
        _, blobJson = get_bucket_and_blob_from_gs_path(storage_client=self.storage_client, gs_path=pathJson)

        cmd_xml = f"gsutil ls gs://{self.portal_bucket.name}/{sample_barcode}/*_NKI-AVL_*_oncoact_wgs_report.xml"
        pathXml = subprocess.check_output(cmd_xml, shell=True, text=True).strip()
        _, blobXml = get_bucket_and_blob_from_gs_path(storage_client=self.storage_client, gs_path=pathXml)

        result.append(blobPdf)
        result.append(blobJson)
        result.append(blobXml)

        return result

    def _share_failure_report(self, report_blobs):
        print(f"Sharing {len(report_blobs)} report files with the GCP bucket")
        for blob in report_blobs:
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.wgs_share_bucket)

    def _share_report(self, report_blobs):
        wgs_blobs = self._get_blobs_from_bucket(bucket=self.portal_bucket, subfolder='RUO')
        wgs_germline_blobs = self._get_blobs_from_bucket(bucket=self.portal_bucket, subfolder='RUO_germline')

        print(f"Sharing ${len(report_blobs)} report files with the GCP bucket")
        for blob in report_blobs:
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.wgs_share_bucket)

        print(f"Sharing a total of '{len(wgs_blobs)}' wgs files with the GCP bucket")
        for blob in wgs_blobs:
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.wgs_share_bucket, target_sub_folder='RUO')

        if wgs_germline_blobs:
            print(f"Sharing a total of '{len(wgs_germline_blobs)}' wgs germline files with the GCP bucket")
            for blob in wgs_germline_blobs:
                self._copy_blob_to_bucket(blob=blob, destination_bucket=self.wgs_share_bucket, target_sub_folder='RUO_germline')

    def _get_blobs_from_bucket(self, bucket, subfolder):
        print(bucket)
        blobs = list(bucket.list_blobs(prefix=self.sample_barcode + f"/{subfolder}/"))
        print(blobs)
        return blobs

    def _copy_blob_to_bucket(self, blob, destination_bucket, target_sub_folder=''):
        if target_sub_folder != '' and target_sub_folder[-1] != '/':
            target_sub_folder += '/'
        file_name = get_file_name_from_blob_name(blob.name)
        source_bucket = blob.bucket
        print(f"{source_bucket}/{blob.name} ---> {destination_bucket.name}/{target_sub_folder}{file_name}")
        source_bucket.copy_blob(blob=blob,
                                destination_bucket=destination_bucket,
                                new_name=f'{self.sample_barcode}/{target_sub_folder}{file_name}')

    def _create_complete_file(self):
        blob = self.wgs_share_bucket.blob(f'{self.sample_barcode}.complete')
        blob.upload_from_string(date.today().strftime('%Y-%m-%d'))

if __name__ == "__main__":
    main()