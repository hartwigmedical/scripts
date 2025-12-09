import argparse
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

        self.panel_pipeline_output_bucket: Bucket = self.storage_client.bucket(f'targeted-pipeline-output-{profile}-1')
        self.oncoact_bucket: Bucket = self.storage_client.bucket(bucket_name="patient-reporter-final-prod-1")
        self.panel_share_bucket: Bucket = self.storage_client.bucket(f'oncoact-panel-files-nki')

        self.report_created_record = self.rest_client.get_report_created(self.sample_barcode)
        self.run_id = self.report_created_record['run_id']
        self.run = self.rest_client.get_run(self.run_id) if self.run_id else None

    def share_report(self):
        """
        Shares the report by updating the remote portal bucket and the remote archive bucket. It also updates the API.

        """
        print("Deleting existing panel files of GCP shared bucket")
        self._delete_old_artifacts_in_gcp_share_bucket()

        print("Copy panel files to GCP shared bucket")
        self._copy_files_to_remote_buckets()

        print("Done!")

    def _delete_old_artifacts_in_gcp_share_bucket(self):
        print(f"deleting old artifacts from bucket '{self.panel_share_bucket.name}'")
        blobs_old_run = list(self.panel_share_bucket.list_blobs(prefix=self.sample_barcode))
        print(*[blob.name for blob in blobs_old_run], sep='\n')
        self.panel_share_bucket.delete_blobs(blobs=blobs_old_run)

    def _copy_files_to_remote_buckets(self):
        report_blobs = self._get_report_blobs()
        if self._is_panel_failure():
            self._share_panel_failure_report(report_blobs)
        else:
            self._share_panel_report(report_blobs)

    def _get_report_blobs(self):
        sample_barcode = self.sample_barcode
        result = []
        cmd_pdf = f"gsutil ls gs://{self.oncoact_bucket.name}/panel-{sample_barcode}/patient-reporter/*_oncoact_panel_result_report.pdf"
        pathPdf = subprocess.check_output(cmd_pdf, shell=True, text=True).strip()
        _, blobPdf = get_bucket_and_blob_from_gs_path(storage_client=self.storage_client, gs_path=pathPdf)

        cmd_json = f"gsutil ls gs://{self.oncoact_bucket.name}/panel-{sample_barcode}/patient-reporter/*_oncoact_panel_result_report.json"
        pathJson = subprocess.check_output(cmd_json, shell=True, text=True).strip()
        _, blobJson = get_bucket_and_blob_from_gs_path(storage_client=self.storage_client, gs_path=pathJson)

        result.append(blobPdf)
        result.append(blobJson)
        return result

    def _share_panel_failure_report(self, report_blobs):
        print(f"Sharing {len(report_blobs)} report files with the GCP bucket")
        for blob in report_blobs:
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.panel_share_bucket)

    def _share_panel_report(self, report_blobs):
        panel_blobs = self._get_blobs_from_bucket(bucket=self.panel_pipeline_output_bucket,
                                                  file_names=self._panel_files())

        print(f"Sharing ${len(report_blobs)} report files with the GCP bucket")
        for blob in report_blobs:
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.panel_share_bucket)
        print(f"Sharing a total of '{len(panel_blobs)}' panel files with the GCP bucket")
        for blob in panel_blobs:
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.panel_share_bucket, target_sub_folder='RUO')

        igv_config = f"gsutil ls gs://{self.oncoact_bucket.name}/panel-{self.sample_barcode}/create-igv-config/*-igv-config.txt"
        pathIgv_config = subprocess.check_output(igv_config, shell=True, text=True).strip()
        _, blob_Igv_config = get_bucket_and_blob_from_gs_path(storage_client=self.storage_client, gs_path=pathIgv_config)

        self._copy_blob_to_bucket(blob=blob_Igv_config, destination_bucket=self.panel_share_bucket, target_sub_folder='RUO')

    def _get_blobs_from_bucket(self, bucket, file_names):
        result = []
        bucket_contents = list(bucket.list_blobs(prefix=self._set_name()))
        for blob in bucket_contents:
            for file_name in file_names:
                if blob.name[-len(file_name):] == file_name:  # this checks if the blob name ends with the file name.
                    result.append(blob)
                    break
        return result

    def _copy_blob_to_bucket(self, blob, destination_bucket, target_sub_folder=''):
        if target_sub_folder != '' and target_sub_folder[-1] != '/':
            target_sub_folder += '/'
        file_name = get_file_name_from_blob_name(blob.name)
        source_bucket = blob.bucket
        print(f"{source_bucket}/{blob.name} ---> {destination_bucket.name}/{target_sub_folder}{file_name}")
        source_bucket.copy_blob(blob=blob,
                                destination_bucket=destination_bucket,
                                new_name=f'{self.sample_barcode}/{target_sub_folder}{file_name}')

    def _set_name(self):
        return self.run['set']['name']

    def _is_panel_failure(self):
        return self.report_created_record['report_type'] in ['panel_result_report_fail',
                                                             'panel_result_report_fail_corrected_internal',
                                                             'panel_result_report_fail_corrected_external']

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
            "orange.pdf",
            ".reconCNV.html",
            ".sage.visualisation.zip",
            ".vchord.prediction.tsv",
            ".gnomad.pon.vcf",
            "-igv-config.txt"
        }

if __name__ == "__main__":
    main()