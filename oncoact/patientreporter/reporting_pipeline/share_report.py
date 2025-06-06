import argparse
from rest_util import RestClient
from google.cloud.storage import Bucket, Client
from gsutil import get_bucket_and_blob_from_gs_path, get_file_name_from_blob_name


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

        self.pipeline_output_bucket: Bucket = self.storage_client.bucket(f"diagnostic-pipeline-output-{profile}-1")

        self.archive_bucket: Bucket = self.storage_client.bucket(f'patient-reporter-final-{profile}-1')
        self.portal_bucket: Bucket = self.storage_client.bucket(f'hmf-customer-portal-report-shared-{profile}')
        self.panel_pipeline_output_bucket: Bucket = self.storage_client.bucket(f'targeted-pipeline-output-{profile}-1')
        self.panel_share_bucket: Bucket = self.storage_client.bucket(f'oncoact-panel-files-nki')

        self.report_created_record = self.rest_client.get_report_created(self.sample_barcode)
        self.run_id = self.report_created_record['run_id']
        self.run = self.rest_client.get_run(self.run_id) if self.run_id else None

    def share_report(self):
        """
        Shares the report by updating the remote portal bucket and the remote archive bucket. It also updates the API.

        """
        if not self.run:
            self._prompt_user_no_run()
        elif self.run['status'] != 'Validated':
            self._prompt_user_non_validated_run()

        self._delete_old_artifacts_in_portal_bucket()
        self._delete_old_artifacts_in_gcp_share_bucket()

        lama = self.rest_client.get_lama_patient_reporter_data(self.report_created_record['barcode'])

        if "reportSettings" not in lama:
            print("key reportSettings is not present")
            exit(1)
        else:
            report_settings = lama["reportSettings"]

        if "isSharedThroughPortal" not in report_settings:
            print("isSharedThroughPortal is not present")
            exit(1)
        else:
            is_shared_through_portal = report_settings['isSharedThroughPortal']

        self._copy_files_to_remote_buckets(publish_to_portal=is_shared_through_portal, report_settings=report_settings, lama=lama)

        print('Updating api')
        response = self.rest_client.post_report_shared(report_created_id=self.report_created_record['id'],
                                                       publish_to_portal=is_shared_through_portal,
                                                       notify_users=is_shared_through_portal)
        print("API response:", response)

        delete_response = self._delete_run_from_reporting_pipeline()
        print("Reporting pipeline response:", delete_response)

        if is_shared_through_portal:
            print("Report is shared through portal")
        else:
            print("Reports should be shared through nextcloud")

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

    def _delete_run_from_reporting_pipeline(self):
        print(f"deleting runs with sample barcode '{self.sample_barcode}' from reporting pipeline")
        return self.rest_client.delete_finished_execution(self.sample_barcode)

    def _copy_files_to_remote_buckets(self, publish_to_portal, report_settings, lama):
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
                self._share_wgs_report(report_blobs, report_settings, lama)

    def _archive_blobs(self, blobs):
        print(f"Copying {len(blobs)} blobs to the archive bucket")
        for blob in blobs:
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.archive_bucket)

    def _share_failure_report(self, report_blobs):
        print(f"Sharing {len(report_blobs)} report files with the portal")
        for blob in report_blobs:
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.portal_bucket)

    def _share_panel_failure_report(self, report_blobs):
        print(f"Sharing {len(report_blobs)} report files with the portal")
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

    def _share_wgs_report(self, report_blobs, report_settings, lama):
        if "shareGermlineData" not in report_settings:
            print("shareGermlineData is not present")
            exit(1)
        else:
            share_germline_data = report_settings['shareGermlineData']

        if "reportingId" not in lama:
            print("reportingId is not present")
            exit(1)
        else:
            reporting_id = lama['reportingId']

        if "hospitalSampleLabel" not in lama:
            print("hospitalSampleLabel is not present")
            exit(1)
        else:
            hospital_sample_label = lama['hospitalSampleLabel']

        if hospital_sample_label is not None:
            converted_reporting_id = f"{reporting_id}-{hospital_sample_label}"
        else:
            converted_reporting_id = f"{reporting_id}"


        molecular_blobs = self._get_blobs_from_bucket(bucket=self.pipeline_output_bucket,
                                                      file_names=self._molecular_files(converted_reporting_id))
        germline_blobs = self._get_blobs_from_bucket(bucket=self.pipeline_output_bucket,
                                                     file_names=self._germline_files(converted_reporting_id))

        print(f"Sharing {len(report_blobs)} report files with the portal")
        for blob in report_blobs:
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.portal_bucket)
        print(f"Sharing a total of '{len(molecular_blobs)}' molecular files with the portal")
        for blob in molecular_blobs:
            self._copy_blob_to_bucket(blob=blob, destination_bucket=self.portal_bucket, target_sub_folder='RUO')

        if share_germline_data:
            print(f"Sharing a total of '{len(germline_blobs)}' germline files with the portal")
            for blob in germline_blobs:
                self._copy_blob_to_bucket(blob=blob, destination_bucket=self.portal_bucket,
                                          target_sub_folder='RUO_germline')
            import subprocess

        # Create + share merged vcf files with id snps
        try:
            set_name = self._set_name()
            print(f"Running combine_vcfs.sh with set_name: {set_name}")
            subprocess.run(["bash", "id_snpcaller_vcfmerge", set_name], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Snp check + vcf combine shell script failed: {e}")
            return

        # Prepare expected VCF file name
        combined_vcf_filename = f"{self.sample_barcode}_merged_final.vcf"
        combined_vcf_files = {combined_vcf_filename}

        # Get blob from the 'wgs-combined-snps-vcfs' bucket
        vcf_bucket = self.storage_client.bucket("wgs-combined-snps-vcfs")
        combined_vcf_blobs = self._get_blobs_from_bucket(bucket=vcf_bucket, file_names=combined_vcf_files)

        # If found, copy it to the portal
        if combined_vcf_blobs:
            print(f"Sharing combined VCF file: {combined_vcf_blobs[0].name}")
            self._copy_blob_to_bucket(blob=combined_vcf_blobs[0], destination_bucket=self.portal_bucket, target_sub_folder="RUO")
        else:
            print(f"No combined VCF file found for {self.sample_barcode}")


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
            "orange.pdf",
            ".reconCNV.html",
            ".sage.visualisation.zip"
        }

    def _molecular_files(self, converted_reporting_id):
        set_name = self._set_name()

        return {
            "purple.driver.catalog.somatic.tsv",
            "linx.driver.catalog.tsv",
            "purple.somatic.vcf.gz",
            "purple.sv.vcf.gz",
            "linx.fusion.tsv",
            f"{set_name}/orange_no_germline/{converted_reporting_id}.orange.pdf",
            f"{set_name}/purple/{converted_reporting_id}.purple.cnv.gene.tsv"
        }

    def _germline_files(self, converted_reporting_id):
        set_name = self._set_name()

        return {
            "purple.germline.vcf.gz",
            f"{set_name}/purple/{converted_reporting_id}.purple.germline.vcf.gz.tbi",
            "purple.germline.deletion.tsv",
            "purple.driver.catalog.germline.tsv",
            "gripss.filtered.germline.vcf.gz",
            f"{set_name}/gripss_germline/{converted_reporting_id}.gripss.filtered.germline.vcf.gz.tbi",
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
        print(f"{source_bucket}/{blob.name} ---> {destination_bucket.name}/{target_sub_folder}{file_name}")
        source_bucket.copy_blob(blob=blob,
                                destination_bucket=destination_bucket,
                                new_name=f'{self.sample_barcode}/{target_sub_folder}{file_name}')

    def _set_name(self):
        return self.run['set']['name']

    def _is_panel(self):
        return self.report_created_record['report_type'] in ['panel_result_report',
                                                             'panel_result_report_corrected_internal',
                                                             'panel_result_report_corrected_external',
                                                             'panel_result_report_fail',
                                                             'panel_result_report_fail_corrected_internal',
                                                             'panel_result_report_fail_corrected_external']

    def _is_failure(self):
        return self.report_created_record['report_type'] in ['wgs_processing_issue',
                                                             'wgs_processing_issue_corrected_internal',
                                                             'wgs_processing_issue_corrected_external',
                                                             'wgs_isolation_fail',
                                                             'wgs_isolation_fail_corrected_internal',
                                                             'wgs_isolation_fail_corrected_external',
                                                             'wgs_tcp_shallow_fail',
                                                             'wgs_tcp_shallow_fail_corrected_internal',
                                                             'wgs_tcp_shallow_fail_corrected_external',
                                                             'wgs_preparation_fail',
                                                             'wgs_preparation_fail_corrected_internal',
                                                             'wgs_preparation_fail_corrected_external',
                                                             'wgs_tumor_processing_issue',
                                                             'wgs_tumor_processing_issue_corrected_internal',
                                                             'wgs_tumor_processing_issue_corrected_external',
                                                             'wgs_pipeline_fail',
                                                             'wgs_pipeline_fail_corrected_internal',
                                                             'wgs_pipeline_fail_corrected_external',
                                                             'wgs_tcp_fail',
                                                             'wgs_tcp_fail_corrected_internal',
                                                             'wgs_tcp_fail_corrected_external']

    def _is_panel_failure(self):
        return self.report_created_record['report_type'] in ['panel_result_report_fail',
                                                             'panel_result_report_fail_corrected_internal',
                                                             'panel_result_report_fail_corrected_external']


if __name__ == "__main__":
    main()
