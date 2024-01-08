import argparse
import pathlib
import shutil

import rest_util
from google.cloud.storage import Client, Bucket, Blob
import subprocess
import os
import gzip
from reports_to_nc import upload_to_nextcloud


def main():
    argument_parser = argparse.ArgumentParser('Generates the panel artifacts for a given tumor_sample_barcode')
    argument_parser.add_argument('tumor_sample_barcode')
    argument_parser.add_argument('--profile', choices=['pilot', 'prod'], default='pilot')

    args = argument_parser.parse_args()

    ArtifactGenerator(profile=args.profile, sample_barcode=args.tumor_sample_barcode).generate_artifacts()


class ArtifactGenerator:

    def __init__(self, profile, sample_barcode):
        self._sample_barcode = sample_barcode
        self._rest_client = rest_util.RestClient(profile=profile)
        self._storage_client = Client()

        self._report_created_record = self._rest_client.get_report_created(sample_barcode=self._sample_barcode)
        self._panel_pipeline_output_bucket: Bucket = self._storage_client.bucket(
            f'targeted-pipeline-output-{profile}-1')

        run_id = self._report_created_record['run_id']
        self.run = self._rest_client.get_run(run_id)
        self.hmf_id = self.run['set']['tumor_sample']
        self.set_name = self.run['set']['name']

    def generate_artifacts(self):
        input_folder, output_folder = self._generate_input_and_output_folders()
        self._download_required_resources(download_to=input_folder)

        self._run_scripts(input_folder=input_folder, output_folder=output_folder)
        # self._generate_vcf(input_folder=input_folder, output_folder=output_folder)
        # upload_to_nextcloud(output_folder)
        # self._copy_output_to_bucket(output_folder=output_folder)


    def _generate_input_and_output_folders(self):
        home_dir = os.path.expanduser('~')
        path = f'{home_dir}/tmp/{self._sample_barcode}/panel_artifacts'

        input_folder = f'{path}/input/'
        output_folder = f'{path}/output/'
        shutil.rmtree(input_folder, ignore_errors=True)
        shutil.rmtree(output_folder, ignore_errors=True)
        pathlib.Path(input_folder).mkdir(parents=True, exist_ok=False)
        pathlib.Path(output_folder).mkdir(parents=True, exist_ok=False)
        return input_folder, output_folder

    def _run_scripts(self, input_folder, output_folder):
        sample_report_script_path = '/data/repos/scripts/panel/createSampleQcReport.R'
        deamination_script_path = '/data/repos/scripts/panel/createSampleQcReport_Deamination.R'

        self._run_r_script_qc(sample_report_script_path, input_folder=input_folder, output_folder=output_folder)
        self._run_r_script_deamination(deamination_script_path, input_folder=input_folder, output_folder=output_folder)
        print(f"The scripts have finished. Output is stored at '{output_folder}'")

    def _run_r_script_qc(self, script_location, input_folder, output_folder):
        sample_name = self._report_created_record['sample_name']
        barcode = self._rest_client.get_tumor_sample_barcode_from_run_id(self._report_created_record['run_id'])
        gb_yield = self._rest_client.get_yield(barcode)
        subprocess.run(['Rscript', script_location, sample_name, input_folder, output_folder, gb_yield], check=False)

    def _run_r_script_deamination(self, script_location, input_folder, output_folder):
        sample_name = self._report_created_record['sample_name']
        subprocess.run(['Rscript', script_location, sample_name, input_folder, output_folder], check=False)

    def _generate_vcf(self, input_folder, output_folder):
        res = []
        with gzip.open(f"{input_folder}purple/{self.hmf_id}.purple.somatic.vcf.gz", 'rt') as file:
            for line in file.readlines():
                if "REPORTED" in line or (len(line) > 0 and line[0] == '#'):
                    res.append(line)

        if res:
            with open(f"{output_folder}{self.hmf_id}.reported.somatic.vcf", 'x') as file:
                file.writelines(res)

    def _download_required_resources(self, download_to):
        required_resources: list[Blob] = self._get_required_resources_as_blobs()

        for blob in required_resources:
            target_file_path = blob.name.split('/', 1)[1]
            dir_name = os.path.dirname(f'{download_to}/{target_file_path}')
            os.makedirs(dir_name, exist_ok=True)

            with open(f'{download_to}/{target_file_path}', 'wb') as file:
                blob.download_to_file(file)

    def _get_required_resources_as_blobs(self):
        required_file_suffixes = {
            "driver.catalog.somatic.tsv",
            "purple.somatic.vcf.gz",
            "purple.purity.tsv",
            "purple.somatic.vcf.gz",
            "purple.cnv.gene.tsv",
            "sage.exon.medians.tsv"
        }
        run_blobs = list(self._panel_pipeline_output_bucket.list_blobs(prefix=self.set_name))

        result = []
        for suffix in required_file_suffixes:
            for run_blob in run_blobs:
                if run_blob.name[-len(suffix):] == suffix:  # this checks if the blob name ends with the suffix.
                    result.append(run_blob)
                    break

        return result

    def _copy_output_to_bucket(self, output_folder):
        output_blob_prefix = f"{self.set_name}/reporting"

        artifacts = dict()
        for dir_path, dir_names, filenames in os.walk(output_folder):
            for filename in filenames:
                artifact_path = os.path.join(dir_path, filename)
                artifacts[artifact_path] = filename

        for path, filename in artifacts.items():
            output_blob_name = f'{output_blob_prefix}/{filename}'
            output_blob = self._panel_pipeline_output_bucket.blob(blob_name=output_blob_name)
            output_blob.upload_from_filename(path)
            print(f"Copied {filename} to gs://{self._panel_pipeline_output_bucket.name}/{output_blob.name}")


def _generate_file(path, content):
    with open(path, 'x') as file:
        file.write(content)



if __name__ == '__main__':
    main()
