import argparse
import pathlib
import rest_util
from google.cloud.storage import Client, Bucket, Blob
import subprocess
import os


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

    def generate_artifacts(self):
        output_folder = self._generate_output_folder()
        self._download_required_resources(output_folder)

        self._run_scripts(output_folder=output_folder)
        print(f"The scripts have finished. Output is stored at '{output_folder}'")

    def _generate_output_folder(self):
        home_dir = os.path.expanduser('~')
        path = f'{home_dir}/tmp/{self._sample_barcode}/panel_artifacts'
        pathlib.Path(path).mkdir(parents=True, exist_ok=True)
        return path

    def _run_scripts(self, output_folder):
        sample_report_script_path = '/data/repos/scripts/panel/createSampleQcReport.R'
        deamination_script_path = '/data/repos/scripts/panel/createSampleQcReport_Deamination.R'

        self._run_r_script(sample_report_script_path, output_folder=output_folder,
                           log_file_name='createSampleQcReport.log')
        self._run_r_script(deamination_script_path, output_folder=output_folder,
                           log_file_name='createSampleQcReport_Deamination.log')

    def _run_r_script(self, script_location, output_folder, log_file_name):
        sample_name = self._report_created_record['sample_name']
        output = subprocess.check_output(['Rscript', script_location, sample_name, output_folder, output_folder])
        _generate_file(f'{output_folder}/{log_file_name}', content=output.decode())

    def _download_required_resources(self, download_folder):
        required_resources: list[Blob] = self._get_required_resources_as_blobs()

        for blob in required_resources:
            target_file_name = blob.name.split('/', 4)[4]
            with open(f'{download_folder}/{target_file_name}', 'wb') as file:
                blob.download_to_file(file)

    def _get_required_resources_as_blobs(self):
        required_file_suffixes = {
            "driver.catalog.somatic.tsv",
            "purple.somatic.vcf.gz",
            "purple.purity.tsv",
            "purple.somatic.vcf.gz",
            "purple.cnv.gene.tsv"
        }
        run_id = self._report_created_record['run_id']
        set_name = self._rest_client.get_run(run_id)['set']['name']
        run_blobs = self._panel_pipeline_output_bucket.list_blobs(prefix=set_name)

        result = []
        for suffix in required_file_suffixes:
            for run_blob in run_blobs:
                if run_blob.name[-len(suffix):] == suffix:  # this checks if the blob name ends with the suffix.
                    result.append(run_blob)

        return result


def _generate_file(path, content):
    with open(path, 'x') as file:
        file.write(content)


if __name__ == '__main__':
    main()
