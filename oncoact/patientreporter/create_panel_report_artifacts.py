import argparse
import pathlib
import rest_util
from google.cloud.storage import Client, Bucket
from gsutil import get_bucket_and_blob_names_from_gs_path
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

    def _download_required_resources(self, destination_folder):
        run_id = self._report_created_record['run_id']
        run_files = self._rest_client.get_run_files(run_id)
        required_files_by_type = {
            'sage.exon.medians.tsv',  # TODO
            'purple.cnv.gene.tsv',  # TODO
            'purple_purity',
            'somatic_variants_purple'
        }
        required_file_paths = [file['filepath'] for file in run_files if file['datatype'] in required_files_by_type]
        for file_path in required_file_paths:
            destination_path = f'{destination_folder}/{file_path.split("/", 4)[4]}'
            bucket_name, blob_name = get_bucket_and_blob_names_from_gs_path(file_path)
            with open(destination_path, 'wb') as file_obj:
                self._storage_client.bucket(bucket_name).get_blob(blob_name).download_to_file(file_obj)


def _generate_file(path, content):
    with open(path, 'x') as file:
        file.write(content)


if __name__ == '__main__':
    main()
