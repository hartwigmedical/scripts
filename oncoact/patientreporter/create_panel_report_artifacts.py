import argparse
import pathlib
import rest_util
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

        self._sample_name = self._rest_client.get_report_created(sample_barcode=self._sample_barcode)['sample_name']

    def generate_artifacts(self):
        output_folder = self.generate_output_folder()
        print(f"Running R scripts.")
        self.run_scripts(output_folder=output_folder)
        print(f"The scripts have finished. Output is stored at '{output_folder}'")

    def generate_output_folder(self):
        home_dir = os.path.expanduser('~')
        path = f'{home_dir}/tmp/{self._sample_barcode}/panel_artifacts'
        pathlib.Path(path).mkdir(parents=True, exist_ok=True)
        return path

    def run_scripts(self, output_folder):
        script_location = '/data/repos/scripts/panel/createSampleQcReport.R'
        script2_location = '/data/repos/scripts/panel/createSampleQcReport_Deamination.R'

        self.run_r_script(script_location, output_folder=output_folder, log_file_name='createSampleQcReport.log')
        self.run_r_script(script2_location, output_folder=output_folder,
                          log_file_name='createSampleQcReport_Deamination.log')

    def run_r_script(self, script_location, output_folder, log_file_name):
        output = subprocess.check_output(['Rscript', script_location, self._sample_name, output_folder, output_folder])
        generate_file(f'{output_folder}/{log_file_name}', content=output.decode())


def generate_file(path, content):
    with open(path, 'x') as file:
        file.write(content)


if __name__ == '__main__':
    main()
