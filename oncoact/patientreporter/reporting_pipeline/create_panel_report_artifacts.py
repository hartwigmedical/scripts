import argparse
import json
import pathlib
import shutil

import rest_util
from google.cloud.storage import Client, Bucket, Blob
import subprocess
import os
import gzip
from reports_to_nc_for_sharing import upload_to_nextcloud


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

        self.run_id = self._report_created_record['run_id']
        self.run = self._rest_client.get_run(self.run_id)
        self.set_name = self.run['set']['name']

    def generate_artifacts(self):
        input_folder, output_folder = self._generate_input_and_output_folders()
        self._download_required_resources(download_to=input_folder)

      #  self._run_scripts(input_folder=input_folder, output_folder=output_folder)
        self._generate_vcf(input_folder=input_folder, output_folder=output_folder)
      #  upload_to_nextcloud(output_folder)
        self._copy_output_to_bucket(output_folder=output_folder)

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

        self._run_r_script(sample_report_script_path, input_folder=input_folder, output_folder=output_folder)
        self._run_r_script(deamination_script_path, input_folder=input_folder, output_folder=output_folder)
        print(f"The scripts have finished. Output is stored at '{output_folder}'")

    def _run_r_script(self, script_location, input_folder, output_folder):
        metaDataFile = open(f"{input_folder}metadata.json", "r")
        pipelineSampleName = json.load(metaDataFile)["tumor"]["sampleName"]

        subprocess.run(['Rscript', script_location, pipelineSampleName, input_folder, output_folder], check=False)



    def _get_germline_reportable_genes(self,fn):
        germline_reportable_genes = []
        fh = open(fn,'r')

        for lines in fh:
            lines = lines.rstrip('\n')
            germline_reportable_genes.append(lines)
        return germline_reportable_genes

    def _get_gene_from_line(self,line):
        arr = line.split('\t')
        arr2 = arr[7].split(';')

        for i in arr2:
            if i[0:6]=='IMPACT':
                arr3 = i.split(',')
                arr4 = arr3[0].split('=')
                return(arr4[1])
        return('')

    def _generate_vcf(self, input_folder, output_folder):
        reported_res = []
        annotated_res = []

        metaDataFile = open(f"{input_folder}metadata.json", "r")
        pipelineSampleName = json.load(metaDataFile)["tumor"]["sampleName"]

        gene_list_path = '/data/resources/reporting-resources/panel/genenames.pass.pon.csv'
        germline_genes = self._get_germline_reportable_genes(gene_list_path)

        with gzip.open(f"{input_folder}purple/{pipelineSampleName}.purple.somatic.vcf.gz", 'rt') as file:
            for line in file.readlines():
                if "REPORTED" in line or (len(line) > 0 and line[0] == '#'):
                    reported_res.append(line)
                if len(line) > 0 and line[0] == '#':
                    annotated_res.append(line)
                else:
                    arr = line.split('\t')
                    if arr[6] == "PASS":
                       annotated_res.append(line)

                    elif (not 'PONArtefact' in arr[6]) and (self._get_gene_from_line(line) in germline_genes):
                       annotated_res.append(line)

        if reported_res:
            with open(f"{output_folder}{pipelineSampleName}.reported.somatic.vcf", 'x') as file:
                file.writelines(reported_res)
        if annotated_res:
            with open(f"{output_folder}{pipelineSampleName}.gnomad.pon.vcf", 'x') as file:
                file.writelines(annotated_res)

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
            "purple.driver.catalog.somatic.tsv",
            "purple.somatic.vcf.gz",
            "purple.purity.tsv",
            "purple.somatic.vcf.gz",
            "purple.cnv.gene.tsv",
            "sage.exon.medians.tsv",
            "genenames.pass.pon.csv",
            "metadata.json"
        }
        run_blobs = list(self._panel_pipeline_output_bucket.list_blobs(prefix=self.set_name))

        result = []
        for suffix in required_file_suffixes:
            for run_blob in run_blobs:

                if run_blob.name[-len(suffix):] == suffix:  # this checks if the blob name ends with the suffix.
                    result.append(run_blob)

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
