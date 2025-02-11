# /path/to/your_script.py
import yaml
import argparse
import subprocess
import os

class MyDumper(yaml.Dumper):
    def increase_indent(self, flow=False, indentless=False):
        return super(MyDumper, self).increase_indent(flow, False)

def CheckExt(choices):
    class Act(argparse.Action):
        def __call__(self, parser, namespace, fname, option_string=None):
            ext = os.path.splitext(fname)[1][1:]
            if ext not in choices:
                option_string = f'({option_string})' if option_string else ''
                parser.error(f"File {option_string} doesn't end with one of {choices}")
            else:
                setattr(namespace, self.dest, fname)
    return Act

def get_file_paths_barcode(fastq_directory, input_file):
    file_path_dict = {}
    
    with open(input_file, 'r') as f:
        sample_entries = [line.strip().split(',') for line in f]
    
    print(f"Debug: Found {len(sample_entries)} sample entries")

    result = subprocess.run(['gsutil', 'ls', '-lR', str(fastq_directory)], stdout=subprocess.PIPE).stdout.decode('utf-8').splitlines()
    
    file_count = 0
    for file in result:
        file_count += 1
        file_path = file.split(" ")[-1]
        file_name = os.path.basename(file_path)
        for entry in sample_entries:
            if len(entry) == 4:
                tumor_barcode, reference_barcode, tumor_sample_id, reference_sample_id = entry
                if (f"{tumor_barcode}") in file_name or (f"{tumor_sample_id}") in file_name:
                    if tumor_sample_id not in file_path_dict:
                        file_path_dict[tumor_sample_id] = {'normal': [], 'tumor': []}
                    file_path_dict[tumor_sample_id]['tumor'].append(file_path)
                elif (f"{reference_barcode}") in file_name or f"{reference_sample_id}") in file_name:
                    if reference_sample_id not in file_path_dict:
                        file_path_dict[reference_sample_id] = {'normal': [], 'tumor': []}
                    file_path_dict[reference_sample_id]['normal'].append(file_path)
            elif len(entry) in [1, 2]:
                sample_id = entry[0]
                if (f"{sample_id}R") in file_name:
                    if sample_id not in file_path_dict:
                        file_path_dict[sample_id] = {'normal': [], 'tumor': []}
                    file_path_dict[sample_id]['normal'].append(file_path)
                elif (f"{sample_id}T") in file_name:
                    if sample_id not in file_path_dict:
                        file_path_dict[sample_id] = {'normal': [], 'tumor': []}
                    file_path_dict[sample_id]['tumor'].append(file_path)

    return file_path_dict

def create_data_dict(list_of_dicts, sample_ids=None):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    options_file = os.path.join(script_dir, 'options.yaml')
    
    with open(options_file, 'r') as file:
        options_data = yaml.safe_load(file)
    
    if sample_ids is not None:
        options_data['sampleIds'] = sample_ids
    else:
        options_data['samples'] = list_of_dicts

    return options_data

def main(input_file, directory, output_file, bam_flag, tumor_only_flag, sample_id_flag):
    list_of_dicts = []
    sample_ids = []
    
    if sample_id_flag:
        with open(input_file, 'r') as f:
            sample_ids = [line.strip().split(',')[0] for line in f]
        data = create_data_dict(list_of_dicts, sample_ids)
    else:
        file_path_dict = get_file_paths_barcode(directory, input_file)
        with open(input_file, "r") as f:
            for line in f:
                entry = line.strip().split(",")
                per_sample_dict = dict()
                
                if len(entry) == 4:
                    tumor_barcode, reference_barcode, tumor_sample_id, reference_sample_id = entry
                    sample_name_processed = tumor_sample_id[:-1]

                    if not tumor_only_flag:
                        if not bam_flag:
                            tumor_fastq_files = sorted(file_path_dict.get(tumor_sample_id, {}).get('tumor', []))
                            tumor_fastq_pairs = []
                            for i in range(0, len(tumor_fastq_files), 2):
                                if i + 1 < len(tumor_fastq_files):
                                    tumor_fastq_pairs.append(dict(read1=tumor_fastq_files[i], read2=tumor_fastq_files[i+1]))

                            normal_fastq_files = sorted(file_path_dict.get(reference_sample_id, {}).get('normal', []))
                            normal_fastq_pairs = []
                            for i in range(0, len(normal_fastq_files), 2):
                                if i + 1 < len(normal_fastq_files):
                                    normal_fastq_pairs.append(dict(read1=normal_fastq_files[i], read2=normal_fastq_files[i+1]))

                            per_sample_dict = dict(
                                name=sample_name_processed,
                                tumors=[dict(
                                    name=tumor_sample_id,
                                    fastq=tumor_fastq_pairs
                                )],
                                normal=dict(
                                    name=reference_sample_id,
                                    fastq=normal_fastq_pairs
                                )
                            )
                        else:
                            per_sample_dict = dict(
                                name=sample_name_processed,
                                tumors=[dict(
                                    name=tumor_sample_id,
                                    bam=file_path_dict[tumor_sample_id].get('tumor', [])
                                )],
                                normal=dict(
                                    name=reference_sample_id,
                                    bam=file_path_dict[reference_sample_id].get('normal', [])
                                )
                            )
                else:
                    sample_id = entry[0]
                    sample_name_processed = sample_id if len(entry) == 1 else "".join(entry[1].split("-"))

                    if tumor_only_flag:
                        if not bam_flag:
                            fastq_files = sorted(file_path_dict.get(sample_id, {}).get('tumor', []))
                            fastq_pairs = []
                            for i in range(0, len(fastq_files), 2):
                                if i + 1 < len(fastq_files):
                                    fastq_pairs.append(dict(read1=fastq_files[i], read2=fastq_files[i+1]))
                            per_sample_dict = dict(
                                name=sample_name_processed,
                                tumors=[dict(
                                    name=sample_name_processed,
                                    fastq=fastq_pairs
                                )]
                            )
                        else:
                            per_sample_dict = dict(
                                name=sample_name_processed,
                                tumors=[dict(
                                    name=sample_name_processed,
                                    bam=file_path_dict[sample_id].get('tumor', [])
                                )]
                            )
                    else:
                        if not bam_flag:
                            tumor_fastq_files = sorted(file_path_dict.get(sample_id, {}).get('tumor', []))
                            tumor_fastq_pairs = []
                            for i in range(0, len(tumor_fastq_files), 2):
                                if i + 1 < len(tumor_fastq_files):
                                    tumor_fastq_pairs.append(dict(read1=tumor_fastq_files[i], read2=tumor_fastq_files[i+1]))

                            normal_fastq_files = sorted(file_path_dict.get(sample_id, {}).get('normal', []))
                            normal_fastq_pairs = []
                            for i in range(0, len(normal_fastq_files), 2):
                                if i + 1 < len(normal_fastq_files):
                                    normal_fastq_pairs.append(dict(read1=normal_fastq_files[i], read2=normal_fastq_files[i+1]))

                            per_sample_dict = dict(
                                name=sample_name_processed,
                                tumors=[dict(
                                    name=sample_name_processed,
                                    fastq=tumor_fastq_pairs
                                )],
                                normal=dict(
                                    name=sample_name_processed,
                                    fastq=normal_fastq_pairs
                                )
                            )
                        else:
                            per_sample_dict = dict(
                                name=sample_name_processed,
                                tumors=[dict(
                                    name=sample_name_processed,
                                    bam=file_path_dict[sample_id].get('tumor', [])
                                )],
                                normal=dict(
                                    name=sample_name_processed,
                                    bam=file_path_dict[sample_id].get('normal', [])
                                )
                            )

                list_of_dicts.append(per_sample_dict)
        data = create_data_dict(list_of_dicts)

    with open(output_file, "w") as outfile:
        yaml.dump(data, outfile, Dumper=MyDumper, sort_keys=False, default_flow_style=False)
    print("Yaml file created!")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',
        help='Name of txt file containing barcodes and/or sample IDs',
        action=CheckExt({'txt'}))
    parser.add_argument('-d', '--directory',
        help='Path to directory containing input files (.bam or .fastq.gz)',
        type=lambda s: s if s.endswith('/') else s + '/',
        required=False)
    parser.add_argument('-o', '--output',
        help='Name of output yaml file',
        action=CheckExt({'yaml'}))
    parser.add_argument('-b', '--bam',
        help='Flag to indicate if running from bam files',
        action='store_true')
    parser.add_argument('-t', '--tumor_only',
        help='Flag to indicate if running in tumor only mode',
        action='store_true')
    parser.add_argument('-s', '--sample_id_mode',
        help='Flag to indicate if running from sample IDs directly without bucket directory',
        action='store_true')
    args = parser.parse_args()
    main(args.input, args.directory, args.output, args.bam, args.tumor_only, args.sample_id_mode)
