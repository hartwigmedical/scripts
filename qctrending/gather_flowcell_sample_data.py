#!/usr/bin/env python3

import argparse
import concurrent.futures
import csv
import datetime
import json
import subprocess
import sys

from typing import Dict, List, NamedTuple, Optional


def main():
    config = parse_args()
    flowcells = get_flowcells()
    flowcells = select_flowcells(flowcells, config.start_date, config.end_date)
    samples = get_samples(flowcells)
    write_complete_samples_info(samples, config.output_file)


def get_flowcells() -> List[Dict]:
   bashCommand = "hmf_api_get flowcells"
   flowcells = subprocess.run(bashCommand.split(), capture_output=True, text=True)
   if flowcells.stderr != '':
       print(f"[ERROR]: {flowcells.stderr}")
       sys.exit()
   return json.loads(flowcells.stdout)

def check_flowcell_time_range(flowcell: Dict, s: datetime.date, e: datetime.date) -> Optional[bool]:
    convert_time = flowcell["convertTime"]
    if convert_time == None:
       return False
    convert_date = datetime.datetime.strptime(convert_time, "%Y-%m-%dT%H:%M:%S").date()
    if convert_date >= s and convert_date <= e:
        return True
    else:
        return False

def select_flowcells(flowcells: List[Dict], s: datetime.date, e: datetime.date) -> List[Dict]:
    selected_flowcells = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        future_to_flowcell = {executor.submit(check_flowcell_time_range, flowcell, s, e): flowcell for flowcell in flowcells}
        for future in concurrent.futures.as_completed(future_to_flowcell):
            flowcell = future_to_flowcell[future]
            try:
                check = future.result()
            except Exception as exc:
                print(f"[WARN]: Flowcell {flowcell['id']} generated an exception: {exc}")
            else:
                if check:
                    selected_flowcells.append(flowcell)
    return selected_flowcells


def get_samples_from_flowcell(flowcell_id: str) -> List[Dict]:
    bashCommand = f"hmf_api_get samples?flowcell_id={flowcell_id}"
    samples = subprocess.run(bashCommand.split(), capture_output=True, text=True)
    return json.loads(samples.stdout)

def get_runs_info(sample_barcode: str) -> List[Dict]:
    bashCommand = f"hmf_api_get runs?barcode={sample_barcode}"
    runs = subprocess.run(bashCommand.split(), capture_output=True, text=True)
    runs = json.loads(runs.stdout)
    return runs

def select_valid_run(runs: List[Dict]) -> Dict:
    if len(runs) == 0:
        return {"context": None, "ini": None, "status": None, "failure": None, "coverage": None, "bqr": None}
    candidate_runs = []
    for run in runs:
        if run["context"] == "DIAGNOSTIC" and run["ini"] == "Somatic.ini" and (run["status"] == "Validated" or run["status"] == "Deleted" or run["status"] == "Failed"):
            candidate_runs.append(run)
    if len(candidate_runs) == 0:
        run = runs[0]
        return run
    elif len(candidate_runs) == 1:
        return candidate_runs[0]
    else:
        return max(candidate_runs, key=lambda k: k["startTime"]) # Most recent pipeline run that is valid

def get_sample_run_qc_info(sample_barcode: str) -> Dict:
    runs_info = get_runs_info(sample_barcode)
    run = select_valid_run(runs_info)
    return run

def get_samples_qc_info(samples: List[Dict]) -> List[Dict]:
    samples_qc = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        future_to_sample = {executor.submit(get_sample_run_qc_info, sample["barcode"]): sample for sample in samples}
        for future in concurrent.futures.as_completed(future_to_sample):
            sample = future_to_sample[future]
            try:
                run = future.result()
            except Exception as exc:
                print(f"[WARN]: Sample {sample['id']} generated an exception: {exc}")
            else:
                samples_qc.append({**sample, **run})
    return samples_qc


def merge_sample_and_flowcell_info(sample: Dict, flowcell: Dict) -> Dict:
    if sample["failure"] == None:
        sample["failure"] = {"type": None}
    return {"flowcell": flowcell["flowcell_id"],
            "sequencer": flowcell["sequencer"],
            "date": flowcell["convertTime"],
            "undetermined_percentage": flowcell["undet_rds_p"],
            "flowcell_yield": flowcell["yld"],
            "flowcell_q30": flowcell["q30"],
            "submission": sample["submission"],
            "sample": sample["barcode"],
            "type": sample["type"],
            "sample_q30": sample["q30"],
            "sample_yield": sample["yld"],
            "context": sample["context"],
            "ini": sample["ini"],
            "status": sample["status"],
            "run_qc_status": sample["failure"]["type"],
            #"version": sample["version"]
            #"coverage": sample["coverage"],
            #"bqr": sample["bqr"]
           }

def merge_samples_and_flowcell_info(samples: List[Dict], flowcell: Dict) -> List[Dict]:
    complete_sample_info = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        future_to_sample = {executor.submit(merge_sample_and_flowcell_info, sample, flowcell): sample for sample in samples}
        for future in concurrent.futures.as_completed(future_to_sample):
            try:
                data = future.result()
            except Exception as exc:
                print(f"[WARN]: Sample generated an exception: {exc}")
            else:
                complete_sample_info.append(data)
    return complete_sample_info


def get_samples(flowcells: List[Dict]):
    samples = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        future_to_flowcell = {executor.submit(get_samples_from_flowcell, flowcell["id"]): flowcell for flowcell in flowcells}
        for future in concurrent.futures.as_completed(future_to_flowcell):
            flowcell = future_to_flowcell[future]
            try:
                samples_data = future.result()
            except Exception as exc:
                print(f"[WARN]: Flowcell {flowcell['id']} generated an exception: {exc}")
            else:
                data = get_samples_qc_info(samples_data)
                data = merge_samples_and_flowcell_info(data, flowcell)
                for sample in data:
                   samples.append(sample)
    return samples


def parse_args() -> NamedTuple:
    parser = argparse.ArgumentParser(prog="gather_data_qc_trends.py", description="Gather flowcell/sample qc data within given time period")
    parser.add_argument("--start_date", "-s", type=valid_date, required=True, help="start date in format YYYY-MM-DD")
    parser.add_argument("--end_date", "-e", type=valid_date, required=True, help="end date in format YYYY-MM-DD")
    parser.add_argument("--output_file", "-o", type=str, required=True, help="path to output file")
    args = parser.parse_args()
    check_dates(args.start_date, args.end_date)
    return args

def valid_date(date: str) -> datetime.date:
    try:
        return datetime.datetime.strptime(date, "%Y-%m-%d").date()
    except ValueError:
        raise argparse.ArgumentTypeError(f"ERROR: not a valid date: {date}")

def check_dates(s: datetime.date, e: datetime.date):
    if s > e:
        print(f"[ERROR]: start_date {s} greater than end_date {e}")
        sys.exit()
    if s >= datetime.date.today() or e >= datetime.date.today():
        print(f"[ERROR]: start_date or end_date greater than current date")
        sys.exit()

def write_complete_samples_info(samples: List[Dict], outfile: str):
    with open(outfile, 'w') as output_file:
        dw = csv.DictWriter(output_file, sorted(samples[0].keys()), delimiter='\t')
        dw.writeheader()
        dw.writerows(samples)


if __name__ == "__main__":
    main()
