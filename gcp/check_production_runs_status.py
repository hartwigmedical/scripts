import json

import requests
import subprocess


def main():
    """
    Goes through all the production runs and will give a summary in different categories which
    will be printed to the stdout.
    """
    print('Retrieving all non-research somatic and CPCT runs from the API...')
    all_non_research_runs = get_all_runs_from_api()

    def get_run_by(status: str, version: str = "5.31") -> list[json]:
        return [run for run in all_non_research_runs if get_run_status(run) == status]  # and run['version'] == version]

    validated_runs = get_run_by(status='Validated')
    fail_runs = get_run_by(status='Failed')
    fini_runs = get_run_by(status='Finished')
    proc_runs = get_run_by(status='Pending') + get_run_by(status='Processing')
    wait_runs = get_run_by(status='Waiting')

    tumor_samples_with_validated_runs = get_sample_names_from_runs(validated_runs)
    tumor_sample_with_finished_runs = get_sample_names_from_runs(fini_runs)
    tumor_samples_with_processing_runs = get_sample_names_from_runs(proc_runs)
    tumor_samples_with_waiting_runs = get_sample_names_from_runs(wait_runs)

    _handle_failed_runs(fail_runs, tumor_samples_with_validated_runs, tumor_sample_with_finished_runs,
                        tumor_samples_with_processing_runs, tumor_samples_with_waiting_runs)
    _handle_finished_runs(fini_runs, tumor_samples_with_validated_runs, tumor_samples_with_processing_runs,
                          tumor_samples_with_waiting_runs)
    _handle_validated_runs(validated_runs)
    _handle_processing_runs(proc_runs)
    _handle_waiting_runs(wait_runs)

    print("All done ʕっ•ᴥ•ʔっ♪♬ !")
    exit(0)


def _handle_failed_runs(failed_runs, tumor_samples_with_validated_runs, tumor_sample_with_finished_runs,
                        tumor_samples_with_processing_runs, tumor_samples_with_waiting_runs):
    failed_but_not_reported = []
    failed_but_already_reported = []
    for run in failed_runs:
        tumor_sample = get_sample_name_from_run(run)
        if not tumor_sample:
            continue
        shared_count = 0

        if not (tumor_sample in tumor_samples_with_validated_runs or
                tumor_sample in tumor_sample_with_finished_runs or
                tumor_sample in tumor_samples_with_processing_runs or
                tumor_sample in tumor_samples_with_waiting_runs):
            if shared_count > 0:
                failed_but_already_reported.append(run)
                continue
            failed_but_not_reported.append(run)

    print("### FAILED RUNS THAT NOT HAVE BEEN REPORTED BEFORE ( please solve! ):")
    for i, run in enumerate(failed_but_not_reported):
        set_name = get_set_name_from_run(run)
        if not set_name:
            continue
        db_status = get_db_status(run)
        research_db_status = get_research_db_status(run)
        # TODO this is not complete yet
        print(f"""
            ---------
            -{i}- {set_name} has failed. healthchecker issue? (see below)
            {get_current_tat(run)}
            db_status:diagnostic={db_status}/research={research_db_status}
            Action - inspect the error above and solve if possible with lab (extra sequencing for coverage and/or second biopt)...
            ...  ONLY if solving is not possible, proceed to the next steps: ..
            Action - check whether the SNP check was done and was ok: perform_snpcheck_run ${set_name}
            Action - patch the api to validated: patch_api_run_validate ${set_name} 
            """
              )

    print("### FAILED RUNS THAT ALREADY HAVE BEEN REPORTED ( no action required ):")
    for i, run in enumerate(failed_but_already_reported):
        set_name = get_set_name_from_run(run)
        db_status = get_db_status(run)
        research_db_status = get_research_db_status(run)
        print(
            f"""
            ---------
            -${i}- ${set_name} has Failed but has been reported before
            # db_status:diagnostic={db_status}/research={research_db_status}"
            {get_reported_info_sample(run)}
            """)


def _handle_finished_runs(finished_runs, tumor_samples_with_validated_runs, tumor_samples_with_processing_runs,
                          tumor_samples_with_waiting_runs):
    finished_but_not_reported = []
    finished_but_already_reported = []
    for run in finished_runs:
        tumor_sample = get_sample_name_from_run(run)
        if not tumor_sample:
            continue
        reporting_id = get_reporting_id(run)
        if not reporting_id:
            continue
        shared_count = get_shared_count(reporting_id)
        if not (tumor_sample in tumor_samples_with_validated_runs or
                tumor_sample in tumor_samples_with_processing_runs or
                tumor_sample in tumor_samples_with_waiting_runs):
            if shared_count > 0:
                finished_but_already_reported.append(run)
                continue
            finished_but_not_reported.append(run)

    print('### FINISHED RUNS THAT HAVE NOT BEEN REPORTED BEFORE ( please solve! ):')
    for i, run in enumerate(finished_but_not_reported):
        set_name = get_set_name_from_run(run)
        get_current_tat(run)
        db_status = get_db_status(run)
        research_db_status = get_research_db_status(run)
        print(  # TODO this is not complete yet does snpcheck have to run?
            f"""
            ---------
            -{i}- {set_name} has been finished but not validated. snpcheck issue (see below)?
            {get_current_tat(run)}
            db_status:diagnostic={db_status}/research={research_db_status}
            """)

    print('### FINISHED RUNS THAT HAVE ALREADY BEEN REPORTED ( no action required ):')
    for i, run in enumerate(finished_but_already_reported):
        set_name = get_set_name_from_run(run)
        db_status = get_db_status(run)
        research_db_status = get_research_db_status(run)
        print(f"""
        ---------
        -{i}- -{set_name}- has Finished but has been reported before:
        # db_status:diagnostic={db_status}/research={research_db_status}"
        {get_reported_info_sample(run)}
        """)


def _handle_validated_runs(validated_runs):
    validated_runs_incorrect_summary = []
    validated_runs_with_warnings = []
    validated_runs_no_report = []
    validated_runs_health_error = []

    for i, run in enumerate(validated_runs):
        if i % 25 == 0:
            print(f'Processing validated runs... ({i}/{len(validated_runs)})')
        set_name = get_set_name_from_run(run)
        if not set_name:
            continue

        set_samples = get_samples_from_run(run)
        if not set_samples or 'tumor' not in set_samples or 'barcode' not in set_samples['tumor']:
            continue
        reporting_id = get_reporting_id(run)
        if not reporting_id:
            continue
        shared_count = get_shared_count(reporting_id)
        validation_ss_run = "_HMFregVAL_" in set_name
        if not (shared_count == 0 and validation_ss_run):
            continue

        sample_name = get_sample_name_from_run(run)
        log_file_path = get_reporting_log_file_for_sample_from_api(sample_name)

        log_contents = subprocess.check_output(['gsutil', 'cat', log_file_path]).decode()
        rose_error = is_rose_error(sample_name, get_sample_name_from_run(run))
        doid_error = get_doid_error(set_name)
        health_error = get_health_error(set_name)

        if 'summary' in log_contents or rose_error:
            validated_runs_incorrect_summary.append(run)
        elif 'WARN' in log_contents or ('The primary tumor location provided is Unknown' in doid_error and
                                        'The primary tumor location provided is Unknown,'
                                        ' but it is a CUP so this is correct' not in doid_error):
            validated_runs_with_warnings.append(run)
        elif 'WARN' in health_error:
            validated_runs_health_error.append(run)
        else:
            validated_runs_no_report.append(run)

    print("### VALIDATED RUNS WITHOUT WARNINGS BUT INCORRECT SUMMARY ( please solve! ):")
    for i, run in enumerate(validated_runs_incorrect_summary):
        set_name = get_set_name_from_run(run)
        db_status = get_db_status(run)
        research_db_status = get_research_db_status(run)
        print(f"""
            -{i}- ${set_name} is validated, there were no information warnings while creating the report but the summary is missing
            # db_status:diagnostic={db_status}/research={research_db_status}
            {get_current_tat(run)}
            """)


def _handle_processing_runs(processing_runs):
    print('### PROCESSING/UPLOADING/DOWNLOADING RUNS ( no action required ):')
    for i, run in enumerate(processing_runs):
        set_name = get_set_name_from_run(run)
        status = get_run_status(run)
        print(f"""
        -{i}- {set_name} is {status}
        {get_current_tat(run)}
        """)


def _handle_waiting_runs(waiting_runs):
    print('### WAITING RUNS ( no action required ):')
    for i, run in enumerate(waiting_runs):
        set_name = get_set_name_from_run(run)
        status = get_run_status(run)
        print(f"""
        -{i}- {set_name} is {status}
        {get_current_tat(run)}
        """)


def get_all_runs_from_api() -> list[json]:
    """
    Gets all the non-research somatic CPCT and Somatic runs from the API
    :return: a json list containing all the relevant runs
    :raises ValueError: If the request returned a non 2XX status code
    """
    all_cpct_runs = requests.get('http://api.prod-1/hmf/v1/runs?ini=CPCT.ini')
    if not all_cpct_runs.ok:
        raise ValueError(f'Unable to get CPCT runs: {all_cpct_runs.status_code}')
    all_somatic_runs = requests.get('http://api.prod-1/hmf/v1/runs?ini=Somatic.ini')
    if not all_somatic_runs.ok:
        raise ValueError(f'Unable to get Somatic runs: {all_somatic_runs.status_code}')

    all_runs = all_cpct_runs.json() + all_somatic_runs.json()
    non_research_runs = [run for run in all_runs if run['context'] != 'RESEARCH']
    return non_research_runs


def get_reporting_log_file_for_sample_from_api(sample_name: str):
    """
    Gets the log file associated with the patient_reporter stage for a given sample
    :param sample_name: the name of the sample
    :return: the path of the log file in the following format: gs://<bucket-name>/<blob-name>
    :raises: ValueError:
            - if the initial request returns a non 2XX status code
            - if the initial request returns more than 1 result
            - if the returned data does not contain the path to the log file
    """
    response = requests.get('http://api.prod-1/hmf/v1/reports/2/created', params={'sample_name': sample_name})
    if not response.ok:
        raise ValueError(f'Could not get created report for {sample_name}')
    json_response = response.json()
    if len(json_response) != 1:
        raise ValueError(f'Request returned too many results: {len(json_response)}')

    report_created_data = response.json()[0]
    report_files = report_created_data['report_files']

    logs = [r for r in report_files if r['datatype'] == 'report_log']
    if len(logs) == 0:
        raise ValueError(f'Log file not found!')
    if len(logs) > 1:
        print(f'Multiple log files detected ("{len(logs)}"), returning the first one')
    return logs[0]


def get_sample_name_from_run(run: json) -> str:
    """
    Gets the tumor sample name from a given run.

    :param run: the run to get the tumor sample from.
    :return: the tumor sample.
    """
    return run['set']['tumor_sample']


def get_sample_names_from_runs(runs: list[json]) -> set[str]:
    """
    Gets the tumor sample name strings from a list of runs.

    :param runs: the runs to get the tumor samples from.
    :return: a set containing all the unique tumor sample strings.
    """
    return {get_sample_name_from_run(r) for r in runs}


def get_set_name_from_run(run: json) -> str:
    """
    Get the set name from a run

    :param run: the run
    :return: the set name
    """
    return run['set']['name']


def get_samples_from_run(run) -> {str: json}:
    """
    Gets the samples from the set information.

    :param run: the run to get the samples for.
    :return: A dictionary where the key is the sample type ('ref', 'tumor', etc.) and the value is the sample json
    if samples are found, `None` otherwise
    """
    set_id = run['set']['id']
    response = requests.get("http://api.prod-1/hmf/v1/sets", params={'sample_id': set_id})
    if not response.ok:
        raise ValueError(f'Request failed: {response.status_code}')
    json_response = response.json()
    if len(json_response) == 0:
        return None

    entry = json_response[-1]
    samples: list[json] = entry['samples']
    return {s['type']: s for s in samples}


def get_shared_count(reporting_id: str) -> int:
    """
    Gets the amount of times a report has been shared before.

    :param reporting_id: the id of the report to check.
    :return: the amount of times the given report has been shared before.
    """
    result = subprocess.check_output(
        ['bash', '-c', f'source locate_reporting_api; extract_object_shared {reporting_id} prod'])
    words = result.split()
    return len(words)


def get_reporting_id(run: json) -> str:
    """
    Get the reporting_id for a given run

    :param run: the run to get the reporting_id for
    :return: the reporting_id if found, `None` otherwise
    """
    isolation_barcode = get_tumor_barcode(run)
    if isolation_barcode:
        reporting_id = subprocess.check_output(
            ['bash', '-c', 'source locate_reporting_api; extract_most_recent_reporting_id_on_barcode',
             isolation_barcode, 'prod']).decode()
        if 'null' in reporting_id:
            return None
        return reporting_id


def get_doid_error(set_name: str) -> str:
    return subprocess.check_output(
        ['doid_check_validated_run', set_name, '2>&1']).decode()


def get_health_error(set_name: str) -> str:
    return subprocess.check_output(['health_check_validated_run', set_name, '2>&1']).decode()


def is_rose_error(set_name: str, sample_name: str):
    """
    Finds whether there is a error with the 'rose' tool for a given set and sample

    :param set_name: the set to check for
    :param sample_name: the sample to check for within the set
    :return: true if a 'rose' error was found, else false
    """
    res = subprocess.check_output(['gsutil', 'cat',
                                   f'gs://diagnostic-pipeline-output-prod-1/{set_name}/purple/{sample_name}.driver.catalog.somatic.tsv']).decode()
    return res.find('AMP') > 0


def get_db_status(run: json) -> str:
    return run['db_status']


def get_run_status(run: json) -> str:
    return run['status']


def get_research_db_status(run):
    """
    Gets the research db status of a given run.

    :param run: the run.
    :return: the research db status if found, else `None`.
    """
    set_name = get_set_name_from_run(run)
    response = requests.get('http://api.prod-1/hmf/v1/runs',
                            params={'set_name': set_name, 'bucket': 'research-pipeline-output-prod-1'})
    if not response.ok:
        raise ValueError(f"Request returned status: {response.status_code}")
    json_response = response.json()
    if len(json_response) == 0:
        return None
    return response.json()[-1]['db_status']


def get_tumor_barcode(run: json) -> str:
    """
    Gets the tumor barcode of a given run
    :param run: the run to get the tumor barcode for
    :return: the tumor barcode
    """
    samples = get_samples_from_run(run)
    if samples and 'tumor' in samples and 'barcode' in samples['tumor']:
        return samples['tumor']['barcode']


def get_current_tat(run: json) -> str:
    """
    Get the current tat for a run (TODO: not sure what this is?)
    :param run: the run to get the tat for
    :return: the tat if it is there, else an empty string
    """
    samples = get_samples_from_run(run)
    if samples and 'tumor' in samples:
        sample_name = samples['tumor']['name']
        res = subprocess.check_output(
            ['bash', '-c', f'../oncoact/patientreporter/ops-vm/oncoact/get_current_TAT.sh {sample_name}']).decode()
        if res:
            return res
    return ''


def get_reported_info_sample(run: json) -> str:
    """
    Gets the reported info for a given run (TODO: not sure what this is?)
    :param run: the run to get the reported info for
    :return: the reported info if any, else an empty string
    """
    sample_name = get_sample_name_from_run(run)
    res = subprocess.check_output(
        ['bash', '-c', f'../oncoact/patientreporter/ops-vm/get_reported_info_sample {sample_name}']).decode()
    if res:
        return res
    return ''


if __name__ == '__main__':
    main()
