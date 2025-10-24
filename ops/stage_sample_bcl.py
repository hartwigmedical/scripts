#!/usr/bin/env python3

import argparse
import re
from dataclasses import dataclass
from typing import List, Union
import math

parser = argparse.ArgumentParser(description="Stages the BCL sample conversion results for the sample sheet excel")
parser.add_argument("file_path", nargs='+')


@dataclass
class SampleInput:
    flowcell_name: str
    submission: str
    isolation_barcode: str
    sample_name: str
    reporting_id: str
    output_type: str
    sequencing_result_status: str
    set_name: str


@dataclass
class SampleParseError:
    data: SampleInput
    error: str


@dataclass
class Sample:
    submission: str
    isolation_barcode: str
    sample_name: str
    reporting_id: str
    output_type: str
    sample_type: str
    sample_status: str
    set_name: str
    q30: float
    yield_cur: float
    yield_total: float
    yield_requirement: float

    def by_seq_in_gb(self):
        return math.ceil(max(self.yield_requirement - self.yield_total, 0))


def main():
    args = parser.parse_args()
    sample_input = [s for fp in args.file_path for s in read_input(fp)]
    processed = [process_input(si) for si in sample_input]

    samples = [s for s in processed if isinstance(s, Sample)]
    err = [e for e in processed if isinstance(e, SampleParseError)]

    print_output_for_excel(samples)
    print_err(err)


def read_input(path: str) -> List[SampleInput]:
    result = []
    with open(path, 'r') as samples:
        for line in samples.readlines():
            seperated = line.split('\t')
            result.append(
                SampleInput(
                    flowcell_name=seperated[0].strip(),
                    submission=seperated[1].strip(),
                    isolation_barcode=seperated[2].strip(),
                    sample_name=seperated[3].strip(),
                    reporting_id=seperated[4].strip(),
                    output_type=seperated[5].strip(),
                    sequencing_result_status=seperated[6].strip(),
                    # index 7 is skipped intentionally since this contains team information which is never used
                    set_name=seperated[8].strip()
                )
            )
    return result


def process_input(sample_input: SampleInput) -> Union[Sample, SampleParseError]:
    try:
        m = re.match(r"q=(\d+(?:\.\d+)?)\s+y=(\d+(?:\.\d)*)\+(\d+(?:\.\d)*)/(?:(\d+(?:\.\d)*)|\?)\s+(.+)", sample_input.sequencing_result_status)

        q30 = float(m.group(1))
        yield_cur = float(m.group(2))
        yield_total = float(m.group(3))
        yield_required = float(m.group(4))

        sample_status, run_status, ini = m.group(5).split("|")

        if re.search(r"T\d*$", sample_input.sample_name):
            sample_type = 'TUMOR'
        elif re.search(r"R\d*$", sample_input.sample_name):
            sample_type = 'REF'
        else:
            sample_type = 'OTHER'

        return Sample(
            submission=sample_input.submission,
            isolation_barcode=sample_input.isolation_barcode,
            sample_name=sample_input.sample_name,
            reporting_id=sample_input.reporting_id,
            output_type=sample_input.output_type,
            sample_type=sample_type,
            sample_status=sample_status,
            set_name=sample_input.set_name,
            q30=q30,
            yield_cur=yield_cur,
            yield_total=yield_total,
            yield_requirement=yield_required
        )
    except Exception as e:
        return SampleParseError(
            data=sample_input,
            error=str(e)
        )

def print_err(err: List[SampleParseError]):
    print('\n\n')
    print('Error parsing the following entries:')
    for e in err:
        print(e.data, e.error)


def print_output_for_excel(samples: List[Sample]):
    results = []
    by_seq = []
    for s in samples:
        if s.by_seq_in_gb() > 0:
            result = compute_result(
                s,
                f'Extra Seq {s.by_seq_in_gb()} GBase',
                'AB'
            )
            by_seq.append(result)
        else:
            if s.sample_type == 'TUMOR' and s.output_type in ['Rna', 'Targeted']:
                result = compute_result(
                    s,
                    'Processing',
                    'BFX'
                )
                results.append(result)
            elif s.sample_type == 'TUMOR' and s.output_type in ['Somatic', 'ShallowSeq']:
                ref = find_ref_or_none(s, samples)
                if ref and ref.by_seq_in_gb() > 0:
                    result = compute_result(
                        s,
                        'Waiting on R',
                        'AB'
                    )
                    results.append(result)
                else:
                    result = compute_result(
                        s,
                        'Processing',
                        'BFX'
                    )
                    results.append(result)
            elif s.output_type == 'FastQ':
                result = compute_result(
                    s,
                    'Done',
                    ''
                )
                results.append(result)
            else:
                # samples that don't belong in sample-overview are skipped (i.e., ref samples)
                continue

    results.sort(key=lambda x: x[2])  # sort by sample name
    results = ['\t'.join(r) for r in results]

    print(*results, sep='\n')
    print("\n\n\n")
    print("The following samples need by-seq:", *['\t'.join(bs) for bs in by_seq], sep="\n")


def compute_result(s: Sample, output_status: str, team: str):
    return (
        s.submission,
        s.isolation_barcode,
        s.sample_name,
        s.reporting_id,
        s.output_type,
        output_status,
        team,
        '',
        '',
        '',
        '',
        '',
        s.set_name
    )


def find_ref_or_none(tumor: Sample, samples: List[Sample]):
    base = re.sub(r"T\d*$", "", tumor.sample_name)
    pattern = f"{base}R\\d*$"
    return next((s for s in samples if re.match(pattern, s.sample_name)), None)


if __name__ == '__main__':
    main()
