#!/usr/bin/env python3

import argparse
import re
from dataclasses import dataclass
from typing import List

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
    yield_cur: int
    yield_total: int
    yield_requirement: int

    def by_seq_in_gb(self):
        return max(self.yield_requirement - self.yield_total, 0)


def main():
    args = parser.parse_args()
    samples = [s for fp in args.file_path for s in read_input(fp)]
    print_output_for_excel([process_input(s) for s in samples])


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


def process_input(sample_input: SampleInput) -> Sample:
    m = re.match(r"q=(\d+(?:\.\d+)?)\s+y=(\d+(?:\.\d)*)\+(\d+(?:\.\d)*)/((\d+(?:\.\d)*)|\?)\s+(.+)", sample_input.sequencing_result_status)

    if not m:
        raise Exception('Error parsing', sample_input.sequencing_result_status)

    q30 = float(m.group(1))
    yield_cur = int(m.group(2))
    yield_total = int(m.group(3))
    yield_required = int(m.group(4)) if m.group(4) != '?' else None
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


def print_output_for_excel(samples: List[Sample]):
    results = []
    by_seq = []
    err = []
    for s in samples:
        if s.yield_requirement is None:
            result = None
            err.append(compute_result(
                s,
                'Unknown yield req in API',
                'AB'
            ))
        elif s.by_seq_in_gb() > 0:
            result = compute_result(
                s,
                f'Extra Seq {s.by_seq_in_gb()} GBase',
                'AB'
            )
            by_seq_result = result + (str(s.yield_total), str(s.yield_requirement))
            by_seq.append(by_seq_result)
        else:
            if s.sample_type == 'TUMOR' and s.output_type in ['Rna', 'Targeted']:
                result = compute_result(
                    s,
                    'Processing',
                    'BFX'
                )
            elif s.sample_type == 'TUMOR' and s.output_type in ['Somatic', 'ShallowSeq']:
                ref = find_ref_or_none(s, samples)
                if ref and ref.by_seq_in_gb() > 0:
                    result = compute_result(
                        s,
                        'Waiting on R',
                        'AB'
                    )
                else:
                    result = compute_result(
                        s,
                        'Processing',
                        'BFX'
                    )
            elif s.output_type == 'FastQ':
                result = compute_result(
                    s,
                    'Done',
                    ''
                )
            else:
                # samples that don't belong in sample-overview are skipped (i.e., ref samples)
                continue
        if result is not None:
            results.append(result)

    results.sort(key=lambda x: x[2])  # sort by sample name
    results = ['\t'.join(r) for r in results]

    print(*results, sep='\n')
    print("\n\n\n")
    print("The following samples need by-seq:", *['\t'.join(bs) for bs in by_seq], sep="\n")

    if err:
        print("\n\n")
        print("The following samples had errors:", *['\t'.join(e) for e in err], sep="\n")


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
