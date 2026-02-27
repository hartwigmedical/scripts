#!/usr/bin/env python3
"""
Generate YAML execution file for the hotspot-coverage workflow.

This script scans a GCP bucket for sample directories and creates a single YAML
file with all samples separated by ---.

Usage:
    python generate_execution_hotspot-coverage.py <bucket_name> [options]

Examples:
    python generate_execution_hotspot-coverage.py my-input-bucket
    python generate_execution_hotspot-coverage.py my-input-bucket --tag volta -o hotspot-coverage.yaml
    python generate_execution_hotspot-coverage.py my-input-bucket --hotspot-file custom_hotspots.vcf.gz
"""

import argparse
import os
import sys
from pathlib import Path

try:
    from google.cloud import storage
except ImportError:
    print("Error: google-cloud-storage is required. Install with:")
    print("  pip install google-cloud-storage")
    sys.exit(1)

import yaml


def get_samples_from_bucket(bucket_name: str, prefix: str = "") -> list[dict]:
    """
    Scan a GCP bucket and extract sample information.

    Assumes bucket structure like:
        gs://bucket/SAMPLE123/SAMPLE123/aligner/SAMPLE123.bam

    Args:
        bucket_name: Name of the GCP bucket
        prefix: Optional prefix to filter objects

    Returns:
        List of dicts containing sample_id and alignment info
    """
    client = storage.Client()
    bucket = client.bucket(bucket_name)

    samples = {}

    blobs = bucket.list_blobs(prefix=prefix)

    for blob in blobs:
        if blob.name.endswith('/'):
            continue

        if blob.name.endswith('.bam') or blob.name.endswith('.cram'):
            parts = blob.name.split('/')

            if len(parts) >= 1:
                sample_id = parts[0].split('-')[0]

                alignment_file = os.path.basename(blob.name)
                alignment_dir = os.path.dirname(blob.name)

                if sample_id not in samples:
                    samples[sample_id] = {
                        'sample_id': sample_id,
                        'alignment_uri': f"gs://{bucket_name}/{alignment_dir}",
                        'alignment_file': alignment_file,
                        'full_path': blob.name
                    }

    return list(samples.values())


def create_yaml_content(
    sample: dict,
    tag: str = "default",
    workflow: str = "hotspot-coverage",
    version: str = "0.0.1",
    hotspot_uri: str = "gs://hmf-crunch-innovation/info/hawe-resources/",
    hotspot_file: str = "KnownHotspots.somatic.38.vcf.gz"
) -> dict:
    """
    Create YAML content for a sample.

    Args:
        sample: Sample info dict with sample_id, alignment_uri, alignment_file
        tag: Tag to append to the name
        workflow: Workflow name
        version: Workflow version
        hotspot_uri: GCS URI to the hotspot VCF directory
        hotspot_file: Hotspot VCF filename

    Returns:
        Dict representing the YAML structure
    """
    return {
        'name': f"{sample['sample_id']}-{tag}",
        'workflow': workflow,
        'version': version,
        'params': {
            'alignment_uri': sample['alignment_uri'],
            'sample_id': sample['sample_id'],
            'alignment_file': sample['alignment_file'],
            'hotspot_uri': hotspot_uri,
            'hotspot_file': hotspot_file
        }
    }


def main():
    parser = argparse.ArgumentParser(
        description='Generate YAML execution file for hotspot-coverage workflow',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s my-bucket
  %(prog)s my-bucket --tag volta
  %(prog)s my-bucket -o hotspot-coverage.yaml --hotspot-file custom.vcf.gz
  %(prog)s my-bucket --prefix project1/ --dry-run
        """
    )

    parser.add_argument(
        'bucket',
        help='GCP bucket name (without gs:// prefix)'
    )

    parser.add_argument(
        '--prefix',
        default='',
        help='Optional prefix to filter objects in the bucket'
    )

    parser.add_argument(
        '--tag',
        default='default',
        help='Tag to include in the YAML name field (default: "default")'
    )

    parser.add_argument(
        '--workflow',
        default='hotspot-coverage',
        help='Workflow name (default: "hotspot-coverage")'
    )

    parser.add_argument(
        '--version',
        default='0.0.1',
        help='Workflow version (default: "0.0.1")'
    )

    parser.add_argument(
        '--hotspot-uri',
        default='gs://hmf-crunch-innovation/info/hawe-resources/',
        help='GCS URI to the hotspot VCF directory (default: "gs://hmf-crunch-innovation/info/hawe-resources/")'
    )

    parser.add_argument(
        '--hotspot-file',
        default='KnownHotspots.somatic.38.vcf.gz',
        help='Hotspot VCF filename (default: "KnownHotspots.somatic.38.vcf.gz")'
    )

    parser.add_argument(
        '--output',
        '-o',
        default='hotspot-coverage-execution.yaml',
        help='Output YAML file (default: hotspot-coverage-execution.yaml)'
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print what would be created without writing files'
    )

    parser.add_argument(
        '--verbose',
        '-v',
        action='store_true',
        help='Print verbose output'
    )

    args = parser.parse_args()

    bucket_name = args.bucket.replace('gs://', '').rstrip('/')

    print(f"Scanning bucket: gs://{bucket_name}")
    if args.prefix:
        print(f"Using prefix: {args.prefix}")

    try:
        samples = get_samples_from_bucket(bucket_name, args.prefix)
    except Exception as e:
        print(f"Error accessing bucket: {e}")
        sys.exit(1)

    if not samples:
        print("No samples found in bucket.")
        sys.exit(0)

    print(f"Found {len(samples)} sample(s)")

    yaml_contents = []

    for sample in sorted(samples, key=lambda x: x['sample_id']):
        content = create_yaml_content(
            sample=sample,
            tag=args.tag,
            workflow=args.workflow,
            version=args.version,
            hotspot_uri=args.hotspot_uri,
            hotspot_file=args.hotspot_file
        )
        yaml_contents.append(content)

        if args.verbose:
            print(f"  - {sample['sample_id']}: {sample['alignment_file']}")

    output_parts = []
    for content in yaml_contents:
        output_parts.append(yaml.dump(content, default_flow_style=False, sort_keys=False))

    output_str = '---\n'.join(output_parts)

    if args.dry_run:
        print("\n--- Dry run output: ---\n")
        print(output_str)
        return

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        f.write(output_str)

    print(f"Written {len(yaml_contents)} sample(s) to: {output_path}")


if __name__ == '__main__':
    main()
