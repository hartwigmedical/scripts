#!/usr/bin/env python3
"""
Generate YAML workflow file for samples in a GCP bucket.

This script scans a GCP bucket for sample directories and creates a single YAML
file with all samples separated by ---.

Usage:
    python generate_sample_yamls.py <bucket_name> [options]

Examples:
    python generate_sample_yamls.py my-input-bucket
    python generate_sample_yamls.py my-input-bucket --tag production -o samples.yaml
    python generate_sample_yamls.py my-input-bucket --workflow my-workflow --version 1.0.0
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
    
    # List all blobs in the bucket
    blobs = bucket.list_blobs(prefix=prefix)
    
    for blob in blobs:
        # Skip directories (they end with /)
        if blob.name.endswith('/'):
            continue
            
        # Look for .bam files
        if blob.name.endswith('.bam'):
            parts = blob.name.split('/')
            
            # Extract sample ID from the path
            # Expected structure: SAMPLE123-xxx/SAMPLE123-xxx/aligner/SAMPLE123-xxx.bam
            # Split on '-' and take only the first part as sample_id
            if len(parts) >= 1:
                sample_id = parts[0].split('-')[0]
                
                # Get the alignment file name
                alignment_file = os.path.basename(blob.name)
                
                # Get the directory containing the alignment
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
    bucket_name: str,
    tag: str = "default",
    workflow: str = "collect-hs-metrics",
    version: str = "0.0.1"
) -> dict:
    """
    Create YAML content for a sample.
    
    Args:
        sample: Sample info dict with sample_id, alignment_uri, alignment_file
        bucket_name: Source bucket name
        tag: Tag to append to the name
        workflow: Workflow name
        version: Workflow version
        
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
            'alignment_file': sample['alignment_file']
        }
    }


def main():
    parser = argparse.ArgumentParser(
        description='Generate YAML workflow file for samples in a GCP bucket',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s my-bucket
  %(prog)s my-bucket --tag production
  %(prog)s my-bucket -o workflows.yaml --workflow my-workflow
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
        default='collect-hs-metrics',
        help='Workflow name (default: "collect-hs-metrics")'
    )
    
    parser.add_argument(
        '--version',
        default='0.0.1',
        help='Workflow version (default: "0.0.1")'
    )
    
    parser.add_argument(
        '--output',
        '-o',
        default='samples.yaml',
        help='Output YAML file (default: samples.yaml)'
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
    
    # Remove gs:// prefix if provided
    bucket_name = args.bucket.replace('gs://', '').rstrip('/')
    
    print(f"Scanning bucket: gs://{bucket_name}")
    if args.prefix:
        print(f"Using prefix: {args.prefix}")
    
    # Get samples from bucket
    try:
        samples = get_samples_from_bucket(bucket_name, args.prefix)
    except Exception as e:
        print(f"Error accessing bucket: {e}")
        sys.exit(1)
    
    if not samples:
        print("No samples found in bucket.")
        sys.exit(0)
    
    print(f"Found {len(samples)} sample(s)")
    
    # Generate YAML content for all samples
    yaml_contents = []
    
    for sample in sorted(samples, key=lambda x: x['sample_id']):
        content = create_yaml_content(
            sample=sample,
            bucket_name=bucket_name,
            tag=args.tag,
            workflow=args.workflow,
            version=args.version
        )
        yaml_contents.append(content)
        
        if args.verbose:
            print(f"  - {sample['sample_id']}: {sample['alignment_file']}")
    
    # Build the output string with --- separators
    output_parts = []
    for content in yaml_contents:
        output_parts.append(yaml.dump(content, default_flow_style=False, sort_keys=False))
    
    output_str = '---\n'.join(output_parts)
    
    if args.dry_run:
        print("\n--- Dry run output: ---\n")
        print(output_str)
        return
    
    # Write to file
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        f.write(output_str)
    
    print(f"Written {len(yaml_contents)} sample(s) to: {output_path}")


if __name__ == '__main__':
    main()