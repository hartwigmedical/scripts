#!/usr/bin/env python3
import argparse, sys, re

def parse_gcs_path(uri: str):
    uri = uri.strip()
    if not uri or uri.startswith("#"):
        return None, None
    if uri.endswith(".crai"):
        return None, None
    if not uri.startswith("gs://"):
        raise ValueError(f"Not a GCS URI: {uri}")
    folder, _, fname = uri.rpartition("/")
    if not folder or not fname:
        raise ValueError(f"Malformed GCS URI: {uri}")
    if not fname.endswith(".cram"):
        raise ValueError(f"Not a .cram file: {uri}")
    sample_id = fname[:-5]  # strip .cram
    return folder + "/", sample_id

def extract_key_part(sample_id: str, max_len: int = 20) -> str:
    """Try to pull meaningful part of ID, fallback to shortening"""
    # look for HD_x or similar pattern
    m = re.search(r"(HD_\d+)", sample_id)
    if m:
        return m.group(1)
    # else: fallback to shortened version
    if len(sample_id) <= max_len:
        return sample_id
    half = (max_len - 3) // 2
    return f"{sample_id[:half]}...{sample_id[-half:]}"

def main():
    ap = argparse.ArgumentParser(description="Convert CRAM GCS URIs to workflow executions YAML.")
    ap.add_argument("input_tsv", help="TSV/line file with one gs://.../SAMPLE.cram per line")
    ap.add_argument("-o", "--output", default="-", help="Output YAML file (default: stdout)")
    ap.add_argument("--workflow", default="snp-check-matching", help="Workflow name")
    ap.add_argument("--version", default="0.0.1-beta1", help="Workflow version")
    ap.add_argument("--max-name-len", type=int, default=20, help="Maximum length for fallback shortening")
    args = ap.parse_args()

    with open(args.input_tsv, "r", encoding="utf-8") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    docs = []
    for ln in lines:
        try:
            folder, sample = parse_gcs_path(ln)
        except ValueError as e:
            print(f"Warning: {e}", file=sys.stderr)
            continue
        if not sample:
            continue

        short_name = extract_key_part(sample, args.max_name_len)
        doc = (
            f"name: {short_name}\n"
            f"workflow: {args.workflow}\n"
            f"version: {args.version}\n"
            f"params:\n"
            f"  cram_folder_uri: {folder}\n"
            f"  sample_id: {sample}\n"
        )
        docs.append(doc)

    out_text = "---\n".join(docs) + ("\n" if docs else "")
    if args.output == "-" or args.output.lower() == "stdout":
        sys.stdout.write(out_text)
    else:
        with open(args.output, "w", encoding="utf-8") as out:
            out.write(out_text)

if __name__ == "__main__":
    main()
