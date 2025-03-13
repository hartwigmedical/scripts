import os
import glob
import gzip
import pandas as pd
import argparse
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows

# Functions to count truth entries per category

def count_germline_deletion(truth_dir):
    pattern = os.path.join(truth_dir, "*", "purple", "*.purple.germline.deletion.tsv")
    files = glob.glob(pattern)
    total = 0
    for f in files:
        try:
            df = pd.read_csv(f, sep="\t")
            # Count rows where the last column is "true" (case insensitive)
            count = df[df.iloc[:, -1].astype(str).str.lower() == "true"].shape[0]
            total += count
        except Exception as e:
            print(f"Error processing file {f}: {e}")
    return total

def count_somatic_variant(truth_dir):
    pattern = os.path.join(truth_dir, "*", "purple", "*.purple.somatic.vcf.gz")
    files = glob.glob(pattern)
    total = 0
    for f in files:
        try:
            with gzip.open(f, 'rt') as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    if "REPORTED" in line:
                        total += 1
        except Exception as e:
            print(f"Error processing file {f}: {e}")
    return total

def count_germline_sv(truth_dir):
    pattern = os.path.join(truth_dir, "*", "linx_germline", "*.linx.germline.driver.catalog.tsv")
    files = glob.glob(pattern)
    total = 0
    for f in files:
        try:
            df = pd.read_csv(f, sep="\t")
            total += df.shape[0]
        except Exception as e:
            print(f"Error processing file {f}: {e}")
    return total

def count_disruption(truth_dir):
    pattern = os.path.join(truth_dir, "*", "linx", "*.linx.driver.catalog.tsv")
    files = glob.glob(pattern)
    total = 0
    for f in files:
        try:
            df = pd.read_csv(f, sep="\t")
            if "driver" in df.columns:
                count = df[df["driver"].astype(str).str.contains("DISRUPTION", na=False)].shape[0]
            else:
                count = 0
            total += count
        except Exception as e:
            print(f"Error processing file {f}: {e}")
    return total

def count_germline_variant(truth_dir):
    pattern = os.path.join(truth_dir, "*", "purple", "purple", "*.purple.germline.vcf.gz")
    files = glob.glob(pattern)
    total = 0
    for f in files:
        try:
            with gzip.open(f, 'rt') as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    if "REPORTED" in line:
                        total += 1
        except Exception as e:
            print(f"Error processing file {f}: {e}")
    return total

def count_fusion(truth_dir):
    pattern = os.path.join(truth_dir, "*", "linx", "*.linx.fusion.tsv")
    files = glob.glob(pattern)
    total = 0
    for f in files:
        try:
            df = pd.read_csv(f, sep="\t")
            # Count rows where the fourth column (index 3) equals "true"
            count = df[df.iloc[:, 3].astype(str).str.lower() == "true"].shape[0]
            total += count
        except Exception as e:
            print(f"Error processing file {f}: {e}")
    return total

def count_driver(truth_dir):
    pattern = os.path.join(truth_dir, "*", "purple", "*.purple.driver.catalog.somatic.tsv")
    files = glob.glob(pattern)
    total = 0
    for f in files:
        try:
            df = pd.read_csv(f, sep="\t")
            if df.shape[1] >= 6:
                # Count only rows where the 6th column (index 5) is one of DEL, AMP, PARTIAL_AMP
                count = df[df.iloc[:, 5].isin(["DEL", "AMP", "PARTIAL_AMP"])].shape[0]
            else:
                count = 0
            total += count
        except Exception as e:
            print(f"Error processing file {f}: {e}")
    return total

# Mapping from category name (as parsed from the compar file) to the corresponding function.
truth_count_functions = {
    "germline_deletion": count_germline_deletion,
    "somatic_variant": count_somatic_variant,
    "germline_sv": count_germline_sv,
    "disruption": count_disruption,
    "germline_variant": count_germline_variant,
    "fusion": count_fusion,
    "driver": count_driver,
}

def categorize_files_to_excel(input_dir, output_file, truth_dir):
    # Define colors for tabs
    COLOR_STANDARD = "ffc000"  # Yellowish
    COLOR_INVALID = "ff0000"   # Red
    COLOR_EMPTY = "92d050"     # Green
    
    # Initialize workbook and remove default sheet
    workbook = Workbook()
    if 'Sheet' in workbook.sheetnames:
        workbook.remove(workbook['Sheet'])
    
    summary_data = []
    
    # Process compar TSV files from the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith(".tsv"):
            file_path = os.path.join(input_dir, file_name)
            # Here we assume filenames are like "sample.compar.<category>.tsv"
            try:
                parts = file_name.split('.')
                # Extract category as the third token (index 2)
                tab_name = parts[2]
                df = pd.read_csv(file_path, sep='\t')
                
                # Determine the tab color based on file content
                if df.empty or (len(df.columns) == 1 and all(df.columns == df.iloc[0])):
                    color = COLOR_EMPTY
                elif "INVALID" in df.to_string(index=False):
                    color = COLOR_INVALID
                else:
                    color = COLOR_STANDARD
                
                # Count occurrences in the "MismatchType" column if present
                if "MismatchType" in df.columns:
                    new_only_count = (df["MismatchType"] == "NEW_ONLY").sum()
                    ref_only_count = (df["MismatchType"] == "REF_ONLY").sum()
                    value_count = (df["MismatchType"] == "VALUE").sum()
                else:
                    new_only_count = ref_only_count = value_count = 0
                
                summary_data.append([tab_name, new_only_count, ref_only_count, value_count])
                
                # Create worksheet and set its tab color
                sheet = workbook.create_sheet(title=tab_name[:31])
                for row in dataframe_to_rows(df, index=False, header=True):
                    sheet.append(row)
                sheet.sheet_properties.tabColor = color
            except Exception as e:
                print(f"Error processing file '{file_name}': {e}")
    
    # Build the summary DataFrame with collected counts
    summary_df = pd.DataFrame(summary_data, columns=["Category", "#NEW_ONLY", "#REF_ONLY", "#VALUE"])
    
    # Compute CONCORDANCE for applicable categories.
    # The concordance is calculated as:
    #   (pipeline_truth_count - REF_ONLY) / (pipeline_truth_count + NEW_ONLY)
    concordance_values = []
    for index, row in summary_df.iterrows():
        category = row["Category"]
        new_only = row["#NEW_ONLY"]
        ref_only = row["#REF_ONLY"]
        if category in truth_count_functions:
            pipeline_count = truth_count_functions[category](truth_dir)
            if (pipeline_count + new_only) > 0:
                concordance = (pipeline_count - ref_only) / (pipeline_count + new_only)
            else:
                concordance = None
        else:
            concordance = None
        concordance_values.append(concordance)
    
    summary_df["CONCORDANCE"] = concordance_values
    
    # Create summary worksheet at the start and add the DataFrame rows
    summary_sheet = workbook.create_sheet(title="Summary", index=0)
    for row in dataframe_to_rows(summary_df, index=False, header=True):
        summary_sheet.append(row)
    
    # Save the workbook to the output file
    workbook.save(output_file)
    print(f"Excel file '{output_file}' created successfully with summary tab.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Organize TSV files into Excel worksheet tabs with specific colors, add a summary tab, and compute concordance for specific categories."
    )
    parser.add_argument(
        "input_dir",
        type=str,
        help="The directory containing the compar TSV files to process."
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="The name of the output Excel file."
    )
    parser.add_argument(
        "truth_dir",
        type=str,
        help="The directory containing the truth pipeline output files."
    )

    args = parser.parse_args()
    categorize_files_to_excel(args.input_dir, args.output_file, args.truth_dir)