#!/usr/bin/env python3

import pandas as pd
import sys
import numpy as np


def main():
    if len(sys.argv) != 2:
        print("Error (incorrect usage), exiting ...")
        sys.exit(1)

    if sys.argv[1] == "Ampliseq":
        bedfile = pd.read_csv("/Users/mbaksi/scratch/CancerHotSpot-v2.dna_manifest.20180509.bed", sep="\t", header=None)
        bedfile.columns = ["Chrom", "Start", "End", "Target"]
        variants = pd.read_csv("/Users/mbaksi/scratch/temp.tsv", sep="\t")
        new_col = []
        for index, row in bedfile.iterrows():
            new_col.append(row["Chrom"][3:])
        bedfile["Chrom"] = new_col
        analysis(bedfile, variants)

    elif sys.argv[1] == "TSO500":
        bedfile = pd.read_csv("/Users/mbaksi/scratch/ActionableCodingPanel.tso500.37.bed", sep="\t", header=None)
        bedfile.columns = ["Chrom", "Start", "End", "Target"]
        variants = pd.read_csv("/Users/mbaksi/scratch/temp.tsv", sep="\t")
        analysis(bedfile, variants)

    elif sys.argv[1] == "Protect":
        df = pd.read_csv("/home/mbaksi/scratch/protect.csv", encoding='cp1252', sep="\t")
        total_dict = {}
        for index, row in df.iterrows():
            sampleId = row["sampleId"]
            onlabel = row['onLabel']
            event = row['event']
            gene = row['gene']
            identifier = sampleId + "_" + str(onlabel) + "_" + str(event) + "_" + str(gene)
            updated = False
            if identifier in total_dict.keys():
                sample_dict = total_dict[identifier]
            else:
                sample_dict = {}
            if 'event' in sample_dict.keys():
                if sample_dict['onLabel'] == row['onLabel']:
                    if sample_dict['event'] == row['event']:
                        if sample_dict['gene'] == row['gene']:
                            if sample_dict['level'] > row['level']:
                                updated = True
                                sample_dict['sampleId'] = row["sampleId"]
                                sample_dict['gene'] = row['gene']
                                sample_dict['event'] = row['event']
                                sample_dict['sourceEvent'] = row['sourceEvent']
                                sample_dict['onLabel'] = row['onLabel']
                                sample_dict['level'] = row['level']
            else:
                sample_dict['sampleId'] = row["sampleId"]
                sample_dict['gene'] = row['gene']
                sample_dict['event'] = row['event']
                sample_dict['sourceEvent'] = row['sourceEvent']
                sample_dict['onLabel'] = row['onLabel']
                sample_dict['level'] = row['level']

            total_dict[identifier] = sample_dict
        df2 = pd.DataFrame(total_dict)
        df2 = df2.T
        col_dict = {}
        for index, row in df2.iterrows():
            dictkey = str(row["sampleId"]) + "_" + str(row["sampleId"]) + "_" + str(row['gene']) + "_" + str(row['event'])
            if dictkey in col_dict.keys():
                sample_dict = col_dict[dictkey]
            else:
                sample_dict = {'event': row['event'], 'sourceEvent': row['sourceEvent'], 'gene': row['gene'], 'sampleId': row["sampleId"]}
            if row['onLabel'] == "1":
                sample_dict["offLabel"] = row['level']
            elif row['onLabel'] == "0":
                sample_dict["onLabel"] = row['level']
            col_dict[dictkey] = sample_dict
        df3 = pd.DataFrame(total_dict)
        df3 = df3.T
        df3.to_csv("/home/mbaksi/scratch/" + sys.argv[1] + "_result.tsv", sep="\t", index=False, header=False)

    elif sys.argv[1] == "MolDx":
        df = pd.read_csv("/Users/mbaksi/scratch/MolDx vdHoeven-Roepman list - SOC assays & tumor types.csv", sep=",")
        moldx(df)
        pass

    else:
        print("Error (incorrect usage), exiting ...")
        sys.exit(1)


def analysis(bedfile, variants):
    # Make a map
    mapping = {}
    for index, row in bedfile.iterrows():
        if row["Chrom"] in mapping.keys():
            value = mapping[row["Chrom"]]
            couple = (row["Start"], row["End"], row["Target"])
            value.append(couple)
            mapping[row["Chrom"]] = value
        else:
            value = []
            couple = (row["Start"], row["End"], row["Target"])
            value.append(couple)
            mapping[row["Chrom"]] = value

    # Get result
    results = {}
    for index, row in variants.iterrows():
        if row["chromosome"] in mapping.keys():
            allowed_range = mapping[row["chromosome"]]
            for pos in allowed_range:
                if pos[0] <= row["position"] <= pos[1]:
                    results[row["id"]] = pos[2]

    # Parse results
    variants_results = []
    test_results = []
    for index, row in variants.iterrows():
        variants_results.append(row["id"])
        if row["id"] in results.keys():
            test_results.append(results[row["id"]])
        else:
            test_results.append("Not found")
    df_results = pd.DataFrame()
    df_results["id"] = variants_results
    df_results[sys.argv[1]] = test_results
    df_results.to_csv("/Users/mbaksi/scratch/" + sys.argv[1] + "_result.tsv", sep="\t", index=False, header=False)


def moldx(df):
    doid_row = df.iloc[0].tolist()
    tumortype_row = df.iloc[1].tolist()
    for i in range(len(doid_row)):
        if i + 1 < len(doid_row):
            if str(doid_row[i + 1]).lower() == "nan":
                doid_row[i + 1] = doid_row[i]
    for i in range(len(tumortype_row)):
        if i + 1 < len(tumortype_row):
            if str(tumortype_row[i + 1]).lower() == "nan":
                tumortype_row[i + 1] = tumortype_row[i]

    total = [["Gene", "Event type", "DOID", "Tumor type", "Routine", "ESCAT"]]
    for index, row in df.iterrows():
        event = row[0]
        for i in range(4, len(doid_row) - 3):
            new = []
            if len(str(event).split()) > 1:
                new.append(str(event).split()[0])
                new.append(" ".join(str(event).split()[1:]))
                new.append(doid_row[i])
                new.append(tumortype_row[i])
                if row[i] in ("L", "R"):
                    new.append(row[i])
                    new.append(np.nan)
                else:
                    new.append(np.nan)
                    new.append(row[i])
                total.append(new)
    df_flattened = pd.DataFrame(total)
    new_header = df_flattened .iloc[0]
    df_flattened = df_flattened[1:]
    df_flattened.columns = new_header
    df_flattened.to_csv("/Users/mbaksi/scratch/" + sys.argv[1] + "_result.tsv", sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
