#!/usr/bin/env python3

import argparse
import numpy as np

def read_adapters_list_from_file(input_file):
    with open(input_file) as f:
        adapters_list = [x for x in f.read().splitlines() if x.strip()]
    return adapters_list

def create_distance_matrix(adapters_list):
    distance_matrix = []
    for adapter_1 in adapters_list:
        distance_row = []
        for adapter_2 in adapters_list:
            distance_row.append((np.array(list(adapter_1)) != np.array(list(adapter_2))).sum())
        distance_matrix.append(distance_row)
    return distance_matrix

def find_conflicting_adapters_in_distance_matrix(max_distance, labels, distances):
    conflicts = []
    for x in range(len(labels)-1):
        for y in range(x+1, len(labels)):
            distance = distances[x][y]
            if distance <= max_distance:
                conflicts.append([f"{x}-{y}", labels[x], labels[y], str(distance)])
    return conflicts

def write_info_to_file(conflicts, output_file):
    with open(output_file, "w") as f:
        f.write("Pos\tAdapter1\tAdapter2\tDistance\n")
        for line in conflicts:
            f.write('\t'.join(line) + "\n")

def parse_args():
    parser = argparse.ArgumentParser(prog="check_adapters_conflict.py", description="Check if list of adapters has conflicts with each other")
    parser.add_argument("--input_file", "-i", type=str, required=True, help="path to input file, containing list of adapters (each line one adapter sequence)")
    parser.add_argument("--output_file", "-o", type=str, required=True, help="path to output file")
    parser.add_argument("---distance", "-d", type=int, required=False, default=2, help="minimum allowable distance between adapters (default 2)")
    args = parser.parse_args()
    return args

def main():
    config = parse_args()
    adapters_list = read_adapters_list_from_file(config.input_file)
    distance_matrix = create_distance_matrix(adapters_list)
    conflicts = find_conflicting_adapters_in_distance_matrix(config.distance, adapters_list, distance_matrix)
    write_info_to_file(conflicts, config.output_file)

if __name__ == "__main__":
    main()