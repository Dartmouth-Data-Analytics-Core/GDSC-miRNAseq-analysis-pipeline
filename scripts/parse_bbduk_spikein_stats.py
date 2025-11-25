#!/usr/bin/env python3

################################################################################
# Script Name: parse_bbduk_spikein_stats.py
#
# Description:
#   This script parses BBDuk spike-in statistics files produced by the spikein quantification.
#   It extracts spike-in read counts (miND-01 to miND-07), combines them across multiple samples, 
#   and outputs a tab-separated table with spike-ins as rows and samples as columns.
#
# Usage:
#   python parse_bbduk_spikein_stats.py --stats <files> --samples "<sample1 sample2 ...>" --output <file>
#
# Required Arguments:
#   --stats           One or more BBDuk stats files (in miND spike-in format)
#   --samples         Space-separated sample names, in the same order as --stats
#   --output          Path to output TSV file
#
# Example:
#   python parse_bbduk_spikein_stats.py \
#       --stats spikein_alignment/{sample1}.stats spikein_alignment/{sample2}.stats spikein_alignment/{sample3}.stats \
#       --samples "S1 S2 S3" \
#       --output spikein_counts.tsv
#
# Input Format (BBDuk .stats file):
#   #File
#   #Total
#   #Matched
#   #Name      Reads     ReadsPct
#   #miND-01   9093      0.25671%
#   #miND-02   5000      0.14281%
#   ...
#
# Output:
#   - <output>.tsv: Table with spike-ins (miND-01..miND-07) as rows and sample
#                   read counts as columns. Missing counts are filled with zero.
#
# Author:
#   Beatriz Bergamo
################################################################################

import sys
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description='Parse BBDuk spike-in stats')
    parser.add_argument('--stats', nargs='+', required=True, help='BBDuk stats files')
    parser.add_argument('--samples', type=str, required=True, help='Space-separated sample names')
    parser.add_argument('--output', type=str, required=True, help='Output file')
    return parser.parse_args()

def parse_bbduk_stats(stats_file):
    """
    Parse BBDuk stats file with miND spike-in format:
    #File
    #Total
    #Matched
    #Name   Reads   ReadsPct
    #miND-01   9093   0.25671%
    """
    data = {}

    with open(stats_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Skip only metadata / header lines
            if line.startswith("#File") or line.startswith("#Total") or line.startswith("#Matched") or line.startswith("#Name"):
                continue

            parts = line.split("\t")
            if len(parts) < 2:
                continue

            spike_id = parts[0]
            reads_str = parts[1]

            # Keep only miND spike-ins
            if spike_id.startswith("#miND-"):
                spike_id = spike_id.lstrip("#")  # remove leading #
                try:
                    reads = int(reads_str)
                except ValueError:
                    continue

                data[spike_id] = reads

    return data


def main():
    args = parse_args()
    samples = args.samples.split()

    spike_data = {}

    # Parse each stats file paired with corresponding sample
    for stats_file, sample in zip(args.stats, samples):
        print(f"Parsing {stats_file} ({sample})...")
        stats = parse_bbduk_stats(stats_file)

        for spike_id, read_count in stats.items():
            if spike_id not in spike_data:
                spike_data[spike_id] = {}
            spike_data[spike_id][sample] = read_count

    # Convert to DataFrame
    df = pd.DataFrame.from_dict(spike_data, orient='index')
    df.index.name = "spikein_ID"
    df = df.fillna(0).astype(int)

    df.to_csv(args.output, sep="\t")
    print(f"\n✅ Spike-in counts saved to: {args.output}")
    print(f"Shape: {df.shape}")
    print(df)


if __name__ == '__main__':
    main()
