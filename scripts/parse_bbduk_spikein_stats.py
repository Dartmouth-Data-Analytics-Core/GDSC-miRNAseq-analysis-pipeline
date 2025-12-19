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
#       --output spikein.readcounts.tsv
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
# Author:
#   Beatriz Bergamo, Owen Wilkins
################################################################################

import sys
import argparse
import re
import pandas as pd

# for debugging 
#sys.argv = ['parse_bbduk_spikein_stats.py', 
#            '--stats', 'spikein_alignment/acaro_5BS54_left_90_spk.stats spikein_alignment/acaro_5BS54_left_90.stats spikein_alignment/nhbcs_E33_spk.stats spikein_alignment/nhbcs_E33.stats spikein_alignment/nhbcs_E4_spk.stats spikein_alignment/nhbcs_E4.stats',
#            '--samples', 'acaro_5BS54_left_90_spk acaro_5BS54_left_90 nhbcs_E33_spk nhbcs_E33 nhbcs_E4_spk nhbcs_E4',
#            '--output', 'test_output.tsv']

# Expected spike-ins
EXPECTED_SPIKEINS = [f"miND-{i:02d}" for i in range(1, 8)]

# define command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Parse BBDuk spike-in stats')
    parser.add_argument('--stats', type=str, required=True, help='BBDuk stats files')
    # Keep samples as a single space-separated string for backward compatibility
    parser.add_argument('--samples', type=str, required=True, help='Space-separated sample names (in the same order as --stats)')
    parser.add_argument('--output', type=str, required=True, help='Output file')
    return parser.parse_args()

### add documentation 

def parse_bbduk_stats(stats_file, verbose=False):
    """
    Script Name: parse_bbduk_stats
    Description:
        Parse a single BBDuk spike-in stats file and extract miND spike-in read counts.
        It returns a dictionary mapping spike-in identifiers (e.g., 'miND-01') to integer
        read counts. Header/metadata lines are ignored and malformed rows are skipped.
    Usage:
        parse_bbduk_stats(stats_file, verbose=False)
    Required Arguments:
        stats_file (str): Path to a BBDuk .stats file.
    Optional Arguments:
        verbose (bool): If True, print warnings for skipped or malformed lines to stderr.
    Example:
        >>> parse_bbduk_stats('spikein_alignment/S1.stats', verbose=True)
        {'miND-01': 9093, 'miND-02': 5000}
    Output:
        dict: Mapping from spike-in ID to integer read count. Example:
              {'miND-01': 9093, 'miND-02': 5000}
    """
    data = {}
    # read file line by line and extract read counts for each spike-in
    with open(stats_file, 'r', encoding='utf-8') as f:
        for lineno, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            # Skip only metadata / header lines
            if line.startswith("#File") or line.startswith("#Total") or line.startswith("#Matched") or line.startswith("#Name"):
                continue
            # Split line on any whitespace (handles tab or space separated files)
            parts = re.split(r"\s+", line)
            if len(parts) < 2:
                if verbose:
                    print(f"Skipping malformed line {lineno} in {stats_file}: less than 2 columns: {line}", file=sys.stderr)
                continue
            spike_id, reads_str = parts[0], parts[1]
            # extract reads, keeping only miND spike-ins
            if spike_id.startswith("#miND-"):
                spike_id = spike_id.lstrip("#")
                try:
                    reads = int(reads_str)
                except ValueError:
                    if verbose:
                        print(f"Skipping line {lineno} in {stats_file} due to invalid reads value '{reads_str}': {line}", file=sys.stderr)
                    continue
                data[spike_id] = reads
    # Ensure all expected spike-ins are present; fill missing with zero
    for sp in EXPECTED_SPIKEINS:
        data.setdefault(sp, 0)
    return data

def main():
    args = parse_args()
    samples = args.samples.split()
    stats_files = args.stats.split()
    # Validate that the provided stats files and sample names map 1:1
    if len(stats_files) != len(samples):
        sys.exit('Error: The number of --stats files must match the number of --samples names')
    spike_data = {}
    # Parse each stats file paired with corresponding sample
    for stats_file, sample in zip(stats_files, samples):
        print(f"Parsing {stats_file} ({sample})...")
        stats = parse_bbduk_stats(stats_file, verbose=True)
        # extract read counts for each spike-in and add to spike_data as nested dicts 
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
