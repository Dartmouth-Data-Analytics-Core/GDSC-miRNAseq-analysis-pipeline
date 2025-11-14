#!/usr/bin/env python3

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
