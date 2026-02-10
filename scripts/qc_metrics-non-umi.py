#!/usr/bin/env python3

# Script:   qc_metrics-non-umi.py
# Purpose:  Parse miRBase alignment logs, and genome alignment logs
#           to generate a combined QC metrics table.
#
# Usage (positional arguments):
#   python qc_metrics-non-umi.py <mirbase_alignments_dir> <genome_alignment_dir>
#
# Output:
#   metrics/mirna_genome_alignment_metrics.tsv
#   metrics/mirna_genome_alignment_metrics.xlsx

import sys
from glob import glob
import pandas as pd
import subprocess
from pathlib import Path

def run_samtools_count(bam, include_flag=None, exclude_flag=None):
    """Run samtools view from Python and return count as int."""
    cmd = ["samtools", "view", bam, "-c"]
    if include_flag is not None:
        cmd.extend(["-f", str(include_flag)])
    if exclude_flag is not None:
        cmd.extend(["-F", str(exclude_flag)])

    result = subprocess.run(cmd, text=True, check=True, capture_output=True)
    return int(result.stdout.strip())

def parse_bowtie1_log(logfile):
    """
    Parse Bowtie1 log file (e.g. *.mature.bowtie1.log) and return:
        total_reads, aligned_reads, unaligned_reads

    Expected Bowtie1 summary lines look like:
        # reads processed: 9263619
        # reads with at least one alignment: 14616 (0.16%)
        # reads that failed to align: 9249003 (99.84%)
        Reported 212080 alignments
        Time searching: 00:03:08
        Overall time: 00:03:08
    """
    total_reads = None
    aligned_reads = None
    unaligned_reads = None

    with open(logfile) as f:
        for line in f:
            line = line.strip()

            if line.startswith("# reads processed:"):
                # "# reads processed: 9263619"
                total_reads = int(line.split(":", 1)[1].strip().split()[0])
            elif line.startswith("# reads with at least one alignment:"):
                # "# reads with at least one alignment: 14616 (0.16%)"
                aligned_reads = int(line.split(":", 1)[1].strip().split()[0])
            elif line.startswith("# reads that failed to align:"):
                # "# reads that failed to align: 9249003 (99.84%)"
                unaligned_reads = int(line.split(":", 1)[1].strip().split()[0])

    if total_reads is None or aligned_reads is None or unaligned_reads is None:
        raise ValueError(f"Bowtie1 log file {logfile} missing expected summary fields.")

    return total_reads, aligned_reads, unaligned_reads


def parse_bowtie2_log(logfile):
    """
    Parse Bowtie2-style log file (e.g. *_genome.log.txt) and return:
        total_reads, aligned_reads, unaligned_reads, multimapped_reads

    Expected Bowtie2 summary lines look like:
        9247601 reads; of these:
          9247601 (100.00%) were unpaired; of these:
            8958899 (96.88%) aligned 0 times
            19851 (0.21%) aligned exactly 1 time
            268851 (2.91%) aligned >1 times
        3.12% overall alignment rate
    """
    total_reads = None
    aligned_0 = None
    aligned_1 = None
    aligned_gt1 = None

    with open(logfile) as f:
        for line in f:
            line = line.strip()

            if line.endswith("reads; of these:"):
                # "9247601 reads; of these:"
                total_reads = int(line.split()[0])
            elif "aligned 0 times" in line:
                aligned_0 = int(line.split()[0])
            elif "aligned exactly 1 time" in line:
                aligned_1 = int(line.split()[0])
            elif "aligned >1 times" in line:
                aligned_gt1 = int(line.split()[0])

    if None in (total_reads, aligned_0, aligned_1, aligned_gt1):
        raise ValueError(f"Bowtie2 log file {logfile} missing expected summary fields.")

    aligned_reads = aligned_1 + aligned_gt1
    unaligned_reads = aligned_0
    multimapped_reads = aligned_gt1

    return total_reads, aligned_reads, unaligned_reads, multimapped_reads

# ------------------------------
# Command-line args
# ------------------------------
mir_dir = Path(sys.argv[1])          # parent directory, e.g. mirbase_alignments
genome_dir = Path(sys.argv[2])

mir_mature_dir = mir_dir / "mature"
mir_hairpin_dir = mir_dir / "hairpin"

# ------------------------------
# Parse miRBase mapping logs (mature and hairpin)
# ------------------------------
mir_data = {}
sample_list = []

for logfile in sorted(glob(str(mir_mature_dir / "*.mature.bowtie1.log"))):
    sample_id = Path(logfile).name.replace(".mature.bowtie1.log", "")
    sample_list.append(sample_id)
    reads, mapped, unmapped = parse_bowtie1_log(logfile)
    mir_data[sample_id] = {
        "reads_mature": reads,
        "mapped_mature": mapped,
        "unaligned_mature": unmapped,
    }

hairpin_data = {}
for logfile in sorted(glob(str(mir_hairpin_dir / "*.hairpin.bowtie1.log"))):
    sample_id = Path(logfile).name.replace(".hairpin.bowtie1.log", "")
    reads, mapped, unmapped = parse_bowtie1_log(logfile)
    hairpin_data[sample_id] = {"reads": reads, "mapped": mapped, "unaligned": unmapped}

# ------------------------------
# Parse genome mapping logs
# ------------------------------
genome_data = {}
for logfile in sorted(glob(str(genome_dir / "*genome.log.txt"))):
    sample_id = Path(logfile).name.replace("_genome.log.txt", "")
    reads, mapped, unmapped, multimap = parse_bowtie2_log(logfile)
    genome_data[sample_id] = {"reads": reads, "mapped": mapped, "unaligned": unmapped, "multimap": multimap}

# ------------------------------
# Build metrics DataFrame
# ------------------------------
sample_list = sorted(sample_list)
metrics = pd.DataFrame(index=[
    "# of reads (mature)",
    "# of reads mapping to mirbase (mature)",
    "% of reads mapping to mirbase (mature)",
    "# of reads mapping to miRBase after filters (mature)",
    "% of reads mapping to miRBase after filters (mature)",
    "# of reads mapping to mirbase (hairpin)",
    "% of reads mapping to mirbase (hairpin)",
    "# of reads mapped genome",
    "% of reads mapped genome"
    ],
    columns=sample_list)

# ------------------------------
# miRBase counts (mature)
# ------------------------------
metrics.loc["# of reads (mature)"] = [mir_data[s]["reads_mature"] for s in sample_list]
metrics.loc["# of reads mapping to mirbase (mature)"] = [mir_data[s]["mapped_mature"] for s in sample_list]
metrics.loc["% of reads mapping to mirbase (mature)"] = (
    metrics.loc["# of reads mapping to mirbase (mature)"].astype(float)
    / metrics.loc["# of reads (mature)"].astype(float) * 100
).round(2)

filtered_bams = sorted(glob(str(mir_mature_dir / "*.mature.srt.bam")))
metrics.loc["# of reads mapping to miRBase after filters (mature)"] = [
    run_samtools_count(bam, exclude_flag=4) for bam in filtered_bams
]
metrics.loc["% of reads mapping to miRBase after filters (mature)"] = (
    metrics.loc["# of reads mapping to miRBase after filters (mature)"].astype(float)
    / metrics.loc["# of reads (mature)"].astype(float) * 100
).round(2)

# ------------------------------
# miRBase counts (hairpin)
# ------------------------------
metrics.loc["# of reads mapping to mirbase (hairpin)"] = [hairpin_data[s]["mapped"] for s in sample_list]
metrics.loc["% of reads mapping to mirbase (hairpin)"] = (
    metrics.loc["# of reads mapping to mirbase (hairpin)"].astype(float)
    / metrics.loc["# of reads (mature)"].astype(float) * 100
).round(2)

# ------------------------------
# Genome counts
# ------------------------------
genome_bams = sorted(glob(str(genome_dir / "*.srt.bam")))
metrics.loc["# of reads mapped genome"] = [
    run_samtools_count(bam, exclude_flag=4) for bam in genome_bams
]
metrics.loc["% of reads mapped genome"] = (
    metrics.loc["# of reads mapped genome"].astype(float)
    / metrics.loc["# of reads (mature)"].astype(float) * 100
).round(2)

# ------------------------------
# Output
# ------------------------------
metrics.to_csv("metrics/mirna_genome_alignment_metrics.tsv", sep="\t", index=True)
metrics.to_excel("metrics/mirna_genome_alignment_metrics.xlsx", index=True)