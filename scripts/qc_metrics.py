#!/usr/bin/env python3
################################################################################
# Script:   qc_metrics.py
# Purpose:  Parse UMI logs, miRBase alignment logs, and genome alignment logs
#           to generate a combined QC metrics table.
#
# Usage (positional arguments):
#   python qc_metrics.py <umi_reads_dir> <mirbase_alignments_dir> <genome_alignment_dir>
#
# Arguments:
#   umi_reads_dir             Directory containing UMI extraction log files.
#   mirbase_alignment_dir     Directory containing miRBase alignment log files.
#   genome_alignment_dir      Directory containing genome alignment log files.
#
# Example:
#   python qc_metrics.py umi_reads mirbase_alignment genome_alignment
#
# Output:
#   metrics/mirna_genome_alignment_metrics.tsv
#   metrics/mirna_genome_alignment_metrics.xlsx
#
# Author: Beatriz Bergamo
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

def parse_log_value(logfile, key):
    """Return the last integer found in a log file line containing key."""
    value = None
    with open(logfile) as f:
        for line in f:
            if key in line:
                value = int(line.strip().split()[-1])
    return value

def parse_bowtie1_log(logfile):
    """
    Parse Bowtie1 log file (e.g. *.hairpin.bowtie1.log) and return:
        total_reads, aligned_reads, unaligned_reads

    Expected Bowtie1 summary lines look like:
        # reads processed: 9249003
        # reads with at least one alignment: 1402 (0.02%)
        # reads that failed to align: 9247601 (99.98%)
        Reported 2849 alignments
        Time searching: 00:03:36
        Overall time: 00:03:36
    """
    total_reads = None
    aligned_reads = None
    unaligned_reads = None

    with open(logfile) as f:
        for line in f:
            line = line.strip()

            if line.startswith("# reads processed:"):
                total_reads = int(line.split(":", 1)[1].strip().split()[0])
            elif line.startswith("# reads with at least one alignment:"):
                aligned_reads = int(line.split(":", 1)[1].strip().split()[0])
            elif line.startswith("# reads that failed to align:"):
                unaligned_reads = int(line.split(":", 1)[1].strip().split()[0])

    if total_reads is None or aligned_reads is None or unaligned_reads is None:
        raise ValueError(f"Bowtie1 log file {logfile} missing expected summary fields.")

    return total_reads, aligned_reads, unaligned_reads


def parse_bowtie2_log(logfile):
    """
    Parse Bowtie2-style log file (e.g. *_mirbase.log.txt, *_genome.log.txt) and return:
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

# set variables for input directories from parser 
umi_dir = Path(sys.argv[1])
mir_dir = Path(sys.argv[2])      # parent directory, e.g. mirbase_alignments
genome_dir = Path(sys.argv[3])

mir_mature_dir = mir_dir / "mature"
mir_hairpin_dir = mir_dir / "hairpin"

# ------------------------------
# Parse UMI log files
# ------------------------------
umi_data = {}
sample_list = []
for logfile in sorted(glob(str(umi_dir / "*.log.txt"))):
    sample_id = Path(logfile).name.split(".")[0]
    sample_list.append(sample_id)
    reads_before = parse_log_value(logfile, "INFO Input Reads:")
    reads_after  = parse_log_value(logfile, "INFO Reads output:")
    reads_removed = reads_before - reads_after
    umi_data[sample_id] = {
        "reads_before": reads_before,
        "reads_removed": reads_removed,
        "reads_after": reads_after,
    }

# ------------------------------
# Parse miRBase mapping logs (mature and hairpin)
# ------------------------------
mir_data = {}
for logfile in sorted(glob(str(mir_mature_dir / "*mirbase.log.txt"))):
    sample_id = Path(logfile).name.replace("_mirbase.log.txt", "")
    _, mapped, _, multimap = parse_bowtie2_log(logfile)
    if sample_id not in mir_data:
        mir_data[sample_id] = {}
    mir_data[sample_id].update({"mapped_mature": mapped, "multimap_mature": multimap})

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
    _, mapped, _, multimap = parse_bowtie2_log(logfile)
    genome_data[sample_id] = {"mapped": mapped, "multimap": multimap}

# ------------------------------
# Build metrics DataFrame
# ------------------------------
sample_list = sorted(sample_list)
metrics = pd.DataFrame(index=[
    "# of reads",
    "# of UMI-containing reads",
    "% of UMI-containing reads",
    "# of reads mapping to mirbase (mature)",
    "% of reads mapping to mirbase (mature)",
    "# of reads multimapping mirbase (mature)",
    "% of reads multimapping mirbase (mature)",
    "# of reads mapping to miRBase after filters (mature)",
    "% of reads mapping to miRBase after filters (mature)",
    "# of reads mapping to mirbase (hairpin)",
    "% of reads mapping to mirbase (hairpin)",
    "# of reads mapped genome",
    "% of reads mapped genome",
    "# of reads multimapping genome",
    "% of reads multimapping genome",
    "# of reads after deduplication",
    "% of reads after deduplication",
    "# of reads assigned in featurecounts",
    "% of reads assigned in featurecounts"
], columns=sample_list)

# ------------------------------
# Fill UMI rows
# ------------------------------
metrics.loc["# of reads"] = [umi_data[s]["reads_before"] for s in sample_list]
metrics.loc["# of UMI-containing reads"] = [umi_data[s]["reads_after"] for s in sample_list]
metrics.loc["% of UMI-containing reads"] = (
    metrics.loc["# of UMI-containing reads"].astype(float)
    / metrics.loc["# of reads"].astype(float) * 100
).round(2)

# ------------------------------
# miRBase counts (mature)
# ------------------------------
metrics.loc["# of reads mapping to mirbase (mature)"] = [mir_data[s]["mapped_mature"] for s in sample_list]
metrics.loc["% of reads mapping to mirbase (mature)"] = (
    metrics.loc["# of reads mapping to mirbase (mature)"].astype(float)
    / metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)

metrics.loc["# of reads multimapping mirbase (mature)"] = [mir_data[s]["multimap_mature"] for s in sample_list]
metrics.loc["% of reads multimapping mirbase (mature)"] = (
    metrics.loc["# of reads multimapping mirbase (mature)"].astype(float) 
    / metrics.loc["# of UMI-containing reads"].astype(float) * 100).round(2)

filtered_bams = sorted(glob(str(mir_mature_dir / "*.srt.bam")))
metrics.loc["# of reads mapping to miRBase after length/gap/MAPQ/mismatch filters (mature)"] = [
    run_samtools_count(bam, exclude_flag=4) for bam in filtered_bams
]

# ------------------------------
# miRBase counts (hairpin)
# ------------------------------
metrics.loc["# of reads mapping to mirbase (hairpin)"] = [hairpin_data[s]["mapped"] for s in sample_list]
metrics.loc["% of reads mapping to mirbase (hairpin)"] = (
    metrics.loc["# of reads mapping to mirbase (hairpin)"].astype(float)
    / metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)

# ------------------------------
# Genome counts
# ------------------------------
# genome BAMs 
genome_bams = sorted(glob(str(genome_dir / "*.srt.bam")))
metrics.loc["# of reads mapped genome"] = [
    run_samtools_count(bam, exclude_flag=4) for bam in genome_bams
]
metrics.loc["# of reads multimapping genome"] = [genome_data[s]["multimap"] for s in sample_list]
metrics.loc["% of reads multimapping genome"] = (
    metrics.loc["# of reads multimapping genome"].astype(float) 
    / metrics.loc["# of UMI-containing reads"].astype(float) * 100).round(2)
# deduplicated BAMs 
dedup_bams = sorted(glob(str(genome_dir / "*.srt.dedup.bam")))
metrics.loc["# of reads after deduplication"] = [
    run_samtools_count(bam, exclude_flag=4) for bam in dedup_bams
]

# ------------------------------
# FeatureCounts assignment
# ------------------------------
# load featurecounts summary file 
fc = pd.read_csv("genome_counts/featurecounts.readcounts.raw.tsv.summary",
                 sep="\t", index_col=0)
# assign valus to metrics 
metrics.loc["# of reads assigned in featurecounts"] = fc.loc["Assigned"].tolist()
metrics.loc["% of reads assigned in featurecounts"] = (
    metrics.loc["# of reads assigned in featurecounts"].astype(float)
    / metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)

# ------------------------------
# additional calculations 
# ------------------------------
# % of reads mapped genome
metrics.loc["% of reads mapped genome"] =  (
    metrics.loc["# of reads mapped genome"].astype(float)
    / metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)
# % of reads after deduplication
metrics.loc["% of reads after deduplication"] =  (
    metrics.loc["# of reads after deduplication"].astype(float)
    / metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)
# % of reads after filter (mature)
metrics.loc["% of reads mapping to miRBase after length/gap/MAPQ/mismatch filters (mature)"] =  (
    metrics.loc["# of reads mapping to miRBase after length/gap/MAPQ/mismatch filters (mature)"].astype(float)
    / metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)

# ------------------------------
# write output files 
# ------------------------------
metrics.to_csv( "metrics/mirna_genome_alignment_metrics.tsv", sep="\t", index=True)
metrics.to_excel("metrics/mirna_genome_alignment_metrics.xlsx", index=True)



