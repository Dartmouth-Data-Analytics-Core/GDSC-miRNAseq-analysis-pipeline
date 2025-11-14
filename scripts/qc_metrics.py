#!/usr/bin/env python3
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


def parse_bowtie_log(logfile):
    """Parse bowtie-style log and return mapped, unmapped."""
    reads = None
    unaligned = None
    
    with open(logfile) as f:
        for line in f:
            if "reads; of these" in line:
                reads = int(line.strip().split()[0])
            if "aligned 0 times" in line:
                unaligned = int(line.strip().split()[0])

    if reads is None or unaligned is None:
        raise ValueError(f"Log file {logfile} missing expected fields.")

    return reads - unaligned, unaligned



umi_dir = sys.argv[1]
mir_dir = sys.argv[2]
genome_dir = sys.argv[3]



# ------------------------------
# Parse UMI log files
# ------------------------------
umi_data = {}
sample_list = []

for logfile in sorted(glob(umi_dir+"/*.log.txt")):
    sample_id = logfile.split(".")[0].split("/")[-1]
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
# Parse miRBase mapping logs
# ------------------------------

mirmap_data_dict = {}
mir_data = {}

for logfile in sorted(glob(mir_dir+"/*mirbase.log.txt")):
    sample_id = Path(logfile).name.replace("_mirbase.log.txt", "")
    mapped, unmapped = parse_bowtie_log(logfile)
    mir_data[sample_id] = {"mapped": mapped, "unmapped": unmapped}

# ------------------------------
# Parse genome mapping logs
# ------------------------------

genome_data = {}
for logfile in sorted(glob(genome_dir+"/*genome.log.txt")):
    sample_id = Path(logfile).name.replace("_genome.log.txt", "")
    mapped, unmapped = parse_bowtie_log(logfile)
    genome_data[sample_id] = {"mapped": mapped, "unmapped": unmapped}



# ------------------------------
# Build metrics DataFrame
# ------------------------------


sample_list = sorted(sample_list)
metrics = pd.DataFrame(index=[
    "# of reads",
    "# of reads missing UMI",
    "% of reads missing UMI",
    "# of UMI-containing reads",
    "% of UMI-containing reads",
    "# of reads mapping to mirbase",
    "% of reads mapping to mirbase",
    "# of reads unmapped to mirbase",
    "% of reads unmapped to mirbase",
    "# of reads mapped genome",
    "% of reads mapped genome",
    "# of reads unmapped genome",
    "% of reads unmapped genome",
    "# of reads after deduplication",
    "% of reads after deduplication",
    "# of reads after filter",
    "% of reads after filter",
    "# of reads assigned in feature counts",
    "% of reads assigned in feature counts"
], columns=sample_list)


# ------------------------------
# Fill UMI rows
# ------------------------------

metrics.loc["# of reads"] = [umi_data[s]["reads_before"] for s in sample_list]
metrics.loc["# of reads missing UMI"] = [umi_data[s]["reads_removed"] for s in sample_list]
metrics.loc["# of UMI-containing reads"] = [umi_data[s]["reads_after"] for s in sample_list]

metrics.loc["% of UMI-containing reads"] = (
    metrics.loc["# of UMI-containing reads"].astype(float)
    / metrics.loc["# of reads"].astype(float) * 100
).round(2)

metrics.loc["% of reads missing UMI"] = (
    metrics.loc["# of reads missing UMI"].astype(float)
    / metrics.loc["# of reads"].astype(float) * 100
).round(2)


# ------------------------------
# miRBase counts
# ------------------------------



metrics.loc["# of reads mapping to mirbase"] = [mir_data[s]["mapped"] for s in sample_list]
metrics.loc["% of reads mapping to mirbase"] = (
    metrics.loc["# of reads mapping to mirbase"].astype(float)
    / metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)

# Count unmapped from BAM
mir_bams = sorted(glob(mir_dir+"/*.aln.sam"))
metrics.loc["# of reads unmapped to mirbase"] = [
    run_samtools_count(bam, include_flag=4) for bam in mir_bams
]



# ------------------------------
# Genome counts
# ------------------------------

genome_bams = sorted(glob(genome_dir+"/*.srt.bam"))
metrics.loc["# of reads mapped genome"] = [
    run_samtools_count(bam, exclude_flag=4) for bam in genome_bams
]
metrics.loc["# of reads unmapped genome"] = [
    run_samtools_count(bam, include_flag=4) for bam in genome_bams
]

dedup_bams = sorted(glob(genome_dir+"/*.srt.dedup.bam"))
metrics.loc["# of reads after deduplication"] = [
    run_samtools_count(bam, exclude_flag=4) for bam in dedup_bams
]

filtered_bams = sorted(glob(genome_dir+"/*.srt.dedup.filt.bam"))
metrics.loc["# of reads after filter"] = [
    run_samtools_count(bam, exclude_flag=4) for bam in filtered_bams
]


# ------------------------------
# FeatureCounts assignment
# ------------------------------

fc = pd.read_csv("genome_counts/featurecounts.readcounts.raw.tsv.summary",
                 sep="\t", index_col=0)

metrics.loc["# of reads assigned in feature counts"] = fc.loc["Assigned"].tolist()

metrics.loc["% of reads assigned in feature counts"] = (
    metrics.loc["# of reads assigned in feature counts"].astype(float)
    / metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)




# ------------------------------
# % Calculations
# ------------------------------
# % of reads unmapped to mirbase
metrics.loc["% of reads unmapped to mirbase"] =  (
    metrics.loc["# of reads unmapped to mirbase"].astype(float)
    / metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)
# % of reads mapped genome
metrics.loc["% of reads mapped genome"] =  (
    metrics.loc["# of reads mapped genome"].astype(float)
    / metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)
# % of reads unmapped genome
metrics.loc["% of reads unmapped genome"] =  (
    metrics.loc["# of reads unmapped genome"].astype(float)
    / metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)
# % of reads after deduplication
metrics.loc["% of reads after deduplication"] =  (
    metrics.loc["# of reads after deduplication"].astype(float)
    / metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)
# % of reads after filter
metrics.loc["% of reads after filter"] =  (
    metrics.loc["# of reads after filter"].astype(float)
    / metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)

# ------------------------------
# Output
# ------------------------------

metrics.to_csv( "metrics/mirna_genome_alignment_metrics", sep="\t", index=True)
metrics.to_excel("metrics/mirna_genome_alignment_metrics.xlsx", index=True)



