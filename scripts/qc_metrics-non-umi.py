#!/usr/bin/env python3

import sys
from glob import glob
import pandas as pd
import subprocess
from pathlib import Path


# -------------------------------------------------------
# Functions
# -------------------------------------------------------

def run_samtools_count(bam, include_flag=None, exclude_flag=None):
    """Run samtools view and return read count."""
    cmd = ["samtools", "view", bam, "-c"]

    if include_flag is not None:
        cmd.extend(["-f", str(include_flag)])

    if exclude_flag is not None:
        cmd.extend(["-F", str(exclude_flag)])

    result = subprocess.run(cmd, text=True, check=True, capture_output=True)
    return int(result.stdout.strip())


def parse_bowtie2_log(logfile):
    """Parse Bowtie2 alignment log."""
    total_reads = mapped = unaligned = multimap = None

    with open(logfile) as f:
        for line in f:

            if "reads; of these" in line:
                total_reads = int(line.strip().split()[0])

            elif "aligned 0 times" in line:
                unaligned = int(line.strip().split()[0])

            elif "aligned >1 times" in line:
                multimap = int(line.strip().split()[0])

    if total_reads is None or unaligned is None or multimap is None:
        raise ValueError(f"Missing info in log {logfile}")

    mapped = total_reads - unaligned

    return total_reads, mapped, unaligned, multimap


def parse_bowtie1_log(logfile):
    """Parse Bowtie1 alignment log (SeqCluster hairpin step)."""
    total = mapped = unmapped = None

    with open(logfile) as f:
        for line in f:

            if "reads processed" in line:
                total = int(line.split(":")[1].strip())

            elif "reads with at least one alignment" in line:
                mapped = int(line.split(":")[1].split()[0])

            elif "failed to align" in line:
                unmapped = int(line.split(":")[1].split()[0])

    if total is None or mapped is None or unmapped is None:
        raise ValueError(f"Missing info in log {logfile}")

    return total, mapped, unmapped


# -------------------------------------------------------
# Input directories
# -------------------------------------------------------

mir_dir = Path(sys.argv[1])
filt_mir_dir = Path(sys.argv[2])
genome_dir = Path(sys.argv[3])
seqcluster_dir = Path(sys.argv[4])


# -------------------------------------------------------
# Sample IDs
# -------------------------------------------------------

sample_list = sorted([
    Path(f).stem.replace(".bowtie2.mature", "")
    for f in glob(str(mir_dir / "*.bowtie2.mature.log"))
])


# -------------------------------------------------------
# Parse miRBase logs
# -------------------------------------------------------

mir_data = {}

for sample in sample_list:

    log_file = mir_dir / f"{sample}.bowtie2.mature.log"

    total, mapped, unaligned, multimap = parse_bowtie2_log(log_file)

    mir_data[sample] = {
        "reads": total,
        "mapped": mapped,
        "unaligned": unaligned,
        "multimap": multimap
    }


# -------------------------------------------------------
# Parse genome logs
# -------------------------------------------------------

genome_data = {}

for sample in sample_list:

    log_file = genome_dir / f"{sample}.bowtie2.genome.log"

    total, mapped, unaligned, multimap = parse_bowtie2_log(log_file)

    genome_data[sample] = {
        "reads": total,
        "mapped": mapped,
        "unaligned": unaligned,
        "multimap": multimap
    }


# -------------------------------------------------------
# Parse SeqCluster hairpin logs (Bowtie1)
# -------------------------------------------------------

hairpin_data = {}

for sample in sample_list:

    log_file = seqcluster_dir / f"{sample}.bowtie1.hairpin.aln.log"

    if log_file.exists():

        total, mapped, unmapped = parse_bowtie1_log(log_file)

        hairpin_data[sample] = {
            "reads": total,
            "mapped": mapped,
            "unmapped": unmapped
        }

    else:

        hairpin_data[sample] = {
            "reads": 0,
            "mapped": 0,
            "unmapped": 0
        }


# -------------------------------------------------------
# Metrics table
# -------------------------------------------------------

metrics = pd.DataFrame(index=[

    "# of reads",
    "# of reads mapping to mirbase",
    "% of reads mapping to mirbase",
    "# of reads multimapping mirbase",
    "% of reads multimapping mirbase",
    "# of reads mapping to miRBase after filters",
    "% of reads mapping to miRBase after filters",
    "# of unique seqcluster clusters",
    "% of clustered sequences aligning to hairpins",
    "# of reads mapped genome",
    "% of reads mapped genome",
    "# of reads multimapping genome",
    "% of reads multimapping genome",
    "# of reads assigned in featurecounts",
    "% of reads assigned in featurecounts"

], columns=sample_list)


# -------------------------------------------------------
# miRBase metrics
# -------------------------------------------------------

metrics.loc["# of reads"] = [mir_data[s]["reads"] for s in sample_list]

metrics.loc["# of reads mapping to mirbase"] = [
    mir_data[s]["mapped"] for s in sample_list
]

metrics.loc["% of reads mapping to mirbase"] = (
    metrics.loc["# of reads mapping to mirbase"].astype(float) /
    metrics.loc["# of reads"].astype(float) * 100
).round(2)

metrics.loc["# of reads multimapping mirbase"] = [
    mir_data[s]["multimap"] for s in sample_list
]

metrics.loc["% of reads multimapping mirbase"] = (
    metrics.loc["# of reads multimapping mirbase"].astype(float) /
    metrics.loc["# of reads"].astype(float) * 100
).round(2)


# -------------------------------------------------------
# Filtered miRBase BAMs
# -------------------------------------------------------

filtered_bams = sorted(glob(str(filt_mir_dir / "*.mature.srt.bam")))

bam_dict = {
    Path(bam).stem.replace(".mature.srt", ""): bam
    for bam in filtered_bams
}

metrics.loc["# of reads mapping to miRBase after filters"] = [
    run_samtools_count(bam_dict[s], exclude_flag=4) if s in bam_dict else 0
    for s in sample_list
]

metrics.loc["% of reads mapping to miRBase after filters"] = (
    metrics.loc["# of reads mapping to miRBase after filters"].astype(float) /
    metrics.loc["# of reads"].astype(float) * 100
).round(2)


# -------------------------------------------------------
# SeqCluster hairpin metrics
# -------------------------------------------------------

metrics.loc["# of unique seqcluster clusters"] = [
    hairpin_data[s]["reads"] for s in sample_list
]

metrics.loc["% of clustered sequences aligning to hairpins"] = (
    pd.Series(
        [hairpin_data[s]["mapped"] for s in sample_list],
        index=sample_list
    ).astype(float)
    /
    metrics.loc["# of unique seqcluster clusters"].astype(float)
    * 100
).round(2)


# -------------------------------------------------------
# Genome metrics (from Bowtie2 logs)
# -------------------------------------------------------

metrics.loc["# of reads mapped genome"] = [
    genome_data[s]["mapped"] for s in sample_list
]

metrics.loc["# of reads multimapping genome"] = [
    genome_data[s]["multimap"] for s in sample_list
]

metrics.loc["% of reads mapped genome"] = (
    metrics.loc["# of reads mapped genome"].astype(float) /
    metrics.loc["# of reads"].astype(float) * 100
).round(2)

metrics.loc["% of reads multimapping genome"] = (
    metrics.loc["# of reads multimapping genome"].astype(float) /
    metrics.loc["# of reads"].astype(float) * 100
).round(2)


# -------------------------------------------------------
# FeatureCounts
# -------------------------------------------------------

fc = pd.read_csv(
    "genome_counts/featurecounts.tsv.summary",
    sep="\t",
    index_col=0
)

fc.columns = [
    Path(c).name.replace(".genome.srt.filt.bam", "")
    for c in fc.columns
]

metrics.loc["# of reads assigned in featurecounts"] = [
    fc.loc["Assigned"].get(sample, 0)
    for sample in sample_list
]

metrics.loc["% of reads assigned in featurecounts"] = (
    metrics.loc["# of reads assigned in featurecounts"].astype(float) /
    metrics.loc["# of reads"].astype(float) * 100
).round(2)


# -------------------------------------------------------
# Write output
# -------------------------------------------------------

Path("metrics").mkdir(exist_ok=True)

metrics.to_csv(
    "metrics/mirna_genome_alignment_metrics.tsv",
    sep="\t"
)

metrics.to_excel(
    "metrics/mirna_genome_alignment_metrics.xlsx"
)
