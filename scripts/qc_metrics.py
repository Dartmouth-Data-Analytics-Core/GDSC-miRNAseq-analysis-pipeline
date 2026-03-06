#!/usr/bin/env python3
import sys
from glob import glob
import pandas as pd
import subprocess
from pathlib import Path

# ----- Functions -----
def run_samtools_count(bam, include_flag=None, exclude_flag=None):
    """Run samtools view from Python and return count as int."""
    cmd = ["samtools", "view", bam, "-c"]
    if include_flag is not None:
        cmd.extend(["-f", str(include_flag)])
    if exclude_flag is not None:
        cmd.extend(["-F", str(exclude_flag)])
    result = subprocess.run(cmd, text=True, check=True, capture_output=True)
    return int(result.stdout.strip())

def parse_bowtie_log(logfile):
    """Parse bowtie-style log and return mapped, unmapped, multimap reads."""
    mapped = multimap = unaligned = None
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
    return total_reads - unaligned, unaligned, multimap

def parse_umi_log(logfile):
    """Parse umi extraction log file."""
    reads_before = reads_after = None
    with open(logfile) as f:
        for line in f:
            if "INFO Input Reads:" in line:
                reads_before = int(line.strip().split()[-1])
            elif "INFO Reads output:" in line:
                reads_after = int(line.strip().split()[-1])
    if reads_before is None or reads_after is None:
        raise ValueError(f"Missing info in UMI log {logfile}")
    reads_removed = reads_before - reads_after
    return reads_before, reads_after, reads_removed

# ----- Directories from command-line -----
umi_dir = Path(sys.argv[1])        # logs/umi_reads
mir_dir = Path(sys.argv[2])        # logs/mirbase_padded_aln
filt_mir_dir = Path(sys.argv[3])   # filtered miRBase BAMs
genome_dir = Path(sys.argv[4])     # genome alignment logs

# ----- Collect sample IDs -----
sample_list = sorted([
    Path(f).stem.replace(".umi", "").replace(".bowtie2", "").replace(".mature", "").replace(".genome", "")
    for f in glob(str(umi_dir / "*.log"))
])

# ----- Parse UMI logs -----
umi_data = {}
for sample in sample_list:
    log_file = umi_dir / f"{sample}.umi.log"
    reads_before, reads_after, reads_removed = parse_umi_log(log_file)
    umi_data[sample] = {
        "reads_before": reads_before,
        "reads_after": reads_after,
        "reads_removed": reads_removed
    }

# ----- Parse miRBase logs -----
mir_data = {}
for sample in sample_list:
    log_file = mir_dir / f"{sample}.bowtie2.mature.log"
    mapped, unaligned, multimap = parse_bowtie_log(log_file)
    mir_data[sample] = {"mapped": mapped, "multimap": multimap}

# ----- Parse genome logs -----
genome_data = {}
for sample in sample_list:
    log_file = genome_dir / f"{sample}.bowtie2.genome.log"
    mapped, unaligned, multimap = parse_bowtie_log(log_file)
    genome_data[sample] = {"mapped": mapped, "multimap": multimap}

# ----- Build metrics DataFrame -----
metrics = pd.DataFrame(index=[
    "# of reads",
    "# of UMI-containing reads",
    "% of UMI-containing reads",
    "# of reads mapping to mirbase",
    "% of reads mapping to mirbase",
    "# of reads multimapping mirbase",
    "% of reads multimapping mirbase",
    "# of reads mapping to miRBase after filters",
    "% of reads mapping to miRBase after filters",
    "# of reads mapped genome",
    "% of reads mapped genome",
    "# of reads multimapping genome",
    "% of reads multimapping genome",
    "# of reads after deduplication",
    "% of reads after deduplication",
    "# of reads assigned in featurecounts",
    "% of reads assigned in featurecounts"
], columns=sample_list)

# ----- Fill UMI counts -----
metrics.loc["# of reads"] = [umi_data[s]["reads_before"] for s in sample_list]
metrics.loc["# of UMI-containing reads"] = [umi_data[s]["reads_after"] for s in sample_list]
metrics.loc["% of UMI-containing reads"] = (
    metrics.loc["# of UMI-containing reads"].astype(float) /
    metrics.loc["# of reads"].astype(float) * 100
).round(2)

# ----- miRBase counts -----
metrics.loc["# of reads mapping to mirbase"] = [mir_data[s]["mapped"] for s in sample_list]
metrics.loc["% of reads mapping to mirbase"] = (
    metrics.loc["# of reads mapping to mirbase"].astype(float) /
    metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)

metrics.loc["# of reads multimapping mirbase"] = [mir_data[s]["multimap"] for s in sample_list]
metrics.loc["% of reads multimapping mirbase"] = (
    metrics.loc["# of reads multimapping mirbase"].astype(float) /
    metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)

# ----- Filtered miRBase BAMs (updated) -----
filtered_bams = sorted(glob(str(filt_mir_dir / "*.mature.srt.bam")))
metrics.loc["# of reads mapping to miRBase after filters"] = [
    run_samtools_count(bam, exclude_flag=4) for bam in filtered_bams
]
metrics.loc["% of reads mapping to miRBase after filters"] = (
    metrics.loc["# of reads mapping to miRBase after filters"].astype(float) /
    metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)

# ----- Genome BAMs -----
genome_bams = sorted(glob(str(genome_dir / "*.genome.srt.filt.bam")))
metrics.loc["# of reads mapped genome"] = [
    run_samtools_count(bam, exclude_flag=4) for bam in genome_bams
]
metrics.loc["# of reads multimapping genome"] = [genome_data[s]["multimap"] for s in sample_list]
metrics.loc["% of reads mapped genome"] = (
    metrics.loc["# of reads mapped genome"].astype(float) /
    metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)
metrics.loc["% of reads multimapping genome"] = (
    metrics.loc["# of reads multimapping genome"].astype(float) /
    metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)

# ----- Deduplicated BAMs -----
dedup_bams = sorted(glob(str(genome_dir / "*.genome.srt.filt.dedup.bam")))
metrics.loc["# of reads after deduplication"] = [
    run_samtools_count(bam, exclude_flag=4) for bam in dedup_bams
]
metrics.loc["% of reads after deduplication"] = (
    metrics.loc["# of reads after deduplication"].astype(float) /
    metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)

# ----- FeatureCounts -----
fc = pd.read_csv("genome_counts/featurecounts.tsv.summary", sep="\t", index_col=0)

# Extract sample names from BAM paths
fc.columns = [Path(c).name.replace(".genome.srt.filt.bam", "") for c in fc.columns]

# Assign counts by matching sample names
metrics.loc["# of reads assigned in featurecounts"] = [
    fc.loc["Assigned"].get(sample, 0) for sample in sample_list
]

metrics.loc["% of reads assigned in featurecounts"] = (
    metrics.loc["# of reads assigned in featurecounts"].astype(float) /
    metrics.loc["# of UMI-containing reads"].astype(float) * 100
).round(2)

# ----- Write outputs -----
Path("metrics").mkdir(exist_ok=True)
metrics.to_csv("metrics/mirna_genome_alignment_metrics.tsv", sep="\t", index=True)
metrics.to_excel("metrics/mirna_genome_alignment_metrics.xlsx", index=True)