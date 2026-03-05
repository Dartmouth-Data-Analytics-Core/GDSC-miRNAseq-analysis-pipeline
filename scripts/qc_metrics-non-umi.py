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
    """Parse bowtie-style log and return total, mapped, unmapped, multimap reads."""
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

# ----- Directories from command-line -----
mir_dir = Path(sys.argv[1])       # mirbase alignment logs
filt_mir_dir = Path(sys.argv[2])  # filtered mirbase BAMs
genome_dir = Path(sys.argv[3])    # genome alignment logs

# ----- Collect sample IDs from miRBase logs -----
sample_list = sorted([
    Path(f).stem.replace(".bowtie2.mature", "") 
    for f in glob(str(mir_dir / "*.bowtie2.mature.log"))
])

# ----- Parse miRBase logs -----
mir_data = {}
for sample in sample_list:
    log_file = mir_dir / f"{sample}.bowtie2.mature.log"
    total, mapped, unaligned, multimap = parse_bowtie_log(log_file)
    mir_data[sample] = {"reads": total, "mapped": mapped, "unaligned": unaligned, "multimap": multimap}

# ----- Parse genome logs -----
genome_data = {}
for sample in sample_list:
    log_file = genome_dir / f"{sample}.bowtie2.genome.log"
    total, mapped, unaligned, multimap = parse_bowtie_log(log_file)
    genome_data[sample] = {"reads": total, "mapped": mapped, "unaligned": unaligned, "multimap": multimap}

# ----- Build metrics DataFrame -----
metrics = pd.DataFrame(index=[
    "# of reads",
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
    "# of reads assigned in featurecounts",
    "% of reads assigned in featurecounts"
], columns=sample_list)

# ----- miRBase counts -----
metrics.loc["# of reads"] = [mir_data[s]["reads"] for s in sample_list]
metrics.loc["# of reads mapping to mirbase"] = [mir_data[s]["mapped"] for s in sample_list]
metrics.loc["% of reads mapping to mirbase"] = (
    metrics.loc["# of reads mapping to mirbase"].astype(float) /
    metrics.loc["# of reads"].astype(float) * 100
).round(2)

metrics.loc["# of reads multimapping mirbase"] = [mir_data[s]["multimap"] for s in sample_list]
metrics.loc["% of reads multimapping mirbase"] = (
    metrics.loc["# of reads multimapping mirbase"].astype(float) /
    metrics.loc["# of reads"].astype(float) * 100
).round(2)

# ----- Filtered miRBase BAMs -----
filtered_bams = sorted(glob(str(filt_mir_dir / "*.mature.srt.bam")))
bam_dict = {Path(bam).stem.replace(".mature.srt",""): bam for bam in filtered_bams}

metrics.loc["# of reads mapping to miRBase after filters"] = [
    run_samtools_count(bam_dict[s], exclude_flag=4) if s in bam_dict else 0
    for s in sample_list
]
metrics.loc["% of reads mapping to miRBase after filters"] = (
    metrics.loc["# of reads mapping to miRBase after filters"].astype(float) /
    metrics.loc["# of reads"].astype(float) * 100
).round(2)

# ----- Genome counts -----
genome_bams = sorted(glob(str(genome_dir / "*.genome.srt.filt.bam")))
bam_dict_genome = {Path(bam).stem.replace(".genome.srt.filt",""): bam for bam in genome_bams}

metrics.loc["# of reads mapped genome"] = [
    run_samtools_count(bam_dict_genome[s], exclude_flag=4) if s in bam_dict_genome else 0
    for s in sample_list
]
metrics.loc["# of reads multimapping genome"] = [genome_data[s]["multimap"] for s in sample_list]
metrics.loc["% of reads mapped genome"] = (
    metrics.loc["# of reads mapped genome"].astype(float) /
    metrics.loc["# of reads"].astype(float) * 100
).round(2)
metrics.loc["% of reads multimapping genome"] = (
    metrics.loc["# of reads multimapping genome"].astype(float) /
    metrics.loc["# of reads"].astype(float) * 100
).round(2)

# ----- FeatureCounts assignment -----
fc = pd.read_csv("genome_counts/featurecounts.tsv.summary", sep="\t", index_col=0)
metrics.loc["# of reads assigned in featurecounts"] = [
    fc.loc["Assigned", s] if s in fc.columns else 0 for s in sample_list
]
metrics.loc["% of reads assigned in featurecounts"] = (
    metrics.loc["# of reads assigned in featurecounts"].astype(float) /
    metrics.loc["# of reads"].astype(float) * 100
).round(2)

# ----- Write output -----
Path("metrics").mkdir(exist_ok=True)
metrics.to_csv("metrics/mirna_genome_alignment_metrics.tsv", sep="\t", index=True)
metrics.to_excel("metrics/mirna_genome_alignment_metrics.xlsx", index=True)