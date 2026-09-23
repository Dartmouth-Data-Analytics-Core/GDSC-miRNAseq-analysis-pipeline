#!/usr/bin/env python3
"""
smRNA_biotype_metrics.py

Compute per-sample small RNA gene biotype metrics from featureCounts output.

Usage:
    python gene_biotype_metrics.py genome_counts/featurecounts.readcounts.biotype.tsv
"""

import sys
from pathlib import Path
import pandas as pd

# ----- Input -----
if len(sys.argv) < 2:
    print("Usage: python smRNA_biotype_metrics.py <featureCounts_biotype_file>")
    sys.exit(1)

fc_file = Path(sys.argv[1])
if not fc_file.exists():
    print(f"Error: {fc_file} not found.")
    sys.exit(1)

# ----- Read featureCounts annotated table -----
# Columns: Geneid, gene_name, gene_biotype, Chr, Start, End, Strand, Length, sample1, sample2, ...
fc_biotype = pd.read_csv(fc_file, sep="\t")

# Identify sample columns (counts start at column 9, 0-indexed)
sample_cols = fc_biotype.columns[8:]

# ----- Filter biotypes to only those ending with "RNA" or "ribozyme" -----
is_smRNA = fc_biotype['gene_biotype'].str.endswith("RNA") | fc_biotype['gene_biotype'].str.endswith("ribozyme")
fc_biotype_filtered = fc_biotype[is_smRNA]

# ----- Sum counts per gene_biotype -----
biotype_sum = fc_biotype_filtered.groupby("gene_biotype")[sample_cols].sum()

# ----- Calculate percentages per sample -----
biotype_pct = (biotype_sum / biotype_sum.sum(axis=0)) * 100
biotype_pct = biotype_pct.round(2)

# ----- Combine counts and percentages -----
multiindex_metrics = pd.concat(
    [biotype_sum, biotype_pct],
    keys=["counts", "percent"],
    names=["metric", "gene_biotype"]
)

# ----- Create output folder -----
Path("metrics").mkdir(exist_ok=True)

# ----- Write standard outputs -----
multiindex_metrics.to_csv("metrics/smRNA_biotype_metrics.tsv", sep="\t")
multiindex_metrics.to_excel("metrics/smRNA_biotype_metrics.xlsx")

print("smRNA biotype metrics written to metrics/smRNA_biotype_metrics.tsv and .xlsx")

# ----- Generate MultiQC custom content file -----
# Transpose data so samples are rows
mqc_data = pd.DataFrame()

# Add count columns
for biotype in biotype_sum.index:
    mqc_data[f"{biotype}_counts"] = biotype_sum.loc[biotype]

# Add percentage columns
for biotype in biotype_pct.index:
    mqc_data[f"{biotype}_percent"] = biotype_pct.loc[biotype]

# Reset index to make sample names a column
mqc_data.index.name = "Sample"
mqc_data = mqc_data.reset_index()

# Generate header configuration
biotypes_list = sorted(biotype_sum.index.tolist())

header_lines = [
    "# id: 'smrna_biotype_metrics'",
    "# section_name: 'Small RNA Biotype Statistics'",
    "# description: 'Distribution of small RNA reads across different gene biotypes'",
    "# plot_type: 'table'",
    "# pconfig:",
    "#     id: 'smrna_biotype_table'",
    "#     title: 'Small RNA Gene Biotype Counts and Percentages'",
    "#     save_file: true",
    "#     col1_header: 'Sample'",
    "# headers:",
]

# Add count headers
for biotype in biotypes_list:
    header_lines.extend([
        f"#     {biotype}_counts:",
        f"#         title: '{biotype}'",
        f"#         description: '{biotype} counts'",
        "#         format: '{:,.0f}'",
        "#         scale: false",
    ])

# Add percentage headers with color scale
for biotype in biotypes_list:
    header_lines.extend([
        f"#     {biotype}_percent:",
        f"#         title: '{biotype} %'",
        f"#         description: '{biotype} percentage'",
        "#         suffix: '%'",
        "#         format: '{:.2f}'",
        "#         scale: 'RdYlGn-rev'",
        "#         max: 100",
    ])

# Write MultiQC file
mqc_file = "metrics/smRNA_biotype_metrics_mqc.tsv"
with open(mqc_file, 'w') as f:
    # Write header configuration
    f.write('\n'.join(header_lines) + '\n')
    
    # Write data table
    mqc_data.to_csv(f, sep='\t', index=False)

print(f"MultiQC custom content written to {mqc_file}")
print("\nTo include in your MultiQC report, run MultiQC in the directory containing this file:")
print("  multiqc .")