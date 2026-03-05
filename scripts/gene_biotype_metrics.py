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

# ----- Write outputs -----
multiindex_metrics.to_csv("metrics/smRNA_biotype_metrics.tsv", sep="\t")
multiindex_metrics.to_excel("metrics/smRNA_biotype_metrics.xlsx")

print("smRNA biotype metrics written to metrics/smRNA_biotype_metrics.tsv and .xlsx")