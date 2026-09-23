#!/usr/bin/env python3

import pandas as pd
import re
import sys

if len(sys.argv) < 4:
    print("Usage: python add_biotype_and_gene_name.py <featureCounts_file> <annotation_gtf> <output_file>")
    sys.exit(1)

fc_file = sys.argv[1]
gtf_file = sys.argv[2]
out_file = sys.argv[3]

# ---- 1. Read featureCounts output ----
# Skip commented lines starting with #
fc = pd.read_csv(fc_file, sep="\t", comment="#", low_memory=False)

# ---- 2. Parse gene_id, gene_biotype, and gene_name from GTF ----
gene_ids = []
gene_biotypes = []
gene_names = []

with open(gtf_file) as gtf:
    for line in gtf:
        if line.startswith("#"):
            continue
        cols = line.strip().split("\t")
        attr = cols[8]

        # Extract gene_id
        m_id = re.search(r'gene_id "([^"]+)"', attr)
        gene_id = m_id.group(1) if m_id else None

        # Extract gene_biotype (blank if missing)
        m_bt = re.search(r'gene_biotype "([^"]+)"', attr)
        gene_biotype = m_bt.group(1) if m_bt else ""

        # Extract gene_name (blank if missing)
        m_name = re.search(r'gene_name "([^"]+)"', attr)
        gene_name = m_name.group(1) if m_name else ""

        if gene_id:
            gene_ids.append(gene_id)
            gene_biotypes.append(gene_biotype)
            gene_names.append(gene_name)

# Create dataframe with unique gene_id -> gene_biotype, gene_name
gtf_df = pd.DataFrame({
    'gene_id': gene_ids,
    'gene_biotype': gene_biotypes,
    'gene_name': gene_names
})
gtf_df = gtf_df.drop_duplicates(subset='gene_id')

# ---- 3. Merge with featureCounts ----
fc = fc.merge(gtf_df, left_on='Geneid', right_on='gene_id', how='left')
fc.drop(columns=['gene_id'], inplace=True)

# Reorder columns: Geneid, gene_name, gene_biotype, rest...
cols = fc.columns.tolist()
cols.insert(1, cols.pop(cols.index('gene_name')))
cols.insert(2, cols.pop(cols.index('gene_biotype')))
fc = fc[cols]

# ---- 4. Write output ----
fc.to_csv(out_file, sep="\t", index=False)
print(f"Annotated counts saved to {out_file}")