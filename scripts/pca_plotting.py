#!/usr/bin/env python3
"""
smRNA PCA & Heatmap Analysis

This script performs:
1. Median-of-ratios normalization on raw read counts.
2. HVG (highly variable genes/miRNAs) detection using variance plateau.
3. PCA plots on all miRNAs and top HVGs.
4. Clustered heatmap of top HVGs.

Input: CSV/TSV file with first column `miRNA` and remaining columns as samples.
Output: Plots saved to the specified output folder.
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.optimize import curve_fit
from plotnine import *
import seaborn as sns
import matplotlib.pyplot as plt

# ----------------------------
# Argument Parsing
# ----------------------------
parser = argparse.ArgumentParser(description="smRNA PCA and heatmap analysis")
parser.add_argument("tsv_file", help="Path to CSV/TSV read counts file (miRNA in first column)")
parser.add_argument("output_path", help="Folder to save plots")
parser.add_argument("-p", "--pca_comp", type=int, default=10, help="Number of PCA components")
args = parser.parse_args()

os.makedirs(args.output_path, exist_ok=True)

# ----------------------------
# Load Data
# ----------------------------
try:
    df = pd.read_csv(args.tsv_file, sep=None, engine="python", index_col=0)
except Exception as e:
    print(f"Failed to load input file: {e}")
    sys.exit(1)

# Clean sample names
sample_names = [s.split("/")[-1].replace(".srt.bam", "") for s in df.columns.tolist()]
df.columns = sample_names

# ----------------------------
# Median-of-Ratios Normalization
# ----------------------------
def median_of_ratios(df_counts):
    df_counts = df_counts.astype(float)
    df_nozero = df_counts[(df_counts != 0).all(axis=1)]
    if df_nozero.shape[0] == 0:
        return df_counts
    gene_geom = np.exp(np.log(df_nozero).mean(axis=1))
    size_factors = np.median(df_nozero.div(gene_geom, axis=0), axis=0)
    return df_counts.div(size_factors, axis=1)

df_norm = median_of_ratios(df)
df_log = np.log2(df_norm + 1)

# ----------------------------
# Variance & HVG Detection
# ----------------------------
gene_var = df_log.var(axis=1).sort_values(ascending=False)
gene_var_df = pd.DataFrame({"Rank": np.arange(1, len(gene_var)+1), "Variance": gene_var.values})

def decay(x, a, b, c):
    return a * np.exp(-b * x) + c

xdata = np.arange(len(gene_var))
ydata = gene_var.values

try:
    params, _ = curve_fit(decay, xdata, ydata, p0=[ydata.max()-ydata.min(), 0.001, ydata.min()])
    a, b, c = params
except Exception:
    c = ydata.min()
    print("Curve fitting failed, using minimum variance as plateau")

tol = 0.01 * (ydata.max() - ydata.min())
plateau_indices = np.where(ydata <= c + tol)[0]
plateau_idx = plateau_indices[0]+1 if len(plateau_indices) > 0 else len(ydata)
print(f"Auto-detected plateau at rank {plateau_idx}, variance ~ {ydata[plateau_idx-1]:.4f}")

# Save HVG log
log_file = os.path.join(args.output_path, "pca_hvg_log.txt")
with open(log_file, "w") as log:
    log.write("PCA automatic HVG thresholding log\n")
    log.write("---------------------------------\n")
    log.write(f"Total miRNAs: {df_log.shape[0]}\n")
    log.write(f"Selected miRNAs (plateau): {plateau_idx}\n")
    log.write(f"Plateau variance: {ydata[plateau_idx-1]:.6f}\n")

# ----------------------------
# Variance Plot
# ----------------------------
p1 = (
    ggplot(gene_var_df, aes(x="Rank", y="Variance")) +
    geom_point(color="red", size=2.5) +
    geom_vline(xintercept=plateau_idx, linetype="dashed", color="blue", size=1) +
    theme_bw(base_size=16) +
    labs(title="miRNA Expression Variance over Samples",
         x="miRNA Rank",
         y="Variance")
)
p1.save(f"{args.output_path}/miRNA_Variance_Plot.png", width=8, height=6, dpi=300)

# ----------------------------
# Select Top miRNAs
# ----------------------------
top_genes = gene_var.index[:plateau_idx]
df_top = df_log.loc[top_genes]

# ----------------------------
# Scaling
# ----------------------------
X_scaled_all = StandardScaler().fit_transform(df_log.T)
X_scaled_top = StandardScaler().fit_transform(df_top.T)

# ----------------------------
# PCA Function
# ----------------------------
def run_pca(X_scaled, sample_names, prefix, ncomp=args.pca_comp):
    ncomp = min(ncomp, X_scaled.shape[0], X_scaled.shape[1])
    pca = PCA(n_components=ncomp)
    pcs = pca.fit_transform(X_scaled)

    pcadf = pd.DataFrame(pcs, columns=[f"PC{i+1}" for i in range(ncomp)])
    pcadf["sample"] = sample_names

    def plot_pc_pair(pc_x, pc_y):
        x = f"PC{pc_x}"
        y = f"PC{pc_y}"
        pcadf["x_offset"] = pcadf[x] + 0.02*(pcadf[x].max() - pcadf[x].min())
        pcadf["y_offset"] = pcadf[y] + 0.02*(pcadf[y].max() - pcadf[y].min())
        p = (
            ggplot(pcadf, aes(x=x, y=y)) +
            geom_point(color="black", size=4.5, fill="red") +
            geom_text(aes(x="x_offset", y="y_offset", label="sample"),
                      ha="left", va="bottom", size=8) +
            theme_bw(base_size=16) +
            labs(x=f"{x}: {pca.explained_variance_ratio_[pc_x-1]*100:.2f}%",
                 y=f"{y}: {pca.explained_variance_ratio_[pc_y-1]*100:.2f}%")
        )
        p.save(f"{args.output_path}/{prefix}_{x}_vs_{y}.png", width=8, height=6, dpi=300)

    plot_pc_pair(1, 2)
    if ncomp >= 3: plot_pc_pair(2, 3)
    if ncomp >= 4: plot_pc_pair(3, 4)

    # Variance bar plot
    pca_var_df = pd.DataFrame({
        "PC": [f"PC{i+1}" for i in range(ncomp)],
        "Variance": pca.explained_variance_ratio_ * 100
    })
    pca_bar = (
        ggplot(pca_var_df, aes(x="PC", y="Variance")) +
        geom_bar(stat="identity", fill="skyblue", color="black") +
        theme_bw(base_size=16) +
        labs(x="Principal Component", y="Percent Variance Explained")
    )
    pca_bar.save(f"{args.output_path}/{prefix}_PCA_variance_bar.png", width=8, height=6, dpi=300)

# ----------------------------
# Run PCA
# ----------------------------
run_pca(X_scaled_all, sample_names, prefix="PCA_all")
run_pca(X_scaled_top, sample_names, prefix="PCA_top")

# ----------------------------
# Heatmap of Top miRNAs
# ----------------------------
row_mean = df_top.mean(axis=1)
row_std = df_top.std(axis=1)
valid_rows = row_std > 0
df_top_filtered = df_top.loc[valid_rows]
row_mean = row_mean[valid_rows]
row_std = row_std[valid_rows]

df_top_scaled = (df_top_filtered - row_mean.values[:, None]) / row_std.values[:, None]
df_top_scaled = df_top_scaled.dropna(axis=0, how='any').dropna(axis=1, how='any')
col_variance = df_top_scaled.var(axis=0)
df_top_scaled = df_top_scaled.loc[:, col_variance > 0]

if df_top_scaled.shape[0] >= 2 and df_top_scaled.shape[1] >= 2:
    sns_clustermap = sns.clustermap(
        df_top_scaled,
        cmap="vlag",
        method="average",
        metric="euclidean",
        col_cluster=True,
        row_cluster=True,
        yticklabels=True,
        figsize=(12, 14)
    )
    sns_clustermap.fig.suptitle(f"Heatmap of Top {len(df_top_scaled)} miRNAs", y=1.05)
    sns_clustermap.savefig(f"{args.output_path}/Top_miRNAs_Heatmap.png", dpi=300, bbox_inches="tight")
    plt.close()
else:
    fig, ax = plt.subplots(figsize=(12, 14))
    sns.heatmap(df_top_scaled if not df_top_scaled.empty else df_top.fillna(0),
                cmap="vlag",
                yticklabels=True,
                ax=ax)
    ax.set_title(f"Heatmap of Top {len(df_top)} miRNAs (unclustered)")
    fig.savefig(f"{args.output_path}/Top_miRNAs_Heatmap.png", dpi=300, bbox_inches="tight")
    plt.close()